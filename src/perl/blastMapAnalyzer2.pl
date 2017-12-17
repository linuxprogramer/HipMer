#!/usr/bin/env perl
#
# blastMapAnalyzer.pl by Jarrod Chapman <jchapman@lbl.gov> Fri Jun 12 12:58:11 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use warnings;
use lib "$ENV{MERACULOUS_ROOT}/lib";
use M_Utility qw(parse_fastq_header);

use Getopt::Std;
my %opts = ();
my $validLine = getopts('b:i:m:I:M:RAS:E:L:F:T:U:v', \%opts);
my @required = ("b","m");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= (($nRequired == @required) && (exists($opts{"i"}) || exists($opts{"I"})));
if (exists($opts{"U"}) && ($opts{"U"} != 3) && ($opts{"U"} != 5)) {
    $validLine = 0;
}
if ($validLine != 1) {
    print "Usage: ./blastMapAnalyzer.pl <-b blastMap> <-m minMatch> ( <-i insertSize:sigma> || <-I testInsertSize:sigma> ) <<-R(everseComplement?)>> <<-A(rtifactRemoval?)>> <<-M minFreqReported>> <<-S scaffoldReportFile>> <<-E endAversionDistance>> <<-L minTestLength>> <<-F fivePrimeWiggleRoom>> <<-T threePrimeWiggleRoom>> <<-U trUncateAlignmentEnd(5|3)>> <<-v(erbose)\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\nblastMapAnalyzer2";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $blastMapFile = $opts{"b"};

my $reverseComplement = 0;
if (exists($opts{"R"})) {
    $reverseComplement = 1;
}

my $innieRemoval = 0;
if (exists($opts{"A"})) {
    $innieRemoval = 1;
}

my $endAversion = 0;
if (exists($opts{"E"})) {
    $endAversion = $opts{"E"};
}

my $minTestLength = 0;
if (exists($opts{"L"})) {
    $minTestLength = $opts{"L"};
}

my $srfFile = undef;
if (exists($opts{"S"})) {
    $srfFile = $opts{"S"};
}

my $minFreqReported = 10;
if (exists($opts{"M"})) {
    $minFreqReported = $opts{"M"};
}

my $fivePrimeWiggleRoom = 5;
if (exists($opts{"F"})) {
    $fivePrimeWiggleRoom = $opts{"F"};
}
my $threePrimeWiggleRoom = 5;
if (exists($opts{"T"})) {
    $threePrimeWiggleRoom = $opts{"T"};
}

my $truncate = 0;
if (exists($opts{"U"})) {
    $truncate = $opts{"U"};
}
my $minMatch = $opts{"m"};

my $verbose = $opts{"v"};

my $insertSize = undef;
my $endDistance = undef;
my $testMode = 0;
my $nFullPairs = 0;
my $nNonStandard = 0;
my %testLocations = ();
my $shortPair = 600;
my $estimatedInsertSize = 0;
my $estimatedSigma = 0;
my $sampleSize = 200000;

if (exists($opts{"I"})) {
    $testMode = 1;
    ($estimatedInsertSize,$estimatedSigma) = $opts{"I"} =~ /^(\d+)\:(\d+)$/;
} else {
    ($insertSize,$insertSigma) = $opts{"i"} =~ /^(\d+)\:(\d+)$/;
    $endDistance = $insertSize + 3*$insertSigma;
}

# If an existing scaffold report exists, load scaffolding information
 
my %contigScaffoldMap = ();  # contigName -> contigOrientation scaffoldName . scaffoldStart . scaffoldEnd 
my %scaffLengths = ();    # scaffoldName -> scaffoldLength
 
if (defined($srfFile)) {
    
    print STDERR "Reading scaffold report file: $srfFile...\n";
    open (S,$srfFile) || die "Couldn't open $srfFile\n";
    while (my $line = <S>) {
	chomp $line;
	my @cols = split(/\t/,$line);
	if ($cols[1] =~ /^CONTIG/) {
	    my ($scaffID,$contigID,$sStart,$sEnd) = @cols[0,2,3,4];
	    
	    if (exists($scaffLengths{$scaffID})) {
		if ($sEnd > $scaffLengths{$scaffID}) {
		    $scaffLengths{$scaffID} = $sEnd;
		}
	    } else {
		$scaffLengths{$scaffID} = $sEnd;
	    }
	    
	    my ($cStrand,$cName) = $contigID =~ /^([\+\-])(.+)$/;
	    
	    if (exists($contigScaffoldMap{$cName})) {
		die "Error in $srfFile: Current implementation limits contigs to a single scaffold location: $cName\n";
	    }
	    
	    $contigScaffoldMap{$cName} = "$cStrand$scaffID.$sStart.$sEnd";
	    
	}
    }
    close S;
    
    print STDERR "Done.\n";
}

my $currentPair = undef;
my $pairInfo = "";
my %pairSizes = ();
my %omitted_hits;
my $sampleCohortTested;

print STDERR "Reading blast map file: $blastMapFile...\n";
open (B,$blastMapFile) || die "Couldn't open $blastMapFile\n";

my $n =0;
while (my $line = <B>) {
    chomp $line;
    next if ( $line =~ /^BLAST_TYPE/ );

    my ($blastType,$query,$qStart,$qStop,$qLength,$subject,$sStart,$sStop,$sLength,
	$strand,$score,$eValue,$identities,$alignLength) = split(/\t/,$line);

    if ($verbose) {print STDERR "$line\n"}

    $n++;
    
    unless (($blastType eq "BLASTN") && ($identities >= $minMatch))  {
	$omitted_hits{"NOT A SUFFICIENT BLAST HIT"}++;
	if ($verbose) {print STDERR "HIT OMITTED: 0 NOT A SUFFICIENT BLAST HIT (identities: $identities  minMatch: $minMatch)\n"}
	next;
    }

    if ($truncate == 5) {
#	if (($qStart > $fivePrimeWiggleRoom) || ($qStop < $minMatch)) {

	# Strict alignment truncation.  Does not play nice with fivePrimeWiggleRoom option.
	# If necessary use the above condition instead.

	if (($qStart > 1) || ($qStop < $minMatch)) {
	    $omitted_hits{"5-TRUNCATED ALIGNMENT"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 1 TRUNCATED ALIGNMENT (qStart: $qStart qstop: $qStop  minMatch: $minMatch)\n"}
	    next;
	}
	my $trimOff = $qStop-$minMatch;
	$qStop = $minMatch;
	$qLength = $minMatch;
	if ($strand eq "Plus") {
	    $sStop -= $trimOff;
	} else {
	    $sStart += $trimOff;
	}
    } elsif ($truncate == 3) {
	# Strict alignment truncation.  Does not play nice with threePrimeWiggleRoom option.

	if (($qLength-$qStart+1 < $minMatch) || ($qStop < $qLength)) {
	    $omitted_hits{"3-TRUNCATED ALIGNMENT"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 1 TRUNCATED ALIGNMENT (qStart: $qStart qstop: $qStop  minMatch: $minMatch)\n"}
	    next;
	}
	my $trimOff = ($qLength-$minMatch+1)-$qStart;
	$qStart = 1;
	$qStop = $minMatch;
	$qLength = $minMatch;
	if ($strand eq "Plus") {
	    $sStart += $trimOff;
	} else {
	    $sStop -= $trimOff;
	}
    }

    if (defined($srfFile)) {
	unless (exists($contigScaffoldMap{$subject})) {
	    $omitted_hits{"NO SCAFFOLD INFO"}++;
#	    print STDERR "UNSCAFFOLDED CONTIG:$line\n";
	    if ($verbose) {print STDERR "HIT OMITTED: 2 NO SCAFFOLD INFO\n"}
	    next;
	} else {
	    my $mapInfo = $contigScaffoldMap{$subject};
	    my ($contigScaffStrand,$scaffID,$contigScaffStart,$contigScaffEnd) = $mapInfo =~ /^([\+\-])(.+)\.(\d+)\.(\d+)$/;
	    unless(exists($scaffLengths{$scaffID})) {
		die "Error: Length information for scaffold $scaffID missing.\n";
	    }
	    my $scaffLen = $scaffLengths{$scaffID};
	    $contigScaffStrand = ($contigScaffStrand eq "+") ? "Plus" : "Minus";
	    my $scaffStrand = ($contigScaffStrand eq $strand) ? "Plus" : "Minus";

	    my $scaffStart = ($contigScaffStrand eq "Plus") ? ($contigScaffStart + $sStart - 1) : ($contigScaffEnd - $sStop + 1);
	    my $scaffStop = ($contigScaffStrand eq "Plus") ? ($contigScaffStart + $sStop - 1) : ($contigScaffEnd - $sStart + 1);

	    $subject = $scaffID;
	    $sLength = $scaffLen;
	    $strand = $scaffStrand;
	    $sStart = $scaffStart;
	    $sStop = $scaffStop;

	    $line = join("\t",$blastType,$query,$qStart,$qStop,$qLength,
			 $subject,$sStart,$sStop,$sLength,$strand,$score,$eValue,$identities,$alignLength);
	    
	}
    }

    if ($endAversion) {   # to completely avoid short insert close to ends
	unless (($sStop > $endAversion) && ($sStart < $sLength-$endAversion)) {
	    $omitted_hits{"END AVERSION FAILED"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 3 END AVERSION FAILED (sStart: $sStart sStop: $sStop  sLength: $sLength  endAversion: $endAversion)\n"}
	    next;
	}
    }

    if ($testMode) {
	my $sLengthCutoff = 2*$estimatedInsertSize+8*$estimatedSigma;
	unless ($sLength > $sLengthCutoff) {
	    $omitted_hits{"SCAFFOLD LENGTH CUTOFF FAILED ($sLengthCutoff)"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 4 SCAFFOLD LENGTH CUTOFF FAILED (length: $sLength cutoff: $sLengthCutoff)\n"}
	    next;
	}
    }

    if ($reverseComplement) {
	$strand = ($strand eq "Plus") ? "Minus" : "Plus";
	$line = join("\t",$blastType,$query,$qStart,$qStop,$qLength,
		     $subject,$sStart,$sStop,$sLength,$strand,$score,$eValue,$identities,$alignLength);
    }

    my $unalignedStart = $qStart-1;
    my $projectedStart = undef;
    if ($strand eq "Plus") {
	$projectedStart = $sStart - $unalignedStart;
    } else {
	$projectedStart = $sStop + $unalignedStart;
    }

    if ($innieRemoval) {
	my $minBackDistance = $shortPair;
	my $backDistance = $projectedStart;
	if ($strand eq "Minus") {
	    $backDistance = $sLength - $projectedStart;
	}
	unless ($backDistance > $minBackDistance) {
	    $omitted_hits{"INNIE CONDITION FAILED"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 5 INNIE CONDITION FAILED (backDistance: $backDistance  minBackDistance: $minBackDistance)\n";}
	    next;
	}
    }

    my $startStatus = undef;

    my $projectedOff = 0;
    if ($projectedStart < 1) {
	$projectedOff = 1-$projectedStart;
    } elsif ($projectedStart > $sLength) {
	$projectedOff = $projectedStart-$sLength;
    }
    my $missingStartBases = $unalignedStart-$projectedOff;

    if ($unalignedStart == 0) {
	$startStatus = "FULL";
    } elsif ( ($projectedOff > 0) && ($missingStartBases < $fivePrimeWiggleRoom) ) {
	$startStatus = "GAP";
    } elsif ($unalignedStart < $fivePrimeWiggleRoom) {
	$startStatus = "INC";
    } else {
	$startStatus = "TRUNC";
    }

    my $unalignedEnd = $qLength-$qStop;
    my $projectedEnd = undef;
    if ($strand eq "Plus") {
	$projectedEnd = $sStop + $unalignedEnd;
    } else {
	$projectedEnd = $sStart - $unalignedEnd;
    }

    my $endStatus = undef;

    $projectedOff = 0;
    if ($projectedEnd < 1) {
	$projectedOff = 1-$projectedEnd;
    } elsif ($projectedEnd > $sLength) {
	$projectedOff = $projectedEnd-$sLength;
    }
    my $missingEndBases = $unalignedEnd-$projectedOff;

    if ($unalignedEnd == 0) {
	$endStatus = "FULL";
    } elsif ( ($projectedOff > 0) && ($missingEndBases < $threePrimeWiggleRoom) ) {
	$endStatus = "GAP";
    } elsif ($unalignedEnd < $threePrimeWiggleRoom) {
	$endStatus = "INC";
    } else {
	$endStatus = "TRUNC";
    }

    my $location = "UNK";

    unless ($testMode) {

	if ($strand eq "Plus") {
	    if ($projectedStart > ($sLength-$endDistance)) {
		$location = "OUT";
	    } elsif ($projectedEnd < $endDistance) {
		$location = "IN";
	    } else {
		$location = "MID";
	    }
	} else {
	    if ($projectedStart < $endDistance) {
		$location = "OUT";
	    } elsif ($projectedEnd > ($sLength-$endDistance)) {
		$location = "IN";
	    } else {
		$location = "MID";
	    }
	}
    }

    my $readStatus = "$startStatus.$endStatus.$location";

    if ($testMode) {
	unless ((($startStatus eq "FULL") || ($startStatus eq "INC") ) &&
		(($endStatus eq "FULL") || ($endStatus eq "INC") ) ) {
	    $omitted_hits{"UNRELIABLE READSTATUS"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 6 UNRELIABLE READSTATUS: $readStatus\n"}
	    next;
	}
    } else {
	unless (($readStatus =~ /OUT/) || ($readStatus =~ /GAP/)) {
	    $omitted_hits{"UNINFORMATIVE READSTATUS"}++;
	    if ($verbose) {print STDERR "HIT OMITTED: 7 UNINFORMATIVE READSTATUS: $readStatus\n"}
	    next;
	}
    }

    my ($null,$pairName,$pairEnd) = parse_fastq_header($query);
    unless ($pairName && $pairEnd)
    { die "Error: couldn't extract pairing info from the read name:  $query\n" }


    if (defined($currentPair)) {
	if ($pairName eq $currentPair) {
	    $pairInfo .= "$readStatus\t$line\n";
	} else {
	    unless ($testMode) { 
		processPair($pairInfo);
	    } else {
		$sampleCohortTested = testProcessPair($pairInfo);
		last if ($sampleCohortTested);
	    }

	    $currentPair = $pairName;
	    $pairInfo = "$readStatus\t$line\n";
	}
    } else {
	$currentPair = $pairName;
	$pairInfo = "$readStatus\t$line\n";
    }
}
close B;

if (defined($currentPair)) {
    if ($testMode) {  # we never reached the sampling cutoff, so test the frequecies  using whatever we've got.
	if (!$sampleCohortTested) {
	    testResults();
	}
    } else {
	processPair($pairInfo);
    }
}
else
{
    print STDERR "No pairs passed the cutoff criteria!\n";
}

print STDERR "TOTAL BLAST MAP LINES EVALUATED:  $n \n";

print STDERR "\nHITS OMITTED SUMMARY:\n";
for (keys %omitted_hits)
{
    print STDERR "omitted due to $_ : $omitted_hits{$_}\n";
}

print STDERR "\nPAIRS WITH NON-STANDARD ORIENTATION OR DISTANCE: $nNonStandard\n";

print STDERR "Done.\n";    

# -------------
# |SUBROUTINES|
# -------------


sub testResults {
    my $meanRedundancy = 0;
    my $nLoci = 0;
    while (my ($locus,$redundancy) = each(%testLocations)) {
	$nLoci++;
	$meanRedundancy += $redundancy;
    }
    $meanRedundancy = sprintf("%.1f",$meanRedundancy/($nLoci||1));

    print "# $nFullPairs sampled canonical pairs ($nNonStandard non-standard pairs discarded); mean redundancy = $meanRedundancy\n";

    my $mean = 0;
    my $stdDev = 0;
    my $meanSquare = 0;
    my $n = 0;

    my @sortedPairSizes = sort {$a <=> $b} keys(%pairSizes);
    foreach my $size (@sortedPairSizes) {

	print STDERR "size: $size  freq: $pairSizes{$size}\n";
	if ($size < $minTestLength) {
	    next;
	}

	my $freq = $pairSizes{$size};
	if ($freq >= $minFreqReported) {
	    $mean += $freq*$size;
	    $meanSquare += $freq*$size*$size;
	    $n += $freq;
	}
    }
    if ($n > 0) {
	$mean = $mean/$n;
	$stdDev = sqrt($meanSquare/$n - $mean*$mean);
    }
    printf("# %.0f +/- %.0f (if this is very different from $estimatedInsertSize +/- $estimatedSigma consider rerun!)\n", $mean, $stdDev); 
    print "# Pair separation distribution:\n";


    # The frequency of pairs mapped to the same contig may bee too low (e.g. long-insert libraries) which would result in nothing passing the -M cutoff. 
    # To mitigate this we apply binning

    my %tally = (); 
    my %binrecords = ();
    my $binsize=$mean/100; 
    if ($binsize < 1) {$binsize = "1"};

    foreach my $size (@sortedPairSizes) {
	if ($size < $minTestLength) {
	    next;
	}

	my $freq = $pairSizes{$size};

	my $bin = sprintf("%.0f", $size/$binsize);
	$binrecords{$bin}{$size} = $freq;   # add this size+frequency to the corresponding bin as a hash element

        if (!exists($tally{$bin})) {  # add to the total frequency tally for this bin
            $tally{$bin} = $freq;
        } else {
            $tally{$bin} += $freq;
        }
    }

    foreach $bin ( sort {$a<=>$b} keys %tally ) {
	my $binned_freq = $tally{$bin};	
	if ($binned_freq >= $minFreqReported) {
	    for $size ( sort {$a<=>$b} keys %{ $binrecords{$bin} } ) {
		print "$size\t$binrecords{$bin}{$size} \n";
	    }
	}
    }
}

sub testProcessPair {

    my ($info) = @_;

    my @alignments = split(/\n/,$info);
    my $nAligns = scalar(@alignments);

    unless ($nAligns == 2) {
	if ($verbose) {
	    print STDERR "TPP: Reject pair - too many / too few alignments\n";
	}
	return;
    }

    my @sPositions = ();

    my ($status1,$blastType1,$query1,$qStart1,$qStop1,$qLength1,$subject1,$sStart1,$sStop1,$sLength1,
	$strand1,$score1,$eValue1,$identities1,$alignLength1) = split(/\t/,$alignments[0]);
    push (@sPositions,$sStart1);
    push (@sPositions,$sStop1);
    my ($status2,$blastType2,$query2,$qStart2,$qStop2,$qLength2,$subject2,$sStart2,$sStop2,$sLength2,
	$strand2,$score2,$eValue2,$identities2,$alignLength2) = split(/\t/,$alignments[1]);
    push (@sPositions,$sStart2);
    push (@sPositions,$sStop2);

    unless ( ($query1 ne $query2) && ($subject1 eq $subject2) ) {
        if ($verbose) {
	    print STDERR "TPP: Reject pair - different contigs\n";
	}
	return;
    }

    if ($strand1 eq $strand2) {
	$nNonStandard++;
        if ($verbose) {
	    print STDERR "TPP: Reject pair - not in expected orientation (same strand) \n";
	}
	return;
    }
    
    my $orientation = "undefined";
    if ( (($strand1 eq "Plus") && ($sStart1 < $sStart2) && ($sStart1 < $sStop2)) || 
	 (($strand2 eq "Plus") && ($sStart2 < $sStart1) && ($sStart2 < $sStop1)) ) {
	$orientation = "canonical";
    } elsif ( (($strand1 eq "Plus") && ($sStart1 > $sStop2)) || 
	      (($strand2 eq "Plus") && ($sStart2 > $sStop1)) ) {
	$orientation = "non-canonical";
    }

    my @sPosSorted = sort {$a <=> $b} @sPositions;
    my $extent = $sPosSorted[3]-$sPosSorted[0]+1;

    unless ($orientation eq "canonical") {
	$nNonStandard++;
	if ($verbose) {
	    print STDERR "TPP: Reject pair - not in expected orientation (different strands) \n";
	    if ($innieRemoval && $orientation eq "non-canonical")
	    {
		print STDERR "INNIE\t$extent\n";

	    }
	}
	return;
    }

    my $safeDistance = $estimatedInsertSize+4*$estimatedSigma;
    my $leftBoundary = $safeDistance;
    my $rightBoundary = $sLength1-$safeDistance;
    unless (($sPosSorted[3] > $leftBoundary) && ($sPosSorted[0] < $rightBoundary)) {
        if ($verbose) { print STDERR "TPP: Reject pair - alignment beyond safe boundaries\n"}
	return;
    }

    my $testLocation = "$subject1.$sPosSorted[0].$sPosSorted[3]";
    if (exists($testLocations{$testLocation})) {
	$testLocations{$testLocation}++;
	return;
    } else {
	$testLocations{$testLocation} = 1;
    }

    if ($innieRemoval) {
	if ($extent < $shortPair) {
	    $nNonStandard++;
	    if ($verbose) { print STDERR "TPP: Reject pair - short\n" }
	    return;
	}
    }

    $pairSizes{$extent}++;
    $nFullPairs++;

    if ($nFullPairs == $sampleSize) {
	testResults();
	return 1;
    }
}

sub processPair {
    my ($info) = @_;
    my @alignments = split(/\n/,$info);
    my $nAligns = scalar(@alignments);

    my @aStatus = ();
    my @aContig = ();
    my @aStrand = ();
    my %reads = ();

    my @dataSlice = (0,2,3,4,5,6,7,8,9,10);

    for (my $a = 0; $a < $nAligns; $a++) {
	my $alignment = $alignments[$a];
	my ($readStatus,$blastType,$query,$qStart,$qStop,$qLength,$subject,$sStart,$sStop,$sLength,
	    $strand,$score,$eValue,$identities,$alignLength) = split(/\t/,$alignment);
	push(@aStatus,$readStatus);
	push(@aContig,$subject);
	push(@aStrand,$strand);

	if (exists($reads{$query})) {
	    push(@{$reads{$query}},$a);
	} else {
	    $reads{$query} = [$a];
	}
    }

    my @reads = keys(%reads);
    my @rStatus = ();

    foreach my $read (@reads) {
	my $nAnchors = 0;
	my $nOutOfGap = 0;
	my $nIntoGap = 0;
	my $anchor = "";
	my $outOfGapSplint = "";
	my $outOfGapNoSplint = "";
	my $intoGap = "";
	my @aIndex = @{$reads{$read}};
	foreach my $a (@aIndex) {
	    my $status = $aStatus[$a];
	    my ($startStat,$stopStat,$dirStat) = $status =~ /^(.+)\.(.+)\.(.+)$/;
	    my $contig = $aContig[$a];
	    my $strand = $aStrand[$a];

	    unless ( ($startStat =~ /GAP/) || ($stopStat =~ /GAP/) ) {
		$nAnchors++;
		if ($strand eq "Plus") {
		    $anchor = "$contig.3.$a";
		} else {
		    $anchor = "$contig.5.$a";
		}

	    } else {

		if ($startStat =~ /GAP/) {
		    $nOutOfGap++;
		    if ($strand eq "Plus") {
			$outOfGapSplint = "$contig.5.$a";
			$outOfGapNoSplint = "$contig.5.$a";
			if ($dirStat eq "OUT") {
			    $outOfGapNoSplint = "$contig.3.$a";
			}
		    } else {
			$outOfGapSplint = "$contig.3.$a";
			$outOfGapNoSplint = "$contig.3.$a";
			if ($dirStat eq "OUT") {
			    $outOfGapNoSplint = "$contig.5.$a";
			}
		    }
		}
		if ($stopStat =~ /GAP/) {
		    $nIntoGap++;
		    if ($strand eq "Plus") {
			$intoGap = "$contig.3.$a";
		    } else {
			$intoGap = "$contig.5.$a";
		    }
		}
	    }
	}

	if (($nAnchors == 1) && ($nOutOfGap == 0) && ($nIntoGap == 0)) {
	    push(@rStatus,"ANCHOR.$anchor");
	} elsif (($nAnchors == 0) && ($nOutOfGap == 1) && ($nIntoGap == 0)) {
	    push(@rStatus,"OUTGAP.$outOfGapNoSplint");
	} elsif (($nAnchors == 0) && ($nOutOfGap == 0) && ($nIntoGap == 1)) {
	    push(@rStatus,"INTGAP.$intoGap");
	} elsif (($nAnchors == 0) && ($nOutOfGap == 1) && ($nIntoGap == 1)) {
	    my ($outOfGapContig) = $outOfGapSplint =~ /^(.+)\.[35]\.\d+$/;
	    my ($intoGapContig) = $intoGap =~ /^(.+)\.[35]\.\d+$/;
	    if ($outOfGapContig ne $intoGapContig) {
		push(@rStatus,"SPLINT.$outOfGapSplint.$intoGap");
	    } else {
		push(@rStatus,"OUTGAP.$outOfGapNoSplint");
	    }
	} else {
	    push(@rStatus,"CONFLICT");
	}
    }

    my $nReads = scalar(@reads);
    if ($nReads == 2) {
	my $read1 = $reads[0];
	my $readStatus1 = $rStatus[0];
	my $read2 = $reads[1];
	my $readStatus2 = $rStatus[1];

	unless (($readStatus1 eq "CONFLICT") || ($readStatus2 eq "CONFLICT")) {
	    my ($type1,$details1) = $readStatus1 =~ /^([A-Z]{6})\.(.+)$/;
	    my ($type2,$details2) = $readStatus2 =~ /^([A-Z]{6})\.(.+)$/;

	    my $printInfo1 = "";
	    if ($type1 eq "SPLINT") {
		my ($c1,$a1,$c2,$a2) = $details1 =~ /^(.+\.[35])\.(\d+)\.(.+\.[35])\.(\d+)$/;
		my @info1 = (split(/\t/,$alignments[$a1]))[@dataSlice];
		my @info2 = (split(/\t/,$alignments[$a2]))[@dataSlice];
		$printInfo1 = "$c1\t[@info1]\t$c2\t[@info2]";
	    } else {
		my ($c,$a) = $details1 =~ /^(.+\.[35])\.(\d+)$/;
		my @info = (split(/\t/,$alignments[$a]))[@dataSlice];
		$printInfo1 = "$c\t[@info]";
	    }

	    my $printInfo2 = "";
	    if ($type2 eq "SPLINT") {
		my ($c1,$a1,$c2,$a2) = $details2 =~ /^(.+\.[35])\.(\d+)\.(.+\.[35])\.(\d+)$/;
		my @info1 = (split(/\t/,$alignments[$a1]))[@dataSlice];
		my @info2 = (split(/\t/,$alignments[$a2]))[@dataSlice];
		$printInfo2 = "$c1\t[@info1]\t$c2\t[@info2]";
	    } else {
		my ($c,$a) = $details2 =~ /^(.+\.[35])\.(\d+)$/;
		my @info = (split(/\t/,$alignments[$a]))[@dataSlice];
		$printInfo2 = "$c\t[@info]";
	    }

	    print "PAIR\t$type1.$type2\t$printInfo1\t$printInfo2\n";

	}

    } elsif ($nReads == 1) {
	my $read = $reads[0];
	my $readStatus = $rStatus[0];
	unless ($readStatus eq "CONFLICT") {
	    my ($type,$details) = $readStatus =~ /^([A-Z]{6})\.(.+)$/;
	    if ($type eq "SPLINT") {
		my ($c1,$a1,$c2,$a2) = $details =~ /^(.+\.[35])\.(\d+)\.(.+\.[35])\.(\d+)$/;
		my @info1 = (split(/\t/,$alignments[$a1]))[@dataSlice];
		my @info2 = (split(/\t/,$alignments[$a2]))[@dataSlice];
		print "SINGLE\t$type\t$c1\t[@info1]\t$c2\t[@info2]\n";
	    } else {
		my ($c,$a) = $details =~ /^(.+\.[35])\.(\d+)$/;
		my @info = (split(/\t/,$alignments[$a]))[@dataSlice];
		print "SINGLE\t$type\t$c\t[@info]\n";
	    }
	}
    }
}
