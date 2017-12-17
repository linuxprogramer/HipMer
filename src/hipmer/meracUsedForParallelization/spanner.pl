#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# spanner.pl by Jarrod Chapman <jchapman@lbl.gov> Fri Jun 12 12:58:11 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
my %opts = ();
my $validLine = getopts('b:m:F:T:i:RIS:U:D:', \%opts);
my @required = ("b","m","i");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print STDERR "Usage: ./spanner.pl <-b blastMapGlob> <-m minMatch> <-i insertSize:sigma> <<-R(everseComplement?)>> <<-I(nnieRemoval?)>> <<-S scaffoldReportFile>> <<-U truncateAlignmentEnd(5|3)>> <<-D minEndSeparation>> <<-F fivePrimeWiggleRoom>> <<-T threePrimeWiggleRoom>>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\nspanner.pl";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $blastMapFileGlob = $opts{"b"};
my @blastMapFiles = glob($blastMapFileGlob);

my $minMatch = $opts{"m"};

my ($insertSize,$insertSigma) = $opts{"i"} =~ /^(\d+)\:(\d+)$/;
my $endDistance = $insertSize + 3*$insertSigma;

my $fivePrimeWiggleRoom = 5;
if (exists($opts{"F"})) {
    $fivePrimeWiggleRoom = $opts{"F"};
}
my $threePrimeWiggleRoom = 5;
if (exists($opts{"T"})) {
    $threePrimeWiggleRoom = $opts{"T"};
}

my $reverseComplement = 0;
if (exists($opts{"R"})) {
    $reverseComplement = 1;
}

my $innieRemoval = 0;
if (exists($opts{"I"})) {
    $innieRemoval = 1;
}

my $minEndSeparation = 0;
if (exists($opts{"D"})) {
    $minEndSeparation = $opts{"D"};
}

my $srfFile = undef;
if (exists($opts{"S"})) {
    $srfFile = $opts{"S"};
}

my $truncate = 0;
if (exists($opts{"U"})) {
    $truncate = $opts{"U"};
}

my %rejected = ();
my @rejectReasons = ("FORMAT","MINLEN","UNINFORMATIVE","5-TRUNCATED","3-TRUNCATED","SINGLETON",
		     "TRUNC5_FAIL","TRUNC3_FAIL", "NO_SCAFFOLD_INFO", "POTENTIAL_INNIE", "POTENTIAL_SHORTY");
foreach my $r (@rejectReasons) {
    $rejected{$r} = 0;
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
my $currentPairAlignments = 0;
my $pairInfo = "";
my %pairSizes = ();
my $totalSpansFound = 0;
my $totalHitsFound = 0;

foreach my $blastMapFile (@blastMapFiles) {

    print STDERR "Reading blast map file: $blastMapFile...\n";
    open (B,$blastMapFile) || die "Couldn't open $blastMapFile\n";

    while (my $line = <B>) {
	chomp $line;
	if ($line =~ /^BLAST_TYPE/) {
	    next;
	}
	
	my ($blastType,$query,$qStart,$qStop,$qLength,$subject,$sStart,$sStop,$sLength,
	    $strand,$score,$eValue,$identities,$alignLength) = split(/\t/,$line);

	# Eugene's code for checking valid query names
	if (($query !~ /\/[12]/) && ($query !~ /\s[12]:/)) {  #known Illumina-style demarkations of fwd/rev reads
	    die "Query name $query is not in valid paired-end format\n";
	}
	# trim any extra junk in the query name following /1 or /2
	$query =~ s/(\S+\/[12]).*/$1/;
	# re-form the line
	$line = join("\t",$blastType,$query,$qStart,$qStop,$qLength,
		     $subject,$sStart,$sStop,$sLength,$strand,$score,$eValue,$identities,$alignLength);
	# End query validation
	
	$totalHitsFound++;

	if ($blastType ne "BLASTN") {
	    $rejected{"FORMAT"}++;
	    next;
	}
	
	if ($identities < $minMatch) {
	    $rejected{"MINLEN"}++;
	    next;
	}

	# Assess alignment for completeness (do this before scaffold coordinate conversion!)
	# Important: Truncations are performed before reverse complementation
	# and apply to the end of the actual read
	my $unalignedStart = $qStart-1;
	my $projectedStart = undef;
	if ($strand eq "Plus") {
	    $projectedStart = $sStart - $unalignedStart;
	} else {
	    $projectedStart = $sStop + $unalignedStart;
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
	} elsif (($unalignedStart < $fivePrimeWiggleRoom) || ($truncate == 5)) {
	    $startStatus = "INC";
	} else {
	    $rejected{"5-TRUNCATED"}++;
	    next;
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
	} elsif (($unalignedEnd < $threePrimeWiggleRoom) || ($truncate == 3)) {
	    $endStatus = "INC";
	} else {
	    $rejected{"3-TRUNCATED"}++;
	    next;
	}

	# Re-orient alignment if requested
	if ($reverseComplement) {
	    $strand = ($strand eq "Plus") ? "Minus" : "Plus";
	    $qStart = $qLength-$qStop+1;
	    $qStop = $qLength-$qStart+1;
	    $line = join("\t",$blastType,$query,$qStart,$qStop,$qLength,
			 $subject,$sStart,$sStop,$sLength,$strand,$score,$eValue,$identities,$alignLength);
	}

	# Convert to scaffold coordinate system if srfFile is specified
	if (defined($srfFile)) {
	    unless (exists($contigScaffoldMap{$subject})) {
		$rejected{"NO_SCAFFOLD_INFO"}++;
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

		$projectedStart = ($strand eq "Plus") ? $sStart - $unalignedStart : $sStop + $unalignedStart;
		$projectedEnd = ($strand eq "Plus") ? $sStop + $unalignedEnd : $sStart - $unalignedEnd;

		$line = join("\t",$blastType,$query,$qStart,$qStop,$qLength,
			     $subject,$sStart,$sStop,$sLength,$strand,$score,$eValue,$identities,$alignLength);
	    }
	}
	
	# Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle)
	my $location = "UNK";
	
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
	
	my $readStatus = "$startStatus.$endStatus.$location";
	
	unless ($location eq "OUT") {
#	    print STDERR "UNINFORMATIVE: $location:$projectedStart:$projectedEnd:$sLength:$endDistance $line\n";
	    $rejected{"UNINFORMATIVE"}++;
	    next;
	}

	# Assumes read names of the form NAME/1, NAME/2, etc.
	my ($pairName,$pairEnd) = $query =~ /^(.+)\/(\d+)$/;

	if (defined($currentPair)) {
	    if ($pairName eq $currentPair) {
		$pairInfo .= "$readStatus\t$line\n";
		$currentPairAlignments++;
	    } else {
		if ($currentPairAlignments > 1) {
		    $totalSpansFound += processPair($pairInfo);
		} else {
#		    print STDERR "SINGLETON: $line\n";
		    $rejected{"SINGLETON"}++;
		}
		$currentPair = $pairName;
		$pairInfo = "$readStatus\t$line\n";
		$currentPairAlignments = 1;
	    }
	} else {
	    $currentPair = $pairName;
	    $pairInfo = "$readStatus\t$line\n";
	    $currentPairAlignments = 1;
	}
    }
    close B;
    print STDERR "Done.\n";

    if (defined($currentPair)) {
	if ($currentPairAlignments > 1) {
	    $totalSpansFound += processPair($pairInfo);
	} else {
	    $rejected{"SINGLETON"}++;
	}
    }
}

print STDERR "Total alignments found: $totalHitsFound\n";
print STDERR "Unused alignments:\n";
foreach my $r (@rejectReasons) {
    my $n = $rejected{$r};
    if ($n) {
	print STDERR "\t$r\t$n\n";
    }
}
print STDERR "Total spans found = $totalSpansFound\n";

$date = `date`;
chomp $date;
print STDERR "Done. $date\n";

# -------------
# |SUBROUTINES|
# -------------

sub processPair {
    my ($info) = @_;
    my @alignments = split(/\n/,$info);
    my $nAligns = scalar(@alignments);

    my %firstAlignments = ();

    my %alignInfo = ();
    for (my $a = 0; $a < $nAligns; $a++) {
	$alignInfo{$a} = [split(/\t/,$alignments[$a])];
	my ($readName,$readStart) = ($alignInfo{$a}->[2],$alignInfo{$a}->[3]);
	if (exists($firstAlignments{$readName})) {
	    if ($readStart < $firstAlignments{$readName}->[1]) {
		$firstAlignments{$readName} = [$a,$readStart];
	    }
	} else {
	    $firstAlignments{$readName} = [$a,$readStart];
	}
    }

    my @readNames = keys(%firstAlignments);
    my $nReads = scalar(@readNames);
    unless ($nReads == 2) {
	$rejected{"SINGLETON"}++;
	return 0;
    }

    my @dataSlice = (0,2,3,4,5,6,7,8,9,10);
#   0=$readStatus,1=$query,2=$qStart,3=$qStop,4=$qLength,5=$subject,6=$sStart,7=$sStop,8=$sLength,9=$strand

    my $alignInfo1 = $alignInfo{$firstAlignments{$readNames[0]}->[0]};
    my @alignInfo1 = @{$alignInfo1}[@dataSlice];
    my ($status1,$contig1,$strand1) = @alignInfo1[0,5,9];
    
    my $alignInfo2 = $alignInfo{$firstAlignments{$readNames[1]}->[0]};
    my @alignInfo2 = @{$alignInfo2}[@dataSlice];
    my ($status2,$contig2,$strand2) = @alignInfo2[0,5,9];

    unless ($contig1 ne $contig2) {
	$rejected{"UNINFORMATIVE"}++;
	return 0;
    }

    # Exclude pairs whose orientation is ambiguous (under innie-artifact model)
    # Note: Using projected start and stop coordinates would be slightly less conservative than this implementation
    if ($innieRemoval) {
	my $innieMaxSeparation = 800;
	my ($sStart1,$sStop1,$sLength1) = @alignInfo1[6,7,8];
	my ($sStart2,$sStop2,$sLength2) = @alignInfo2[6,7,8];

	my $s1Innie = ($strand1 eq "Plus") ? $sStop1 : ($sLength1-$sStart1+1); 
	my $s2Innie = ($strand2 eq "Plus") ? $sStop2 : ($sLength2-$sStart2+1); 
	my $innieSep = $s1Innie+$s2Innie;
	
	if ($innieSep < $innieMaxSeparation) {
	    $rejected{"POTENTIAL_INNIE"}++;
	    return 0;
	}
    }

    # Exclude pairs whose separation is ambiguous (under short-modality model)
    # Note: Using projected start and stop coordinates would be slightly less conservative than this implementation
    if ($minEndSeparation) {
	my ($sStart1,$sStop1,$sLength1) = @alignInfo1[6,7,8];
	my ($sStart2,$sStop2,$sLength2) = @alignInfo2[6,7,8];

	my $s1Outie = ($strand1 eq "Plus") ? ($sLength1-$sStart1+1) : $sStop1;
	my $s2Outie = ($strand2 eq "Plus") ? ($sLength2-$sStart2+1) : $sStop2;
	my $outieSep = $s1Outie+$s2Outie;

	if ($outieSep < $minEndSeparation) {
	    $rejected{"POTENTIAL_SHORTY"}++;
	    return 0;
	}
    }

    my ($startStat1,$stopStat1,$dirStat1) = $status1 =~ /^(.+)\.(.+)\.(.+)$/;
    my ($startStat2,$stopStat2,$dirStat2) = $status2 =~ /^(.+)\.(.+)\.(.+)$/;

    my $type1 = "ANCHOR";
    if ($startStat1 eq "GAP") {
	$type1 = "OUTGAP";
    } elsif ($stopStat1 eq "GAP") {
	$type1 = "INTGAP";
    }
    my $type2 = "ANCHOR";
    if ($startStat2 eq "GAP") {
	$type2 = "OUTGAP";
    } elsif ($stopStat2 eq "GAP") {
	$type2 = "INTGAP";
    }
    my $contigEnd1 = ($strand1 eq "Plus") ? "3" : "5";
    my $contigEnd2 = ($strand2 eq "Plus") ? "3" : "5";

    print "PAIR\t$type1.$type2\t$contig1.$contigEnd1\t[@alignInfo1]\t$contig2.$contigEnd2\t[@alignInfo2]\n";

    return 1;
}

