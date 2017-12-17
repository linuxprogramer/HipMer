#! /usr/bin/env perl
#
# gapPlacer.pl by Jarrod Chapman <jchapman@lbl.gov> Fri Jun 12 12:58:11 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use warnings;

#use lib "$ENV{MERACULOUS_ROOT}/lib";
use lib "perllib";
use M_Utility;
use Getopt::Std;
#use Data::Dumper;

#my $ALN_CUTOFF = 1000000;
#my $PRINT_READ_DATA = 1;

my %opts = ();
my $validLine = getopts('b:m:F:T:i:s:f:c:GP', \%opts);
my @required = ("b","m","i","s","f","c");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print STDERR "Usage: ./gapPlacer.pl <-b blastMapGlob> <-m merSize> <-i insertSize:sigma> <-s scaffoldReportFile> <-f fastqFileList> <-c contigFastaFile> <<-G(zipped?)>> <<-F fivePrimeWiggleRoom>> <<-T threePrimeWiggleRoom>> <<-P(noPairProjection)>>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\ngapPlacer.pl";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $blastMapFileGlob = $opts{"b"};
my @blastMapFiles = glob($blastMapFileGlob);

my $merSize = $opts{"m"};

my ($insertSize,$insertSigma) = $opts{"i"} =~ /^(\d+)\:(\d+)$/;
my $projectZ = 2;
my $maxProject = $insertSize+$projectZ*$insertSigma;
my $minProject = $insertSize-$projectZ*$insertSigma;

my $srfFile = $opts{"s"};
my $contigFastaFile = $opts{"c"};
my $fastqFileList = $opts{"f"};
my $gzipped = 0;
if (exists($opts{"G"})) {
    $gzipped = 1;
}

my $pairProjection = 1;
if (exists($opts{"P"})) {
    $pairProjection = 0;
}

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

my $truncate = 0;
if (exists($opts{"U"})) {
    $truncate = $opts{"U"};
}

my %rejected = ();
my @rejectReasons = ("FORMAT","MINLEN","UNINFORMATIVE","5-TRUNCATED","3-TRUNCATED","SINGLETON","NO_ASSOCIATED_GAP",
                     "TRUNC5_FAIL","TRUNC3_FAIL", "NO_SCAFFOLD_INFO", "POTENTIAL_INNIE", "POTENTIAL_SHORTY");
foreach my $r (@rejectReasons) {
    $rejected{$r} = 0;
}

my @gapInfo = ();  # gaps indexed along srf file [scaff, prevContigEnd, nextContigEnd, scaffStart:scaffEnd:uncertainty, (read1, read2, ...)]
my %contigGapMap = ();  # map from contig to flanking gap indices [5'-gapID,3'-gapID]
my %readData = ();  # read name -> sequence:quality
my %scaffInfo = (); # scaffName -> [nContigs,scaffLen]

print STDERR "Reading scaffold report file: $srfFile...\n";
open (S,$srfFile) || die "Couldn't open $srfFile\n";
my $gapIndex = 0;
my $currentScaff = "NULL";
my $prevContigEnd = undef;
my $scaffCoord = 0;
my $openGap = 0;
while (my $line = <S>) {
    chomp $line;
    my @cols = split(/\t/,$line);
    my $scaffID = $cols[0];
    
    unless (exists($scaffInfo{$scaffID})) {
        $scaffInfo{$scaffID} = [0,0];
    }

    if ($cols[1] =~ /^CONTIG/) {   # a contig
        my ($contigID,$sStart,$sEnd) = @cols[2,3,4];
        my ($cStrand,$cName) = $contigID =~ /^([\+\-])(.+)$/;

        $prevContigEnd = ($cStrand eq "+") ? "$cName.3" : "$cName.5";
        $scaffCoord = $sEnd;
        $scaffInfo{$scaffID}->[1] = $sEnd;
        $scaffInfo{$scaffID}->[0]++;
        $contigGapMap{$cName} = [undef,undef];
        
        if ($openGap) {
            if ($cStrand eq "+") {
                $contigGapMap{$cName}->[0] = $gapIndex-1;
                $gapInfo[$gapIndex-1]->[2] = "$cName.5";
            } else {
                $contigGapMap{$cName}->[1] = $gapIndex-1;
                $gapInfo[$gapIndex-1]->[2] = "$cName.3";
            }	    
        }
        $openGap = 0;

    } else {  # a gap
        my ($gapSize,$uncertainty) = @cols[2,3];
        my $nextScaffCoord = $scaffCoord+$gapSize+1;
        $gapInfo[$gapIndex] = [$scaffID,$prevContigEnd,"?","$scaffCoord:$nextScaffCoord:$uncertainty"];

        my ($cName,$cEnd) = $prevContigEnd =~ /^(.+)\.([35])$/;
        if ($cEnd eq "5") {
            $contigGapMap{$cName}->[0] = $gapIndex;
        } else {
            $contigGapMap{$cName}->[1] = $gapIndex;
        }

        $gapIndex++;
        $openGap = 1;
    }
    
}
close S;
print STDERR "Done.\n";

# while (($key, $val) = each (%contigGapMap)) {
#     $v1 = -1;
#     if (defined($val->[0])) {
#         $v1 = $val->[0];
#     } 
#     $v2 = -1;
#     if (defined($val->[1])) {
#         $v2 = $val->[1];
#     } 
#     printf("$key $v1 $v2\n");
# }

# exit;

my $totalGaps = $gapIndex;

my $currentPair = undef;
my $currentPairAlignments = 0;
my $pairInfo = "";
my %pairSizes = ();
my $alignedGapReads = 0;
my $projectedGapReads = 0;

foreach my $blastMapFile (@blastMapFiles) {

    print STDERR "Reading blast map file: $blastMapFile...\n";
    open (B,$blastMapFile) || die "Couldn't open $blastMapFile\n";
    
    my $nalns = 0;

    while (my $line = <B>) {
        chomp $line;	
        my ($blastType,$query,$qStart,$qStop,$qLength,$subject,$sStart,$sStop,$sLength,
            $strand,$score,$eValue,$identities,$alignLength) = split(/\t/,$line);
        next if ($blastType eq "BLAST_TYPE");


        #if ($blastType ne "BLASTN") {
        #    $rejected{"FORMAT"}++;
        #    next;
        #}
        
        #if ($identities < $merSize) {
        #    $rejected{"MINLEN"}++;
        #    next;
        #}

        # If there's no gap associated with the contig, skip the alignment
        unless (defined($contigGapMap{$subject}->[0]) || defined($contigGapMap{$subject}->[1])) {
            $rejected{"NO_ASSOCIATED_GAP"}++;
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
        

#	print STDERR "DEBUG: $query $qStart $sStart  ::::  projected start in $subject: $projectedStart\n";

        my $startStatus = undef;
        
        my $projectedOff = 0;
        if ($projectedStart < 1) {
            $projectedOff = 1-$projectedStart;
        } elsif ($projectedStart > $sLength) {
            $projectedOff = $projectedStart-$sLength;
        }
        my $missingStartBases = $unalignedStart-$projectedOff;
        
        # classify this alignment
        if ($unalignedStart == 0) {
            $startStatus = "FUL";
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
            $endStatus = "FUL";
        } elsif ( ($projectedOff > 0) && ($missingEndBases < $threePrimeWiggleRoom) ) {
            $endStatus = "GAP";
        } elsif (($unalignedEnd < $threePrimeWiggleRoom) || ($truncate == 3)) {
            $endStatus = "INC";
        } else {
            $rejected{"3-TRUNCATED"}++;
            next;
        }

        my $readStatus = "$startStatus.$endStatus";

        # Re-orient alignment if requested
        if ($reverseComplement) {
            $strand = ($strand eq "Plus") ? "Minus" : "Plus";
            $qStart = $qLength-$qStop+1;
            $qStop = $qLength-$qStart+1;
        }

        # Form the alignment string
        my $alnStr = join("\t",$blastType,$query,$qStart,$qStop,$qLength,
                          $subject,$sStart,$sStop,$sLength,$strand,$score,$eValue,$identities,$alignLength);

        # Parse pair name and end from the query line 
        my ($null,$pairName,$pairEnd) = M_Utility::parse_fastq_header($query);
        unless (defined($pairName) && defined($pairEnd))
        { die "Error: couldn't extract pairing info from the read name:  $query\n" }

        # Determine this read's placement wrt gaps

        if (defined($currentPair)) {
            if ($pairName eq $currentPair) {
                $pairInfo .= "$readStatus\t$alnStr\n";
                $currentPairAlignments++;
            } else {
                #print "$pairInfo\n";
                my $result = processPair($pairInfo); 
                my ($nAligned,$nProjected) = $result =~ /^(\d+):(\d+)$/;
                $alignedGapReads += $nAligned;
                $projectedGapReads += $nProjected;
                
                $currentPair = $pairName;
                $pairInfo = "$readStatus\t$alnStr\n";
                $currentPairAlignments = 1; 
           }
        } else {
            $currentPair = $pairName;
            $pairInfo = "$readStatus\t$alnStr\n";
            $currentPairAlignments = 1;
        }

#        print "$nalns\n\t$currentPair $currentPairAlignments\n$pairInfo";

        $nalns++;
        if (defined $ALN_CUTOFF && $nalns > $ALN_CUTOFF) {last;}
    }
    close B;
    print STDERR "Done.\n";

    if (defined($currentPair)) {
        my $result = processPair($pairInfo);
        my ($nAligned,$nProjected) = $result =~ /^(\d+):(\d+)$/;
        $alignedGapReads += $nAligned;
        $projectedGapReads += $nProjected;
    }

    if (defined $PRINT_READ_DATA) {
        my $tot_nreads = 0;
        my $max_nreads = 0;
        for (my $g=0;$g<$totalGaps;$g++) {
            my ($scaff,$prevContigEnd,$nextContigEnd,$coords,@reads) = @{$gapInfo[$g]};
            my $nReads = scalar(@reads);
            if ($nReads > $max_nreads) {
                $max_nreads = $nReads;
            }
            $tot_nreads += $nReads;
            #print "gap $g, nreads $nReads\n";
            foreach my $r (sort @reads) {
                unless (exists($readData{$r})) {
                    warn "Warning: no sequence found for read $r. Omitting.\n";
                    next;
                }
                #print "\t$r $readData{$r}\n";
                print "$g $r $readData{$r}\n";
            }
        }
        #printf("max nreads $max_nreads, average %d\n", $tot_nreads / $totalGaps);
    }
}

print STDERR "Unused alignments:\n";
foreach my $r (@rejectReasons) {
    my $n = $rejected{$r};
    if ($n) {
        print STDERR "\t$r\t$n\n";
    }
}
print STDERR "Total reads placed in gaps = $alignedGapReads (aligned) + $projectedGapReads (projected)\n";

open (FOF,"<$fastqFileList");
while (<FOF>)
{
    my @FOFline = split /\s+/, $_;
    my $seqFile = $FOFline[0];
    print STDERR "Reading sequence file $seqFile...\n";

    if ($gzipped) {
        open (S,"gunzip -c $seqFile |") || die "Couldn't open pipe to $seqFile\n";
    } else {
        open (S,$seqFile) || die "Couldn't open $seqFile\n";
    }

    while (my $line = <S>) {
        my $readName = undef;
        my $nts = undef;
        my $quals = undef;

        chomp $line;
        unless ($line =~ s/^@//) {
            die "Invalid fastQ $line\n";
        }
        $readName = $line;

        $line = <S>;
        chomp $line;
        $nts = $line;
        $line = <S>;
        unless ($line =~ /^\+/) {
            die "Invalid fastQ $line\n";
        }
        $line = <S>;
        chomp $line;
        $quals = $line;

        if (exists($readData{$readName})) {
            my $orient = $readData{$readName};
            if ($orient eq "-") {
                $nts = reverse($nts);
                $nts =~ tr/ACGTacgt/TGCAtgca/;
                $quals = reverse($quals);
            }
            $readData{$readName} = "$nts:$quals";
            #print "$readName $readData{$readName}\n";
        }
    }
    close S;
    print STDERR "Done.\n";
}
close FOF;


# Read in contig sequence

print STDERR "Reading $contigFastaFile...\n";

my %contigSequence = ();
my $currentEntry = undef;
my $seq = "";
open (F,$contigFastaFile) || die "Couldn't open $contigFastaFile\n";
while (my $line = <F>) {
    chomp $line;
    if ($line =~ /^>/) {
        if (defined($currentEntry)) {
            $contigSequence{$currentEntry} = $seq;
        }
        $seq = "";
        ($currentEntry) = $line =~ /^>(.+)$/;
    } else {
        $seq .= $line;
    }
}
close F;

if (defined($currentEntry)) {
    $contigSequence{$currentEntry} = $seq;
}
print STDERR "Done.\n";

for (my $g=0;$g<$totalGaps;$g++) {
    my ($scaff,$prevContigEnd,$nextContigEnd,$coords,@reads) = @{$gapInfo[$g]};
    my $nReads = scalar(@reads);
    unless ($nReads) {
#	print STDERR "No reads are associated with gap $g !\n";
        next;
    }

    my ($prevContig,$prevEnd) = $prevContigEnd =~ /^(.+)\.([35])$/;
    unless (exists($contigSequence{$prevContig})) {
        die "Sequence not found for contig $prevContig\n";
    }
    my $prevSeq = $contigSequence{$prevContig};
    my $primer1 = ($prevEnd eq "5") ? substr($prevSeq,0,$merSize) : substr($prevSeq,-$merSize);
    if ($prevEnd eq "5") {
        $primer1 = reverse($primer1);
        $primer1 =~ tr/ACGTacgt/TGCAtgca/;
    }

    my ($nextContig,$nextEnd) = $nextContigEnd =~ /^(.+)\.([35])$/;
    unless (exists($contigSequence{$nextContig})) {
        die "Sequence not found for contig $nextContig\n";
    }
    my $nextSeq = $contigSequence{$nextContig};
    my $primer2 = ($nextEnd eq "5") ? substr($nextSeq,0,$merSize) : substr($nextSeq,-$merSize);
    if ($nextEnd eq "3") {
        $primer2 = reverse($primer2);
        $primer2 =~ tr/ACGTacgt/TGCAtgca/;
    }

    my ($leftScaffCoord,$rightScaffCoord,$uncertainty) = split(/:/,$coords);
    my $gapSize = $rightScaffCoord-$leftScaffCoord-1;
    print "$scaff\t$prevContigEnd\t$primer1\t$nextContigEnd\t$primer2\t$gapSize\t$uncertainty";
    foreach my $r (@reads) {
        unless (exists($readData{$r})) {
            warn "Warning: no sequence found for read $r. Omitting.\n";
            next;
        }
        print "\t$readData{$r}";
    }
    print "\n";
}


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

    # Alignment fields:
    # 0=readStatus,1=BLASTN,2=$query,3=$qStart,4=$qStop,5=$qLength,6=$subject,7=$sStart,8=$sStop,9=$sLength,10=$strand

    my $nPlaced = 0;
    my $nProjected = 0;

    my %readGapPairs = ();
    my %readsToAlignments = ();

    my %orientationMap = ();
    my @combos = ("LPlus5","LPlus3","LMinus5","LMinus3",
                  "RPlus5","RPlus3","RMinus5","RMinus3");

    # e.g. "LPlus5" means the read is aligned on the Plus strand to the contig to the Left of the gap
    # whose 5' end is adjacent to the gap

    my @orients = ("-","+","+","-",
                   "+","-","-","+");
    @orientationMap{@combos} = @orients;

    my $maxAlignedScaffLength = 0;

#    print STDERR "DEBUG: processPair():  evaluating alignments: @alignments \n";

    for (my $a = 0; $a < $nAligns; $a++) {

        my @alignInfo = split(/\t/,$alignments[$a]);
        my ($read5status,$read3status) = $alignInfo[0] =~ /^(...)\.(...)$/;
        my ($readName,$contigName,$strand) = @alignInfo[2,6,10];

        # If there's no gap associated with the contig, skip it
        unless (defined($contigGapMap{$contigName}->[0]) || defined($contigGapMap{$contigName}->[1])) {
#	    print STDERR "DEBUG: No gap is associated with contig $contigName.. Skippping..\n";
            next;
        }
        
        # If the read aligns into a gap, add a read<->gap pair
        if ($read5status eq "GAP") {
            my $gapID = ($strand eq "Plus") ? $contigGapMap{$contigName}->[0] : $contigGapMap{$contigName}->[1];
            if (defined($gapID)) {
                $readGapPairs{"$readName:$gapID"} = "$strand:$contigName";
                #print "XX $readName:$gapID $strand:$contigName\n";
            }
        }
        if ($read3status eq "GAP") {
            my $gapID = ($strand eq "Plus") ? $contigGapMap{$contigName}->[1] : $contigGapMap{$contigName}->[0];
            if (defined($gapID)) {
                $readGapPairs{"$readName:$gapID"} = "$strand:$contigName";
                #print "XX $readName:$gapID $strand:$contigName\n";
            }
        }

        # one alignment is recorded for each read to be used to project the read's mate if needed
        # the alignment to the largest scaffold is retained
        # each read is assigned to a single contig/scaffold for this purpose
        # mate pairs may only be projected into gap(s) on a single scaffold by this construction
        # (although direct read alignments may place reads in gaps on multiple scaffolds)

        my $gapID = ($strand eq "Plus") ? $contigGapMap{$contigName}->[1] : $contigGapMap{$contigName}->[0];
        if (defined($gapID)) {
            my $scaffID = $gapInfo[$gapID]->[0];
            my $scaffLength = $scaffInfo{$scaffID}->[1];
            if (exists $readsToAlignments{$readName}) {
                if ($scaffLength > $maxAlignedScaffLength) {
                    $readsToAlignments{$readName} = $a;
                    $maxAlignedScaffLength = $scaffLength;
                }
            } else {
                $readsToAlignments{$readName} = $a;
                $maxAlignedScaffLength = 0;
            }
        }
    }

    # For reads aligning directly into a gap, record the orientation of the read wrt to the gap

    while (my ($rg,$sc) = each(%readGapPairs)) {
        my ($r,$g) = $rg =~ /^(.+):(\d+)$/;
        my ($s,$c) = $sc =~ /^(.+):(.+)$/;

        my $leftContig = $gapInfo[$g]->[1];
        my $rightContig = $gapInfo[$g]->[2];
        my ($lName,$lEnd) = $leftContig =~ /^(.+)\.([35])$/;
        my ($rName,$rEnd) = $rightContig =~ /^(.+)\.([35])$/;
        my $side = undef;
        my $end = undef;
        if ($c eq $lName) {
            $side = "L";
            $end = $lEnd;
        } elsif ($c eq $rName) {
            $side = "R";
            $end = $rEnd;
        } else {
            die "Error: Neither [$leftContig] nor [$rightContig] match $c\n";
        }
        my $combo = "$side$s$end";
        my $orient = $orientationMap{$combo};

#	print STDERR "DEBUG: Gap alignment:  $r\t$g\t$s\t$c\t$combo\t$orient\n";

        push(@{$gapInfo[$g]},$r);
        # FIXME: this is a hack. Surely there shouldn't be duplicate mappings? (shofmeyr@lbl.gov)
        if (exists $readData{$r}) {
            if ($orient eq "+") {
                $readData{$r} = $orient;
            }
        } else {
            $readData{$r} = $orient;
        }
        $nPlaced++;
    }

    unless ($pairProjection == 1) {
        return "$nPlaced:$nProjected";
    }


    # Projected pairs

    while (my ($r,$a) = each(%readsToAlignments)) {
#        print "$r $a\n";

        my $mateName = M_Utility::mate_name($r);
        unless (defined $mateName) {
            die "Error: Could not determine mate name for read $r\n";
        }

        # Attempt to project the unplaced mate into a gap
        # assumes unplaced mate is same length as the placed one
        unless (exists($readsToAlignments{$mateName})) {         #If the mate has an alignment, don't project it

            my @alignInfo = split(/\t/,$alignments[$a]);
            my ($rStart,$rStop,$rLen,$cName,$cStart,$cStop,$cLen,$strand) = @alignInfo[3..10];
            my $unalignedStart = $rStart-1;
            my $farProjection = $maxProject-$unalignedStart;
            my $nearProjection = $minProject-$rLen-$unalignedStart;

            # in the coordinate system of the aligned contig, the mate is projected to lie between nearProjection and farProjection
            my $startGap = undef;
            if ($strand eq "Plus") {
                $farProjection = $cStart+$farProjection;
                $nearProjection = $cStart+$nearProjection;
                $startGap = $contigGapMap{$cName}->[1];  # the 3'-contig-end gap		
            } else {
                $farProjection = $cStop-$farProjection;
                $nearProjection = $cStop-$nearProjection;
                $startGap = $contigGapMap{$cName}->[0];  # the 5'-contig-end gap
            }

            if (defined($startGap)) {
                my $scaffID = $gapInfo[$startGap]->[0];

                my $testGap = $startGap;
                my ($testScaff,$leftContig,$rightContig,$gapCoords) = @{$gapInfo[$testGap]};

#		print STDERR "DEBUG:  Projects into the gap @{$gapInfo[$testGap]} \n";

                my ($lcName,$lcEnd) = $leftContig =~ /^(.+)\.([35])$/;
                my ($rcName,$rcEnd) = $rightContig =~ /^(.+)\.([35])$/;
                my ($leftGapCoord,$rightGapCoord) = $gapCoords =~ /^(\d+):(\d+):/;
                my $iterator = undef;
                my $contigOrigin = undef;

                # find anchor-contig origin position in scaffold coords, assign direction of scan via iterator,
                # and transform projection bounds to scaffold coordinate system

                my $leftProjection = undef;
                my $rightProjection = undef;
                my $orient = undef;
                
                if ($lcName eq $cName) {  # the contig to the left of the gap is aligned to the projector
                    $iterator = 1;
                    $orient = $orientationMap{"L$strand$lcEnd"};  
                    if ($strand eq "Plus") {
                        $contigOrigin = $leftGapCoord-$cLen+1;
                        $leftProjection = $contigOrigin+$nearProjection;
                        $rightProjection = $contigOrigin+$farProjection;
                    } else {
                        $contigOrigin = $leftGapCoord;
                        $leftProjection = $contigOrigin-$nearProjection;
                        $rightProjection = $contigOrigin-$farProjection;
                    }
                } elsif ($rcName eq $cName) { # the contig to the right of the gap is aligned to the projector
                    $iterator = -1;
                    $orient = $orientationMap{"R$strand$rcEnd"}; 
                    if ($strand eq "Plus") {
                        $contigOrigin = $rightGapCoord+$cLen-1;
                        $leftProjection = $contigOrigin-$farProjection;
                        $rightProjection = $contigOrigin-$nearProjection;
                    } else {
                        $contigOrigin = $rightGapCoord;
                        $leftProjection = $contigOrigin+$farProjection;
                        $rightProjection = $contigOrigin+$nearProjection;
                    }
                } else {
                    die "Error: Flanking contigs ($lcName [$testGap] $rcName) do not match anchor contig $cName.\n"; 
                }

                # The orientation of the projected read is opposite the projector
                if (defined($orient)) {
                    $orient = ($orient eq "+") ? "-" : "+";
                }

                #print "$mateName $iterator $leftProjection $rightProjection $orient $unalignedStart $farProjection $nearProjection $startGap $scaffID\n";
                #next;

                # Iterate over gaps projecting mate into each gap that is within range

                while (1) {
                    my $placed = 0;
                    my ($l,$r) = ($leftGapCoord < $rightGapCoord) ? ($leftGapCoord,$rightGapCoord) : ($rightGapCoord,$leftGapCoord);
                    if (($rightProjection > $l) && ($leftProjection < $r)) {
                        $placed = 1;
#			print STDERR "DEBUG: Projection placed:  rightProj: $rightProjection   leftProj: $leftProjection \n";
                    } elsif ( (($iterator == 1) && ($l > $rightProjection)) ||  
                              (($iterator == -1) && ($r < $leftProjection)) )  {
                        last;
                    }

                    if ($placed) {
                        push(@{$gapInfo[$testGap]},$mateName);
                        $readData{$mateName} = $orient;
                        $nProjected++;
                    }

                    $testGap += $iterator;
                    if (($testGap < 0) || ($testGap == $totalGaps)) {
                        last;
                    }

                    ($testScaff,$leftContig,$rightContig,$gapCoords) = @{$gapInfo[$testGap]};
                    unless ($testScaff eq $scaffID) {
                        last;
                    }
                    ($leftGapCoord,$rightGapCoord) = $gapCoords =~ /^(\d+):(\d+):/;
                }
            }
        }
    }

    return "$nPlaced:$nProjected";
}

