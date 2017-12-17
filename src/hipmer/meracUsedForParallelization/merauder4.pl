#!/usr/bin/perl
#
# merauder.pl by Jarrod Chapman <jchapman@lbl.gov> Tue Jun  2 07:48:59 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
use Time::HiRes qw(gettimeofday tv_interval);

#my $SOLO_GAP = 62;

#my $PRINT_QUALS = 1;
#my $PRINT_WALKS = 1;

$beginTime = [gettimeofday];

my %opts = ();
my $validLine = getopts('c:i:s:m:D:g:PVR:Q:NA', \%opts);
my @required = ("c","i","s","m","g");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print "Usage: ./merauder.pl <-i maxInsertSize> <-s scaffoldReportFile> <-m merSize> <-c contigFastaFile> <-g gapDataFile> <<-D minDepth (default 2)>> <<-P(olymorphicMode?)>> <<-V(erbose?)>> <<-R repeatCopyCount>> <<-Q qualOffset (default 33)>> <<-N(ibbleMode?)>> <<-A(ggressiveClosures?)>>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\nmerauder4";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $maxInsertSize = $opts{"i"};
my $scaffReportFile = $opts{"s"};
my $merSize = $opts{"m"};
my $contigFastaFile = $opts{"c"};
my $minDepth = 2;
if (exists($opts{"D"})) {
    $minDepth = $opts{"D"};
}

my $gapDataFile = $opts{"g"};

my $nibbleMode = 0;
if (exists($opts{"N"})) {
    $nibbleMode = 1;
}

my $aggression = 0;
if (exists($opts{"A"})) {
    $aggression = 1;
}

my $qualOffset = 33;
if (exists($opts{"Q"})) {
    $qualOffset = $opts{"Q"};
}

my $polymorphicMode = 0;
if (exists($opts{"P"})) {
    $polymorphicMode = 1;
}

my $verboseMode = 0;
if (exists($opts{"V"})) {
    $verboseMode = 1;
}

my $excludeRepeats = 0;
if (exists($opts{"R"})) {
    $excludeRepeats = $opts{"R"};
}

my $depthInfoAvailable = 0;
my %scaffoldDepth = ();
my %scaffoldLength = ();

my %contigsOfInterest = ();
my %scaffoldsOfInterest = ();


# Decide which contigs/scaffolds to store in memory

$startTime = [gettimeofday];
open (G,$gapDataFile) || die "Couldn't open $gapDataFile\n";
while (my $line = <G>) {
    chomp $line;

    my ($scaffold,$contig1,$primer1,$contig2,$primer2,$gapSizeEstimate,$gapUncertainty,@fillerSeq) = split(/\t/,$line);
    my ($c1ID) = $contig1 =~ /^(.+)\.[35]$/;
    my ($c2ID) = $contig2 =~ /^(.+)\.[35]$/;
    $contigsOfInterest{$c1ID} = 1;
    $contigsOfInterest{$c2ID} = 1;
    $scaffoldsOfInterest{$scaffold} = 1;
}
close G;
$elapsedTime = tv_interval($startTime);
printf STDERR ("Done (%.4f s).\n", $elapsedTime);
printf STDERR ("contigs of interest: %d\n", scalar(keys %contigsOfInterest));
printf STDERR ("scaffolds of interest: %d\n", scalar(keys %scaffoldsOfInterest));

# Read in contig sequence

print STDERR "Reading $contigFastaFile...\n";
$startTime = [gettimeofday];

my %contigSequence = ();
my $currentEntry = undef;
my $seq = "";
my $contigsFound = 0;
open (F,$contigFastaFile) || die "Couldn't open $contigFastaFile\n";
while (my $line = <F>) {
    chomp $line;
    if ($line =~ /^>/) {
        if (defined($currentEntry)) {
            if (exists($contigsOfInterest{$currentEntry})) { 
                $contigSequence{$currentEntry} = $seq;
                $len = length $seq;
                #print "put contig $currentEntry into table, seq $len:\n$seq\n";
                $contigsFound++;
            }
        }
        $seq = "";
        ($currentEntry) = $line =~ />(.+)/;
    } else {
        $seq .= $line;
    }
}
close F;

if (defined($currentEntry)) {
    if (exists($contigsOfInterest{$currentEntry})) { 
        $contigSequence{$currentEntry} = $seq;
    }
}

$elapsedTime = tv_interval($startTime);
printf STDERR ("Done (%.2f s). Found %d contigs\n", $elapsedTime, $contigsFound);

# Read in scaffold report

my @scaffoldIDs = ();
my %scaffoldInfo = ();
my %contigScaffoldMap = ();

print STDERR "Reading $scaffReportFile...\n";

$startTime = [gettimeofday];

open (S,$scaffReportFile) || die "Couldn't open $scaffReportFile\n";
my $contigOfInterest = 0;
while (my $line = <S>) {
    chomp $line;
    my @cols = split(/\t/,$line);
    my $scaffID = $cols[0];
    next unless exists($scaffoldsOfInterest{$scaffID});
    
    unless (exists ($scaffoldInfo{$scaffID})) {
        push(@scaffoldIDs,$scaffID);
        $scaffoldInfo{$scaffID} = [];
        $scaffoldDepth{$scaffID} = [0,0];
        $scaffoldLength{$scaffID} = 0;
    }
    if ($cols[1] =~ /^CONTIG/) {
        $contigOfInterest = 0;

        my ($cStrand,$cName) = $cols[2] =~ /^([\+\-])(.+)$/;
        if (exists($contigsOfInterest{$cName})) { 
            push(@{$scaffoldInfo{$scaffID}},"$cols[2],$cols[3],$cols[4]");
            $contigOfInterest = 1;

            if (exists($contigScaffoldMap{$cName})) {
                die "Current version only supports unique contig->scaffold mappings ($cName)\n";
            } else {
                $contigScaffoldMap{$cName} = "$cStrand,$scaffID,$cols[3],$cols[4]";
            }
        }

        $scaffoldLength{$scaffID} = $cols[4];

        if ($#cols > 4) {
            $depthInfoAvailable = 1;
            my $depth = $cols[5];
            my $weight = $cols[4]-$cols[3]+1;
            $scaffoldDepth{$scaffID}->[0] += $depth*$weight;
            $scaffoldDepth{$scaffID}->[1] += $weight;
        } elsif ($excludeRepeats) {
            warn "Warning: Depth information not available. Repeat exclusion turned off.\n";
            $excludeRepeats = 0;
        }
    } elsif ($cols[1] =~ /^GAP/) {
        if ($contigOfInterest == 1) { 
            push(@{$scaffoldInfo{$scaffID}},"$cols[2],$cols[3]");
        }
    }
}
close S;

my $peakDepth = 0;
if ($depthInfoAvailable) {
    my %depthHist = ();
    foreach my $scaffID (@scaffoldIDs) {
        my $scaffWeight = $scaffoldDepth{$scaffID}->[1];
        my $meanDepth = $scaffoldDepth{$scaffID}->[0] / $scaffWeight;
        $scaffoldDepth{$scaffID} = $meanDepth;
        my $depthBin = sprintf("%.0f",$meanDepth);
        $depthHist{$depthBin} += $scaffWeight;
    }
    my $maxWeight = 0;
    while (my ($depthBin,$binWeight) = each(%depthHist)) {
        if ($binWeight > $maxWeight) {
            $peakDepth = $depthBin;
            $maxWeight = $binWeight;
        }
    }
    print STDERR "Modal scaffold depth $peakDepth\n";
    foreach my $scaffID (@scaffoldIDs) {
        $scaffoldDepth{$scaffID} /= $peakDepth;
    }
}

$elapsedTime = tv_interval($startTime);
printf STDERR ("Done (%.2f s).\n", $elapsedTime);

# Build scaffold data structures

print STDERR "Building data structures ...\n";
$startTime = [gettimeofday];

my %gapInfo = ();
my @gapOrder = ();
my @gapBoundaries = ();
my %scaffGapRange = ();
my $nGaps = 0;

foreach my $scaffID (@scaffoldIDs) {
    my @scaffSegments = @{$scaffoldInfo{$scaffID}};

    my $nScaffSegments = scalar(@scaffSegments);

    if ($nScaffSegments > 1) {

        $scaffGapRange{$scaffID} = "$nGaps.";

        for (my $cIndex = 0; $cIndex < $nScaffSegments-2; $cIndex += 2) {
            my $contigInfo1 = $scaffSegments[$cIndex];
            my ($cID1,$sStart1,$sEnd1) = split(/\,/,$contigInfo1);
            my ($cStrand1,$cName1) = $cID1 =~ /^([\+\-])(.+)$/;
            my $cSeq1 = $contigSequence{$cName1};
            my $contig1 = $cName1;
            if ($cStrand1 eq "-") {
                my $rc = reverse($cSeq1);
                $rc =~ tr/acgtACGT/tgcaTGCA/;
                $cSeq1 = $rc;
                $contig1 .= ".5";
            } else {
                $contig1 .= ".3";
            }
            my $primer1 = substr($cSeq1,-$merSize);
            
            my $gapInfo = $scaffSegments[$cIndex+1];
            my ($gapSize,$gapUncertainty) = split(/\,/,$gapInfo);

            my $contigInfo2 = $scaffSegments[$cIndex+2];
            my ($cID2,$sStart2,$sEnd2) = split(/\,/,$contigInfo2);
            my ($cStrand2,$cName2) = $cID2 =~ /^([\+\-])(.+)$/;
            my $cSeq2 = $contigSequence{$cName2};
            my $contig2 = $cName2;
            if ($cStrand2 eq "-") {
                my $rc = reverse($cSeq2);
                $rc =~ tr/acgtACGT/tgcaTGCA/;
                $cSeq2 = $rc;
                $contig2 .= ".3";
            } else {
                $contig2 .= ".5";
            }
            my $primer2 = substr($cSeq2,0,$merSize);

            $gapInfo{$contig1} = [$nGaps,$scaffID,$contig1,$primer1,$gapSize,$contig2,$primer2,$gapUncertainty];
            $gapInfo{$contig2} = [$nGaps,$scaffID,$contig1,$primer1,$gapSize,$contig2,$primer2,$gapUncertainty];
            push(@gapOrder,$contig1);

            push(@gapBoundaries,"$sEnd1:$sStart2");
            $nGaps++;
        }

        $scaffGapRange{$scaffID} .= "$nGaps";

    }
}

$elapsedTime = tv_interval($startTime);
printf STDERR ("Done (%.4f s). Found $nGaps (potentially) closable gaps.\n", 
               $elapsedTime, $nGaps);

my $nMappedReads = 0;

my %projectors = ();
my %fillers = ();


my %qBins = ();
$qBins{"A0"} = 0;
$qBins{"A1"} = 1;
$qBins{"A2"} = 2;
$qBins{"A3"} = 3;
$qBins{"C0"} = 4;
$qBins{"C1"} = 5;
$qBins{"C2"} = 6;
$qBins{"C3"} = 7;
$qBins{"G0"} = 8;
$qBins{"G1"} = 9;
$qBins{"G2"} = 10;
$qBins{"G3"} = 11;
$qBins{"T0"} = 12;
$qBins{"T1"} = 13;
$qBins{"T2"} = 14;
$qBins{"T3"} = 15;

# Close gaps!

my $mygap = 0;

my $success = 0;
my $nSuccess = 0;
my $nFailure = 0;
my $maxGapZ = 3;
my $gapUncertaintyWiggle = 5.5;

my $maxFillReads = 5000;

my $gi = -1;
open (G,$gapDataFile) || die "Couldn't open $gapDataFile\n";
while (my $line = <G>) {    
    chomp $line;
    $gi++;
    $mygap = $gi;
    if (defined $SOLO_GAP && $mygap != $SOLO_GAP) {
        next;
    }

    my $reportLine = "";
    my $reportStatus = "";
    my $reportNote = "";
    
    $success = 0;
    my ($scaffold,$contig1,$primer1,$contig2,$primer2,$gapSizeEstimate,$gapUncertainty,@fillerSeq) = split(/\t/,$line);
    $reportLine = "REPORT:\t$scaffold\t$contig1\t$primer1\t$contig2\t$primer2\t$gapSizeEstimate\t$gapUncertainty";
    
    my $nFillReads = scalar(@fillerSeq);
    my $gapInfo = "[$nFillReads\t$primer1\t$gapSizeEstimate\t$primer2]";
    
    if ($nFillReads > $maxFillReads) {
        print STDERR "*************\nGap excluded - excessive filler reads. (nFillReads = $nFillReads) $scaffold:$contig1<-[$gapSizeEstimate +/- $gapUncertainty]->$contig2\n";
        $reportStatus = "FAILED\tExcessReads=$nFillReads\n";
        print STDERR "$reportLine\t$reportStatus";
        next;
    }
    
    if ($excludeRepeats) {
        my $scaffDepth = $scaffoldDepth{$scaffold};
        if ($scaffDepth > $excludeRepeats) {
            print STDERR ("*************\nRepeat scaffold gap excluded. (depth = %.6f) $scaffold:$contig1<-[$gapSizeEstimate +/- $gapUncertainty]->$contig2\n", $scaffDepth);
            print STDERR ("$reportLine\tFAILED\tscaffDepth=%.6f\n", $scaffDepth);
            next;
        }
    }
    
    
    print STDERR "*************\nAttempt to close $gi: $scaffold:$contig1<-[$gapSizeEstimate +/- $gapUncertainty]->$contig2\n";

    $startTime = [gettimeofday];

    $gapUncertainty = $maxGapZ*$gapUncertainty + $gapUncertaintyWiggle;
    
    my ($c1ID) = $contig1 =~ /^(.+)\.[35]$/;
    unless (exists($contigSequence{$c1ID})) {
        die "[$c1ID] no seq\n";
    }
    unless (exists($contigScaffoldMap{$c1ID})) {
        die "[$c1ID] no scaffMapInfo\n";
    }
    my $c1Seq = $contigSequence{$c1ID};
    my $c1Info = $contigScaffoldMap{$c1ID};
    my ($c1Strand) = $c1Info =~ /^([\+\-])/;
    if ($c1Strand eq "-") {
        $c1Seq = reverse($c1Seq);
        $c1Seq =~ tr/acgtACGT/tgcaTGCA/;
    }
    
    my ($c2ID) = $contig2 =~ /^(.+)\.[35]$/; 
    unless (exists($contigSequence{$c2ID})) {
        die "[$c2ID] no seq\n";
    }
    unless (exists($contigScaffoldMap{$c2ID})) {
        die "[$c2ID] no scaffMapInfo\n";
    }
    my $c2Seq = $contigSequence{$c2ID};
    my $c2Info = $contigScaffoldMap{$c2ID};
    my ($c2Strand) = $c2Info =~ /^([\+\-])/;
    if ($c2Strand eq "-") {
        $c2Seq = reverse($c2Seq);
        $c2Seq =~ tr/acgtACGT/tgcaTGCA/;
    }
    
    unless (($c1Seq =~ /$primer1$/) && ($c2Seq =~ /^$primer2/)) {
        die "Contig sequence doesn't match primer ($contig1:$primer1 || $contig2:$primer2)\n";
    }

    my $closure = "";
    
    #  Attempt to cross the gap using splinting reads
    
    my $spanClosure = span($primer1,$primer2,\@fillerSeq);
    my $spanCheck = 0;
    
    if ($spanClosure) {
        my $spanCheck = checkClosure($spanClosure,$primer1,$primer2,$gapSizeEstimate,$gapUncertainty);
        $reportNote .= "spanClosure=$spanClosure;spanCheck=$spanCheck;";
        unless ($spanCheck) {
            $closure = $spanClosure;
        }
    }
    
    #  If splints fail try a mer-walk
    my $bridgeClosureRight = 0;
    my $rightWalk = "";
    my $rightFail = "";
    my $rightCheck = 0;
    
    my $bridgeClosureLeft = 0;
    my $leftWalk = "";
    my $leftFail = "";
    my $leftCheck = 0;
    
    my $patchCheck = 0;
    my $patchClosure = 0;
    
    my $nibbleClosure = 0;

    unless ($closure) {
        my $maxLeftWalk = "";
        my $maxRightWalk = bridgeIter($c1Seq,$c2Seq,\@fillerSeq,"right");
        print STDERR "MAXRIGHT [$maxRightWalk]\n";
        ($rightWalk,$rightFail) = $maxRightWalk =~ /^([ACGT]*)(.*)$/;

        unless ($rightFail) {   # Rightward walk succeeded
            $bridgeClosureRight = $rightWalk;
            $rightCheck = checkClosure($rightWalk,$primer1,$primer2,$gapSizeEstimate,$gapUncertainty);
            $reportNote .= "rightClosure=$rightWalk;rightCheck=$rightCheck;";
            unless ($rightCheck) {
                $closure = $bridgeClosureRight;
            }
        }

        unless ($closure) {
            $maxLeftWalk = bridgeIter($c1Seq,$c2Seq,\@fillerSeq,"left");
            print STDERR "MAXLEFT [$maxLeftWalk]\n";
            ($leftWalk,$leftFail) = $maxLeftWalk =~ /^([ACGT]*)(.*)$/;
            $leftWalk = reverse($leftWalk);
            $leftWalk =~ tr/ACGT/TGCA/;
            $leftFail =~ tr/ACGT/TGCA/;
            
            unless ($leftFail) {   # Leftward walk succeeded
                $bridgeClosureLeft = $leftWalk;
                $leftCheck = checkClosure($leftWalk,$primer1,$primer2,$gapSizeEstimate,$gapUncertainty);
                $reportNote .= "leftClosure=$leftWalk;leftCheck=$leftCheck;";
                unless ($leftCheck) {
                    $closure = $bridgeClosureLeft;
                }
            }
        }
    }

    # If walks failed to close, try to patch between left and right walk
    unless ($closure || $bridgeClosureRight || $bridgeClosureLeft) {
        $patchClosure = patch($rightWalk,$leftWalk,$gapSizeEstimate,$gapUncertainty);
        if ($patchClosure) {
            $patchCheck = checkClosure($patchClosure,$primer1,$primer2,$gapSizeEstimate,$gapUncertainty);
            $reportNote .= "patchClosure=$patchClosure;patchCheck=$patchCheck;";
            unless ($patchCheck) {
                $closure = $patchClosure;
            }
        }
    }
    
    if ($nibbleMode && (!$closure)) {
        $nibbleClosure = nibblePatch($rightWalk,$leftWalk,$gapSizeEstimate);
        if ($nibbleClosure) {
            my $nibbleCheck = checkClosure($nibbleClosure,$primer1,$primer2,$gapSizeEstimate,$gapUncertainty);
            $reportNote .= "nibbleClosure=$nibbleClosure;nibbleCheck=$nibbleCheck;";
            unless ($nibbleCheck) {
                $closure = $nibbleClosure;
            }
        }
    }
    
    if ($closure) {
        my $closedGapSize = length($closure) - 2*$merSize;
        if ($nibbleClosure) { 
            $nFailure++;
            print STDERR "$gapInfo failed to close gap $mygap ($spanCheck;$rightCheck:$rightFail;$leftCheck:$leftFail) (gap nibbled)\n";
            $reportStatus = "FAILED\t$reportNote\n";
        } else {
            $success = 1;
            $nSuccess++;
            print STDERR "$gapInfo successfully closed gap $mygap: $closedGapSize $closure\n";
            $reportStatus = "SUCCESS\t$reportNote\n";
        }

        #  Lower case as little of the primer sequences as possible:
        my $minGapMask = 2*5;
        if ($closure =~ /$primer1.{$minGapMask,}$primer2/) {
            my ($p1,$g,$p2) = $closure =~ /($primer1)(.+)($primer2)/;
            $g = lc($g);
            $closure = $p1.$g.$p2;
        } else {
            my $l = length($closure);
            if ($l % 2 == 1) {
                $minGapMask++;
            }
            my $d = ($l-$minGapMask)/2;
            my ($p1,$g,$p2) = $closure =~ /^(.{$d})(.{$minGapMask})(.{$d})/;
            $g = lc($g);
            $closure = $p1.$g.$p2;
        }
        print "$scaffold\t$contig1\t$primer1\t$contig2\t$primer2\t$closure\n";
    } else {
        $nFailure++;
        print STDERR "$gapInfo failed to close gap $mygap ($spanCheck;$rightCheck:$rightFail;$leftCheck:$leftFail)\n";
        $reportStatus = "FAILED\t$reportNote\n";
    }

#    unless($success) {
#        print STDERR "OPEN:$line\n";
#    }

    $elapsedTime = tv_interval($startTime);

    print STDERR "$reportLine\n";
    print STDERR "$reportStatus\n";
    #printf STDERR ("Closing gap $gi, size $gapSizeEstimate took %.6f s\n", $elapsedTime);

}
close G;
print STDERR "Done.  Successfully closed $nSuccess gaps. ($nFailure failed to close)\n";

$elapsedTime = tv_interval($beginTime);
printf STDERR ("Elapsed time %.2f s\n", $elapsedTime);

#
# SUBROUTINES
#

sub printFasta {
    my ($name,$seq) = @_;
    my $seqLen = length($seq);

    my $bpl = 50;
    my $nLines = sprintf("%d",$seqLen/$bpl);
    if ($seqLen%$bpl != 0) {
        $nLines++;
    }

    print ">$name\n";
    for (my $j = 0;$j<$nLines;$j++) {
        my $seqLine = substr($seq,$j*$bpl,$bpl);
        my $text = $seqLine;
        $text.= "\n";
        print $text;
    }
}

sub span {
    my ($primer1,$primer2,$seqRef) = @_;
    
    my @reads = @{$seqRef};
    my $nReads = scalar(@reads);

    if ($verboseMode) {
        print STDERR "Attempting span: $primer1 -> $primer2 ($nReads read(s) available)\n";
    }

    my %pureSpans = ();
    foreach my $read (@reads) {
        my ($nts,$quals) = $read =~ /^([^:]+)\:(.+)$/;
        if (($nts =~ /$primer1/) && ($nts =~ /$primer2/)) {
            my ($tmpSpan) = $nts =~ /($primer1.*)$/;
            my ($span) = $tmpSpan =~ /^(.*$primer2)/; 
            
            if ($span) {
                $pureSpans{$span}++;
            }
        }
    }
    my @spans = keys(%pureSpans);
    my $nSpans = scalar(@spans);

    if ($verboseMode) {
        print STDERR "$nSpans distinct spanning sequence(s) found\n";
    }

    if (($nSpans == 1) && ($pureSpans{$spans[0]} > 1)) {
        print STDERR "Unique spanning sequence found: $pureSpans{$spans[0]}/$nReads reads span the gap.\n";
        return $spans[0];

    } elsif ($polymorphicMode) {
        my $maxSpanIndex = 0;
        my $maxSpanFreq = 0;
        for (my $i=0; $i<$nSpans; $i++) {
            my $spanFreq = $pureSpans{$spans[$i]};
            if ($spanFreq > $maxSpanFreq) {
                $maxSpanIndex = $i;
                $maxSpanFreq = $spanFreq;
            } elsif ($spanFreq == $maxSpanFreq) {
                my $newLen = length($spans[$i]);
                my $prevLen = length($spans[$maxSpanIndex]);
                if ($newLen > $prevLen) {
                    # always prefer longer span - good for consistency checks
                    $maxSpanIndex = $i;
                } elsif ($newLen == $prevLen && ($spans[$i] cmp $spans[$maxSpanIndex]) < 0) {
                    # same length span, same frequency. Prefer lexicographically smallest
                    $maxSpanIndex = $i;
                }
            }
        }

        if ($maxSpanFreq > 1) {
            print STDERR "Maximum frequency spanning sequence found: $pureSpans{$spans[$maxSpanIndex]}/$nReads reads span the gap.\n";
            return $spans[$maxSpanIndex];
        } else {
            if ($verboseMode) {
                print STDERR "No spanning sequence with frequency > 1 found\n";
            } 
            return 0;
        }

    } else {
        if ($verboseMode) {
            print STDERR "No unique spanning sequence with frequency > 1 found\n";
        } 
        return 0;
    }
}


# Iterative bridging allows variable k-mer size

sub bridgeIter {
    my ($c1seq,$c2seq,$seqRef,$direction) = @_;
    
    my $minQ = 19;
    my $minMerLen = 13;
    my $maxReadLength = 0;

    my @reads = @{$seqRef};
    my $nReads = scalar(@reads);

    if ($verboseMode) {
        print STDERR "Attempting $direction bridge: ($nReads read(s) available)\n";
    }

    my $merLen = $merSize;

    my $iterating = 1;
    my $downshift = 0;
    my $upshift = 0;


    if ($direction eq "left") {
        $c1seq = reverse($c1seq);
        $c1seq =~ tr/ACGT/TGCA/;
        $c2seq = reverse($c2seq);
        $c2seq =~ tr/ACGT/TGCA/;
        ($c1seq,$c2seq) = ($c2seq,$c1seq);

        for (my $r=0; $r < $nReads; $r++) {
            my $read = $reads[$r];
            my ($nts,$quals) = $read =~ /^([^:]+)\:(.+)$/;
            $quals = reverse($quals);
            $nts = reverse($nts);
            $nts =~ tr/ACGT/TGCA/;
            $reads[$r] = "$nts:$quals";
        }

    }

#    for ($i = 0; $i < $nReads; $i++) {
#        printf("%d %s\n", $i, $reads[$i]);
#    }

#    return "";


    my $maxWalk = undef;
    my $maxWalkLen = 0;
    my ($truePrimer1) = $c1seq =~ /(.{$merSize})$/;
    my ($truePrimer2) = $c2seq =~ /^(.{$merSize})/;

    while ($iterating) {

        my ($primer1) = $c1seq =~ /(.{$merLen})$/;
        my ($primer2) = $c2seq =~ /^(.{$merLen})/;

        # Allow k-mers less than merSize to be used 

        if ($merLen < $merSize) {
            my $theRest = $merSize-$merLen;
            ($primer1) = $c1seq =~ /(.{$merLen}).{$theRest}$/;
            ($primer2) = $c2seq =~ /^.{$theRest}(.{$merLen})/;
        }

        unless ($primer1 && $primer2) {
            if ($verboseMode) {
                print STDERR "Unable to find k=$merLen seeds\n";
            }
            if (defined($maxWalk)) {
                return $maxWalk;
            } else {
                return "X";
            }
        }

        my %mers = ();
        my %fullInfo = ();

        my $idx = 0;
        foreach my $read (@reads) {

            my ($nts,$quals) = $read =~ /^([^:]+)\:(.+)$/;
            my $readLength = length($nts);
            if ($readLength > $maxReadLength) {
                $maxReadLength = $readLength;
            }
            if ($merLen >= $readLength) {
                next;
            }

            my $merPlus = $merLen+1;
            while ($nts =~ /([ACGT]{$merPlus})/g) {
                my $merExt = $1;

                my ($mer,$extension) = $merExt =~ /^([^:]+)(.)$/;

                my $q = substr($quals,pos($nts)-1,1);
                my $Q = ord($q)-$qualOffset;

                my $qBin = sprintf("%d",$Q/10);
                if ($qBin > 3) {
                    $qBin = 3;
                }

                $qBin = $extension.$qBin;
                unless (exists($fullInfo{$mer})) {
                    $fullInfo{$mer} = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
                }
                $fullInfo{$mer}->[$qBins{$qBin}]++;

                pos($nts) -= $merLen;
            }
            $idx++;
        }

        if (defined $PRINT_QUALS && $PRINT_QUALS) {
            while (($key, $val) = each (%fullInfo)) {
                printf("$merLen $mygap $key ");
                for ($j = 0; $j < 16; $j++) {
                    printf("%d ", $val->[$j]);
                }
                print "\n";
            }
        }

# Full analysis of quality/frequency profile

        while (my ($mer,$extMatrix) = each(%fullInfo)) {
            my @extensions = @{$extMatrix};

            my %ext = ();
            $ext{"A"} = [@extensions[0..3]];
            $ext{"C"} = [@extensions[4..7]];
            $ext{"G"} = [@extensions[8..11]];
            $ext{"T"} = [@extensions[12..15]];
            
            my %nExt = ();
            my %nHiQ = ();
            my %ratings = ();
            my $nTotal = 0;
            my @bases = ("A","C","G","T");
            
            foreach my $base (@bases) {
#    Ignore q<10 bases
#        $nExt{$base} = $ext{$base}->[0]+$ext{$base}->[1]+$ext{$base}->[2]+$ext{$base}->[3];

                $nExt{$base} = $ext{$base}->[1]+$ext{$base}->[2]+$ext{$base}->[3];
                $nHiQ{$base} = $ext{$base}->[2]+$ext{$base}->[3];
                $nTotal += $nExt{$base};
                $ratings{$base} = categorizeExtension(@{$ext{$base}});
            }

            my @sortedBases = sort {($ratings{$b} <=> $ratings{$a}) || 
                                        ($nHiQ{$b} <=> $nHiQ{$a}) || 
                                        ($nExt{$b} <=> $nExt{$a})} (@bases);

            my @sortedRatings = @ratings{@sortedBases};
            

#            print "$mer @sortedBases | @sortedRatings\n";

# Rules for choosing next base
# ratings:
# 0 = No votes
# 1 = One vote
# 2 = nVotes < minViable
# 3 = minDepth > nVotes >= minViable, nHiQ < minViable
# 4 = minDepth > nVotes >= minViable ; nHiQ >= minViable
# 5 = nVotes >= minDepth ; nHiQ < minViable
# 6 = nVotes >= minDepth ; minViable <= nHiQ < minDepth 
# 7 = nHiQ >= minDepth 

            my $topRating = $sortedRatings[0];
            my $runnerUp = $sortedRatings[1];
            my $topRatedBase = $sortedBases[0];

            if ($topRating < 3) {         # must have at least minViable bases
                $mers{$mer} = "X";
            } elsif ($topRating == 3) {    # must be uncontested   
                if ($runnerUp == 0) {
                    $mers{$mer} = $topRatedBase;
                } else {
                    $mers{$mer} = "X";
                }
            } elsif ($topRating < 6) {
                if ($runnerUp < 3) {
                    $mers{$mer} = $topRatedBase;
                } else {
                    $mers{$mer} = "X";
                }   
            } elsif ($topRating == 6) {  # viable and fair hiQ support
                if ($runnerUp < 4) {
                    $mers{$mer} = $topRatedBase;
                } else {
                    $mers{$mer} = "X";
                }
            } else {                     # strongest rating trumps
                if ($runnerUp < 7) {       
                    $mers{$mer} = $topRatedBase;
                } else {
                    my $fork = "F";
                    for (my $b = 0; $b < 4; $b++) {
                        my $base = $sortedBases[$b];
                        my $rating = $sortedRatings[$b];
                        my $n = $nExt{$base};
                        if ($rating == 7) {
                            $fork .= "$base$n";
                        } else {
                            last;
                        }
                    }
                    $mers{$mer} = $fork;
                }
            }
        }

#        while (my ($key, $val) = each(%mers)) {
#            printf(">>%d '%s' '%s'\n", $mygap, $key, $val);
#        }
        
        my $walk = $primer1;
        unless (defined($maxWalk)) {
            $maxWalk = $walk;
            $maxWalkLen = length($maxWalk);
        }
        my $step = $primer1;
        my %loopCheck = ();
        my $success = 0;
        my $nForks = 0;
        my $nSteps = 0;
        while () {
            if (exists($loopCheck{$step})) {
                $walk .= "R";
                last;
            } else {
                $loopCheck{$step} = 1;
            }

            if (exists($mers{$step})) {
                my $next = $mers{$step};

                $nSteps++;
                if ($verboseMode) {
                    print STDERR "$nSteps : $step->$next\t";
                    print STDERR "[@{$fullInfo{$step}}[0..3]][@{$fullInfo{$step}}[4..7]]";
                    print STDERR "[@{$fullInfo{$step}}[8..11]][@{$fullInfo{$step}}[12..15]]\n";
                }

                if ($next =~ /^([ACGT])$/) {
                    $step = substr($step,1);
                    $step .= $next;
                    $walk .= $next;

                    if ($walk =~ /$truePrimer2$/) {
                        $success = 1;
                        last;
                    }

                } elsif (($polymorphicMode == 1) && ($nForks == 0) && ($next =~ /^F/)) {
                    #  Biallelic positions only (for now) maximum vote path is taken
                    if ($next =~ /^F[ACGT]\d+[ACGT]\d+$/) {
                        my ($o1,$v1,$o2,$v2) = $next =~ /^F([ACGT])(\d+)([ACGT])(\d+)$/;
                        $next = ($v1 > $v2) ? $o1 : $o2;
                        if ($verboseMode) {
                            print STDERR "Polymorphic conditions met .. attempting max frequency resolution.\n";
                        }
                    } else {
                        $walk .= $next;
                        last;
                    }

                    $step = substr($step,1);
                    $step .= $next;
                    $walk .= $next;
                    $nForks++;

                    if ($walk =~ /$truePrimer2$/) {
                        $success = 1;
                        last;
                    }
                } else {
                    $walk .= $next;
                    last;
                }
            } else {
                $walk .= "X";
                last;
            }
        }

        if (defined $PRINT_WALKS && $PRINT_WALKS) {
            print "$mygap $walk\n";
        }

        # Trim off extra lead bases if upshifted
        my $additionalBases = $merLen-$merSize;
        if ($additionalBases > 0) {
            my $reducedWalk = $walk;
            ($walk) = $reducedWalk =~ /^.{$additionalBases}(.+)$/;
        }
        
        my ($walkBases) = $walk =~ /^([ACGT]*)/;
        my $walkLen = length($walkBases);

        if ($success || ($walkLen > $maxWalkLen)) {
            $maxWalk = $walk;
            $maxWalkLen = $walkLen;
        }
        
        if (($walk =~ /F/) || ($walk =~ /R/)) {
            $merLen += 2;
            $upshift = 1;

            if (($downshift==1) || ($merLen >= $maxReadLength)) {
                return($maxWalk);
            }
            if ($verboseMode) {
                print STDERR "Degeneracy encountered; upshifting (k->$merLen)\n";
            }
        } elsif ($walk =~ /X/) {
            $merLen -= 2;
            $downshift = 1;

            if (($upshift==1) || ($merLen < $minMerLen)) {
                return($maxWalk);
            }
            if ($verboseMode) {
                print STDERR "Termination encountered; downshifting (k->$merLen)\n";
            }
        } else {
            return($maxWalk);
        }
    }
}

sub categorizeExtension {
    my (@ext) = @_;

    my $minViable = 3;
    if ($minViable > $minDepth) {
        warn "Warning: in categorizeExtension minViable reset to match minDepth ($minDepth)\n";
        $minViable = $minDepth;
    }

# 0 = No votes
# 1 = One vote
# 2 = nVotes < minViable
# 3 = minDepth > nVotes >= minViable, nHiQ < minViable
# 4 = minDepth > nVotes >= minViable ; nHiQ >= minViable
# 5 = nVotes >= minDepth ; nHiQ < minViable
# 6 = nVotes >= minDepth ; minViable < nHiQ < minDepth 
# 7 = nHiQ >= minDepth 

#    Ignore q<10 bases
#    my $n = $ext[0]+$ext[1]+$ext[2]+$ext[3];
    my $n = $ext[1]+$ext[2]+$ext[3];
    my $nHiQ = $ext[2]+$ext[3];
    my $nLoQ = $n-$nHiQ;

    my $category = undef;
    if ($n == 0) {
        $category = 0;
    } elsif ($n == 1) {
        $category = 1;
    } elsif ($n < $minViable) {
        $category = 2;
    } else {
        if (($n < $minDepth) || ($n == $minViable)) {
            if ($nHiQ < $minViable) {
                $category = 3;
            } else {
                $category = 4;
            }
        } else {
            if ($nHiQ < $minViable) {
                $category = 5;
            } elsif ($nHiQ < $minDepth) {
                $category = 6;
            } else {
                $category = 7;
            }
        }
    }

    unless (defined($category)) {
        die "Undefined extension category for [@ext]\n";
    }
    return $category;
};

sub checkClosure {
    my ($closureToCheck,$primer1,$primer2,$gapSizeEstimate,$gapUncertainty) = @_;

    my $badClosure = 0;

    if (! (($closureToCheck =~ /^$primer1/) && ($closureToCheck =~ /$primer2$/))) {
        $badClosure += 1;
        if ($verboseMode) {
            print STDERR "closure [$closureToCheck] rejected because it disagrees with primers\n";
        }
    }
    
    my $closedGapSize = length($closureToCheck) - 2*$merSize;
    my $gapDiff = $closedGapSize - $gapSizeEstimate;
    
    if (abs($gapDiff) > $gapUncertainty) {
        $badClosure += 2;
        if ($verboseMode) {
            printf STDERR ("closure [$closureToCheck] rejected due to gap estimate differential: |$gapDiff| > %.4f\n", 
                           $gapUncertainty);
        }
    }

    # In aggressive mode ignore size differences between closure and estimate
    if ($aggression && ($badClosure == 2) && ($gapSizeEstimate < 2*$maxInsertSize))
    {
        $badClosure = 0;
        if ($verboseMode) {
            print STDERR "closure [$closureToCheck] allowed via AGRESSIVE mode\n";
        }
    }

    return $badClosure;

}

sub nibblePatch {
    my ($rightWalk,$leftWalk,$gapSizeEstimate) = @_;

    my $rightLen = length($rightWalk);
    my $leftLen = length($leftWalk);

    my $minGapSize = 10;
    my $gapPad = "N"x$minGapSize;

    my $idealLength = 2*$merSize+$gapSizeEstimate;

    if ($rightLen + $leftLen + $minGapSize <=  $idealLength) {

        my $extraPadLength = $idealLength - ($rightLen+$leftLen+$minGapSize);
        if ($extraPadLength > 0) {
            $gapPad .= "N"x$extraPadLength;
        }

        return "$rightWalk$gapPad$leftWalk";

    } else {

        my $extraLength = ($rightLen+$leftLen+$minGapSize)-$idealLength;
        if ($rightLen > $leftLen) {
            my $lengthDiff = $rightLen-$leftLen;
            if ($lengthDiff >= $extraLength) {
                my ($trimmedRightWalk) = $rightWalk =~ /^(.+).{$extraLength}$/;
                $rightWalk = $trimmedRightWalk;

                return "$rightWalk$gapPad$leftWalk";

            } else {
                my ($trimmedRightWalk) = $rightWalk =~ /^(.+).{$lengthDiff}$/;
                $rightWalk = $trimmedRightWalk;
            }
        } elsif ($leftLen > $rightLen) {
            my $lengthDiff = $leftLen-$rightLen;
            if ($lengthDiff >= $extraLength) {
                my ($trimmedLeftWalk) = $leftWalk =~ /^.{$extraLength}(.+)$/;
                $leftWalk = $trimmedLeftWalk;

                return "$rightWalk$gapPad$leftWalk";

            } else {
                my ($trimmedLeftWalk) = $leftWalk =~ /^.{$lengthDiff}(.+)$/;
                $leftWalk = $trimmedLeftWalk;
            }
        }

        $rightLen = length($rightWalk);
        $leftLen = length($leftWalk);

        my $rightTurn = 1;
        while ($rightLen + $leftLen + $minGapSize >  $idealLength) {
            if ($rightTurn) {
                chop $rightWalk;
                $rightLen--;
                $rightTurn = 0;
            } else {
                $leftWalk = substr($leftWalk,1);
                $leftLen--;
                $rightTurn = 1;
            }
        }

        return "$rightWalk$gapPad$leftWalk";
    }
}




sub patch {
    my ($rightWalk,$leftWalk,$gapSizeEstimate,$gapUncertainty) = @_;

    my $minAcceptableOverlap = 10;

    my $rightLen = length($rightWalk);
    my $leftLen = length($leftWalk);

    my $minOffset = sprintf("%.0f",-$gapUncertainty);
    my $maxOffset = sprintf("%.0f",$gapUncertainty);

    my $r0 = 1;
    my $r1 = $rightLen;
    my $l1 = 2*$merSize+$gapSizeEstimate;
    my $l0 = $l1-$leftLen+1;

    if ($verboseMode) {
        print STDERR "Attempting patch [$rightLen][$leftLen] (gap: $gapSizeEstimate +/- $gapUncertainty)\n";
    }
    
    my %overlaps = ();
    for (my $o = $minOffset; $o <= $maxOffset; $o++) {

        my $testRight = $rightWalk;
        my $testLeft = $leftWalk;

        my $l0shifted = $l0 + $o;
        my $l1shifted = $l1 + $o;
        
        my $chopLeft = $r0 - $l0;
        my $chopRight = $r1 - $l1;

        if ($chopLeft > 0) {
            my ($chop,$keep) = $testLeft =~ /^(.{$chopLeft})(.*)$/;
            $testLeft = $keep;
        }
        my $tlLen = length($testLeft);

        if ($chopRight > 0) {
            my ($keep,$chop) = $testRight =~ /^(.*)(.{$chopRight})$/;
            $testRight = $keep;
        }
        my $trLen = length($testRight);

#        if ($o == $minOffset) {
#            printf("%d test_right %s test_left %s chop_left %d chop_right %d tr_len %d tl_len %d gu %f mino %d maxo %d mg %d\n", 
#                   $mygap, $testRight, $testLeft, $chopLeft, $chopRight, $trLen, $tlLen, $gapUncertainty, $minOffset, $maxOffset, $l1);
#        }

        my $overlap = $trLen+$tlLen-(2*$merSize+$gapSizeEstimate+$o);

        unless (($overlap >= $minAcceptableOverlap) && ($overlap <= $trLen) && ($overlap <= $tlLen)) {
            next;
        }

        my ($p1Prefix,$p1Suffix) = $testRight =~ /^(.*)(.{$overlap})$/;
        my ($p2Prefix,$p2Suffix) = $testLeft =~ /^(.{$overlap})(.*)$/;
        if ($p1Suffix eq $p2Prefix) {
            $overlaps{$o} = $testRight.$p2Suffix;
#            printf("%d %d %s\n", $mygap, $o, $overlaps{$o})
        }
    }

    my $bestGuessSeq = "";
    my $bestGuessDelta = $gapUncertainty+1;
    my $nBestGuesses = 0;
    my $nGuesses = 0;
    while (my ($offset,$guessSeq) = each(%overlaps)) {
        my $delta = abs($offset);
        $nGuesses++;
        if ($delta < $bestGuessDelta) {
            $bestGuessSeq = $guessSeq;
            $bestGuessDelta = $delta;
            $nBestGuesses = 1;
        } elsif ($delta == $bestGuessDelta) {
            if (length($guessSeq) < length($bestGuessSeq)) {
                $bestGuessSeq = $guessSeq;
                $bestGuessDelta = $delta;
            }
            $nBestGuesses++;
        }
    }

    if ($nGuesses) {
        print STDERR "$nGuesses potential patches identified.  $nBestGuesses best guess(es) differ from gap estimate by $bestGuessDelta\n";
        return $bestGuessSeq;
    } else {
        print STDERR "No valid patches found.\n";
        return 0;
    }
}


