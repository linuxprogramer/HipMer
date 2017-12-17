#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# oNo.pl by Jarrod Chapman <jchapman@lbl.gov> Tue May 26 09:02:36 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
my %opts = ();
my $validLine = getopts('m:l:p:s:c:B:', \%opts);
my @required = ("m","l");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);

$validLine &= (exists($opts{"s"}) || exists($opts{"c"}));

if ($validLine != 1) {
    print "Usage: ./oNo4.pl <-m merSize> <-l linkDataFile> (<-s scaffoldReportFile> || <-c contigReportFile>)  <<-p pairThreshold>> <<-B blessedLinkFile>>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\noNo4";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $merSize = $opts{"m"};
my $fastaFile = $opts{"f"};
my $linkDataFile = $opts{"l"};
my $ceaFile = $opts{"e"};
my $cmdFile = $opts{"d"};

my $srfFile = undef;
my $crfFile = undef;
if (exists($opts{"s"})) {
    $srfFile = $opts{"s"};
} elsif (exists($opts{"c"})) {
    $crfFile = $opts{"c"};
}

my $pairThreshold = 2;
if (exists($opts{"p"})) {
    $pairThreshold = $opts{"p"};
}

my $blessedLinkFile = undef;
my %blessedLinks = ();
if (exists($opts{"B"})) {
    $blessedLinkFile = $opts{"B"};
    open (B,$blessedLinkFile) || die "Couldn't open $blessedLinkFile\n";
    while (my $line = <B>) {
	chomp $line;
	my ($s1,$s2) = split(/\s+/,$line);
	my $bless = ($s1 lt $s2) ? "$s1.$s2" : "$s2.$s1";
	$blessedLinks{$bless} = 1;
    }
    close B;
}

# Read linkData file and store library information

my %libInfo = ();
my @libs = ();
open (L,$linkDataFile) || die "Couldn't open $linkDataFile\n";
while (my $line = <L>) {
    chomp $line;
    my ($libs,$insertSizes,$stdDevs,$linkFile) = split(/\t/,$line);
    my @l = split(/\,/,$libs);
    my @i = split(/\,/,$insertSizes);
    my @s = split(/\,/,$stdDevs);

    my $nLibs = scalar(@l);
    for (my $l = 0; $l < $nLibs; $l++) {
	$libInfo{$l[$l]} = [$i[$l],$s[$l],$linkFile];
	push(@libs,$l[$l]);
    }
}
close L;

my @sortedLibs = sort {$libInfo{$a}->[0] <=> $libInfo{$b}->[0]} @libs;

my %linkFiles = ();
my @linkFileOrder = ();
foreach my $lib (@sortedLibs) {
    my $linkFile = $libInfo{$lib}->[2];
    unless (exists($linkFiles{$linkFile})) {
	push(@linkFileOrder,$linkFile);
	$linkFiles{$linkFile} = 1;
    }
}

my $nLinkFiles = scalar(@linkFileOrder);

my %onoObjectLengths = ();
my %contigDepths = ();
my $depthInfoAvailable = 0;
my $peakDepth = undef;

my $maxInsertSize = $libInfo{$sortedLibs[-1]}->[0];
print STDERR "Maximum insert size: $maxInsertSize\n";
my $suspendable = $maxInsertSize/2;


if (defined($crfFile)) {

    print STDERR "Reading $crfFile...\n";

    open (C,$crfFile) || die "Couldn't open $crfFile\n";
    while (my $line = <C>) {
	chomp $line;
	my @cols = split(/\t/,$line);

	my $contig= $cols[0];
	my $length = $cols[1];
	
	$onoObjectLengths{$contig} = $length;

	if ($#cols > 1) {
	    my $depth = $cols[2];
	    $contigDepths{$contig} = $depth;
	    $depthInfoAvailable = 1;
	}

    }

    close C;

    print STDERR "Done.\n";

    if ($depthInfoAvailable) {
	my %depthHist = ();
	while (my ($c,$d) = each(%contigDepths)) {
	    my $roundedDepth = sprintf("%.0f",$d);
	    $depthHist{$roundedDepth} += $onoObjectLengths{$c};
	}
	my $maxBases = 0;
	while (my ($d,$b) = each(%depthHist)) {
	    if ($b > $maxBases) {
		$peakDepth = $d;
		$maxBases = $b;
	    }
	}

	print STDERR "Depth information available.  Modal depth is $peakDepth.\n";
    }
}

# If an existing scaffold report exists, load scaffolding information
 
my %scaffReport = ();  # scaffoldName -> srfFileLines 

if (defined($srfFile)) {
    
    print STDERR "Reading scaffold report file: $srfFile...\n";
    open (S,$srfFile) || die "Couldn't open $srfFile\n";

    my %scaffRealBases = ();
    my %scaffDepth = ();
    my %depthHist = ();

    while (my $line = <S>) {
        chomp $line;
        my @cols = split(/\t/,$line);
	
	if (exists($scaffReport{$cols[0]})) {
	    $scaffReport{$cols[0]} .= "$line\n";
	} else {
	    $scaffReport{$cols[0]} = "$line\n";
	}

        if ($cols[1] =~ /^CONTIG/) {
            my ($scaffID,$contigID,$sStart,$sEnd) = @cols[0,2,3,4];

	    my ($contigName) = $contigID =~ /^[\+\-](.+)$/;

	    if ($#cols > 4) {
		my $depth = $cols[5];
		my $cLen = $sEnd-$sStart+1;
		$contigDepths{$contigName} = $depth;
		$depthInfoAvailable = 1;
		$scaffRealBases{$scaffID} += $cLen;
		$scaffDepth{$scaffID} += $cLen*$depth;
		my $roundedDepth = sprintf("%.0f",$depth);
		$depthHist{$roundedDepth} += $cLen;
	    }

            if (exists($onoObjectLengths{$scaffID})) {
                if ($sEnd > $onoObjectLengths{$scaffID}) {
                    $onoObjectLengths{$scaffID} = $sEnd;
                }
            } else {
                $onoObjectLengths{$scaffID} = $sEnd;
            }
	}
    }
    close S;
    
    print STDERR "Done.\n";

    if ($depthInfoAvailable) {
	my $maxBases = 0;
	while (my ($d,$b) = each(%depthHist)) {
	    if ($b > $maxBases) {
		$peakDepth = $d;
		$maxBases = $b;
	    }
	}

	print STDERR "Depth information available.  Modal depth is $peakDepth.\n";

	while (my ($s,$d) = each(%scaffDepth)) {
	    $contigDepths{$s} = $d/$scaffRealBases{$s};
	}
    }


}


# Build up contig linkage table using linkFiles

my %links = ();
for (my $lf = 0; $lf < $nLinkFiles; $lf++) {
    my $linkFile = $linkFileOrder[$lf];
    print STDERR "Reading link file $linkFile...";

    open (L,$linkFile) || die "Couldn't open $linkFile\n";
    while (my $line = <L>) {
	chomp $line;
	my @cols = split(/\t/,$line);
	my $linkType = $cols[0];
	my $ends = $cols[1];

	my @linkData = @cols[2..$#cols];
	my $linkInfo = join(":",$lf,$linkType,@linkData);
	if (exists($links{$ends})) {
	    push(@{$links{$ends}},$linkInfo);
	} else {
	    $links{$ends} = [$linkInfo];
	}
    }
    close L;

    print STDERR "Done.\n";
}


# Build ties between contigs consolidating splint/span links

my %endTies = ();

foreach my $linkedEnds (keys(%links)) {

    my @linkData = @{$links{$linkedEnds}};
    my @splints = ();
    my @spans = ();

    my $maxLikelihoodSplint = undef;
    my $maxSplintCount = 0;
    my $minUncertaintySpan = undef;
    my $minSpanUncertainty = 100000;

    my $blessedLink = 0;
    if (defined($blessedLinkFile)) {
	my ($s1,$s2) = $linkedEnds =~ /^(.+)\.[35]<=>(.+).[35]$/;
	my $bless = ($s1 lt $s2) ? "$s1.$s2" : "$s2.$s1";
	if (exists($blessedLinks{$bless})) {
	    $blessedLink = 1;
	    print STDERR "BLESSED: $linkedEnds\n";
	}
    }

    foreach my $link (@linkData) {
	my ($linkFile,$linkType,@linkInfo) = split(/:/,$link);
	if ($linkType eq "SPLINT") {
	    
	    my ($splintCounts,$gapEstimate) = @linkInfo;
	    my ($nUsedSplints,$nAnomalousSplints,$nAllSplints) = $splintCounts =~ /^(\d+)\|(\d+)\|(\d+)$/;
	    if ( ($nUsedSplints >= $pairThreshold) || $blessedLink) {
		my $splint = "$linkFile:$nUsedSplints:$gapEstimate";
		push(@splints,$splint);
		unless (defined($maxLikelihoodSplint)) {
		    $maxLikelihoodSplint = $splint;
		    $maxSplintCount = $nUsedSplints;
		}
		if ($nUsedSplints > $maxSplintCount) {
		    $maxLikelihoodSplint = $splint;
		    $maxSplintCount = $nUsedSplints;
		}
	    }

	} elsif ($linkType eq "SPAN") {

	    my ($spanCounts,$gapEstimate,$gapUncertainty) = @linkInfo;
	    my ($nAnomalousSpans,$nUsedSpans) = $spanCounts =~ /^(\d+)\|(\d+)$/;
	    if (($nUsedSpans >= $pairThreshold) || $blessedLink) {
		my $span = "$linkFile:$nUsedSpans:$gapEstimate:$gapUncertainty";
		push(@spans,$span);
		unless (defined($minUncertaintySpan)) {
		    $minUncertaintySpan = $span;
		    $minSpanUncertainty = $gapUncertainty;
		}
		if ($gapUncertainty < $minSpanUncertainty) {
		    $minUncertaintySpan = $span;
		    $minSpanUncertainty = $gapUncertainty;
		}
	    }

	} else {
	    die "Invalid linkType [$linkType] in link: $link\n";
	}
    }

    my $minimumGapSize = -($merSize-2);

    my $nSplintGap = undef;
    my $splintGapEstimate = undef;
    my $maxSplintError = 2;
    my $splintAnomaly = 0;
    if (defined($maxLikelihoodSplint)) {
	my ($l0,$n0,$g0) = split(/:/,$maxLikelihoodSplint);
	$splintGapEstimate = ($g0 < $minimumGapSize) ? $minimumGapSize : $g0;
	$nSplintGap = $n0;
	foreach my $splint (@splints) {
	    my ($l,$n,$g) = split(/:/,$splint);
	    my $splintGapDiff = abs($g-$splintGapEstimate);
	    if ($splintGapDiff > $maxSplintError) {
		$splintAnomaly = 1;
		warn "splint anomaly for $linkedEnds ($l,$n,$g) not consistent with ($l0,$n0,$g0) [$splintGapEstimate]\n";
	    }
	}
    }

    my $nSpanGap = undef;
    my $spanGapEstimate = undef;
    my $spanGapUncertainty = undef;
    my $maxSpanZ = 3;
    my $spanAnomaly = 0;
    if (defined($minUncertaintySpan)) {
	my ($l0,$n0,$g0,$u0) = split(/:/,$minUncertaintySpan);
	$spanGapEstimate = ($g0 < $minimumGapSize) ? $minimumGapSize : $g0;
	$spanGapUncertainty = ($u0 < 1) ? 1 : $u0;
	$nSpanGap = $n0;
	foreach my $span (@spans) {
	    my ($l,$n,$g,$u) = split(/:/,$span);
	    my $spanGapZ = abs($g-$spanGapEstimate)/($u+1);
	    if ($spanGapZ > $maxSpanZ) {
		$spanAnomaly = 1;
		warn "span anomaly for $linkedEnds ($l,$n,$g,$u) not consistent with ($l0,$n0,$g0,$u0) [$spanGapEstimate:$spanGapUncertainty]\n";
	    }
	}
    }

    my $gapEstimate = undef;
    my $gapUncertainty = undef;
    my $nGapLinks = undef;

    if (defined($splintGapEstimate)) {
	$gapEstimate = $splintGapEstimate;
	$gapUncertainty = 1;
	$nGapLinks = $nSplintGap;
    } elsif (defined($spanGapEstimate)) {
	$gapEstimate = $spanGapEstimate;
	$gapUncertainty = $spanGapUncertainty;
	$nGapLinks = $nSpanGap;
    } else {
	next;
    }

    if (defined($spanGapEstimate) && defined($splintGapEstimate)) {
	my $gapZ = ($splintGapEstimate-$spanGapEstimate)/$spanGapUncertainty;
	if (abs($gapZ) > $maxSpanZ) {
	    warn "splint/span anomaly for $linkedEnds ($nSplintGap:$splintGapEstimate) not consistent with ($nSpanGap:$spanGapEstimate:$spanGapUncertainty)\n";
	}
    }

    my ($end1,$end2) = $linkedEnds =~ /(.+)\<\=\>(.+)/;
    
    if (exists($endTies{$end1})) {
	push (@{$endTies{$end1}},"$nGapLinks.$gapEstimate:$gapUncertainty.$end2");
    } else {
	$endTies{$end1} = ["$nGapLinks.$gapEstimate:$gapUncertainty.$end2"];
    }
    
    if (exists($endTies{$end2})) {
	push (@{$endTies{$end2}},"$nGapLinks.$gapEstimate:$gapUncertainty.$end1");
    } else {
	$endTies{$end2} = ["$nGapLinks.$gapEstimate:$gapUncertainty.$end1"];
    }
}

my @sortedByLen = sort {$onoObjectLengths{$b} <=> $onoObjectLengths{$a}} (keys(%onoObjectLengths));

my %endLabels = ();
$endLabels{"UNMARKED"} = 0;
$endLabels{"SELF_TIE"} = 1;
$endLabels{"TIE_COLLISION"} = 2;
$endLabels{"DEPTH_DROP"} = 3;
$endLabels{"NO_TIES"} = 4;
$endLabels{"MULTIPLE_BEST_TIE"} = 5;
my %endMarks = ();

print STDERR "Marking ends...\n";

foreach my $piece (@sortedByLen) {
    my $end = "$piece.5";
    my $endMark = markEnd($end);
    $endMarks{$end} = $endLabels{$endMark};
    
    $end = "$piece.3";
    $endMark = markEnd($end);
    $endMarks{$end} = $endLabels{$endMark};
}

print STDERR "Done marking ends.\n";

my %suspended = ();
my %bestTie = ();
my %bestTiedBy = ();

print STDERR "Finding best ties...\n";

foreach my $piece (@sortedByLen) {
    if (exists($suspended{$piece})) {
	print STDERR "BESTTIE: (suspended)\t$piece\t$suspended{$piece}\n";
	next;
    }

    my $end = "$piece.5";
    my $endMark = $endMarks{$end};
    my $bestTie = $endMark;
    unless ($endMark) {
	$bestTie = bestTie($end);
	$bestTie{$end} = $bestTie;

	if ($bestTie) {
	    if (exists($bestTiedBy{$bestTie})) {
		$bestTiedBy{$bestTie} .= ",$end";
	    } else {
		$bestTiedBy{$bestTie} .= $end;
	    }
	}

    }
    my $notice = "BESTTIE: [$bestTie]\t$piece\t";

    $end = "$piece.3";
    $endMark = $endMarks{$end};
    $bestTie = $endMark;
    unless ($endMark) {
	$bestTie = bestTie($end);
	$bestTie{$end} = $bestTie;

	if ($bestTie) {
	    if (exists($bestTiedBy{$bestTie})) {
		$bestTiedBy{$bestTie} .= ",$end";
	    } else {
		$bestTiedBy{$bestTie} .= $end;
	    }
	}
    }
    $notice .= "[$bestTie]\n";

    print STDERR $notice;

}

print STDERR "Done finding best ties\n";

# Place suspended where possible 

print STDERR "Inserting suspensions...\n";

my $nSuspendedPieces = keys(%suspended);

my $nInsertedSuspensions = 0;
foreach my $piece (@sortedByLen) {
    unless (exists($suspended{$piece})) {
	next;
    }

    my $end5 = "$piece.5";
    my $endMark5 = $endMarks{$end5};
    my $end3 = "$piece.3";
    my $endMark3 = $endMarks{$end3};

    if ($endMark3 || $endMark5) {
	print STDERR "UNINSERTED: ($endMark5) $piece ($endMark3)\n";
	next;
    }

    my @ties5 = ();
    if (exists($endTies{$end5})) {
	my @tieInfo = @{$endTies{$end5}};
	foreach my $t (@tieInfo) {
	    my ($nLinks,$gapEstimate,$gapUncertainty,$tiedEnd) = $t =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;
	    push(@ties5,$tiedEnd);
	}   
    } else {
	print STDERR "UNINSERTED: (NO_5_END_TIES) $piece ()\n";
	next;
    }
    my %ties3 = ();
    if (exists($endTies{$end3})) {
	my @tieInfo = @{$endTies{$end3}};
	foreach my $t (@tieInfo) {
	    my ($nLinks,$gapEstimate,$gapUncertainty,$tiedEnd) = $t =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;
	    $ties3{$tiedEnd} = 1;
	}   
    } else {
	print STDERR "UNINSERTED: () $piece (NO_3_END_TIES)\n";
	next;
    }

    my %check = ();
    foreach my $t5 (@ties5) {
	if (exists($bestTie{$t5}) && exists($ties3{$bestTie{$t5}})) {
	    $check{$t5} = $bestTie{$t5};
	}
    }

    my @v5 = ();
    my @v3 = ();
    my $nValidSuspensions = 0;
    while (my ($t5,$t3) = each(%check)) {
	if (mutualUniqueBest($t5,$t3)) {
	    push(@v5,$t5);
	    push(@v3,$t3);
	    $nValidSuspensions++;
	}
    }

    if ($nValidSuspensions == 1) {

	my $tie5 = $v5[0];
	my $tie3 = $v3[0];

	delete($suspended{$piece});
	$bestTie{$tie5} = "$piece.5";
	$bestTie{"$piece.5"} = $tie5;
	$bestTiedBy{$tie5} = "$piece.5";
	$bestTiedBy{"$piece.5"} = $tie5;
	
	$bestTie{$tie3} = "$piece.3";
	$bestTie{"$piece.3"} = $tie3;
	$bestTiedBy{$tie3} = "$piece.3";
	$bestTiedBy{"$piece.3"} = $tie3;
	
	$nInsertedSuspensions++;
	print STDERR "SUSPENDED:  $tie5\t$piece\t$tie3\n";
    
    } else {
	print STDERR "UNINSERTED: (N_VALID_SUSP$nValidSuspensions) $piece ()\n";
    }
}

print STDERR "Done. Inserted $nInsertedSuspensions suspensions. ($nSuspendedPieces total suspended pieces)\n";

# Lock ends together with no competing ties

my %endLocks = ();

while (my ($end,$ties) = each (%bestTiedBy)) {

    if (exists($endLocks{$end})) {
	next;
    } 

    my ($piece) = $end =~ /^(.+)\.[35]$/;
    if (exists($suspended{$piece})) {
	next;
    }

    my @ties = split(/\,/,$ties);
    my $nTies = scalar(@ties);

    my $bestTie = "NONE";
    if (exists($bestTie{$end})) {
	$bestTie = $bestTie{$end};
    } else {
	next;
    }

    if (($nTies == 1) && ($ties[0] eq $bestTie)) {
	my $tieInfo = getTieInfo($end,$bestTie);
	my ($nLinks,$gapEstimate,$gapUncertainty,$nextEnd) = $tieInfo =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;
	my $gap = sprintf("%.0f",$gapEstimate);
	$endLocks{$end} = "$bestTie><$gap:$gapUncertainty";
	$endLocks{$bestTie} = "$end><$gap:$gapUncertainty";

    } else {

	$endLocks{$end} = "DIVERGENCE";

    }
}

foreach my $end (keys(%endLocks)) {
    print STDERR "LOCK: $end $endLocks{$end}\n";
}

# Traverse end locks to build scaffolds

my $scaffoldId = 1;
my %visitedContigs = ();

while (my ($contig,$contigLen) = each (%onoObjectLengths)) {

    my $startContig = $contig;
    my $startEnd = 5;

    if (exists($visitedContigs{$contig})) {
	next;
    }

    my $preState = "TERMINATION";
    my %loopCheck = ();
    while () {
	if (exists($loopCheck{$startContig})) {
	    last;
	}
	$loopCheck{$startContig} = 1;
	
	my $end = "$startContig.$startEnd";
	if (exists($endLocks{$end})) {
	    my $next = $endLocks{$end};
	    if (($next eq "CONVERGENCE") || ($next eq "DIVERGENCE")) {
		$preState = $next;
		last;
	    } else {
		($startContig,$startEnd) = $next =~ /(.+)\.([35])\>\</;
		$startEnd = ($startEnd == 3) ? 5 : 3;
		$preState = "EXTENSION";
	    }
	} else {
	    $preState = "TERMINATION";
	    last;
	}
    }
    
    %loopCheck = ();
    
    print STDERR "Scaffold$scaffoldId: $preState $startContig ($startEnd) ";

    my $inScaffold = 0;
    my $nextContig = $startContig;
    my $nextEnd = ($startEnd == 3) ? 5 : 3;
    my $nextGap = 0;
    my $nextGapUncertainty = 0;
    my $prevContig = $nextContig;
    my $prevEnd = $nextEnd;
    
    my $scaffReport = "";
    my $scaffCoord = 1;
    my $scaffContigIndex = 1;

    while () {
	if (exists($loopCheck{$nextContig})) {
	    print STDERR "$nextContig LOOP\n";
	    last;
	}
	$loopCheck{$nextContig} = 1;	    
	$visitedContigs{$nextContig} = 1;

	my $nextContigLen = $onoObjectLengths{$nextContig};
	
	my $contigOri = "+";
	if ($nextEnd == 5) {
	    $contigOri = "-";
	}
	
	if ($inScaffold) {
	    $scaffReport .= "Scaffold$scaffoldId\tGAP$scaffContigIndex\t$nextGap\t$nextGapUncertainty\n";
	    $scaffCoord += $nextGap;
	    $scaffContigIndex++;
	}

	if (defined($srfFile)) {
	    my $oldReport = $scaffReport{$nextContig};
	    my @reportLines = split(/\n/,$oldReport);
	    if ($contigOri eq "-") {
		@reportLines = reverse(@reportLines);
	    }
	    foreach my $line (@reportLines) {
		my @cols = split(/\t/,$line);
		if ($cols[1] =~ /^CONTIG/) {
		    my ($originalContig,$p0,$p1) = @cols[2,3,4];
		    my $originalContigLength = $p1-$p0+1;
		    my ($originalContigOri,$originalContigName) = $originalContig =~ /^([\+\-])(.+)$/;
		    if ($contigOri eq $originalContigOri) {
			$originalContigOri = "+";
		    } else {
			$originalContigOri = "-";
		    }

		    my $contigStartCoord = $scaffCoord;
		    my $contigEndCoord = $scaffCoord + $originalContigLength - 1;
		    if ($depthInfoAvailable) {
			my $depth = $contigDepths{$originalContigName};
			$scaffReport .= "Scaffold$scaffoldId\tCONTIG$scaffContigIndex\t$originalContigOri$originalContigName\t$contigStartCoord\t$contigEndCoord\t$depth\n";
		    } else {
			$scaffReport .= "Scaffold$scaffoldId\tCONTIG$scaffContigIndex\t$originalContigOri$originalContigName\t$contigStartCoord\t$contigEndCoord\n";
		    }
		    $scaffCoord += $originalContigLength;

		} else {

		    my ($originalGapSize,$originalGapUncertainty) = @cols[2,3];
		    $scaffReport .= "Scaffold$scaffoldId\tGAP$scaffContigIndex\t$originalGapSize\t$originalGapUncertainty\n";
		    $scaffCoord += $originalGapSize;
		    $scaffContigIndex++;
		}
	    }

	} else {
	    my $contigStartCoord = $scaffCoord;
	    my $contigEndCoord = $scaffCoord + $nextContigLen - 1;

	    if ($depthInfoAvailable) {
		my $depth = $contigDepths{$nextContig};
		$scaffReport .= "Scaffold$scaffoldId\tCONTIG$scaffContigIndex\t$contigOri$nextContig\t$contigStartCoord\t$contigEndCoord\t$depth\n";
	    } else {
		$scaffReport .= "Scaffold$scaffoldId\tCONTIG$scaffContigIndex\t$contigOri$nextContig\t$contigStartCoord\t$contigEndCoord\n";
	    }
	    $scaffCoord += $nextContigLen;
	}

	$inScaffold = 1;
	
	my $end = "$nextContig.$nextEnd";
	if (exists($endLocks{$end})) {
	    my $next = $endLocks{$end};
	    if (($next eq "CONVERGENCE") || ($next eq "DIVERGENCE")) {
		print STDERR "$nextContig ($nextEnd) $next\n";
		last;
	    } else {
		print STDERR "$nextContig ($nextEnd) ";
		($prevContig,$prevEnd) = ($nextContig,$nextEnd);
		($nextContig,$nextEnd,$nextGap,$nextGapUncertainty) = $next =~ /(.+)\.([35])\>\<(.+)\:(.+)/;
		print STDERR "[$nextGap +/- $nextGapUncertainty] $nextContig ($nextEnd) ";
		$nextEnd = ($nextEnd == 3) ? 5 : 3;
	    }
	} else {
	    print STDERR "$nextContig ($nextEnd) TERMINATION\n";
	    last;
	}
    }
    
    print $scaffReport;
    $scaffoldId++;
    
}

# ---------------
# | SUBROUTINES |
# ---------------

sub mutualUniqueBest {
    my ($end1,$end2) = @_;

    my ($piece1) = $end1 =~ /^(.+)\.[35]$/;
    my ($piece2) = $end2 =~ /^(.+)\.[35]$/;
    if (exists($suspended{$piece1}) || exists($suspended{$piece2})) {
	return 0;
    }

    unless (exists($bestTie{$end1}) && exists($bestTie{$end2})) {
	return 0;
    }

    my $bestTie1 = $bestTie{$end1};
    my $bestTie2 = $bestTie{$end2};

    unless (($bestTie1 eq $end2) && ($bestTie2 eq $end1)) {
	return 0;
    }

    unless (exists($bestTiedBy{$end1}) && exists($bestTiedBy{$end2})) {
	return 0;
    }

    my @ties1 = split(/\,/,$bestTiedBy{$end1});
    my $nTies1 = scalar(@ties1);
    my @ties2 = split(/\,/,$bestTiedBy{$end2});
    my $nTies2 = scalar(@ties1);

    unless ( ($nTies1 == 1) && ($ties1[0] eq $end2) && ($nTies2 == 1) && ($ties2[0] eq $end1)) {
	return 0;
    }

    return 1;
}



sub bestTie {
    my ($end) = @_;
    
    my @ties = ();
    if (exists($endTies{$end})) {
	@ties = @{$endTies{$end}};
    } else {
	return 0;
    }

    my ($testContig) = $end =~ /^(.+)\.[35]$/; 
    my $testContigLen = $onoObjectLengths{$testContig};

    my $largeObject = 0;
    if ($testContigLen > $suspendable) {
	$largeObject = 1;
    }

    my $nTies = scalar(@ties);
    my @tiedEnds = ();
    my @gapEstimates = ();
    my @gapUncertainties = ();
    my @contigLens = ();
    
    for (my $i = 0; $i < $nTies; $i++) {
	my $tie_i = $ties[$i];
	my ($nLinks_i,$gapEstimate_i,$gapUncertainty_i,$tiedEnd_i) = $tie_i =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;

	my $endMark_i = $endMarks{$tiedEnd_i};
	# Skip ties to previously marked ends
	if ($endMark_i) { 
	    next;
	}

	my ($contig_i) = $tiedEnd_i =~ /^(.+)\.[35]$/; 
	my $contigLen_i = $onoObjectLengths{$contig_i};

	push(@tiedEnds,$tiedEnd_i);
	push(@gapEstimates,$gapEstimate_i);
	push(@gapUncertainties,$gapUncertainty_i);
	push(@contigLens,$contigLen_i);
    }

    my $nGoodTies = scalar(@tiedEnds);

    unless($nGoodTies > 0) {
	print STDERR "BT: NoGoodTies\n";
	return 0;
    }

    my @sortedTieIndex = sort {($gapEstimates[$a] <=> $gapEstimates[$b]) || ($contigLens[$a] <=> $contigLens[$b])} (0..$#tiedEnds);

    # If a large object, return the closest large object (if one exists)
    # small objects may be suspended between large objects

    if ($largeObject) {

	my $closestLarge = undef;
	my @suspended = ();
	for (my $i=0;$i<$nGoodTies;$i++) {
	    my $ti = $sortedTieIndex[$i];
	    my $tiedEnd_i = $tiedEnds[$ti];
	    my $contigLen_i = $contigLens[$ti];
	    if ($contigLen_i > $suspendable) {
		$closestLarge = $tiedEnd_i;
		last;
	    } else {
		my ($sContig) = $tiedEnd_i =~ /^(.+)\.[35]$/; 
		push(@suspended,$sContig);
	    }
	}
	
	if (defined($closestLarge)) {
	    my $nSus = scalar(@suspended);
	    foreach my $sContig (@suspended) {
		$suspended{$sContig} =  1;
	    }
	    
	    print STDERR "BT: closestLarge ($nSus)\n";
	    return $closestLarge;
	    
	}
    }

    # Return the closest extendable object (if one exists)
    # unextendable objects may be suspended

    my $closestExtendable = undef;
    my @suspended = ();
    for (my $i=0;$i<$nGoodTies;$i++) {
	my $ti = $sortedTieIndex[$i];
	my $tiedEnd_i = $tiedEnds[$ti];
	my ($contig_i,$end_i) = $tiedEnd_i =~ /^(.+)\.([35])$/; 
	my $otherEnd_i = ($end_i eq "3") ? "5" : "3";
	my $testEnd_i = "$contig_i.$otherEnd_i";
	my $endMark_i = $endMarks{$testEnd_i};
	if ($endMark_i == 0) { 
	    $closestExtendable = $tiedEnd_i;
	    last;
	} else {
	    my ($sContig) = $tiedEnd_i =~ /^(.+)\.[35]$/; 
	    push(@suspended,$sContig);
	}
    }
    
    if (defined($closestExtendable)) {
	my $nSus = scalar(@suspended);
	foreach my $sContig (@suspended) {
	    $suspended{$sContig} =  1;
	}
	
	print STDERR "BT: closestExtendable ($nSus)\n";
	return $closestExtendable;	
    }
    
    # Just return the closest object

    print STDERR "BT: closest\n";

    return $tiedEnds[$sortedTieIndex[0]];

}

sub getTieInfo {
    my ($end1,$end2) = @_;
    
    my @ties = ();
    if (exists($endTies{$end1})) {
	@ties = @{$endTies{$end1}};
    } else {
	return 0;
    }

    my $nTies = scalar(@ties);

    for (my $i = 0; $i < $nTies; $i++) {
	my $tie_i = $ties[$i];
	my ($nLinks_i,$gapEstimate_i,$gapUncertainty_i,$tiedEnd_i) = $tie_i =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;

	if ($tiedEnd_i eq $end2) {
	    return $tie_i;
	}
    }

    return 0;
}


sub isLinkedTo {
    my ($end1,$end2) = @_;
    
    my @ties = ();
    if (exists($endTies{$end1})) {
	@ties = @{$endTies{$end1}};
    } else {
	return 0;
    }

    my $nTies = scalar(@ties);

    for (my $i = 0; $i < $nTies; $i++) {
	my $tie_i = $ties[$i];
	my ($nLinks_i,$gapEstimate_i,$gapUncertainty_i,$tiedEnd_i) = $tie_i =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;

	if ($tiedEnd_i eq $end2) {
	    return 1;
	}
    }

    return 0;
}



sub markEnd {
    my ($end) = @_;
    my ($testContig) = $end =~ /^(.+)\.[35]$/; 
    
    my @ties = ();
    if (exists($endTies{$end})) {
	@ties = @{$endTies{$end}};
    } else {
	return "NO_TIES";
    }

    my $testContigDepth = undef;
    if ($depthInfoAvailable) {
	$testContigDepth = $contigDepths{$testContig};
    }
    my $nTies = scalar(@ties);
    my @tiedEnds = ();
    my @gapEstimates = ();
    my @gapUncertainties = ();
    my @contigLens = ();
    my $nCollisions = 0;
    my $maxStrain = 3;
    my $maxDepthDropoff = 0.8;

    for (my $i = 0; $i < $nTies; $i++) {
	my $tie_i = $ties[$i];
	my ($nLinks_i,$gapEstimate_i,$gapUncertainty_i,$tiedEnd_i) = $tie_i =~ /(\d+)\.(\-?\d+)\:(\d+)\.(.+)/;
	my ($contig_i) = $tiedEnd_i =~ /^(.+)\.[35]$/; 

	if ($contig_i eq $testContig) {
	    return "SELF_TIE";
	}

	if ($depthInfoAvailable) {
	    my $contigDepth_i = $contigDepths{$contig_i};

	    if ( (($testContigDepth - $contigDepth_i)/$peakDepth) > $maxDepthDropoff) {
		return "DEPTH_DROP";
	    }
	}

	my $contigLen_i = $onoObjectLengths{$contig_i};

	push(@tiedEnds,$tiedEnd_i);
	push(@gapEstimates,$gapEstimate_i);
	push(@gapUncertainties,$gapUncertainty_i);
	push(@contigLens,$contigLen_i);
    }

   
    my @sortedTieIndex = sort {($gapEstimates[$a] <=> $gapEstimates[$b]) || ($contigLens[$a] <=> $contigLens[$b])} (0..$#ties);
    for (my $i=0;$i<$nTies-1;$i++) {
	my $ti = $sortedTieIndex[$i];
	my $start_i = $gapEstimates[$ti];
	my $end_i = $start_i+$contigLens[$ti]-1;
	my $uncertainty_i = $gapUncertainties[$ti];

	my $tj = $sortedTieIndex[$i+1];
	my $start_j = $gapEstimates[$tj];
	my $end_j = $start_j+$contigLens[$tj]-1;
	my $uncertainty_j = $gapUncertainties[$tj];

	my $overlap = $end_i - $start_j + 1;
	if ($overlap > $merSize-2) {
	    my $excessOverlap = $overlap-($merSize-2);
	    my $strain = $excessOverlap / ($uncertainty_i+$uncertainty_j);
	    if ($strain > $maxStrain) {
		return "TIE_COLLISION";
	    }
	}
    }

    return "UNMARKED";

}



