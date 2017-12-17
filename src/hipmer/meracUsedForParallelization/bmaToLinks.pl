#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# bmaToLinks.pl by Jarrod Chapman <jchapman@lbl.gov> Tue May 26 09:02:36 PDT 2009
#

use Getopt::Std;
my %opts = ();
my $validLine = getopts('b:m:I:D:P:', \%opts);
my @required = ("b","m",);
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print "Usage: ./bmaToLinks.pl <-m merSize> <-b bmaDataFile> <<-I innieList>> <<-D minNetEndDistance>> <<-P previousMaxInsertSize>>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\nbmaToLinks";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $bmaDataFile = $opts{"b"};
my $merSize = $opts{"m"};

my %innieLibs = ();
if (exists($opts{"I"})) {
    my $innies = $opts{"I"};
    @innieList = split(/\,/,$innies);
    foreach my $lib (@innieList) {
	$innieLibs{$lib} = 1;
    }
}

my $minNetEndDistance = 0;
if (exists($opts{"D"})) {
    $minNetEndDistance = $opts{"D"};
}

my $previousMaxInsertSize = 0;
if (exists($opts{"P"})) {
    $previousMaxInsertSize = $opts{"P"};
}

my $pi = 3.14159265359;
my $sqrtPi = sqrt($pi);
my $sqrtTwo = sqrt(2);

# Input is reads blastMapped to contigs
# Output is reduced contig linkage information

my $nMaps = 0;
my %mapInfo = ();
my @libs = ();
open (B,$bmaDataFile) || die "Couldn't open $bmaDataFile\n";
while (my $line = <B>) {
    chomp $line;
    my ($lib,$insertSize,$stdDev,$bmaFileGlob) = split(/\t/,$line);
    $mapInfo{$lib} = [$insertSize,$stdDev,$bmaFileGlob,undef];
    push(@libs,$lib);
    $nMaps++;
}
close B;

my %links = ();
my %redundancyCheck = ();
my %libRedundancy = ();
my %libPairSummary = ();
my %pairCategories = ();
my %objectLengths = ();
$pairCategories{"TRUNC"} = 0;
$pairCategories{"SMALL"} = 1;
$pairCategories{"INORM"} = 2;
$pairCategories{"ISHRT"} = 3;
$pairCategories{"INNIE"} = 4;
$pairCategories{"SELFL"} = 5;
$pairCategories{"EDIST"} = 6;
$pairCategories{"REDUN"} = 7;
$pairCategories{"ACCPT"} = 8;


# Read in BlastMapAnalyzer files

foreach my $lib (@libs) {

    my $insertSize = $mapInfo{$lib}->[0];
    my $insertStdDev = $mapInfo{$lib}->[1];
    my $bmaFileGlob = $mapInfo{$lib}->[2];
    $libRedundancy{$lib} = 0;
    $libPairSummary{$lib} = [0,0,0,0,0,0,0,0,0];

    print STDERR "Library $lib ($insertSize +/- $insertStdDev)\n";

    my @bmaFiles = glob($bmaFileGlob);

    foreach my $bmaFile (@bmaFiles) {

	print STDERR "Reading blast map analysis file $bmaFile.\n";

	open (B,$bmaFile) || die "Couldn't open $bmaFile\n";

	while (my $line = <B>) {
	    chomp $line;
	    my @cols = split(/\t/,$line);
	    
	    my ($pairStatus,$linkType) = ($cols[0],$cols[1]);
	    
	    my $nEnds = ( scalar(@cols) - 2 ) / 2;
	    my @ends = ();
	    my @aligns = ();
	    for (my $e = 0; $e < $nEnds; $e++) {
		push(@ends,$cols[2*$e+2]);
		push(@aligns,$cols[2*$e+3]);
	    }
	    
	    my ($end1,$end2) = @ends;
	    my ($align1,$align2) = @aligns;
	    if ($pairStatus eq "SINGLE" && $linkType eq "SPLINT") {
		processSplint($end1,$align1,$end2,$align2,$lib);
	    } elsif ($pairStatus eq "PAIR") {
		processPair($end1,$align1,$end2,$align2,$lib);
	    }
	}
	close B;
    }
}
print STDERR "Done reading blast map analysis files.\n";

#foreach my $lib (@libs) {
#    print STDERR "Discounted $libRedundancy{$lib} redundant spanning links for lib $lib.\n";
#}

foreach my $link (keys(%links)) {
    my @data = @{$links{$link}};
    my $nLinks = scalar(@data);

    my $nSplints = 0;
    my $nSpans = 0;
    my %splints = ();
    my %spans = ();
    foreach my $lib (@libs) {
	$spans{$lib} = [0,0,0];
    }

    my $minimumGap = -($merSize-2);
    my $splintMaxDev = 2;
    my $spanMaxZ = 3;

    my ($object1,$object2) = $link =~ /^(.+)\.[35]<=>(.+)\.[35]$/;

    foreach my $datum (@data) {

	# For SPLINTs, endSep is a direct measurement of gap (or overlap) size
	# For SPANs, endSep uses insertSize to estimate gap (or overlap) size 

	if ($datum =~ /^SPLINT/) {
	    my ($lib,$endSep) = $datum =~ /^SPLINT\.([A-Z]+)\.(\-?\d+)$/;

	    # Discard anomalously negative splints
	    unless ($endSep < $minimumGap-$splintMaxDev) {
		$nSplints++;
		$splints{$endSep}++;
	    }
	} else {
	    my ($lib,$endSep) = $datum =~ /^([A-Z]+)\.(\-?\d+)$/;

	    # Discard anomalously negative spans
	    my $spanAnomaly = 0;
	    if ($endSep < $minimumGap) {
		my $insertStdDev = $mapInfo{$lib}->[1];
		my $spanZ = abs($minimumGap-$endSep)/$insertStdDev;
		if ($spanZ > $spanMaxZ) {
		    $spanAnomaly = 1;
		}
	    }

	    unless ($spanAnomaly) {
		$nSpans++;
		$spans{$lib}->[0] += 1;
		$spans{$lib}->[1] += $endSep;
	    }
	}
    }

    # Gap size as estimated from splints is taken to be the maximum frequency splint
    my $splintGapEstimate = undef;
    my $splintMaxFreq = 0;
    while (my ($splint,$count) = each(%splints)) {
	unless (defined($splintGapEstimate)) {
	    $splintGapEstimate = $splint;
	    $splintMaxFreq = $count;
	}
	if ($count > $splintMaxFreq) {
	    $splintGapEstimate = $splint;
	    $splintMaxFreq = $count;
	}	    
    }

    # Gap size as estimated from spans is taken to be the weighted mean of spans
    foreach my $lib (@libs) {
	my $nLibSpans = $spans{$lib}->[0]; 
	unless ($nLibSpans) {
	    next;
	}
	my $insertSize = $mapInfo{$lib}->[0];
	my $insertStdDev = $mapInfo{$lib}->[1];
	my $readLength = $mapInfo{$lib}->[3];
	my $netSeparation = $spans{$lib}->[1];
	my $meanGapEstimate = ($netSeparation/$nLibSpans);
	my $meanOffset = $insertSize - $meanGapEstimate;

	unless (exists($objectLengths{$object1}) && exists($objectLengths{$object2})) {
	    die "Error: length of $object1 or $object2 not found.\n";
	}

	my $l1 = $objectLengths{$object1};
	my $l2 = $objectLengths{$object2};
	my ($objLen1,$objLen2) = ($l1 < $l2) ? ($l1, $l2) : ($l2,$l1);

	my $gapEstimate = estimateGapSize ($meanOffset,$merSize,$readLength,$objLen1,$objLen2,$insertSize,$insertStdDev);

	$spans{$lib}->[2] = $gapEstimate;
    }

    my $spanGapEstimate = 0;
    my $sumOfWeights = 0;
    my $spanGapUncertainty = undef;
   
    foreach my $lib (@libs) {
	my $insertStdDev = $mapInfo{$lib}->[1];
	my $libWeight = ($spans{$lib}->[0])/($insertStdDev*$insertStdDev);
	$sumOfWeights += $libWeight;
	$spanGapEstimate += ($spans{$lib}->[2])*$libWeight;
    }
    if ($sumOfWeights) {
	$spanGapEstimate /= $sumOfWeights;
	$spanGapUncertainty = sqrt(1/$sumOfWeights);
    }

    my $anomalousSplints = 0;
    my $anomalousSpans = 0;
    
    # Revisit data checking for consistency with estimate
    foreach my $datum (@data) {

	if ($datum =~ /^SPLINT/) {
	    unless (defined($splintGapEstimate)) {
		next;
	    }

	    my ($lib,$endSep) = $datum =~ /^SPLINT\.([A-Z]+)\.(\-?\d+)$/;

	    my $splintDev = abs($endSep-$splintGapEstimate);
	    if ($splintDev > $splintMaxDev) {
		$anomalousSplints++;
	    }

	} else {
	    unless(defined($spanGapUncertainty)) {
		next;
	    }

	    my ($lib,$endSep) = $datum =~ /^([A-Z]+)\.(\-?\d+)$/;
	    my $insertStdDev = $mapInfo{$lib}->[1];
	    my $spanZ = abs($endSep-$spanGapEstimate)/$insertStdDev;
	    if ($spanZ > $spanMaxZ) {
		$anomalousSpans++;
	    }
	}
    }

    if (defined($splintGapEstimate)) {
	print "SPLINT\t$link\t$splintMaxFreq|$anomalousSplints|$nSplints\t$splintGapEstimate\n";
    }

    if (defined($spanGapUncertainty)) {
	$gapEstimate = sprintf("%.0f",$spanGapEstimate);
	$gapUncertainty = sprintf("%.0f",$spanGapUncertainty);
	print "SPAN\t$link\t$anomalousSpans|$nSpans\t$gapEstimate\t$gapUncertainty\n";
    }

}

# Report some summary statistics for each library

my @pairCats = ("TRUNC","SMALL","INORM","ISHRT","INNIE","SELFL","EDIST","REDUN","ACCPT");
print STDERR "LIB";
foreach my $cat (@pairCats) {
    print STDERR "\t$cat";
}
print STDERR "\tTOTAL\n";
foreach my $lib (@libs) {
    print STDERR "$lib";
    my $total = 0;
    foreach my $cat (@pairCats) {
	my $catCount = $libPairSummary{$lib}->[$pairCategories{$cat}];
	$total += $catCount;
    } 
    foreach my $cat (@pairCats) {
	my $catCount = $libPairSummary{$lib}->[$pairCategories{$cat}];
	my $percent = $total ? sprintf("%.3f",$catCount/$total) : 0;
	print STDERR "\t$percent";
    }
    print STDERR "\t$total\n";
}

$date = `date`;
chomp $date;
print STDERR "Done. $date\n";


# ---------------
# | SUBROUTINES |
# ---------------

sub processPair {
    my ($end1,$a1,$end2,$a2,$lib) = @_;

    my $insertSize = $mapInfo{$lib}->[0];
    my $insertStdDev = $mapInfo{$lib}->[1];
    my $endDistance = $insertSize + 3*$insertStdDev;

    my $innieMaxSep = 1000;
    my $innieLib = 0;
    if (exists($innieLibs{$lib})) {
	$innieLib = 1;
    }

    my ($contig1) = $end1 =~ /^(.+)\.[35]$/;
    my ($contig2) = $end2 =~ /^(.+)\.[35]$/;

    my ($align) = $a1 =~ /^\[(.+)\]$/;
    my ($rStat,$r,$r0,$r1,$rL,$c,$c0,$c1,$cL,$strand) = split(/\s+/,$align);
    my ($strand1,$sStart1,$sStop1) = ($strand,$c0,$c1);

# Assumes all reads same length from a given library
    unless (defined($mapInfo{$lib}->[3])) {
	$mapInfo{$lib}->[3] = $rL;
    }
    unless (exists($objectLengths{$c})) {
	$objectLengths{$c} = $cL;
    }

    # Omit truncated alignments as potential sources of error
    if ($rStat =~ /TRUNC/) {
	$libPairSummary{$lib}->[$pairCategories{"TRUNC"}] += 1;
	return;
    }

    # Omit alignments to very short (potentially "popped-out") scaffolds
    if ($cL < $previousMaxInsertSize/2) {
	$libPairSummary{$lib}->[$pairCategories{"SMALL"}] += 1;
	return;
    }

    my $d1 = undef;
    my $leftend1 = $c0;
    my $rightend1 = $c1;
    if ($strand eq "Plus") {
	$d1 = ($cL-$c0+1)+($r0-1);
	$leftend1 -= ($r0-1);
	$rightend1 += ($rL-$r1);
    } else {
	$d1 = $c1+($r0-1);
	$leftend1 -= ($rL-$r1);
	$rightend1 += ($r0-1);
    }
    ($align) = $a2 =~ /^\[(.+)\]$/;
    ($rStat,$r,$r0,$r1,$rL,$c,$c0,$c1,$cL,$strand) = split(/\s+/,$align);
    my ($strand2,$sStart2,$sStop2) = ($strand,$c0,$c1);

    unless (exists($objectLengths{$c})) {
	$objectLengths{$c} = $cL;
    }

    if ($rStat =~ /TRUNC/) {
	$libPairSummary{$lib}->[$pairCategories{"TRUNC"}] += 1;
	return;
    }

    # Omit alignments to very short (potentially "popped-out") scaffolds
    if ($cL < $previousMaxInsertSize/2) {
	$libPairSummary{$lib}->[$pairCategories{"SMALL"}] += 1;
	return;
    }

    my $d2 = undef;
    my $leftend2 = $c0;
    my $rightend2 = $c1;
    if ($strand eq "Plus") {
	$d2 = ($cL-$c0+1)+($r0-1);
	$leftend2 -= ($r0-1);
	$rightend2 += ($rL-$r1);
    } else {
	$d2 = $c1+($r0-1);
	$leftend2 -= ($rL-$r1);
	$rightend2 += ($r0-1);
    }
    
    # check for inappropriate self-linkage

    if ($contig1 eq $contig2) {
	
	my $orientation = "undefined";
	my $separation = 0;
	my $sepZ = 0;
	my $maxSepZ = 5;

	if ($strand1 ne $strand2) {
#	    if ( (($strand1 eq "Plus") && ($sStart1 < $sStart2) && ($sStart1 < $sStop2)) || 
#		 (($strand2 eq "Plus") && ($sStart2 < $sStart1) && ($sStart2 < $sStop1)) ) {
	    if ( (($strand1 eq "Plus") && ($sStart1 <= $sStart2) && ($sStop1 <= $sStop2)) || 
		 (($strand2 eq "Plus") && ($sStart2 <= $sStart1) && ($sStop2 <= $sStop1)) ) {
		$orientation = "convergent";
		$separation = ($strand1 eq "Plus") ? ($rightend2 - $leftend1 + 1) : ($rightend1 - $leftend2 + 1);
		$sepZ = abs(($separation-$insertSize)/$insertStdDev);
	    } elsif ( (($strand1 eq "Plus") && ($sStart1 > $sStart2)) || 
		      (($strand2 eq "Plus") && ($sStart2 > $sStart1)) ) {
		$orientation = "divergent";
		$separation = ($strand1 eq "Plus") ? ($rightend1 - $leftend2 + 1) : ($rightend2 - $leftend1 + 1);
	    }
	}

	if (($orientation eq "convergent") && ($sepZ < $maxSepZ)) {
	    # A normal-looking internal pair
	    $libPairSummary{$lib}->[$pairCategories{"INORM"}] += 1;
	    return;

	} elsif (($orientation eq "convergent") && ($separation < $minNetEndDistance)) {
	    # A "shorty"
	    $libPairSummary{$lib}->[$pairCategories{"ISHRT"}] += 1;
	    return;

	} else {

	    if ($innieLib && ($orientation eq "divergent") && ($separation < $innieMaxSep)) {
		# An "innie"
		$libPairSummary{$lib}->[$pairCategories{"INNIE"}] += 1;
		return;
	    }

	    # Let abnormal self-links go through to try to prevent O&O of repeats
	    $libPairSummary{$lib}->[$pairCategories{"SELFL"}] += 1;
	    print STDERR "SELF-LINK\t$end1\t$strand1\t$sStart1\t$sStop1\t$end2\t$strand2\t$sStart2\t$sStop2\t$orientation\t$sepZ\n";
	}
    }


    unless (defined($d1) && defined($d2)) {
	print STDERR "processPair: Error comprehending: [$a1] [$a2]\n";
	return;
    } else {
	
#	unless ($d1 < $endDistance && $d2 < $endDistance) {
	unless (($d1 < $endDistance) && ($d2 < $endDistance) && ($d1+$d2 > $minNetEndDistance)) {
	    $libPairSummary{$lib}->[$pairCategories{"EDIST"}] += 1;
	    return;
	}

	my $endSeparation = $insertSize - ($d1+$d2);

	my $link = ($end1 lt $end2) ? "$end1<=>$end2" : "$end2<=>$end1";
        my $span = ($end1 lt $end2) ? "$d1.$d2" : "$d2.$d1";

        if (exists($redundancyCheck{"$link.$span"})) {
            $redundancyCheck{"$link.$span"}++;
	    $libRedundancy{$lib}++;
	    $libPairSummary{$lib}->[$pairCategories{"REDUN"}] += 1;
            return;
        } else {
            $redundancyCheck{"$link.$span"} = 1;
        }

	$libPairSummary{$lib}->[$pairCategories{"ACCPT"}] += 1;
	if (exists($links{$link})) {
	    push(@{$links{$link}},"$lib.$endSeparation");
	} else {
	    $links{$link} = ["$lib.$endSeparation"];
	}	

    }
}


sub processSplint {
    my ($end1,$a1,$end2,$a2,$lib) = @_;

    my ($contig1) = $end1 =~ /^(.+)\.[35]$/;
    my ($contig2) = $end2 =~ /^(.+)\.[35]$/;

    my ($align) = $a1 =~ /^\[(.+)\]$/;
    my ($rStat,$r,$r0,$r1,$rL,$c,$c0,$c1,$cL,$strand) = split(/\s+/,$align);

# Assumes all reads same length from a given library
    unless (defined($mapInfo{$lib}->[3])) {
	$mapInfo{$lib}->[3] = $rL;
    }

    if ($rStat =~ /TRUNC/) {
	return "";
    }

    my $rightCoord = undef;
    my $leftCoord = undef;

    my ($leftStat,$rightStat,$locStat) = split(/\./,$rStat);
    unless ($rightStat eq "GAP") {
	warn "warning: processSplint: unexpected rStat:$rStat for SPLINT alignment1 $a1 Omitting.\n";
	return "";
    }

    if ($strand eq "Plus") {
	$leftCoord = $r1+($cL-$c1);
    } else {
	$leftCoord = $r1+($c0-1); 
    }
    
    ($align) = $a2 =~ /^\[(.+)\]$/;
    ($rStat,$r,$r0,$r1,$rL,$c,$c0,$c1,$cL,$strand) = split(/\s+/,$align);

    if ($rStat =~ /TRUNC/) {
	return "";
    }

    ($leftStat,$rightStat,$locStat) = split(/\./,$rStat);
    unless ($leftStat eq "GAP") {
	warn "warning: processSplint: unexpected rStat:$rStat for SPLINT alignment2 $a2 Omitting.\n";
	return "";
    }

    if ($strand eq "Plus") {
	$rightCoord = $r0-($c0-1);
    } else {
	$rightCoord = $r0-($cL-$c1); 
    }
    
    unless (defined($rightCoord) && defined($leftCoord)) {
	warn "warning: processSplint: Error comprehending alignments: $a1 $a2 Omitting\n";
	return "";
    } else {
	my $endSeparation = $rightCoord - $leftCoord - 1;
	my $link = ($end1 lt $end2) ? "$end1<=>$end2" : "$end2<=>$end1";

	# Note: no redundancy check for SPLINTs due to limited possible splinting configurations

	if (exists($links{$link})) {
	    push(@{$links{$link}},"SPLINT.$lib.$endSeparation");
	} else {
	    $links{$link} = ["SPLINT.$lib.$endSeparation"];
	}
    }

    return "$end1=>$end2";

}

# Subroutines for evaluating the gap correction analytically
# assuming Gaussian insert size distribution

sub estimateGapSize {
    my ($meanAnchor,$k,$l,$c1,$c2,$mu,$sigma) = @_;

    my $gMax = $mu+3*$sigma-2*$k;
    my $gMin = -($k-2);
    my $gMid = $mu-$meanAnchor;

    if ($gMid < $gMin) {
        $gMid = $gMin+1;
    } elsif ($gMid > $gMax) {
        $gMid = $gMax-1;
    }

    my $aMax = meanSpanningClone($gMax,$k,$l,$c1,$c2,$mu,$sigma) - $gMax;
    my $aMin = meanSpanningClone($gMin,$k,$l,$c1,$c2,$mu,$sigma) - $gMin;
    my $aMid = meanSpanningClone($gMid,$k,$l,$c1,$c2,$mu,$sigma) - $gMid;

    my $deltaG = $gMax-$gMin;

    my $iterations = 0;

    while ($deltaG > 10) {
        $iterations++;
        if ($meanAnchor > $aMid) {
            $gMax = $gMid;
            $aMax = $aMid;
            $gMid = ($gMid+$gMin)/2;
            $aMid = meanSpanningClone($gMid,$k,$l,$c1,$c2,$mu,$sigma) - $gMid;
        } elsif ($meanAnchor < $aMid) {
            $gMin = $gMid;
            $aMin = $aMid;
            $gMid = ($gMid+$gMax)/2;
            $aMid = meanSpanningClone($gMid,$k,$l,$c1,$c2,$mu,$sigma) - $gMid;
        } else {
            last;
        }
        $deltaG = $gMax-$gMin;
    }

#    print STDERR "#ITER $iterations\n";
    return $gMid;

}

sub meanSpanningClone {
    my ($g,$k,$l,$c1,$c2,$mu,$sigma) = @_;

    my $x1 = $g+2*$k-1;
    my $x2 = $g+$c1+$l;
    my $alpha = $x2-$x1;
    my $x3 = $g+$c2+$l;
    my $x4 = $x3+$alpha;

    my $num = 0;
    my $den = 0;

    my $N1 = G2($x1,$x2,$mu,$sigma)-$x1*G1($x1,$x2,$mu,$sigma);
    my $N2 = ($x2-$x1)*G1($x2,$x3,$mu,$sigma);
    my $N3 = $x4*G1($x3,$x4,$mu,$sigma)-G2($x3,$x4,$mu,$sigma);
    
    my $D1 = G1($x1,$x2,$mu,$sigma)-$x1*G0($x1,$x2,$mu,$sigma);
    my $D2 = ($x2-$x1)*G0($x2,$x3,$mu,$sigma);
    my $D3 = $x4*G0($x3,$x4,$mu,$sigma)-G1($x3,$x4,$mu,$sigma);
    
    $num = $N1+$N2+$N3;
    $den = $D1+$D2+$D3;

    if ($den) {
        return $num/$den;
    } else {
        print STDERR "Warning: meanSpanningClone failed for ($g,$k,$l,$c1,$c2,$mu,$sigma)\n";
        return 0;
    }
}



sub erf {
    my ($x) = @_;
    
    my $absX = ($x < 0) ? -$x : $x;
    my $t = 1/(1+0.5*$absX);

    my $t2 = $t*$t;
    my $t3 = $t*$t2;
    my $t4 = $t*$t3;
    my $t5 = $t*$t4;
    my $t6 = $t*$t5;
    my $t7 = $t*$t6;
    my $t8 = $t*$t7;
    my $t9 = $t*$t8;

    my $a0 = -1.26551223;
    my $a1 = 1.00002368;
    my $a2 = 0.37409196;
    my $a3 = 0.09678418;
    my $a4 = -0.18628806;
    my $a5 = 0.27886807;
    my $a6 = -1.13520398;
    my $a7 = 1.48851587;
    my $a8 = -0.82215223;
    my $a9 = 0.17087277;

    my $tau = $t*exp(-$x*$x + $a0 + $a1*$t + $a2*$t2 + $a3*$t3 + $a4*$t4 +
                     $a5*$t5 + $a6*$t6 + $a7*$t7 + $a8*$t8 + $a9*$t9);

    if ($x < 0) {
        return ($tau-1);
    } else {
        return (1-$tau);
    }
}

sub G0 {
    my ($a,$b,$mu,$sigma) = @_;

    my $rt2sig = $sqrtTwo*$sigma;
    my $erfa = erf(($a-$mu)/$rt2sig);
    my $erfb = erf(($b-$mu)/$rt2sig);

    return ($sqrtPi/$sqrtTwo)*$sigma*($erfb-$erfa);
}

sub G1 {
    my ($a,$b,$mu,$sigma) = @_;

    my $za = ($a-$mu)/$sigma;
    my $zb = ($b-$mu)/$sigma;
    
    my $expa = exp(-0.5*$za*$za);
    my $expb = exp(-0.5*$zb*$zb);

    my $g0 = G0($a,$b,$mu,$sigma);

    return ($sigma*$sigma*($expa-$expb) + $mu*$g0);
}

sub G2 {
    my ($a,$b,$mu,$sigma) = @_;

    my $za = ($a-$mu)/$sigma;
    my $zb = ($b-$mu)/$sigma;
    
    my $expa = exp(-0.5*$za*$za);
    my $expb = exp(-0.5*$zb*$zb);

    my $g0 = G0($a,$b,$mu,$sigma);
    my $g1 = G1($a,$b,$mu,$sigma);

    my $sigma2 = $sigma*$sigma;

    return ($sigma2*$g0 + $mu*$g1 + $sigma2*($a*$expa - $b*$expb));
}
