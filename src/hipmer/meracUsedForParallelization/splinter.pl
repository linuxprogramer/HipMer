#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# splinter.pl by Jarrod Chapman <jchapman@lbl.gov> Fri Jun 12 12:58:11 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
my %opts = ();
my $validLine = getopts('b:m:F:T:', \%opts);
my @required = ("b","m");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print STDERR "Usage: ./splinter.pl <-b blastMapGlob> <-m minMatch> <<-F fivePrimeWiggleRoom>> <<-T threePrimeWiggleRoom>>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\nsplinter.pl";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

my $blastMapFileGlob = $opts{"b"};
my @blastMapFiles = glob($blastMapFileGlob);

my $minMatch = $opts{"m"};

my $fivePrimeWiggleRoom = 5;
if (exists($opts{"F"})) {
    $fivePrimeWiggleRoom = $opts{"F"};
}
my $threePrimeWiggleRoom = 5;
if (exists($opts{"T"})) {
    $threePrimeWiggleRoom = $opts{"T"};
}

my $currentRead = undef;
my $currentReadAlignments = 0;
my $readInfo = "";

my %rejected = ();
my @rejectReasons = ("FORMAT","MINLEN","UNINFORMATIVE","TRUNCATED","SINGLETON","OVERLAP");
foreach my $r (@rejectReasons) {
    $rejected{$r} = 0;
}

my $totalSplintsFound = 0;

foreach my $blastMapFile (@blastMapFiles) {

    print STDERR "Reading blast map file: $blastMapFile...\n";
    open (B,$blastMapFile) || die "Couldn't open $blastMapFile\n";
    
    while (my $line = <B>) {
	chomp $line;
	
	my ($blastType,$query,$qStart,$qStop,$qLength,$subject,$sStart,$sStop,$sLength,
	    $strand,$score,$eValue,$identities,$alignLength) = split(/\t/,$line);
	
	if ($blastType ne "BLASTN") {
	    $rejected{"FORMAT"}++;
	    next;
	}
	
	if ($identities < $minMatch) {
	    $rejected{"MINLEN"}++;
	    next;
	}
	
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
	my $endDistance = $qLength;
	
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
	
	if ($readStatus =~ /TRUNC/) {
	    $rejected{"TRUNCATED"}++;
	    next;
	}

	unless ($readStatus =~ /GAP/) {
	    $rejected{"UNINFORMATIVE"}++;
	    next;
	}
	
	if (defined($currentRead)) {
	    if ($query eq $currentRead) {
		$readInfo .= "$readStatus\t$line\n";
		$currentReadAlignments++;
	    } else {
		if ($currentReadAlignments > 1) {
		    $totalSplintsFound += processRead($readInfo);
		} else {
		    $rejected{"SINGLETON"}++;
		}
		$currentRead = $query;
		$readInfo = "$readStatus\t$line\n";
		$currentReadAlignments = 1;
	    }
	} else {
	    $currentRead = $query;
	    $readInfo = "$readStatus\t$line\n";
	    $currentReadAlignments = 1;
	}
    }
    close B;
    print STDERR "Done.\n";

    if (defined($currentRead)) {
	if ($currentReadAlignments > 1) {
	    $totalSplintsFound += processRead($readInfo);
	} else {
	    $rejected{"SINGLETON"}++;
	}
    }
}

print STDERR "Unused alignments:\n";
foreach my $r (@rejectReasons) {
    my $n = $rejected{$r};
    if ($n) {
	print STDERR "\t$r\t$n\n";
    }
}
print STDERR "Total splints found = $totalSplintsFound\n";

$date = `date`;
chomp $date;
print STDERR "Done. $date\n";


# -------------
# |SUBROUTINES|
# -------------

sub processRead {
    my ($info) = @_;
    my @alignments = split(/\n/,$info);
    my $nAligns = scalar(@alignments);
    my $nSplintsFound = 0;

    my %alignInfo = ();
    for (my $a = 0; $a < $nAligns; $a++) {
	$alignInfo{$a} = [split(/\t/,$alignments[$a])];
    }

    # Alignments are sorted by starting coordinate of the alignment on the read,
    # ties are broken by ending coordinate of alignement on the read
    my @sortedAlignIndex = sort {
	($alignInfo{$a}->[3] <=> $alignInfo{$b}->[3]) || ($alignInfo{$a}->[4] <=> $alignInfo{$b}->[4])
    } keys(%alignInfo);

    my $printSplints = "";

    my @dataSlice = (0,2,3,4,5,6,7,8,9,10);
#   $readStatus,$query,$qStart,$qStop,$qLength,$subject,$sStart,$sStop,$sLength,$strand

    for (my $a = 0; $a < $nAligns-1; $a++) {
	my $thisOne = $sortedAlignIndex[$a];
	my $thisInfo = $alignInfo{$thisOne};
	my @thisInfo = @{$thisInfo}[@dataSlice];
	my ($thisStatus,$thisContig,$thisStrand) = @thisInfo[0,5,9];
	my ($thisStartStat,$thisStopStat,$thisDirStat) = $thisStatus =~ /^(.+)\.(.+)\.(.+)$/;

	my $nextOne = $sortedAlignIndex[$a+1];
	my $nextInfo = $alignInfo{$nextOne};
	my @nextInfo = @{$nextInfo}[@dataSlice];
	my ($nextStatus,$nextContig,$nextStrand) = @nextInfo[0,5,9];
	my ($nextStartStat,$nextStopStat,$nextDirStat) = $nextStatus =~ /^(.+)\.(.+)\.(.+)$/;

	my $thisStop = $thisInfo[3];
	my $nextStart = $nextInfo[2];
	my $overlap = $thisStop-$nextStart+1;
	if ($overlap >= $minMatch) {
#	    warn "@thisInfo\n@nextInfo\n$overlap\n";
	    $rejected{"OVERLAP"} += $nAligns;
	    return 0;
	}

	unless ($thisContig ne $nextContig) {
	    next;
	}

	if (($thisStopStat eq "GAP") && ($nextStartStat eq "GAP")) {
	    my $thisContigEnd = ($thisStrand eq "Plus") ? "3" : "5";
	    my $nextContigEnd = ($nextStrand eq "Plus") ? "5" : "3";
	    $printSplints .= "SINGLE\tSPLINT\t$thisContig.$thisContigEnd\t[@thisInfo]\t$nextContig.$nextContigEnd\t[@nextInfo]\n";
	    $nSplintsFound++;
	}

#	print "[$thisStatus $thisContig $thisStrand] -> [$nextStatus $nextContig $nextStrand]\n";
#	print "[@thisInfo] -> [@nextInfo]\n";

    }

    if ($printSplints) {
	print $printSplints;
    }

    return $nSplintsFound;
}

