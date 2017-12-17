#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# merCounterS.pl by Jarrod Chapman <jchapman@lbl.gov> Mon Jun  1 19:41:59 PDT 2009
# Copyright 2009 Jarrod Chapman. All rights reserved.
#

use Proc::ProcessTable;
use Getopt::Std;
my %opts = ();
my $validLine = getopts('i:m:P:M:Y:RC:', \%opts);
my @required = ("i","m");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print "Usage: ./merCounterS.pl <-m merSize> <-i inputDataFile> <<-P prefixes>> <<-R(evComp?)>> <<-M minReported>> <<-Y maxMemory (Mb) default is 8Gb>> <<-C countFile>>\n";
    exit;
}

my $merSize = $opts{"m"};
my $inputDataFile = $opts{"i"};
# inputDataFile contains one line per sequence file, tab-delimited of the form
# fileName\tF|S\tG|U\tqualOffset
# F|S => fastq or scarf format; G|U => gzipped or unzipped; qualOffset is typically 33 or 64

my @inputFileInfo = ();
open (I,$inputDataFile) || die "Couldn't open $inputDataFile.\n";
while (my $line = <I>) {
    chomp $line;
    my @cols = split(/\t/,$line);
    my ($fileName,$fileFormat,$zipStatus,$qualOffset) = @cols;
    unless (-e $fileName) {
	die "Input file $fileName not found.\n";
    }
    unless ($fileFormat =~ /^[FS]$/) {
	die "File foramt ($fileFormat) should be F (fastq) or S (scarf).\n";
    }
    unless ($zipStatus =~ /^[GU]$/) { 
	die "File compression status ($zipStatus) should be U (unzipped) or G (gzipped).\n";
    }
    if ($qualOffset =~ /\D/) {
	die "Quality offset ($qualOffset) should be integral (usually 33 or 64).\n";
    }
    push(@inputFileInfo,$line);
}
close I;

my $revComp = 0;
if (exists($opts{"R"})) {
    $revComp = 1;
}

my $maxMemory = 8000;
if (exists($opts{"Y"})) {
    $maxMemory = $opts{"Y"};
}

my $minReportedCounts = 0;
if (exists($opts{"M"})) {
    $minReportedCounts = $opts{"M"};
}

my @matchStrings = ("[ACGT]{$merSize}");
if (exists($opts{"P"})) {
    @matchStrings = ();
    my @prefixes = split(/\+/,$opts{"P"});
    
    foreach my $prefix (@prefixes) {
	my $matchString = "";
	my $matchLen = $merSize-length($prefix);
	if ($matchLen > 0) {
	    $matchString = $prefix . "[ACGT]{$matchLen}";
	} else {
	    die "Bad options: prefix ($prefix) must be less than mer-size ($merSize)\n";
	}
	push(@matchStrings,$matchString);
    }
}
my @compiledMatchStrings = map qr/$_/, @matchStrings;

my %mers = ();

if (exists($opts{"C"})) {
    my $countFile = $opts{"C"};

    print STDERR "Reading countFile $countFile...\n";
    open (C,$countFile) || die "Couldn't open $countFile.\n";
    my $keptCounts = 0;
    while (my $line = <C>) {
	chomp $line;
	my ($mer,$count) = split(/\t/,$line);

	foreach my $compiledMatchString (@compiledMatchStrings) {
	    if ($mer =~ /$compiledMatchString/) {
		$mers{$mer} = $count;
		$keptCounts++;
		last;
	    }
	}
    }
    close C;
    my $memUsed = memoryUsed();
    if ($memUsed > $maxMemory) {
	die "Memory used ($memUsed MB) exceeded maximum ($maxMemory MB)\n";
    }
    print STDERR "$keptCounts counts kept.  $memUsed MB in use.\n";
}

my $readCount = 0;
my $checkMemFreq = 500000;

my $readBufferSize = 50000;
my $readBuffer = "";
my $inBuffer = 0;

my $nSeqFiles = scalar(@inputFileInfo);
my $seqFileIndex = 0;
foreach my $inputFile (@inputFileInfo) {

    my @cols = split(/\t/,$inputFile);
    my ($seqFile,$fileFormat,$zipStatus,$qualOffset) = @cols;

    my $q2sym = chr(2+$qualOffset);
    
    my $fastQ = 0;
    if ($fileFormat eq "F") {
	$fastQ = 1;
    }
    
    my $gzipped = 0;
    if ($zipStatus eq "G") {
	$gzipped = 1;
    }

    if ($gzipped) {
	open (S,"gunzip -c $seqFile |") || die "Couldn't open pipe to $seqFile (Error code: ($?) ($!))\n";
    } else {
	open (S,$seqFile) || die "Couldn't open $seqFile (Error code: ($?) ($!))\n";
    }
    $seqFileIndex++;
    print STDERR "Reading $seqFile ($seqFileIndex/$nSeqFiles)\n";

    while (my $line = <S>) {
	chomp $line;
	$readCount++;
	
	if ($readCount%$checkMemFreq == 0) {
	    my $memUsed = memoryUsed();
	    if ($memUsed > $maxMemory) {
		die "Memory used ($memUsed MB) exceeded maximum ($maxMemory MB)\n";
	    }
	    print STDERR "$readCount reads counted.  $memUsed MB in use.\n";
	}

	my $nts = undef;
	my $qString = undef;

	# Parse as fastq
	if ($fastQ) {
	    $line = <S>;
	    chomp $line;
	    $nts = $line;
	    $line = <S>;
	    $line = <S>;
	    chomp $line;
	    $qString = $line;

	} else {  #Parse as scarf
	    
	    my @cols = split(/\:/,$line,7);  # Modified Jan. 15 2013 to accomodate qual sym ":"
	    $nts = $cols[5];
	    $qString = $cols[6];
	}

# Q2 quality trimming standardized Jan 15 2013
#	if ($trimQ2) {

	if ($qString =~ /$q2sym/) {
	    my ($trimmedQ) = $qString =~ /^([^$q2sym]*)$q2sym/;
	    my $trimLen = length($trimmedQ);
	    $nts = substr($nts,0,$trimLen);
	}

#	}

	$readBuffer .= "X$nts";
	$inBuffer++;

	if ($inBuffer >= $readBufferSize) {

	    foreach my $compiledMatchString (@compiledMatchStrings) {

		while ($readBuffer =~ /($compiledMatchString)/g) {
		    $mers{$1}++;
		    pos($readBuffer) -= ($merSize-1);
		}
	    
		if ($revComp == 1) {
		    my $rc = reverse($readBuffer);
		    $rc =~ tr/ACGT/TGCA/;
		    while ($rc =~ /($compiledMatchString)/g) {
			$mers{$1}++;
			pos($rc) -= ($merSize-1);
		    }
		}
	    }

	    $readBuffer = "";
	    $inBuffer = 0;

	}
    }

    if ($inBuffer > 0) {
	foreach my $compiledMatchString (@compiledMatchStrings) {

	    while ($readBuffer =~ /($compiledMatchString)/g) {
		$mers{$1}++;
		pos($readBuffer) -= ($merSize-1);
	    }
	    
	    if ($revComp == 1) {
		my $rc = reverse($readBuffer);
		$rc =~ tr/ACGT/TGCA/;
		while ($rc =~ /($compiledMatchString)/g) {
		    $mers{$1}++;
		    pos($rc) -= ($merSize-1);
		}
	    }
	}

	$readBuffer = "";
	$inBuffer = 0;
	
    }

    close S;
}

my $nReported = 0;
while (my ($mer,$count) = each(%mers)) {
    if ($count >= $minReportedCounts) {
	print "$mer\t$count\n";
	$nReported++;
    }
}
print STDERR "Done. Reported $nReported mers.\n";

sub memoryUsed {
  my $procTable = new Proc::ProcessTable;
  foreach my $proc ( @{$procTable->table} ) {
    unless ($proc->pid eq $$) {
	next;
    }
    return ($proc->size)/(1024*1024);
  }
}

