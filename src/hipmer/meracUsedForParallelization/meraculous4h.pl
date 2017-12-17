#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# meraculous.pl by Jarrod Chapman <jchapman@lbl.gov> Wed Jun 13 11:42:09 PDT 2007
# Copyright 2007 Jarrod Chapman. All rights reserved.
#

#use Devel::DProf;
use B;
use strict;
use Proc::ProcessTable;
use Getopt::Std;
use POSIX qw (ceil);

my %opts = ();
my $validLine = getopts('i:m:c:D:M:U:P:Y:', \%opts);
$validLine &= ( (exists($opts{"i"}) && exists($opts{"m"})) || exists($opts{"M"}) || exists($opts{"U"}));

my $usageString = "Usage: ./meraculous.pl [ <-U UFXHashListFile> || <-M merGraphFile> || (<-m merCountFile> <-i inputDataFile> <<-P prefixes>>) ] <<-c minContigSize>> <<-D significantDepth>> <<-Y maxMemory (Gb) default is 100Gb>>\n"; 

if ($validLine != 1) {
    print STDERR $usageString;
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\nmeraculous4";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";

#my $DEBUGMER = undef;
#my $DEBUGRC = undef;
#if (exists($opts{"d"})) {
#    $DEBUGMER = $opts{"d"};
#    $DEBUGRC = reverse($DEBUGMER);
#    $DEBUGRC =~ tr/ACGT/TGCA/;
#}

my $maxMemory = 100;
if (exists($opts{"Y"})) {
    $maxMemory = $opts{"Y"};
}

my $merCountFile = undef; 
if (exists($opts{"m"})) {
    $merCountFile = $opts{"m"};
}
my $merGraphFile = undef;
if (exists($opts{"M"})) {
    $merGraphFile = $opts{"M"};
}

my $inputDataFile = undef;
my @inputFileInfo = ();
if (exists($opts{"i"})) {
    $inputDataFile = $opts{"i"};

# inputDataFile contains one line per sequence file, tab-delimited of the form
# fileName\tF|S\tG|U\tqualOffset
# F|S => fastq or scarf format; G|U => gzipped or unzipped; qualOffset is typically 33 or 64

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
}

my $UFXHashListFile = undef;
if (exists($opts{"U"})) {
    $UFXHashListFile = $opts{"U"};
}

my $minContigSize = 500;
if (exists($opts{"c"})) {
    $minContigSize = $opts{"c"};
}

my $significantDepth = 2;
if (exists($opts{"D"})) {
    $significantDepth = $opts{"D"};
}

my %baseCodes = ("A",0,"C",1,"G",2,"T",3,"N",4,"X",5);
my @bases = ("A","C","G","T");
my %mers = ();
my $merSize = undef;

my @ufxHashArray = ();
my $ufxWordSize = undef;
my @ufxHashSizes = ();
my %hashBaseCodes = ("A",0b00,"C",0b01,"G",0b10,"T",0b11);
@hashBaseCodes{values %hashBaseCodes} = keys (%hashBaseCodes);

my %hashBaseCodes2 = ();
for (my $i = 0; $i < 4**4; $i++) { 
    my $mer = ""; 
    my $k = $i; 
    for (my $j = 0; $j < 4; $j++) { 
	$mer .= "$bases[$k%4]"; 
	$k = sprintf("%d",$k/4); 
    } 
    $hashBaseCodes2{$mer} = $i;
}
@hashBaseCodes2{values %hashBaseCodes2} = keys (%hashBaseCodes2);


my %hashExtCodes = ("A",0b0100,"C",0b0101,"G",0b0110,"T",0b0111,"F",0b0001,"X",0b0010);
@hashExtCodes{values %hashExtCodes} = keys (%hashExtCodes);
my $doNotVisitCode = 0b1000;
$hashExtCodes{$doNotVisitCode} = "";
my $wordsPerHashElement = undef;
my @hashPrefixSearchOrder = ();
my %hashFilePrefixToIndex = ();
my $merNibbles = undef;


# Input set of k-mers to use and read seq/quals
# Generate mergraph

if (defined($merCountFile) && defined($inputDataFile)) {

    my @prefixSet = ();
    if (exists($opts{"P"})) {
	@prefixSet = split(/\+/,$opts{"P"});
    }

    print STDERR "Reading $merCountFile...\n";
    my $nMers = 0;
    open (M,$merCountFile) || die "Couldn't open $merCountFile\n";
    while (my $line = <M>) {
	chomp $line;

	my $mer = undef;
	my $count = undef;

	if ($line =~ /^([acgtACGT]+)\s+(\d+)$/) {
	    ($mer,$count) = $line =~ /^([acgtACGT]+)\s+(\d+)$/;
	} elsif ($line =~ /^([acgtACGT]+)\s+(\d+)\s+(\d+)$/) {
	    my ($m,$count1,$count2) = $line =~ /^([acgtACGT]+)\s+(\d+)\s+(\d+)$/;
	    $mer = $m;
	    $count = $count1+$count2;
	} else {
	    die "Error:  line [$line] not a valid mercount format.\n";
	}

	$mer = uc($mer);
	
	if (exists($opts{"P"})) {

	    my $goodMer = 0;
	    foreach my $p (@prefixSet) {
		if ($mer =~ /^$p/) {
		    $goodMer = 1;
		    last;
		}
	    }
	    unless($goodMer) {
		next;
	    }

	}

	unless (defined($merSize)) {
	    $merSize = length($mer);
	}
	unless (exists($mers{$mer}) || ($count < $significantDepth)) {
	    $mers{$mer} = [0,0,0,0,0,0,0,0,0,0,0,0,0];
	    $nMers++;
	}
    }
    close M;
    print STDERR "Done. Found $nMers $merSize-mers with at least $significantDepth counts.\n";

    my $nReads = 0;
    my $totalReadLength = 0;
    my $nValidMers = 0;
    my $minQuality = 19;
    my $reportFreq = 1000000;

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

    my $readBufferSize = 50000;
    my $readBuffer = "X";
    my $qualBuffer = "X";
    my $inBuffer = 0;

    my $nSeqFiles = scalar(@inputFileInfo);
    my $seqFileIndex = 0;
    foreach my $inputFile (@inputFileInfo) {

	my ($sequenceFile,$fileFormat,$zipStatus,$qualOffset) = split(/\t/,$inputFile);

	my $fastQFormat = 0;
	if ($fileFormat eq "F") {
	    $fastQFormat = 1;
	}
    
	my $gzipped = 0;
	if ($zipStatus eq "G") {
	    $gzipped = 1;
	}

	if ($gzipped) {
	    open (S,"gunzip -c $sequenceFile |") || die "Couldn't open pipe to $sequenceFile\n";
	} else {
	    open (S,$sequenceFile) || die "Couldn't open $sequenceFile\n";
	}
#	print STDERR "Reading in $sequenceFile...\n";
	$seqFileIndex++;
	print STDERR "Reading in $sequenceFile ($seqFileIndex/$nSeqFiles)\n";

	while (my $line = <S>) {
	    chomp $line;

	    my $nts = undef;
	    my $quals = undef;

	    if ($fastQFormat) {
		$line = <S>;
		chomp $line;
		$nts = $line;
		$line = <S>;
		$line = <S>;
		chomp $line;
		$quals = $line;
	    } else {
#		my @cols = split(/\:/,$line);
		my @cols = split(/\:/,$line,7);  # Modified Feb. 24 2013 to accomodate qual sym ":"
		my $readName = "$cols[0].$cols[1].$cols[2].$cols[3].$cols[4]";
		$nts = $cols[5];
		$quals = $cols[6]; 
	    }

	    my $seqLen = length($nts);
	    $totalReadLength += $seqLen;
	    $nReads++;
	    if ($nReads%$reportFreq == 0) {
		print STDERR ".";
	    }

	    $readBuffer .= $nts."X";
	    $qualBuffer .= $quals."X";
	    $inBuffer++;

	    if ($inBuffer >= $readBufferSize) {

		foreach my $compiledMatchString (@compiledMatchStrings) {

		    while ($readBuffer =~ /(.)($compiledMatchString)(.)/g) {

			my ($pre,$core,$post) = ($1,$2,$3);
			pos($readBuffer) -= ($merSize+1);
			my $sCoord = pos($readBuffer)-1;
		
			if (exists($mers{$core})) {
		    
			    $nValidMers++;
		    
			    my $preCode = $baseCodes{$pre};
			    my $postCode = 6+$baseCodes{$post};
		    
			    my $qpre = substr($qualBuffer,$sCoord,1);
			    my $qpost = substr($qualBuffer,$sCoord+$merSize+1,1);
		    
			    my $Qpre = ord($qpre) - $qualOffset;
			    my $Qpost = ord($qpost) - $qualOffset;
		    
			    if ($Qpre > $minQuality) {
				$mers{$core}->[$preCode] += 1;
			    } 
			    if ($Qpost > $minQuality) {
				$mers{$core}->[$postCode] += 1;
			    }
			}
		    }
		    
		    $readBuffer = revComp($readBuffer);
		    $qualBuffer = reverse($qualBuffer);
	    
		    while ($readBuffer =~ /(.)($compiledMatchString)(.)/g) {

			my ($pre,$core,$post) = ($1,$2,$3);
			pos($readBuffer) -= ($merSize+1);
			my $sCoord = pos($readBuffer)-1;
		
			if (exists($mers{$core})) {
		    
			    my $preCode = $baseCodes{$pre};
			    my $postCode = 6+$baseCodes{$post};
		    
			    my $qpre = substr($qualBuffer,$sCoord,1);
			    my $qpost = substr($qualBuffer,$sCoord+$merSize+1,1);
		    
			    my $Qpre = ord($qpre) - $qualOffset;
			    my $Qpost = ord($qpost) - $qualOffset;
		    
			    if ($Qpre > $minQuality) {
				$mers{$core}->[$preCode] += 1;
			    }
			    if ($Qpost > $minQuality) {
				$mers{$core}->[$postCode] += 1;	    
			    }
			}
		    }
		}

		$readBuffer = "X";
		$qualBuffer = "X";
		$inBuffer = 0;

	    }
	}

	if ($inBuffer > 0) {
	    foreach my $compiledMatchString (@compiledMatchStrings) {

		while ($readBuffer =~ /(.)($compiledMatchString)(.)/g) {

		    my ($pre,$core,$post) = ($1,$2,$3);
		    pos($readBuffer) -= ($merSize+1);
		    my $sCoord = pos($readBuffer)-1;
		
		    if (exists($mers{$core})) {
		    
			$nValidMers++;
		    
			my $preCode = $baseCodes{$pre};
			my $postCode = 6+$baseCodes{$post};
		    
			my $qpre = substr($qualBuffer,$sCoord,1);
			my $qpost = substr($qualBuffer,$sCoord+$merSize+1,1);
		    
			my $Qpre = ord($qpre) - $qualOffset;
			my $Qpost = ord($qpost) - $qualOffset;
		    
			if ($Qpre > $minQuality) {
			    $mers{$core}->[$preCode] += 1;
			} 
			if ($Qpost > $minQuality) {
			    $mers{$core}->[$postCode] += 1;
			}
		    }
		}

		$readBuffer = revComp($readBuffer);
		$qualBuffer = reverse($qualBuffer);
	    
		while ($readBuffer =~ /(.)($compiledMatchString)(.)/g) {

		    my ($pre,$core,$post) = ($1,$2,$3);
		    pos($readBuffer) -= ($merSize+1);
		    my $sCoord = pos($readBuffer)-1;
		
		    if (exists($mers{$core})) {
		    
			my $preCode = $baseCodes{$pre};
			my $postCode = 6+$baseCodes{$post};
		    
			my $qpre = substr($qualBuffer,$sCoord,1);
			my $qpost = substr($qualBuffer,$sCoord+$merSize+1,1);
			
			my $Qpre = ord($qpre) - $qualOffset;
			my $Qpost = ord($qpost) - $qualOffset;
			
			if ($Qpre > $minQuality) {
			    $mers{$core}->[$preCode] += 1;
			}
			if ($Qpost > $minQuality) {
			    $mers{$core}->[$postCode] += 1;	    
			}
		    }
		}
	    }

	    $readBuffer = "X";
	    $qualBuffer = "X";
	    $inBuffer = 0;

	}

	print STDERR "Done.\n";
	close S;
    }
    
    print STDERR "Found $nReads total reads spanning $totalReadLength total bases and including $nValidMers valid $merSize-mers.\n";
    print STDERR "Dumping $merSize-mer count graph...\n";

    while (my ($mer,$countRef) = each(%mers)) {
	my @counts = @{$countRef};
	my @ins = @counts[0..5];
	my $in = 0;
	map {$in += $_} @ins;
	my @outs = @counts[6..11];
	my $out = 0;
	map {$out += $_} @outs;
	if ($in || $out) {
	    print "$mer\t@counts\n";
	}
    }
    print STDERR "Done.\n";

}

# Input mergraph(s)
# Generate UFX

if (defined($merGraphFile)) {
    print STDERR "Reading merGraph from $merGraphFile...\n";

    my $distinctMers = 0;
    my %merTypes = ();
    open (G,$merGraphFile) || die "Couldn't open $merGraphFile\n";
    while (my $line = <G>) {
	chomp $line;
	my ($mer,@counts) = split(/\s+/,$line);
	my $nCounts = scalar(@counts);
	unless ($nCounts == 13) {
	    die "Invalid merGraphFile line: $line\n";
	}
	$mer = uc($mer);

	# Report canonical form only
	my $rcMer = reverse($mer);
	$rcMer =~ tr/ACGT/TGCA/;
	unless ($mer lt $rcMer) {
	    next;
	}

	unless (defined($merSize)) {
	    $merSize = length($mer);
	}

	my @fPaths = @counts[6..11];
	my @rPaths = @counts[0..5];

	my $fPaths = "";	
	my $rPaths = "";	
	for (my $p = 0; $p < 4; $p++) {
	    if ($fPaths[$p] >= $significantDepth) {
		$fPaths .= $bases[$p];
	    }
	    if ($rPaths[$p] >= $significantDepth) {
		$rPaths .= $bases[$p];
	    }
	}
	my $nfPaths = length($fPaths);
	my $nrPaths = length($rPaths);

	my $merType = "";
	if ($nrPaths == 1) {
	    $merType .= $rPaths;
	} elsif ($nrPaths > 1) {
	    $merType .= "F";
	} else {
	    $merType .= "X";
	}
	if ($nfPaths == 1) {
	    $merType .= $fPaths;
	} elsif ($nfPaths > 1) {
	    $merType .= "F";
	} else {
	    $merType .= "X";
	}

	print "$mer\t$merType\n";
	$distinctMers++;
	$merTypes{$merType}++;

    }

    close G;
    print STDERR "Done.\n";
    print STDERR "Found $distinctMers $merSize-mers.\n";

    my @sortedMerTypes = sort {$merTypes{$b} <=> $merTypes{$a}} keys(%merTypes);
    foreach my $merType (@sortedMerTypes) {
	print STDERR "$merType\t$merTypes{$merType}\n";
    }

}

# Input UFX-graph
# Output contigs

if (defined($UFXHashListFile)) {

    my @hashFiles = ();
    open (U,$UFXHashListFile) || die "Couldn't open $UFXHashListFile\n";
    my $hashFileIndex = 0; 
    my %prefixLengths = ();
    while (my $line = <U>) {
	chomp $line;
	my ($file,$prefixSet) = split(/\s+/,$line);
	my @prefixes = split(/\+/,$prefixSet);
	foreach my $p (@prefixes) {
	    my $pLen = length($p);
	    $prefixLengths{$pLen}++;
	    $hashFilePrefixToIndex{$p} = $hashFileIndex;
	}
	push(@hashFiles,$file);
	$hashFileIndex++;
    }
    @hashPrefixSearchOrder = sort {($prefixLengths{$b} <=> $prefixLengths{$a}) || ($a <=> $b)} keys(%prefixLengths);


    my $nHashFiles = scalar(@hashFiles);
    $hashFileIndex = 0; 
    my $merCount = 0;
    foreach my $hashFile (@hashFiles) {
	$hashFileIndex++;
	print STDERR "Reading UFX hash $hashFile ($hashFileIndex/$nHashFiles)\n";
	open (H,$hashFile) || die "Couldn't open $hashFile\n";
	binmode(H);

	my $header = <H>;
	chomp $header;

	my ($hashSize,$hashWordSize,$hashMerSize) = split(/\t/,$header);
	if (defined($merSize)) {
	    unless ($hashMerSize == $merSize) {
		die "hashMerSize inconsistency ($hashMerSize != $merSize)\n";
	    }
	} else {
	    $merSize = $hashMerSize;
	    $merNibbles = 1+ceil($merSize*2/8);
	}
	if (defined($ufxWordSize)) {
	    unless ($hashWordSize == $ufxWordSize) {
		die "hashWordSize inconsistency ($hashWordSize != $ufxWordSize)\n";
	    }
	} else {
	    $ufxWordSize = $hashWordSize;
	}
	push(@ufxHashSizes,$hashSize);
	$merCount += $hashSize/2;  #  Assumes alpha = 2 filling factor

	unless (defined($wordsPerHashElement)) {
	    my $merBits = 2*$hashMerSize;
	    my $extBits = 4;
	    my $storeBits = 0;
	    while ($storeBits < $merBits+2*$extBits) {
		$storeBits += $hashWordSize;
		$wordsPerHashElement++;     
	    }
	}

	my $totalBytes = $hashSize*$wordsPerHashElement*$hashWordSize/8;

	my $hashVector = "";
	my $bytesRead = read(H,$hashVector,$totalBytes);
	unless ($bytesRead == $totalBytes) {
	    die "Error reading hashFile.  Byte count mismatch $totalBytes != $bytesRead\n";
	}

	close H;

	push(@ufxHashArray,$hashVector);

	my $memUsed = memoryUsed();
	if ($memUsed > $maxMemory) {
	    die "Memory used ($memUsed GB) exceeded maximum ($maxMemory GB)\n";
	}
	my $printMem = sprintf("%.2f",$memUsed);
	print STDERR "$merCount mers read.  $printMem GB in use.\n";

    }

    # Traverse U-U k-mers generating contigs
    # Validation step is folded into walk as of v.4h

    print STDERR "Generating contigs...\n";

    my $contigId = 0;
    my $totalShortContigBases = 0;
    my $nShortContigs = 0;
    my $totalLongContigBases = 0;
    my $nLongContigs = 0;

    for (my $hashIndex = 0; $hashIndex < $nHashFiles; $hashIndex++) {	
	my $hashSize = $ufxHashSizes[$hashIndex];
	
	my $printIndex = $hashIndex+1;
	print STDERR "Seeding from hash index $printIndex / $nHashFiles\n";

	for (my $elementIndex = 0; $elementIndex < $hashSize; $elementIndex++) {
	    my $hashElement = getElement($hashIndex,$elementIndex);
	    
	    # a completely empty element is an unused hash entry -Rule 0
	    unless($hashElement) {
		next;
	    }

	    my ($mer,$merCode,$visited) = $hashElement =~ /^(.{$merSize})(.?.?)(.)$/;  #CHECK

#	    my $debug = 0;
#
#	    if ((defined($DEBUGMER)) && (($mer eq $DEBUGMER) || ($mer eq $DEBUGRC))) {
#		print STDERR "DEBUG MER encountered as seed ($mer).  merCode = [$merCode]\n";
#		$debug = 1;
#	    }

	    # if visited is true, do not use as a seed - Rule 3
	    if ($visited) {
		next;
	    }
	    
	    # if non-UU do not use as seed - Rule 1
	    unless ($merCode =~ /^[ACGT][ACGT]$/) {
#  No longer hide elements unless they are visited
#		hideElement($hashIndex,$elementIndex);
		next;
	    }
	    
	    my ($prevBase,$nextBase) = $merCode =~ /^(.)(.)$/;
	    my $takeForwardWalk = 1;
	    my $takeReverseWalk = 1;

	    my $nextMer = substr($mer.$nextBase,1);
	    my $prevMer = $prevBase.$mer;
	    chop $prevMer;

#CHECK 
#  Check that seeding the current mer doesn't preclude the next or previous mer due to palindrome condition -Rule 4* (optional?)
#  If so, discard the current mer. 
#	    my $rcMer = reverse($mer);
#	    $rcMer =~ tr/ACGT/TGCA/;
#	    if ( ($rcMer eq $nextMer) || ($rcMer eq $prevMer)) {
#		hideElement($hashIndex,$elementIndex);
#		next;
#	    }
#ENDCHECK

	    my $nextMerInfo = lookUpMer($nextMer);

	    my $nextMerCode = undef;
	    my $nextHashIndex = undef;
	    my $nextElementIndex = undef;
	    my $nextVisited = undef;
	    if ($nextMerInfo) {
		($nextMerCode,$nextHashIndex,$nextElementIndex,$nextVisited) = split(/\./,$nextMerInfo); #CHECK
		if ($nextMerCode =~ /^[ACGT][ACGT]$/) {
		    my ($nextPrevBase,$nextNextBase) = $nextMerCode =~ /^(.)(.)$/;
		    my $firstBase = substr($mer,0,1);
		    if ($nextPrevBase ne $firstBase) { 	# Non-reciprocal forward U-U linkage - Rule 2a
#  No longer hide elements unless they are visited
#			hideElement($hashIndex,$elementIndex);

#			if ($debug) {
#			    print STDERR "DEBUG non-reciprocal UU: $nextPrevBase $firstBase\n";
#			}

			next;
		    }
		} else {   # Next mer is not U-U
		    $takeForwardWalk = 0;
		}
	    } else {  # Next mer does not exist
		$takeForwardWalk = 0;
	    }

#	    if ($debug) {
#		print STDERR "DEBUG takeForwardWalk = $takeForwardWalk\n";
#	    }

	    my $prevMerInfo = lookUpMer($prevMer);

	    my $prevMerCode = undef;
	    my $prevHashIndex = undef;
	    my $prevElementIndex = undef;
	    my $prevVisited = undef;
	    if ($prevMerInfo) {
		($prevMerCode,$prevHashIndex,$prevElementIndex,$prevVisited) = split(/\./,$prevMerInfo); #CHECK
		if ($prevMerCode =~ /^[ACGT][ACGT]$/) {
		    my ($prevPrevBase,$prevNextBase) = $prevMerCode =~ /^(.)(.)$/;
		    my $lastBase = substr($mer,-1);
		    if ($prevNextBase ne $lastBase) { 	# Non-reciprocal reverse U-U linkage - Rule 2b
#  No longer hide elements unless they are visited
#			hideElement($hashIndex,$elementIndex);

#			if ($debug) {
#			    print STDERR "DEBUG non-reciprocal UU: $prevNextBase $lastBase\n";
#			}

			next;
		    }
		} else {   # Prev mer is not U-U
		    $takeReverseWalk = 0;
		}
	    } else {  # Prev mer does not exist 
		$takeReverseWalk = 0;
	    }

#	    if ($debug) {
#		print STDERR "DEBUG takeReverseWalk = $takeReverseWalk\n";
#	    }

	    # hide the seed - it is now to be visited
	    hideElement($hashIndex,$elementIndex);
	    
	    my $forwardWalk = "";
	    if ($takeForwardWalk) {
#		$forwardWalk = walk($nextMer,$nextMerInfo,$debug);
		$forwardWalk = walk($nextMer,$nextMerInfo);
	    }

#	    if ($debug) {
#		print STDERR "DEBUG forwardWalk returned [$forwardWalk]\n";
#	    }

	    my $reverseWalk = "";
	    if ($takeReverseWalk) {
		my $rcPrevMer = reverse($prevMer);
		$rcPrevMer =~ tr/ACGT/TGCA/;
		my $rcPrevCode = reverse($prevMerCode);
		$rcPrevCode =~ tr/ACGT/TGCA/;
		my $rcPrevMerInfo = "$rcPrevCode.$prevHashIndex.$prevElementIndex.$prevVisited"; #CHECK
#		$reverseWalk = walk($rcPrevMer,$rcPrevMerInfo,$debug);
		$reverseWalk = walk($rcPrevMer,$rcPrevMerInfo);
	    }

	    if ($reverseWalk) {
		$reverseWalk = reverse($reverseWalk);
		$reverseWalk =~ tr/ACGT/TGCA/;
	    }

#	    if ($debug) {
#		print STDERR "DEBUG reverseWalk returned [$reverseWalk]\n";
#	    }

	    my $result = $reverseWalk . $mer . $forwardWalk;
	    my $contigLength = length($result);
	    
	    if ($contigLength > $minContigSize-1) {
		printFasta("Contig$contigId",$result);
#		my $fwl = length($forwardWalk);
#		my $rwl = length($reverseWalk);
#		print STDERR "Contig$contigId generated from seed $mer has length $contigLength [$rwl][$merSize][$fwl]\n";
		$contigId++;
		$totalLongContigBases += $contigLength;
		$nLongContigs++;
	    } else {
		$totalShortContigBases += $contigLength;
		$nShortContigs++;
	    }
	}
    }
    print STDERR "Done.\n";
    print STDERR "Generated $nLongContigs contigs >= $minContigSize bases totalling $totalLongContigBases bases\n";
    print STDERR "Discarded $nShortContigs contigs < $minContigSize bases totalling $totalShortContigBases bases\n";
}


# --- SUBROUTINES ---

sub walk {
#    my ($seed,$seedInfo,$debug) = @_;
    my ($seed,$seedInfo) = @_;

    my $currentMer = $seed;
    my ($currentMerCode,$currentHashIndex,$currentElementIndex,$currentVisited) = split(/\./,$seedInfo);  #CHECK

    my $addBase = substr($seed,-1);
    my $walk = "";

#    my $localDebug = $debug;
    
    while () {
#	if (defined($DEBUGMER) && (($currentMer eq $DEBUGMER) || ($currentMer eq $DEBUGRC))) {
#	    print STDERR "walk: DEBUGMER found ($currentMer)\n";
#	    $localDebug = 1;
#	} else {
#	    $localDebug = $debug;
#	}

#CHECK
	# Check the up-to-date visitation status of the current mer
	$currentVisited = checkVisited($currentHashIndex,$currentElementIndex);

#	print STDERR "$currentHashIndex.$currentElementIndex $currentVisited\n";

	# If we encounter a previously visited mer terminate walk and do not add the base - Rule 3
	if ($currentVisited) {
	    last;
	}
#ENDCHECK

	my ($prevBase,$nextBase) = $currentMerCode =~ /^(.)(.)$/;
	my $nextMer = substr($currentMer.$nextBase,1);

#CHECK
#  Check that adding the current mer doesn't preclude the next mer due to palindrome condition - Rule 4* optional
#  If so, terminate the walk and do not add the base. 
#	my $currentRCmer = reverse($currentMer);
#	$currentRCmer=~ tr/ACGT/TGCA/;
#	if ($nextMer eq $currentRCmer) {
#	    hideElement($currentHashIndex,$currentElementIndex);
#	    last;
#	}
#ENDCHECK

	my $nextMerInfo = lookUpMer($nextMer);

#	if ($localDebug) {
#	    print STDERR "[$currentMer] [$prevBase][$nextBase] [$nextMer] [$nextMerInfo]\n";
#	}

	my $nextMerCode = undef;
	my $nextHashIndex = undef;
	my $nextElementIndex = undef;
	my $nextVisited = undef;
	
	if ($nextMerInfo) {
	    ($nextMerCode,$nextHashIndex,$nextElementIndex,$nextVisited) = split(/\./,$nextMerInfo);
	    if ($nextMerCode =~ /^[ACGT][ACGT]$/) {     # Next mer is U-U
		my ($nextPrevBase,$nextNextBase) = $nextMerCode =~ /^(.)(.)$/;
		my $firstBase = substr($currentMer,0,1);

		if ($nextPrevBase eq $firstBase) { 	# step is reciprocally validated, visit and iterate
		    hideElement($currentHashIndex,$currentElementIndex);
		    $walk .= $addBase;
		    $addBase = $nextBase;
		    $currentMer = $nextMer;
		    ($currentMerCode,$currentHashIndex,$currentElementIndex,$currentVisited) = ($nextMerCode,$nextHashIndex,$nextElementIndex,$nextVisited); #CHECK

#		    if ($localDebug) {
#			print STDERR "VALID\n";
#		    }

		} else {       # Non-reciprocal U-U linkage, terminate without adding last base - Rule 2a
#  No longer hide elements unless they are visited
#		    hideElement($currentHashIndex,$currentElementIndex);

#		    if ($localDebug) {
#			print STDERR "NONRECIP [$nextPrevBase] [$firstBase]\n";
#		    }

		    last;
		}

	    } else {   # Next mer is not U-U, terminate adding last base
		hideElement($currentHashIndex,$currentElementIndex);
#  No longer hide elements unless they are visited
#		hideElement($nextHashIndex,$nextElementIndex);
		$walk .= $addBase;
		
#		if ($localDebug) {
#		    print STDERR "next mer non-UU [$nextMerCode]\n";
#		}

		last;
	    }
	} else {   # Next mer does not exist; terminate walk, and visit/hide current mer 
	    hideElement($currentHashIndex,$currentElementIndex);
	    $walk .= $addBase;

#	    if ($localDebug) {
#		print STDERR "next mer does not exist [$nextMerInfo]\n";
#	    }

	    last;
	}
    }

    return $walk;

}

sub revComp {
    my ($seq) = @_;
    my $rc = reverse($seq);
    $rc =~ tr/ACGT/TGCA/;
    return $rc;
}

sub printFasta {
    my ($name,$seq) = @_;
    my $bpl = 50;
    my $nLines = sprintf("%d",length($seq)/$bpl);
    if ((length($seq) % $bpl) != 0) {
        $nLines++;
    }

    print ">$name\n";
    for (my $j = 0;$j<$nLines;$j++) {
        my $seqLine = substr($seq,$j*$bpl,$bpl);
        my $text = ($seqLine . "\n");
        print $text;
    }
}

sub memoryUsed {
  my $procTable = new Proc::ProcessTable;
  foreach my $proc ( @{$procTable->table} ) {
    unless ($proc->pid eq $$) {
        next;
    }
    return ($proc->size)/(1024*1024*1024);
  }
}

sub hideElement {
    my ($hashIndex,$elementIndex) = @_;

#    reminder: this is a global
#    $doNotVisitCode = 0b1000;

    my $wordOffset = $wordsPerHashElement*$elementIndex;
    my $getWord = vec($ufxHashArray[$hashIndex],$wordOffset,$ufxWordSize);

    # Do not hide vacant elements
    if ($getWord == 0) {
        return 0;
    }

    my $newBits = "";
    vec($newBits,0,$ufxWordSize) = $getWord;
    vec($newBits,0,4) = (vec($newBits,0,4) | 0b1000);   #    $doNotVisitCode = 0b1000;  CHECK
    vec($newBits,1,4) = (vec($newBits,1,4) | 0b1000);   #    $doNotVisitCode = 0b1000;  CHECK
    my $newWord = vec($newBits,0,$ufxWordSize);

    vec($ufxHashArray[$hashIndex],$wordOffset,$ufxWordSize) = $newWord;

    return 1;

}

# Efficient way to check visitation status, will die if a non-existent element is queried
sub checkVisited {

    my ($hashIndex,$elementIndex) = @_;
    
    my $doesNotExist = 0;
    my $wordOffset = $wordsPerHashElement*$elementIndex;
    my $getWord = vec($ufxHashArray[$hashIndex],$wordOffset,$ufxWordSize);
    
    if ($getWord == 0) {
	die "Error: attempted to check visited property of a non-existent element ($hashIndex,$elementIndex)\n";
    }
    
    my $elementBits = "";
    vec($elementBits,0,$ufxWordSize) = $getWord;

    my $v = vec($elementBits,0,4);
    my $visited = $v & 0b1000;

    return $visited;

}

sub getElement {
    my ($hashIndex,$elementIndex) = @_;
    
    my $doesNotExist = 0;
    my $wordOffset = $wordsPerHashElement*$elementIndex;
    my $getWord = vec($ufxHashArray[$hashIndex],$wordOffset,$ufxWordSize);
    
    if ($getWord == 0) {
	return $doesNotExist;
    }
    
    my $elementBits = "";
    vec($elementBits,0,$ufxWordSize) = $getWord;
    for (my $o = 1; $o <= $wordsPerHashElement; $o++) {
	vec($elementBits,$o,$ufxWordSize) = vec($ufxHashArray[$hashIndex],$wordOffset+$o,$ufxWordSize);
    }

    #CHECK
    my $v0 = vec($elementBits,0,4);
    my $v1 = vec($elementBits,1,4);

    my $ext1 = $hashExtCodes{$v0 & 0b0111};
    my $ext2 = $hashExtCodes{$v1 & 0b0111};

    my $visited1 = $v0 & 0b1000;
    my $visited2 = $v1 & 0b1000;

    unless ($visited1 == $visited2) {
	die "That's a bug... visited status should be the same for both extensions: $visited1 != $visited2\n";
    }
    #ENDCHECK

    my $mer = "";
    for (my $b = 1; $b < $merNibbles; $b++) {
	$mer .= $hashBaseCodes2{vec($elementBits,$b,8)};
    }
    $mer = substr($mer,0,$merSize);
    return "$mer$ext1$ext2$visited1"; #CHECK
}

sub lookUpMer {
    my ($mer) = @_;

    # Retrieve extension data for a mer from the hash
    # Return value is extension.hashIndex.elementIndex.visited
    # returns 0 if the mer is not found in the hash 

    my $rcMer = reverse($mer);
    $rcMer =~ tr/ACGT/TGCA/;

    my ($canonicalMer,$rcFlag) = ($mer lt $rcMer) ? ($mer,0) : ($rcMer,1);

    my $hashIndex = undef;

    foreach my $prefixLen (@hashPrefixSearchOrder) {
	my $checkPrefix = substr($canonicalMer,0,$prefixLen);
	if (exists($hashFilePrefixToIndex{$checkPrefix})) {
	    $hashIndex = $hashFilePrefixToIndex{$checkPrefix};
	    last;
	}
    }
    
    unless (defined($hashIndex)) {
	return 0;
#	die "lookUpMer: No hash prefix found for $mer ($canonicalMer)\n";
    }

    my $hashSize = $ufxHashSizes[$hashIndex];
    my $elementIndex = hex(B::hash($canonicalMer)) % $hashSize;

    while (1) {

	my $checkElement = getElement($hashIndex,$elementIndex);
        if (!$checkElement) {
            return 0;
	} else {
	    my ($checkMer,$checkExtensions,$visited) = $checkElement =~ /^(.{$merSize})(.?.?)(.)$/;  #CHECK
	    if ($checkMer eq $canonicalMer) {
		
		my $returnValue = $checkExtensions;
		if ($rcFlag) {
		    $returnValue = reverse($checkExtensions);
		    $returnValue =~ tr/ACGT/TGCA/;
		}
		$returnValue .= ".$hashIndex.$elementIndex.$visited";  #CHECK
		return $returnValue;
	    } else {
		$elementIndex++;
		if ($elementIndex == $hashSize) {
		    $elementIndex = 0;
		}
	    }
	}
    }
}
