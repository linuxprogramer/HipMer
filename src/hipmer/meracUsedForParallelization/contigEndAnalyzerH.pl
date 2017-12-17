#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# contigEndAnalyzer.pl by Jarrod Chapman <jchapman@lbl.gov> Thu Dec  2 13:38:11 PST 2010
# Copyright 2010 Jarrod Chapman. All rights reserved.
#

use B;
use strict;
use Getopt::Std;
use POSIX qw (ceil);

my %opts = ();
my $validLine = getopts('c:u:m:', \%opts);
my @required = ("c","u","m");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print "Usage: ./contigEndAnalyzer.pl <-u ufxHashListFile> <-c contigFile> <-m merSize>\n";
    exit;
}

my $UFXHashListFile = $opts{"u"};
my $contigFile = $opts{"c"};
my $merSize = $opts{"m"};
my $merNibbles = 1+ceil($merSize*2/8);

open (U,$UFXHashListFile) || die "Couldn't open $UFXHashListFile\n";
my $hashFileIndex = 0; 
my %prefixLengths = ();
my %hashFilePrefixToIndex = ();
my @hashPrefixSearchOrder = ();
my @hashFiles = ();
my @mersFromHash = ();
my %hashExtCodes = ("A",0b0100,"C",0b0101,"G",0b0110,"T",0b0111,"F",0b0001,"X",0b0010);
@hashExtCodes{values %hashExtCodes} = keys (%hashExtCodes);
my @bases = ("A","C","G","T");

my %hashBaseCodes = ();
for (my $i = 0; $i < 4**4; $i++) { 
    my $mer = ""; 
    my $k = $i; 
    for (my $j = 0; $j < 4; $j++) { 
        $mer .= "$bases[$k%4]"; 
        $k = sprintf("%d",$k/4); 
    } 
    $hashBaseCodes{$mer} = $i;
}
@hashBaseCodes{values %hashBaseCodes} = keys (%hashBaseCodes);

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
    $mersFromHash[$hashFileIndex] = [];
    $hashFileIndex++;
}
@hashPrefixSearchOrder = sort {($prefixLengths{$b} <=> $prefixLengths{$a}) || ($a <=> $b)} keys(%prefixLengths);

my %contigInfo = ();
my @contigs = ();
my $currentEntry = undef;
my $seq = "";
print STDERR "Reading $contigFile...\n";
open (C,$contigFile) || die "Couldn't open $contigFile\n";
while (my $line = <C>) {
    chomp $line;
    if ($line =~ /^>/) {
        if (defined($currentEntry)) {
            processSequence($currentEntry,$seq);
        }
        $seq = "";
        ($currentEntry) = $line =~ /^>(.+)$/;
    } else {
        $seq .= $line;
    }
}
close C;

if (defined($currentEntry)) {
    processSequence($currentEntry,$seq);
}
my $nContigs = scalar(@contigs);
print STDERR "Done.  Read $nContigs contigs\n";
#for (my $hashIndex = 0; $hashIndex < $hashFileIndex; $hashIndex++) {
#    my $nEnds = scalar(@{$mersFromHash[$hashIndex]});
#    print STDERR "Found $nEnds from hashIndex $hashIndex\n";
#}

my $nHashFiles = scalar(@hashFiles);
$hashFileIndex = 0; 
my $ufxWordSize = undef;
my $wordsPerHashElement = undef;
my @ufxHashSizes = ();
my $hashVector = "";

# First pass through hashes: find extension code of last UU-mer of each contig end

foreach my $hashFile (@hashFiles) {
    my $printIndex = $hashFileIndex+1;
    print STDERR "Loading UFX hash $hashFile ($printIndex/$nHashFiles)\n";
    open (H,$hashFile) || die "Couldn't open $hashFile\n";
    binmode(H);

    my $header = <H>;
    chomp $header;

    my ($hashSize,$hashWordSize,$hashMerSize) = split(/\t/,$header);
    unless ($hashMerSize == $merSize) {
	die "hashMerSize inconsistency ($hashMerSize != $merSize)\n";
    }
    if (defined($ufxWordSize)) {
	unless ($hashWordSize == $ufxWordSize) {
	    die "hashWordSize inconsistency ($hashWordSize != $ufxWordSize)\n";
	}
    } else {
	$ufxWordSize = $hashWordSize;
    }
    push(@ufxHashSizes,$hashSize);

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

    $hashVector = "";
    my $bytesRead = read(H,$hashVector,$totalBytes);
    unless ($bytesRead == $totalBytes) {
	die "Error reading hashFile.  Byte count mismatch $totalBytes != $bytesRead\n";
    }
    close H;

    my $merHashRef = $mersFromHash[$hashFileIndex];
    my $nMersInHash = scalar(@{$merHashRef});

    for (my $m = 0; $m < $nMersInHash; $m++) {

	my $endMer = $merHashRef->[$m];

	my ($contig,$end,$mer) = split(/:/,$endMer);
	my $merInfo = lookUpMer($mer,$hashFileIndex);
	unless ($merInfo) {
	    die "Error: No extension information found for $mer\n";
	}
	my ($ext0,$ext1) = $merInfo =~ /^([ACGT])([ACGT])\./;
	unless (defined($ext0) && defined($ext1)) {
	    die "Error: non-UU extensions for $mer ($merInfo)\n";
	}

	if ($end == 0) {
	    $endMer = "$ext0:$contig:$end:$mer";
	} elsif ($end == 1) {
	    $endMer = "$ext1:$contig:$end:$mer";
	} else {
	    die "Error: invalid hashEntry for at $m ($endMer)\n";
	}

	$merHashRef->[$m] = $endMer;

    }

    my $terminator = 0;
    push(@{$merHashRef},$terminator);

    $hashFileIndex++;

}

# Transition to extensions

for (my $h = 0; $h < $hashFileIndex; $h++) {

    my $merHashRef = $mersFromHash[$h];
    my $currentMer = shift(@{$merHashRef});
    while ($currentMer) {
	my ($ext,$contig,$end,$mer) = split(/:/,$currentMer);

	my $newMer = $mer;
	if ($end == 0) {
	    $newMer = $ext.$newMer;
	    chop $newMer;
	} else {
	    $newMer .= $ext;
	    $newMer = substr($newMer,1);
	}

	my $newIndex = getHashIndex($newMer);
	push(@{$mersFromHash[$newIndex]},$currentMer);

	$currentMer = shift(@{$merHashRef});
    }
}

# Second pass through hashes: find the next extension code

$hashFileIndex = 0; 
foreach my $hashFile (@hashFiles) {
    my $printIndex = $hashFileIndex+1;
    print STDERR "Loading UFX hash $hashFile ($printIndex/$nHashFiles)\n";
    open (H,$hashFile) || die "Couldn't open $hashFile\n";
    binmode(H);

    my $header = <H>;
    chomp $header;

    my ($hashSize,$hashWordSize,$hashMerSize) = split(/\t/,$header);
    my $totalBytes = $hashSize*$wordsPerHashElement*$hashWordSize/8;

    $hashVector = "";
    my $bytesRead = read(H,$hashVector,$totalBytes);
    unless ($bytesRead == $totalBytes) {
	die "Error reading hashFile.  Byte count mismatch $totalBytes != $bytesRead\n";
    }
    close H;

    my $merHashRef = $mersFromHash[$hashFileIndex];
    my $nMersInHash = scalar(@{$merHashRef});

    for (my $m = 0; $m < $nMersInHash; $m++) {

	my $extEndMer = $merHashRef->[$m];
	my ($ext,$contig,$end,$mer) = split(/:/,$extEndMer);

	my $newMer = $mer;
	if ($end == 0) {
	    $newMer = $ext.$newMer;
	    chop $newMer;
	} else {
	    $newMer .= $ext;
	    $newMer = substr($newMer,1);
	}	

	my $merInfo = lookUpMer($newMer,$hashFileIndex);
	$contigInfo{$contig} .= ":$merInfo:$extEndMer";
    }
    $hashFileIndex++;
}

foreach my $contig (@contigs) {

    unless (exists($contigInfo{$contig})) {
	die "Error: information missing for $contig\n";
    }    
    my $contigInfo = $contigInfo{$contig};

    my ($length,$code1,$base1,$c1,$end1,$mer1,$code2,$base2,$c2,$end2,$mer2) = split(/:/,$contigInfo);

    unless ($c1 eq $c2) {
	die "Error: contig ID mismatch $c1 != $c2\n";
    }
    unless ( (($end1 == 0) && ($end2 == 1)) || (($end1 == 1) && ($end2 == 0))) {
	die "Error: Invalid contig ends $end1,$end2\n";
    }
    
    my ($prevBase,$nextBase) = ($end1 == 0) ? ($base1,$base2) : ($base2,$base1);
    my ($prevCode,$nextCode) = ($end1 == 0) ? ($code1,$code2) : ($code2,$code1);
    my ($firstMer,$lastMer) = ($end1 == 0) ? ($mer1,$mer2) : ($mer2,$mer1);

    if ($prevCode =~ /^..\.\d+\.\d+$/) {
	$prevCode = substr($prevCode,0,2);
    } else {
	$prevCode = "00";
    }
    if ($nextCode =~ /^..\.\d+\.\d+$/) {
	$nextCode = substr($nextCode,0,2);
    } else {
	$nextCode = "00";
    }

    print "$contig\t[$prevCode]\t($prevBase)\t$firstMer\t$length\t$lastMer\t($nextBase)\t[$nextCode]\n";
}

sub processSequence {
    my ($name,$seq) = @_;

    my $length = length($seq);

    push(@contigs,$name);
    $contigInfo{$name} = $length;

    my $firstMer = substr($seq,0,$merSize);
    my $lastMer = substr($seq,-$merSize);

    my $hashIndex = getHashIndex($firstMer);
    push(@{$mersFromHash[$hashIndex]},"$name:0:$firstMer");
    $hashIndex = getHashIndex($lastMer);
    push(@{$mersFromHash[$hashIndex]},"$name:1:$lastMer");
}

sub getHashIndex {
    my ($mer) = @_;

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
	die "getHashIndex: No hash prefix found for $mer ($canonicalMer)\n";
    }

    return $hashIndex;
}
    
sub getElement {
#    my ($hashIndex,$elementIndex) = @_;
    my ($elementIndex) = @_;
    
    my $doesNotExist = 0;
    my $wordOffset = $wordsPerHashElement*$elementIndex;
#    my $getWord = vec($ufxHashArray[$hashIndex],$wordOffset,$ufxWordSize);
    my $getWord = vec($hashVector,$wordOffset,$ufxWordSize);
    
    if ($getWord == 0) {
        return $doesNotExist;
    }
    
    my $elementBits = "";
    vec($elementBits,0,$ufxWordSize) = $getWord;
    for (my $o = 1; $o <= $wordsPerHashElement; $o++) {
#        vec($elementBits,$o,$ufxWordSize) = vec($ufxHashArray[$hashIndex],$wordOffset+$o,$ufxWordSize);
        vec($elementBits,$o,$ufxWordSize) = vec($hashVector,$wordOffset+$o,$ufxWordSize);
    }
    
    my $ext1 = $hashExtCodes{vec($elementBits,0,4)};
    my $ext2 = $hashExtCodes{vec($elementBits,1,4)};    
    my $mer = "";
    for (my $b = 1; $b < $merNibbles; $b++) {
        $mer .= $hashBaseCodes{vec($elementBits,$b,8)};
    }
    $mer = substr($mer,0,$merSize);
    return "$mer$ext1$ext2";
}


sub lookUpMer {
    my ($mer,$hashIndex) = @_;

    # Retrieve extension data for a mer from the hash
    # Return value is extension.hashIndex.elementIndex
    # returns 0 if the mer is not found in the hash 

    my $rcMer = reverse($mer);
    $rcMer =~ tr/ACGT/TGCA/;

    my ($canonicalMer,$rcFlag) = ($mer lt $rcMer) ? ($mer,0) : ($rcMer,1);

    my $hashSize = $ufxHashSizes[$hashIndex];
    my $elementIndex = hex(B::hash($canonicalMer)) % $hashSize;

    while (1) {

        my $checkElement = getElement($elementIndex);
        if (!$checkElement) {
            return 0;
        } else {
            my ($checkMer,$checkExtensions) = $checkElement =~ /^(.{$merSize})(.?.?)$/;
            if ($checkMer eq $canonicalMer) {
                
                my $returnValue = $checkExtensions;
                if ($rcFlag) {
                    $returnValue = reverse($checkExtensions);
                    $returnValue =~ tr/ACGT/TGCA/;
                }
                $returnValue .= ".$hashIndex.$elementIndex";
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

