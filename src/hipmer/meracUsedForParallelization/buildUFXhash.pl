#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# buildUFXhash.pl by Jarrod Chapman <jchapman@lbl.gov> Sun Apr 22 08:16:56 PDT 2012
# Copyright 2012 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
use B;

my %opts = ();
my $validLine = getopts('u:n:o:P:U', \%opts);
my @required = ("u","n","o");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print STDERR "Usage: ./buildUFXhash.pl <-u ufxFile> <-n numberOfmers> <-o outputFile> <<-(p)U(rgeNonCanonical?)>> <<-P prefixSet>>\n";
    exit;
}

my $ufxFile = $opts{"u"};
my $nToHash = $opts{"n"};
my $outputFile = $opts{"o"};

my $purgeMode = 0;
if (exists($opts{"U"})) {
    $purgeMode = 1;
}

my $prefixMode = 0;
my @prefixSet = ();
if (exists($opts{"P"})) {
    @prefixSet = split(/\+/,$opts{"P"});
    $prefixMode = 1;
}

my $bufferSize = 1000000;
my @merBuffer = ();
my $inBuffer = 0;
my $alpha = 2;
my $hashSize = $nToHash*$alpha;
my $hashVector = "";
my $loadCollisions = 0;
#my $getCollisions = 0;
my $merSize = undef;
my $storeBits = 0;
my $merBits = undef;
my $extBits = 4;
my $wordSize = 32;
my $wordsPerElement = 0;

my %baseCodes = ("A",0b00,"C",0b01,"G",0b10,"T",0b11);
@baseCodes{values %baseCodes} = keys (%baseCodes);

my %extCodes = ("A",0b0100,"C",0b0101,"G",0b0110,"T",0b0111,"F",0b0001,"X",0b0010);
@extCodes{values %extCodes} = keys (%extCodes);

open (U,$ufxFile) || die "Couldn't open $ufxFile.\n";
print STDERR "Reading $ufxFile...\n";
my $nHashed = 0;
while (my $line = <U>) {
    chomp $line;
    my @cols = split(/\s+/,$line);

    # Discard mers not matching prefix set
    if ($prefixMode) {
	my $goodMer = 0;
	my ($mer,$code) = @cols;
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

    # Discard non-canonical mers

    if ($purgeMode) {
	my ($mer,$code) = @cols;
	my $rcMer = reverse($mer);
	$rcMer =~ tr/ACGT/TGCA/;
	unless ($mer lt $rcMer) {
	    next;
	}
    }

    push(@merBuffer,$cols[0].$cols[1]);
    unless (defined($merSize)) {
	$merSize = length($cols[0]);
	$merBits = 2*$merSize;
	while ($storeBits < $merBits+2*$extBits) {
	    $storeBits += $wordSize;
	    $wordsPerElement++;	    
	}
	vec($hashVector,($hashSize*$wordsPerElement)-1,$wordSize) = 0;
	print STDERR "Building hash using $storeBits ($wordsPerElement x $wordSize words) bits per key/value pair\n";
    }
    $inBuffer++;
    if ($inBuffer >= $bufferSize) {
	foreach my $kv (@merBuffer) {
	    my ($k,$v) = $kv =~ /^(.{$merSize})(..)$/;
	    $loadCollisions += addToHash($k,$v);
	    $nHashed++;
	    if ($nHashed > $nToHash) {
		die "number of hashed elements ($nHashed) exceeds expected ($nToHash)\n";
	    } 

#	    my $check = "";
#	    my $checkCollisions = getFromHash($k,\$check);
#
#	    if (($checkCollisions == -1) || ($check ne $kv)) {
#		die "[$checkCollisions] [$check] [$kv]\n";
#	    } else {
#		$getCollisions += $checkCollisions;
#	    }
	    
	}
	@merBuffer = ();
	$inBuffer = 0;
    }
}
close U;

if ($inBuffer) {
    foreach my $kv (@merBuffer) {
	my ($k,$v) = $kv =~ /^(.{$merSize})(..)$/;
	$loadCollisions += addToHash($k,$v);
	$nHashed++;
	if ($nHashed > $nToHash) {
	    die "number of hashed elements ($nHashed) exceeds expected ($nToHash)\n";
	} 

#	my $check = "";
#	my $checkCollisions = getFromHash($k,\$check);
#	
#	if (($checkCollisions == -1) || ($check ne $kv)) {
#	    die "[$checkCollisions] [$check] [$kv]\n";
#	} else {
#	    $getCollisions += $checkCollisions;
#	}
#


    }
    @merBuffer = ();
    $inBuffer = 0;
}

my $hashLength = length($hashVector);
print STDERR "loaded $nHashed elements into hash of length $hashLength.  Encountered $loadCollisions collisions.\n";

open (OUT, ">$outputFile") || die "Couldn't open $outputFile\n";    
binmode(OUT);
print OUT "$hashSize\t$wordSize\t$merSize\n$hashVector";
close OUT;
print STDERR "Done.\n";

# --- subroutines ----

sub bitify {
    my ($in) = @_;
    my @in = split(//,$in);
    my $out = "";
    vec($out,0,4) = $extCodes{$in[$merSize]};
    vec($out,1,4) = $extCodes{$in[$merSize+1]};
    for (my $b = 0; $b < $merSize; $b++) {
	vec($out,($b+4),2) = $baseCodes{$in[$b]};
   }
    
    return $out;
}

sub unbitify {
    my ($bits)  = @_;
    my $ext1 = $extCodes{vec($bits,0,4)};
    my $ext2 = $extCodes{vec($bits,1,4)};

    my $mer = "";
    for (my $b = 0; $b < $merSize; $b++) {
	$mer .= $baseCodes{vec($bits,($b+4),2)};
    }
    return "$mer$ext1$ext2";
}

sub loadElementToVector {
    my ($data,$elementOffset,$startWord,$endWord) = @_;

    my $wordOffset = $wordsPerElement*$elementOffset;
    for (my $o = $startWord; $o <= $endWord; $o++) {
	my $loadWord = vec($data,$o,$wordSize);
	vec($hashVector,$wordOffset+$o,$wordSize) = $loadWord;
    }
}

sub getElementFromVector {
    my ($elementOffset,$startWord,$endWord) = @_;

    my $wordOffset = $wordsPerElement*$elementOffset;
    my $data = "";
    for (my $o = $startWord; $o <= $endWord; $o++) {
	my $getWord = vec($hashVector,$wordOffset+$o,$wordSize);
	vec($data,$o,$wordSize) = $getWord;
    }

    return $data;

}

sub addToHash {
    my ($key,$value) = @_;

    my $storeVector = bitify("$key$value");

    my $hashVal = hex(B::hash($key)) % $hashSize;
    my $collisions = 0;

    my $isOccupied = 1;
    while ($isOccupied) {
	my $checkWord = getElementFromVector($hashVal,0,0);
	my $checkSum = unpack("%$wordSize"."b*",$checkWord);
	if ($checkSum) {
	    $hashVal++;
	    $collisions++;
	    if ($hashVal == $hashSize) {
		$hashVal = 0;
	    }
	} else {
	    $isOccupied = 0;
	}
    }

    loadElementToVector($storeVector,$hashVal,0,$wordsPerElement-1);

    return $collisions;
}

sub getFromHash {
    my ($key,$valueRef) = @_;

    my $hashVal = hex(B::hash($key)) % $hashSize;
    my $collisions = 0;
    my $doesNotExist = -1;

    my $foundKey = 0;
    while ($foundKey == 0) {
	my $checkWord = getElementFromVector($hashVal,0,0);
	my $checkSum = unpack("%$wordSize"."b*",$checkWord);

	if ($checkSum == 0) {
	    return $doesNotExist;
	}

	my $checkElementBits = getElementFromVector($hashVal,0,$wordsPerElement-1);
	my $checkElement = unbitify($checkElementBits);
	my ($checkMer,$checkExtensions) = $checkElement =~ /^(.{$merSize})(..)$/;
	
	if ($checkMer eq $key) {
	    $$valueRef = $checkElement;
	    return $collisions;
	} else {
	    $hashVal++;
	    $collisions++;
	    if ($hashVal == $hashSize) {
		$hashVal = 0;
	    }
	}
    }
}
    
