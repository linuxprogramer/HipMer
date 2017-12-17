#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# contigMerDepth.pl by Jarrod Chapman <jchapman@lbl.gov> Wed Jan 19 20:48:45 PST 2011
# Copyright 2011 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
my %opts = ();
my $validLine = getopts('m:f:P:', \%opts);
my @required = ("m","f",);
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print "Usage: ./contigMerDepth.pl <<-P prefixes>> <-m merCountFile> <-f fastaFile>\n";
    exit;
}

my $date = `date`;
chomp $date;
print STDERR "$date\ncontigMerDepth";
while (my ($opt,$val) = each(%opts)) {
    print STDERR " -$opt=$val";
}
print STDERR "\n";


my $merCountFile = $opts{"m"};
my $fastaFile = $opts{"f"};
my @prefixes = ();
my $nPrefixes = 0;
if (exists($opts{"P"})) {
    @prefixes = split(/\+/,$opts{"P"});
    $nPrefixes = scalar(@prefixes);
}

my %mers = ();
my $merSize = undef;

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
        
    if ($nPrefixes > 0) {
	my $goodMer = 0;
	foreach my $p (@prefixes) {
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
    if (exists($mers{$mer})) {
	warn "$mer occurs more than once in $merCountFile.  Using last occurrence.\n";
    }
    $mers{$mer} = $count;
    $nMers++;
}
close M;
print STDERR "Done. Found $nMers $merSize-mers.\n";

my @matchStrings = ("[ACGT]{$merSize}");
if (exists($opts{"P"})) {
    @matchStrings = ();

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

my $currentEntry = undef;
my $seq = "";
open (F,$fastaFile) || die "Couldn't open $fastaFile\n";
while (my $line = <F>) {
    chomp $line;
    if ($line =~ /^>/) {
        if (defined($currentEntry)) {
            processSequence($currentEntry,$seq);
        }
        $seq = "";
        ($currentEntry) = $line =~ />(\S+)/;
    } else {
        $seq .= $line;
    }
}
close F;

if (defined($currentEntry)) {
    processSequence($currentEntry,$seq);
}

$date = `date`;
chomp $date;
print STDERR "Done. $date\n";

sub processSequence {
    my ($name,$seq) = @_;
    my $seqLen = length($seq);
    
    unless ($seqLen >= $merSize) {
        return;
    }

    my %counts = ();
    my $n = 0;
    my $mean = 0;
    my $nDistinct = 0;

    foreach my $compiledMatchString (@compiledMatchStrings) {
	while ($seq =~ /($compiledMatchString)/g) {
	    my $mer = $1;
	    pos($seq) -= ($merSize-1);
	    if (exists($mers{$mer})) {
		my $count = $mers{$mer};
		$mean += $count;
		$n++;
		if (exists($counts{$count})) {
		    $counts{$count} += 1;
		} else {
		    $counts{$count} = 1;
		    $nDistinct++;
		}
	    }
	}
    }

    if ($n > 0) {
	$mean = sprintf("%.2f",$mean/$n);
	my @sortedDistinctCounts = sort {$a <=> $b} keys(%counts);
	my $midPoint = sprintf("%d",$n/2);
	my $soFar = 0;
	my $countIndex = 0;
	while ($soFar < $midPoint && $countIndex < $nDistinct-1) {
	    $soFar += $counts{$sortedDistinctCounts[$countIndex]};
	    $countIndex++;
	}
	my $median = $sortedDistinctCounts[$countIndex];
	print "$name\t$n\t$mean\t$median\n";
    }
}

    
