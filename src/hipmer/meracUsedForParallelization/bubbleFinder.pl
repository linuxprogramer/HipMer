#!/usr/common/usg/languages/perl/5.16.0/bin/perl -w
#
# bubbleFinder.pl by Jarrod Chapman <jchapman@lbl.gov> Sat Jan 29 20:12:18 PST 2011
# Copyright 2011 Jarrod Chapman. All rights reserved.
#

use Getopt::Std;
my %opts = ();
my $validLine = getopts('c:f:d:', \%opts);
my @required = ("c","f","d");
my $nRequired = 0;
map {$nRequired += exists($opts{$_})} @required;
$validLine &= ($nRequired == @required);
if ($validLine != 1) {
    print "Usage: ./bubbleFinder.pl <-c ceaFile> <-f fastaFile> <-d depthFile>\n";
    exit;
}

my $ceaFile = $opts{"c"};
my $fastaFile = $opts{"f"};
my $cmdFile = $opts{"d"};

my $merLength = undef;
my $hairLength = undef;

my %contigInfo = ();  # map from contig to cea info
my %tipInfo = ();     # map from mer-tip pair to all contigs they flank

my %bubbleMap = ();  # map from contigs in bubbles to their bubble IDs 
my %bubbletigs = (); # map from bubble ID to mer-tips defining the bubble 
my %linkertigs = (); # map from non-bubble contigs to their mer-tips


print STDERR "Reading contig end data from $ceaFile...\n";

open (F,$ceaFile) || die "Couldn't open $ceaFile\n";
while (my $line = <F>) {
    chomp $line;
    my ($contig,$prevCode,$prevBase,$firstMer,$contigLength,$lastMer,$nextBase,$nextCode) = split(/\s+/,$line);
    unless (defined($merLength)) {
	$merLength = length($firstMer);
	$hairLength = 2*$merLength - 1;
    }

    ($prevCode) = $prevCode =~ /\[(..)\]/;
    ($nextCode) = $nextCode =~ /\[(..)\]/;
    ($prevBase) = $prevBase =~ /\((.)\)/;
    ($nextBase) = $nextBase =~ /\((.)\)/;

    my $prevMer = $prevBase . $firstMer;
    $prevMer = substr($prevMer,0,$merLength);

    my $nextMer = $lastMer . $nextBase; 
    $nextMer = substr($nextMer,1,$merLength);
    
    # hair and sub-hair must be FU-UF or UF-FU

    if ($contigLength <= $hairLength) {
	unless ( ( ($prevCode =~ /F[ACGT]/) && ($nextCode =~ /[ACGT]F/) ) || ( ($prevCode =~ /[ACGT]F/) && ($nextCode =~ /F[ACGT]/) ) ) {
	    next;
	}
    } 

    my $statusField = ($contigLength > $hairLength) ? 1 : 0;
    $contigInfo{$contig} = "$prevCode:$prevBase:$firstMer:$contigLength:$lastMer:$nextBase:$nextCode:$statusField";

    if (($prevCode =~ /F[ACGT]/) || ($nextCode =~ /[ACGT]F/)) {

	my $extensions = "";
	if ($prevCode =~ /F[ACGT]/) {
	    $extensions = "$prevMer.";
	} else {
# I think this is a bug
#	    $extensions = "0";    
# should be
	    $extensions = "0.";    
	}
	if ($nextCode =~ /[ACGT]F/) {
	    $extensions .= "$nextMer";
	} else {
	    $extensions .= "0";
	}

	$linkertigs{$contig} = $extensions;

   } elsif (($prevCode =~ /[ACGT]F/) && ($nextCode =~ /F[ACGT]/)) {

	my $tips = $prevMer . "." . $nextMer;
	my $rcTips = reverse($tips);
	$rcTips =~ tr/ACGT/TGCA/;

	my $rcFlag = "+";
	my $tipKey = $tips;
	if ($rcTips lt $tips) {
	    $tipKey = $rcTips;
	    $rcFlag = "-";
	}
    
	if (exists($tipInfo{$tipKey})) {
	    $tipInfo{$tipKey} .= ",$rcFlag$contig";
	} else {
	    $tipInfo{$tipKey} = "$rcFlag$contig";
	}

	$bubbleMap{$contig} = 0;

    } else {

    }

}
close F;

print STDERR "Done.\n";

# store contig sequence for later use 

print STDERR "Reading contig sequence from $fastaFile...\n";

my %contigSequences = ();
my $currentEntry = undef;
my $currentSeq = "";
open (F,$fastaFile) || die "Couldn't open $fastaFile\n";
while (my $line = <F>) {
    chomp $line;
    if ($line =~ /^>/) {
        if (defined($currentEntry) && exists($contigInfo{$currentEntry})) {
	    $contigSequences{$currentEntry} = $currentSeq;
        }
        $currentSeq = "";
        ($currentEntry) = $line =~ /^>(.+)$/;
    } else {
        $currentSeq .= $line;
    }
}
if (defined($currentEntry) && exists($contigInfo{$currentEntry})) {
    $contigSequences{$currentEntry} = $currentSeq;
}
close F;

print STDERR "Done.\n";

# store contig mean depth for later use 

print STDERR "Reading contig depth info from $cmdFile...\n";

open (F,$cmdFile) || die "Couldn't open $cmdFile\n";
while (my $line = <F>) {
    chomp $line;

    my ($contigID,$nMers,$meanDepth) = split(/\s+/,$line);

    if (exists($contigInfo{$contigID})) {
	$contigInfo{$contigID} .= ":$nMers:$meanDepth";
    }
}
close F;

print STDERR "Done.\n";


my $bubbleID = 0;

# find bubbles from convergent tips

print STDERR "Finding bubbles...\n";

while (my ($tipKey,$contigs) = each(%tipInfo)) {
    my @contigs = split(/\,/,$contigs);
    my $nContigs = scalar(@contigs);

    if ($nContigs > 1) {        # a bubble

	$bubbleID++;
	$bubbletigs{$bubbleID} = $tipKey;

	print STDERR "BUBBLE $bubbleID $tipKey [@contigs]\n";

	foreach my $contigID (@contigs) {
	    my ($rcFlag,$contig) = $contigID =~ /^([\+\-])(.+)/;
	   
	    unless(exists($bubbleMap{$contig})) {
		die "Inconsistency in bubbleMap for contig $contig.\n";
	    } else {
		$bubbleMap{$contig} = $bubbleID;
	    }
	}
    }
	
}

print STDERR "Done.\n";

# Build the bubble-contig graph (%pointsTo contains the edge info)

print STDERR "Building bubble-contig graph...\n";

my %pointsTo = ();

while (my ($bubbleID,$links) = each(%bubbletigs)) {
    my ($in,$out) = $links =~ /^(.+)\.(.+)$/;
    my $rcLinks = reverse($links);
    $rcLinks =~ tr/ACGT/TGCA/;
    my ($rcin,$rcout) = $rcLinks =~ /^(.+)\.(.+)$/;

    if (exists($pointsTo{$in})) {
	$pointsTo{$in} .= ",B+$bubbleID";
    } else {
	$pointsTo{$in} = "B+$bubbleID";
    }
    
    if (exists($pointsTo{$rcin})) {
	$pointsTo{$rcin} .= ",B-$bubbleID";
    } else {
	$pointsTo{$rcin} = "B-$bubbleID";
    }
}

while (my ($contigID,$links) = each(%linkertigs)) {
    my ($in,$out) = $links =~ /^(.+)\.(.+)$/;
    my $rcLinks = reverse($links);
    $rcLinks =~ tr/ACGT/TGCA/;
    my ($rcin,$rcout) = $rcLinks =~ /^(.+)\.(.+)$/;

    if ($in) {
	if (exists($pointsTo{$in})) {
	    $pointsTo{$in} .= ",C+$contigID";
	} else {
	    $pointsTo{$in} = "C+$contigID";
	}
    }
    
    if ($rcin) {
	if (exists($pointsTo{$rcin})) {
	    $pointsTo{$rcin} .= ",C-$contigID";
	} else {
	    $pointsTo{$rcin} = "C-$contigID";
	}
    }
}

print STDERR "Done.\n";

# traverse the bubble-contig graph, starting from each contig

print STDERR "Generating diplotigs...\n";

my %visited = ();
my $diplotigID = 0;

while (my ($contigID,$links) = each(%linkertigs)) {

    my $current = "C+$contigID";
    my $currentID = "C$contigID";

    if (exists($visited{$currentID})) {
	next;
    } else {
	$visited{$currentID} = 1;
    }

    my @downstream = ();
    my $next = getNext($current);

    while ($next) {

	my ($type,$strand,$id) = $next =~ /^([CB])([\+\-])(.+)$/;

	my $nextID = $type . $id;

	if (exists($visited{$nextID})) {
	    warn "warning: Inconsistency encountered in forward walk from $contigID at $nextID.  Truncating.\n";
	    warn "warning: [$type][$strand][$id] [$nextID] [$visited{$nextID}]\n";
	    $next = 0;
	} else {
	    $visited{$nextID} = 1;
	    push (@downstream,$next);
	    $current = $next;
	    $next = getNext($current);
	}
    }

    $current = "C-$contigID";

    my @upstream = ();
    $next = getNext($current);

    while ($next) {

	my ($type,$strand,$id) = $next =~ /^([CB])([\+\-])(.+)$/;

	my $nextID = $type . $id;
	my $rStrand = "-";
	if ($strand eq "-") {
	    $rStrand = "+";
	}

	if (exists($visited{$nextID})) {
	    warn "warning: Inconsistency encountered in reverse walk from $contigID at $nextID.  Truncating.\n";
	    warn "warning: [$type][$strand][$id] [$nextID] [$visited{$nextID}]\n";
	    $next = 0;
	} else {
	    $visited{$nextID} = 1;
	    unshift(@upstream,"$type$rStrand$id");
	    $current = $next;
	    $next = getNext($current);
	}
    }

    # remove leading and trailing bubbles

    my $nUpstream = scalar(@upstream);
    if ($nUpstream > 0) {
	my $first = $upstream[0];
	my ($type,$strand,$id) = $first =~ /^([CB])([\+\-])(.+)$/;
	if ($type eq "B") {
	    shift(@upstream);
	}
    }
    my $nDownstream = scalar(@downstream);
    if ($nDownstream > 0) {
	my $last = $downstream[-1];
	my ($type,$strand,$id) = $last =~ /^([CB])([\+\-])(.+)$/;
	if ($type eq "B") {
	    pop(@downstream);
	}
    }

    my $nTigs = 1 + scalar(@upstream) + scalar(@downstream);

    my $printMode = "maxDepth";  # maxLength, minLength, maxDepth, or all

    if ($nTigs > 1) {       # a bubble-contig string

	my $finalSequence = "";
	$diplotigID++;
	my $nContigMers = 0;
	my $meanContigDepth = 0;
	my $nBubbles = 0; 
	my $nContigs = 0;

	my $warningString = "";

	foreach my $tig (@upstream,"C+$contigID",@downstream) {

	    print STDERR "$tig->";

	    my ($type,$strand,$id) = $tig =~ /^([CB])([\+\-])(.+)$/;

	    if ($type eq "B") {

		$nBubbles++;
		my $links = $bubbletigs{$id};
		my $btigs = $tipInfo{$links};
		my @btigs = split(/\,/,$btigs);

		unless (scalar(@btigs) == 2) {
		    $warningString .= "warning: Only two-path bubbles are currently supported. ($tig = [@btigs])\n";
		}
		
		my %btigSeqs = ();
		my %btigDepths = ();

		foreach my $btig (@btigs) {
		    my ($s,$cid) = $btig =~ /^([\+\-])(.+)$/;

		    labelContig($cid,2);
		    
		    my $depth = getContigInfo($cid,9);
		    $btigDepths{$btig} = $depth;

		    my $seq = $contigSequences{$cid};
		    if ($s ne $strand) {
			my $rcSeq = reverse($seq);
			$rcSeq =~ tr/ACGT/TGCA/;
			$seq = $rcSeq;
		    }
		    $btigSeqs{$btig} = $seq;
		}
		
		my @sortedBtigs = sort {length($btigSeqs{$a}) <=> length($btigSeqs{$b})} @btigs;
		my @btigLens = map {length($btigSeqs{$_})} @sortedBtigs;

		my $minLen = $btigLens[0];
		my $minClipLen = $minLen - (2*$merLength - 4);
		my $tailSeq = "";
		my @paths = ();

		if ($minClipLen < 0) {
		    $tailSeq = substr($finalSequence,$minClipLen);
		    for (my $n = 0; $n < -$minClipLen; $n++) {
			chop $finalSequence;
		    }
		}

		my $maxDepthIndex = 0;
		my $maxDepth = 0;
		my $btigIndex = 0;

		foreach my $btig (@sortedBtigs) {

		    my $seq = $btigSeqs{$btig};
		    my $btigLen = length($seq);
		    my $clipLen = $btigLen - (2*$merLength - 4);
		    if ($clipLen < 0) {
			my $diff = $clipLen - $minClipLen;
			my $addBack = substr($tailSeq,0,$diff);
			push(@paths,$addBack);
		    } else {
			my $clipSeq = substr($seq,$merLength-2,$clipLen);
			push(@paths,$tailSeq.$clipSeq);
		    }

		    my $depth = $btigDepths{$btig};
		    if ($depth >= $maxDepth) {
			$maxDepth = $depth;
			$maxDepthIndex = $btigIndex;
		    }
		    $btigIndex++;

		}

		if ($printMode eq "maxLength") {
		    
		    $finalSequence .= $paths[-1];

		} elsif ($printMode eq "minLength") {

		    $finalSequence .= $paths[0];

		} elsif ($printMode eq "maxDepth") {
		    
		    $finalSequence .= $paths[$maxDepthIndex];

		} elsif ($printMode eq "all") {
		    
		    my $paths = join(",",@paths);
		    $finalSequence .= "[$paths]";

		}

	    } elsif ($type eq "C") {

		$nContigs++;
		labelContig($id,2);

		my $depth = getContigInfo($id,9);
		my $nMers = getContigInfo($id,8);
		$nContigMers += $nMers;
		$meanContigDepth += $nMers*$depth;
		
		my $seq = $contigSequences{$id};
		if ($strand eq "-") {
		    my $rcSeq = reverse($seq);
		    $rcSeq =~ tr/ACGT/TGCA/;
		    $seq = $rcSeq;
		}

		$finalSequence .= "$seq";

	    }
	}
	print STDERR "\n";
	if ($warningString) {
	    warn $warningString;
	}

	my $name = "diplotig$diplotigID";
	my $length = length($finalSequence);
	$meanContigDepth /= $nContigMers;
	print STDERR "$name\t$length\t$nContigs\t$nBubbles\t$meanContigDepth\n";
	printFasta($name,$finalSequence);
    }
}

my $isotigID = 0;
foreach my $cid (keys(%contigInfo)) {

    my $contigStatus = getContigInfo($cid,7);
    if ($contigStatus == 1) {
	$isotigID++;
	my $seq = $contigSequences{$cid};
	my $name = "isotig$isotigID";
	my $length = length($seq);
	my $depth = getContigInfo($cid,9);	
	print STDERR "$name\t$length\t1\t0\t$depth\n";
	printFasta($name,$seq);
    }
}

print STDERR "Done.\n";

# utilities

sub getContigInfo {
    my ($contig, $field) = @_;

    my $contigInfo = $contigInfo{$contig};
    my @contigInfo = split(/:/,$contigInfo);

    my $nFields = scalar(@contigInfo);
    unless ($field < $nFields) {
	die "getContigInfo ( $contig , $field ) failed.  Not enough fields in $contigInfo.\n";
    }

    my $info = $contigInfo[$field];
    return $info;
}


sub labelContig {
    my ($contig,$label) = @_;

    my $contigInfo = $contigInfo{$contig};
    my @contigInfo = split(/:/,$contigInfo);
    $contigInfo[7] = $label;
    $contigInfo = join(":",@contigInfo);
    $contigInfo{$contig} = $contigInfo;    
}


sub getNext {
    my ($tigID) = @_;

    my ($type,$strand,$id) = $tigID =~ /^([BC])([\+\-])(.+)$/;

    my $links = "";
    if ($type eq "C") {
	unless (exists($linkertigs{$id})) {
	    die "LINKERTIGS NOT DEFINED for [$id]\n";
	}
	$links = $linkertigs{$id};
    } else {
	unless (exists($bubbletigs{$id})) {
	    die "BUBBLETIGS NOT DEFINED for [$id]\n";
	}

	$links = $bubbletigs{$id};
    }
    my ($in,$out) = $links =~ /^(.+)\.(.+)$/;
    my $rcLinks = reverse($links);
    $rcLinks =~ tr/ACGT/TGCA/;
    my ($rcin,$rcout) = $rcLinks =~ /^(.+)\.(.+)$/;

    my $next = 0;
    if ($strand eq "-") {
	if (($rcout) && (exists($pointsTo{$rcout}))) {
	    $next = $pointsTo{$rcout};
	}
    } else {
	if (($out) && exists($pointsTo{$out})) {
	    $next = $pointsTo{$out};
	}
    }

    my @next = split(/\,/,$next);
    if ($#next > 0) {
	warn "warning: getNext finds multiple connections for $tigID (@next). Terminating.\n";
	$next = 0;
    }

    return $next;

}

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
        print "$seqLine\n";
    }
}

