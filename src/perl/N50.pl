#!/usr/bin/env perl

use warnings;

$n_args = @ARGV;
if ($n_args != 2)  {
    print "Usage: N50.pl <input file> <0 = fasta, 1 = name/size>\n";
    exit
}

open(F,$ARGV[0]) || die "Couldn't open $ARGV[0]\n";;
if (($ARGV[1] != 0) && ($ARGV[1] != 1)) {
    die "Usage: N50.pl <input file> <0 = fasta, 1 = name/size>\n";
}

my $input_style = $ARGV[1];

my $sequence = "";
my $id = "NO_CURRENT_ID";
my $n_bases = 0;
my $total_bases = 0;
my %name2size = ();

# fasta input
if ($input_style == 0) {

    while (my $i = <F>) {
	chomp $i;
	
	if ($i =~ /^>/) {
	    
	    if ($id ne "NO_CURRENT_ID") {
		$n_bases = length($sequence);
#	    print "$id\t$n_bases\n";
		if (exists($name2size{$id})) {
		    warn "Warning $id multiply defined - results inaccurate\n";
		}
		$name2size{$id} = $n_bases;
		$total_bases += $n_bases;
		
		$sequence = "";
	    }
	    
	    ($id) = $i =~ /^>(\S+)/;
	    
	} else {
	    
	    $sequence .= $i;
	}
    }

    $n_bases = length($sequence);

#print "$id\t$n_bases\n";
    if (exists($name2size{$id})) {
	warn "Warning $id multiply defined - results inaccurate\n";
    }
    $name2size{$id} = $n_bases;
    $total_bases += $n_bases;

} else { 

    while (my $i = <F>) {
	chomp $i;

	($id, $n_bases) = $i =~ /^(\S+)\s+(\d+)/;

#	    print "$id\t$n_bases\n";
	if (exists($name2size{$id})) {
	    warn "Warning $id multiply defined - results inaccurate\n";
	}
	$name2size{$id} = $n_bases;
	$total_bases += $n_bases;
    }
}

close F;

my @size_sorted_names = sort {$name2size{$b} <=> $name2size{$a}} keys(%name2size);

my $running_total = 0;
my $n_seqs = 0;
foreach(@size_sorted_names) {
    my ($name,$size) = ($_,$name2size{$_});
    $running_total += $size;
    $n_seqs++;
    my $fraction = $running_total/$total_bases;
    print "$n_seqs\t$name\t$size\t$running_total\t$fraction\n";
}
