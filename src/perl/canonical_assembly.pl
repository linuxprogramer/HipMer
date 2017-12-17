#!/usr/bin/perl

use strict;
use warnings;

sub revcomp {
  my($seq) = @_;
#  return $seq;
  $seq =~ tr/AaCcGgTt/TtGgCcAa/;
  return reverse($seq);
}

my @contigs;
my $fasta = undef;;

my $minLength = 0;
my $USAGE = "$0 [minLength=0] contig1.fa [...]";
if( scalar(@ARGV)>0 && ! -f $ARGV[0] ) {
  $minLength = shift;
}

my $t_start = time;
while (<>) {
  if (/^>/ || /:/) { push (@contigs, $fasta) if (defined $fasta); $fasta = undef; next; }
  chomp;
  $fasta .= $_;
}
push(@contigs, $fasta) if (defined $fasta);
my $elapsed_t = time - $t_start;
#print STDERR "File scanning time: $elapsed_t\n";

$t_start = time;
foreach my $fa (@contigs) {
  my $rc = revcomp($fa);

  if (($fa cmp $rc) < 0) {
  } else { 
    $fa = $rc;
  }
}
$elapsed_t = time - $t_start;
#print STDERR "revcomp time: $elapsed_t\n";

$t_start = time;
my $contig_id = 0;
for my $contig ( sort( { my $x = length($a) <=> length($b); if ($x == 0) { return $a cmp $b; } else { return $x; } } @contigs)) {
   if (length($contig) >= $minLength) {
       print ">$contig_id-" . length($contig) . "\n$contig\n";
       $contig_id++;
   }
}
$elapsed_t = time - $t_start;
#print STDERR "Sort time: $elapsed_t\n";




