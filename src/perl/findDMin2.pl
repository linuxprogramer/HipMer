#!/usr/bin/env perl

use strict;
use warnings;

# Purpose:  Finds the minimal bin in a histogram.  Accepts as input histograms in these two formats:  
#     a) two columns: bin and count
#     b) 11 columns:  standard histogram2.pl output where column 2 is the bin and column 7 is the count in that bin
# Input: histogram output of mercounts

if (! @ARGV || @ARGV > 2) {
    print STDERR "Usage:  $0 <histogram> [peaks_to_skip] \n";
    print STDERR "Purpose:  Finds the minimal bin in a histogram, optionally skipping first n peaks.  
                  Accepts as input histograms in either the histogram2.pl or a simple two-column fomrat \n";
    exit;
}


open( F, $ARGV[0] ) || die "couldn't open $ARGV[0]\n";

my $prevCnt = -1;
my $peakCnt = -1;
my $prevPeakCnt = -1;
my $peakM = -1;
my $minPeakM = 5;  # don't look for peak inside of this point
my $bMinimaFound = 0;
my $prevM = -1;
my $peaksToSkip = $ARGV[1] || 0 ;
my $peaksSkipped = 0;


#determine if histogram is in the histogram2.pl format or just 2 columns
my $fields = <F>;
my ($bin, $count);

my @fields = split(/\s+/, $fields);
if (@fields  == 2) { 
    $bin =0;
    $count =1;
}
elsif (@fields == 9) {

    $bin =1;
    $count =6;
}
else {
    die "unrecognized input histogram format\n";
}


# find the peak
while( <F> )
{
	if ( $_ !~ /^\#/ )
	{
	    my @F = split;
	    my ( $m, $cnt ) = ( $F[$bin], $F[$count] );
	 
	    # once the minima has been found, then start tracking the peakCnt
	    if ( $bMinimaFound && ( $cnt > $peakCnt && $m > $minPeakM  ) )
	    { 
		$peakCnt = $cnt;
		$peakM = $m;
	    }
	    # if we've reached the bottom, set the minima
	    if ( !$bMinimaFound && $cnt > $prevCnt && $prevCnt != -1 ) 
	    { 
		$bMinimaFound = 1; 
		$minPeakM = $prevM; 
	    }
	    # if skipping peaks, reset the min values when we've gone over a peak until we skip the specified number
	    if ( ($peaksSkipped < $peaksToSkip) && $bMinimaFound && ( $cnt < $peakCnt && $m > $minPeakM  ))
	    {
		$prevPeakCnt = $peakCnt;
		$peakCnt = -1; 
		$peakM = -1;
		$bMinimaFound = 0;
		$peaksSkipped++;
	    }
	    $prevCnt = $cnt; $prevM = $m;
	}
}
if ( !$bMinimaFound || $peakM > 1000 ) { die "no minima found / peakM suspiciously high at $peakM\n" }

# when dealing with multiple peaks, we're suspicious about peaks that are too differnt in size
if ( $peaksToSkip && $bMinimaFound &&  ($peakCnt < ($prevPeakCnt / 10) || $peakCnt > ($prevPeakCnt * 10))) { die "the peaks around the local minimum ($minPeakM) are too different in amplitude ($peakCnt, $prevPeakCnt)\n" }

# If there are multiple peaks, report dmin as the minimum after skipping the desired n peaks
if ($peaksToSkip)
{
    print STDERR "Skipped first $peaksSkipped peaks\n";
    print STDERR "New local peak at $peakM\n";
    print STDERR "New local minimum at $minPeakM\n";
    print STDERR "D-min cutoff picked at: $minPeakM\n";
    print $minPeakM;  # d-min found here
    exit 0;
}
else 
{
    print STDERR "Peak: $peakM\n";
    print STDERR "Minimum: $minPeakM\n";

}

# Otherwise pick dmin based on count at dmax
# look only from [ 0, peak ]
sysseek( F, 0, 0 );
while( <F> )
{
	if ( $_ !~ /^\#/ )
	{
		my @F = split;
		my ( $m, $cnt ) = ( $F[$bin], $F[$count] );

		# d-min should be less than the trough
		if ( $m > $minPeakM ) 
                {
                   die "Couldn't find suitable d-min"; 
                }
		if ( $cnt < $peakCnt )
		{
		    print STDERR "D-min cutoff picked as: $m\n";
		    print $m;  # d-min found here
		    exit 0;
		}
	}
}
die "d-min could not be chosen";
