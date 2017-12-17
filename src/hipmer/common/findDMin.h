#ifndef FIND_DMIN_H
#define FIND_DMIN_H

/*
# Purpose:  Finds the minimal bin in a histogram. 
#     a) two columns: bin and count
*/

static int findDMin(int64_t *histo, int size, int peaksToSkip) {
	int prevCnt = -1;
	int revCnt = -1;
	int peakCnt = -1;
	int prevPeakCnt = -1;
	int peakM = -1;
	int minPeakM = 5;  // don't look for peak inside of this point
	int bMinimaFound = 0;
	int prevM = -1;
	int peaksSkipped = 0;
	int bin, m, cnt;

	for(bin = 0; bin < size; bin++) {
		m = bin+1;
		cnt = histo[bin];

		// once the minima has been found, then start tracking the peakCnt
		if ( bMinimaFound > 0 && ( cnt > peakCnt && m > minPeakM  ) ) { 
			peakCnt = cnt;
			peakM = m;
		}

		// if we've reached the bottom, set the minima
		if ( !bMinimaFound && cnt > prevCnt && prevCnt != -1 ) { 
			bMinimaFound = 1; 
			minPeakM = prevM; 
		}

		// if skipping peaks, reset the min values when we've gone over a peak until we skip the specified number
		if ( (peaksSkipped < peaksToSkip) && bMinimaFound && ( cnt < peakCnt && m > minPeakM  )) {
			prevPeakCnt = peakCnt;
			peakCnt = -1; 
			peakM = -1;
			bMinimaFound = 0;
			peaksSkipped++;
		}
		prevCnt = cnt; prevM = m;
	}

	if ( !bMinimaFound || peakM > 1000 ) {
		fprintf(stderr, "no minima found / peakM suspiciously high at %d\n", peakM);
		return -1;
	}

	// when dealing with multiple peaks, we're suspicious about peaks that are too differnt in size
	if ( peaksToSkip && bMinimaFound &&  (peakCnt < (prevPeakCnt / 10) || peakCnt > (prevPeakCnt * 10))) {
		fprintf(stderr, "the peaks around the local minimum (%d) are too different in amplitude (%d, %d)\n", minPeakM, peakCnt, prevPeakCnt);
		return -1;
	} 

	// If there are multiple peaks, report dmin as the minimum after skipping the desired n peaks
	if (peaksToSkip) {
		fprintf(stderr, "Skipped first %d peaks\n", peaksSkipped);
		fprintf(stderr, "New local peak at %d\n", peakM);
		fprintf(stderr, "New local minimum at %d\n", minPeakM);
		fprintf(stderr, "D-min cutoff picked at: %d\n", minPeakM);
		return minPeakM;  // d-min found here
	} else {
		fprintf(stderr, "Peak: %d\n", peakM);
		fprintf(stderr, "Minimum: %d\n", minPeakM);
	}

	// Otherwise pick dmin based on count at dmax
	// look only from [ 0, peak ]
	for(bin = 0; bin < size; bin++) {
		m = bin+1;
		cnt = histo[bin];

		// d-min should be less than the trough
		if ( m > minPeakM ) 
                {
			fprintf(stderr, "Couldn't find suitable d-min"); 
			return -1;
                }
		if ( cnt < peakCnt )
		{
		    fprintf(stderr, "D-min cutoff picked as: %d\n", m);
		    return m;
		}
	}
	fprintf(stderr, "d-min could not be chosen");
	return -1;
}

#endif
