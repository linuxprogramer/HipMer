#include <math.h>
#include <assert.h>
#include <string.h>

unsigned int packedCodeToFourMer[256];

void init_LookupTable()
{
	// Work with 4-mers for the moment to have small lookup tables
	int merLen = 4, i, slot, valInSlot;
	unsigned char mer[4];

	for ( i = 0; i < 256; i++ ) {
		// convert a packedcode to a 4-mer
		int remainder = i;
		int pos = 0;
		for( slot = merLen-1; slot >= 0; slot-- ) {
			valInSlot =  (int)  (remainder / pow( 4, slot ) );
			char base;
			
			if( valInSlot == 0 ) { base = 'A'; }
			else if( valInSlot == 1 ) { base = 'C'; }
			else if( valInSlot == 2 ) { base = 'G'; }
			else if( valInSlot == 3 ) { base = 'T'; }
			else{ assert( 0 ); }
			
			mer[pos] = base;
			pos++;
			remainder -= valInSlot * pow( 4, slot );
		}
		unsigned int *merAsUInt = (unsigned int*) mer;
		packedCodeToFourMer[i] = (unsigned int) (*merAsUInt);
	}
}

unsigned char convertFourMerToPackedCode(unsigned char *fourMer)
{
	int retval = 0;
	int code, i;
	int pow = 64;

	for ( i=0; i < 4; i++) {
		char base = fourMer[i];
		switch ( base ) {
			case 'A':
				code = 0;	
				break;
			case 'C':
				code = 1;
				break;
			case 'G':
				code = 2;
				break;
			case 'T':
				code = 3;
				break;
		}
		retval += code * pow;
		pow /= 4;
	}
	return ((unsigned char) retval);
}

unsigned char convertExtensionsToPackedCode(unsigned char *ext)
{
	/* Packs the k-mer extensions to a single byte */
	int retval = 0, r_code, l_code;
	char left = ext[0];
	char right = ext[1];
	switch ( right ) {
		case 'A':
			r_code = 0;	
			break;
		case 'C':
			r_code = 1;
			break;
		case 'G':
			r_code = 2;
			break;
		case 'T':
			r_code = 3;
			break;
		case 'F':
			r_code = 4;
			break;
		case 'X':
			r_code = 5;
			break;
	}
	
	switch ( left ) {
		case 'A':
			l_code = 0;	
			break;
		case 'C':
			l_code = 1;
			break;
		case 'G':
			l_code = 2;
			break;
		case 'T':
			l_code = 3;
			break;
		case 'F':
			l_code = 4;
			break;
		case 'X':
			l_code = 5;
			break;
	}
	return ((unsigned char) (r_code+l_code*6));
}

void convertPackedCodeToExtension(unsigned char packed_ext, char *left_ext, char *right_ext)
{
	/* Unpacks the k-mer extensions */	
	int l_code = packed_ext / 6;
	int r_code = packed_ext % 6;
	
	switch ( r_code ) {
		case 0 :
			*right_ext = 'A';	
			break;
		case 1:
			*right_ext = 'C';
			break;
		case 2:
			*right_ext = 'G';
			break;
		case 3:
			*right_ext = 'T';
			break;
		case 4:
			*right_ext = 'F';
			break;
		case 5:
			*right_ext = 'X';
			break;
	}
	
	switch ( l_code ) {
		case 0 :
			*left_ext = 'A';	
			break;
		case 1:
			*left_ext = 'C';
			break;
		case 2:
			*left_ext = 'G';
			break;
		case 3:
			*left_ext = 'T';
			break;
		case 4:
			*left_ext = 'F';
			break;
		case 5:
			*left_ext = 'X';
			break;
	}
}

void packSequence(const unsigned char *seq_to_pack, unsigned char *m_data, int m_len)
{
	/* The pointer to m_data points to the result of the packing */
	
	int ind, j = 0;     // coordinate along unpacked string ( matches with m_len )
    int i = 0;     // coordinate along packed string
		
    // do the leading seq in blocks of 4
    for ( ; j <= m_len - 4; i++, j+=4 ) {
		m_data[i] = convertFourMerToPackedCode( ( unsigned char * ) ( seq_to_pack + j )) ;
    }
    
    // last block is special case if m_len % 4 != 0: append "A"s as filler        
    int remainder = m_len % 4;
    unsigned char blockSeq[4] = "AAA\0";
    for(ind = 0; ind < remainder; ind++) {
		blockSeq[ind] = seq_to_pack[j + ind];
    }
    m_data[i] = convertFourMerToPackedCode(blockSeq);
}

void unpackSequence(const unsigned char *seq_to_unpack, unsigned char *unpacked_seq, int kmer_len)
{
	/* Result string is pointer unpacked_seq */
	int i = 0, j = 0;
	int packed_len = ceil(kmer_len/4.0);
	for( ; i < packed_len ; i++, j += 4 )
	{
		*( ( unsigned int * )( unpacked_seq + j ) ) = packedCodeToFourMer[ seq_to_unpack[i] ];
	}	
	*(unpacked_seq + kmer_len) = '\0';
}

inline int comparePackedSeq(const unsigned char *seq1, const unsigned char *seq2, int seq_len)
{
	return memcmp(seq1, seq2, seq_len);
}

/* TODO: Add function that computes hash using the packed form */

