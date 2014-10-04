#include "localalign.h"

#define MININT -9999999

//const int GAPOPEN   = 4;
//const int GAPEXTEND = 1;
//const int MATCH     = 2;
//const int MISMATCH  = -1;

const int GAPOPEN   = 10;
const int GAPEXTEND = 10;
const int MATCH     = 1;
const int MISMATCH  = -10;

int local_align(char* seq1,
                const int seq1len,
                char* seq2,
                const int seq2len,
                const int indx1,
                const int indx2,
                int* const psi,
                int* const psj,
                int* const pei,
                int* const pej,
                int* const S)
{
    // what am I doing here?
//    fprintf(stderr, "Seq1 : %s\n", seq1);
//    fprintf(stderr, "Seq2 : %s\n", seq2);
//    fprintf(stderr, "Diagonals : %d - %d\n", indx1, indx2);
    forceassert(strlen(seq1) > 0);
    forceassert(strlen(seq2) > 0);

// if A,B are null-terminated strings; then
// A is seq1 - 1
// B is seq2 - 1
// M is strlen(A+1)
// N is strlen(B+1)
// low and up are indices


    char* A  = seq1 - 1;
    char* B  = seq2 - 1;
    int  M  = seq1len;
    int  N  = seq2len;
    int low = indx1;
    int up  = indx2;
    int  G  = GAPOPEN;
    int  H  = GAPEXTEND;

    // Webb's code   	
    int band;
	int i, j, si, ei;
	int c, d, e, t, m;
	int leftd, rightd;
	int best_score, starti = 0, startj = 0, endi = 0, endj = 0;
	int *wa, curd;
	int ib;
	char flag;

    int W[128][128];
    for (i = 0; i < 128; i++)
        for (j = 0; j < 128; j++)
            if (i == j)
                W[i][j] = MATCH;
            else
                W[i][j] = MISMATCH;

	m = G+H;
	low = MAX(-M, low);
	up = MIN(N, up);

	band = up-low+1;
	if (band < 1) {
		printf("low > up is unacceptable! low:%d up:%d\n", low, up);
		exit(1);
	}

    int* CC = ckalloc((band+2) * sizeof(int));
    int* DD = ckalloc((band+2) * sizeof(int));

	if (low > 0) leftd = 1;
	else if (up < 0) leftd = band;
	else leftd = 1-low;
	rightd = band;
	si = MAX(0,-up);
	ei = MIN(M,N-low);
	CC[leftd] = 0;
	for (j = leftd+1; j <= rightd; j++) {
		CC[j] = 0;
		DD[j] = -G;
	}
	CC[rightd+1] = MININT;
	DD[rightd+1] = MININT;
	best_score = 0;
	endi = si;
	endj = si+low;
	CC[leftd-1] = MININT;
	DD[leftd] = -G;
	for (i = si+1; i <= ei; i++) {
	    if (i > N-up) rightd--;
	    if (leftd > 1) leftd--;
	    wa = W[(int)A[i]];
	    if ((c = CC[leftd+1]-m) > (d = DD[leftd+1]-H)) d = c;
	    if ((ib = leftd+low-1+i) > 0) c = CC[leftd]+wa[(int)B[ib]];
	    if (d > c) c = d;
	    if (c < 0) c = 0;
	    e = c-G;
	    DD[leftd] = d;
	    CC[leftd] = c;
	    if (c > best_score) {
		    best_score = c;
		    endi = i;
		    endj = ib;
	    }
	    for (curd=leftd+1; curd <= rightd; curd++) {
		if ((c = c-m) > (e = e-H)) e = c;
		if ((c = CC[curd+1]-m) > (d = DD[curd+1]-H)) d = c;
		c = CC[curd] + wa[(int)B[curd+low-1+i]];
		if (e > c) c = e;
		if (d > c) c = d;
		if (c < 0) c = 0;
		CC[curd] = c;
		DD[curd] = d;
		if (c > best_score) {
		   best_score = c;
		   endi = i;
		   endj = curd+low-1+i;
		}
	    }
	}
	leftd = MAX(1,-endi-low+1);
	rightd = band-(up-(endj-endi));
	CC[rightd] = 0;
	t = -G;
	for (j = rightd-1; j >= leftd; j--) {
		CC[j] = t = t-H;
		DD[j] = t-G;
	}
	for (j = rightd+1; j <= band; ++j) CC[j] = MININT;
	CC[leftd-1] = DD[leftd-1] = MININT;
	DD[rightd] = -G;
	flag = 0;
	for (i = endi; i >= 1; i--) {
	    if (i+low <= 0) leftd++;
	    if (rightd < band) rightd++;
	    wa = W[(int)A[i]];
	    if ((c = CC[rightd-1]-m) > (d = DD[rightd-1]-H)) d = c;
	    if ((ib = rightd+low-1+i) <= N) c = CC[rightd]+wa[(int)B[ib]];
	    if (d > c) c = d;
	    e = c-G;
	    DD[rightd] = d;
	    CC[rightd] = c;
	    if (c == best_score) {
		starti = i;
	 	startj = ib;
	  	flag = 1;
		break;
	    }
	    for (curd=rightd-1; curd >= leftd; curd--) {
		if ((c = c-m) > (e = e-H)) e = c;
		if ((c = CC[curd-1]-m) > (d = DD[curd-1]-H)) d = c;
		c = CC[curd] + wa[(int)B[curd+low-1+i]];
		if (e > c) c = e;
		if (d > c) c = d;
		CC[curd] = c;
		DD[curd] = d;
		if (c == best_score) {
		   starti = i;
		   startj = curd+low-1+i;
		   flag = 1;
		   break;
		}
	    }
	    if (flag == 1) break;
	}

	free(CC);
	free(DD);
	if (starti < 0 || starti > M || startj < 0 || startj > N) {
        return 0;
//		printf("starti=%d, startj=%d\n",starti,startj);
//		*psi = *psj = *pei = *pej;
//		exit(1);
	}
	*psi = starti;
	*psj = startj;
	*pei = endi;
	*pej = endj;
    
    if(((endi - starti) == 0) || ((endj - startj) == 0)){
        return 0;
    }

  	return ALIGN(A+starti-1,B+startj-1,endi-starti+1,endj-startj+1,low-(startj-starti),up-(startj-starti),W,G,H,S);
}
