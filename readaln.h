#ifndef READALN_H
#define READALN_H

#include <inttypes.h>

#include "bam.h"

#include "slinklist.h"

#define BAM_CPMATCH 7
#define BAM_CMMATCH 8

typedef struct readseg_st
{
    struct readseg_st* next;
    char* sequence;
    uint32_t oplen:28, op:4;
    int32_t start; // 0-based leftmost coordinate for this segment
    int32_t end;   // 0-based halfopen right most coordinate for this segment
}readseg;
// for a deleted segment the [start,end) and for inserted segment the inserted
// sequence is after the position denoted by start.

typedef struct readaln_st
{
    char* qname;
    int32_t tid;   // ID of the chromosome this read aligns to
    char strand;   // did I have to reverse complement the read?
    char index;    // index 1/2 for the read
    uint8_t qual;  // mapping quality of the read
    readseg* segments;
}readaln;

// convert 4-bit integer representing the base to nucleotide
char bit2char(const int encodedbase);

// create a new read segment from this information
readseg* new_readseg(const char* const sequence,
                     const uint32_t cigar,
                     int* const prefindx,
                     int* const preadindx);

// glean the information from this BAM alignment
readaln* new_readaln(const bam1_t* const alignment);

// create a read from this BAM alignment
readaln* new_unaligned_readaln(const bam1_t* const alignment, 
                               const uint8_t mmq);

// for debugging only: print the read structure
void print_readaln(const readaln* const rln);

// free the resources used by this
void free_readaln(readaln** prln);

// free the resources held by these readsegs
void free_readsegs(readseg** prs);

// modify the segments for this readaln, as per the given cigarstring
void update_readsegs(readaln* const rln, 
                     const int r1,
                     const uint32_t* cigarstring1,
                     const int numcigarops1,
                     const int index,
                     const int q2,
                     const int r2,
                     const uint32_t* cigarstring2,
                     const int numcigarops2);

// duplicate this read segment. Duplicate all members except the next pointer.
readseg* duplicate_readseg(const readseg* const rs);

// combine all these readseg to one [begin, end)
readseg* combine_readsegs(const readseg* const begin, const readseg* const end);
#endif
