#ifndef BAMOPERATIONS_H
#define BAMOPERATIONS_H

#include <math.h>

#include "bam.h"

#include "hashtable.h"

// simple stucture to store the sum of coverage and the number of bases covered
// on a contig/scaffold/chromosome
typedef struct covsum_st
{
    uint64_t covsum;
    uint32_t numcov;
}covsum;

// print the names of the readgroups and the insert ranges
void print_insertlength_ranges(hashtable* const hash);

// calculate the insert length range for all the read groups in this
// file:chromosome(s)
void estimate_insertlengths(const char* const bam_name,
                            hashtable* const insertlengths,
                            const int chromid);

// observed coverage on the chromosome(s)
void estimate_average_coverage(const char* const bam_name,
                               const int chromid,
                               uint* const averagecoverage);

#endif
