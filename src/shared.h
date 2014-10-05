#ifndef SHARED_H
#define SHARED_H

#include <ctype.h>
#include <math.h>

#include "bam.h"

#include "sequences.h"
#include "hashtable.h"

// all random stuff to be shared between indelminerSR and indelminerPE goes in
// here

// wrapper to transfer data back and forth when looking at coverage using pileup
typedef struct covdata_st
{
    uint32_t tid;
    uint32_t begin;
    uint32_t end;    
    uint32_t* coverage;
}covdata;

// a simple structure to hold the sum of squares of MQ, number of reads with
// MQ=0 and the number of reads covering the location
typedef struct mapqcov_st
{
    uint64_t mqsqsum;
    uint32_t numcov;
    uint32_t numcov0;
    uint coefficent;
}mapqcov;

// read the configuration file which has the information about the BAM file and
// the read group insert length ranges, and average coverages
void read_configuration(const char* const filename,
                        hashtable* const id2chroms,
                        hashtable* insertlengths,
                        uint* const averagecoverage);

// read this reference sequence and return a array of pointers to them
char** read_reference(const char* const reference_fasta_name,
                      const bam_header_t* const hin,
                      const char* const chromosomes);

// print the preamble for the VCF file
void print_vcf_preamble(const float majorversion);

// return the average coverage in this interval
uint calculate_cov_params(const char* const bam_name,
                          const int32_t tid,
                          const int32_t start,
                          const int32_t stop);

void calculate_mq_params(bamFile* const pfp,
                         bam_index_t* const fp_index,
                         const int32_t tid, 
                         const int32_t start,
                         const int32_t stop, 
                         uint* const mqrms, 
                         uint* const mq0);
#endif
