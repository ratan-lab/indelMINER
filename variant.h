#ifndef VARIANT_H
#define VARIANT_H

#include <inttypes.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include "bam.h"

#include "strings.h"
#include "resources.h"
#include "files.h"
#include "slinklist.h"
#include "hashtable.h"
#include "evidence.h"
#include "shared.h"

typedef struct variant_st
{
    struct variant_st* next;
    varianttype type;
    int32_t tid;
    uint start;
    uint stop;
    uint support;
    uint lw;
    uint rw;
    evidence** evidence;
    evidencetype evdnctype;
    bool isused;
} variant;

typedef struct knownvariant_st
{
    struct knownvariant_st* next;
    int32_t tid;
    uint start;
    char* reference;
    char* alternate;
    varianttype type;
    evidencetype evdnctype;
    uint support;
    uint stop;
    uint bpstop;
    char* addntlinfo;
    bool diffsample_support;
    char* sequence;
} knownvariant;

// sort the variants based on the chromosomal coordinates
variant* sort_variants(variant* vs);

// print the variants
void print_variants(variant** const pvariants, 
                    char** const sequences,
                    const char* const bam_name,
                    const bam_header_t* const hin,
                    const char* const outputformat,
                    const uint minsupport,
                    const uint maxdiffsallowed);

// read the variants from the VCF file
knownvariant* read_variants(const char* const vcfname, 
                            char** const sequences,
                            const int tid,
                            const char* const chromname);

// merge the SR variants if the exact position of the SV is unclear since bases
// at the edge of one read-half could equally well be appended to the other
// read-half. Then merge the SR and PE variants if they share the same
// breakpoints
void merge_variants(variant** pvs, 
                    const char* const reference, 
                    const bool join_sr_pe);

// just print the variant
void print_vcf_line(const knownvariant* const kiter,
                    const bam_header_t* const hin);

// print the variants that we already know about along with the tag for this 
// sample (in case this variant was also found in this sample) 
knownvariant* print_knownvariants(knownvariant* const kvariants, 
                                  variant* const variants, 
                                  const char* const sample_name,
                                  const char* const bam_name,
                                  const char* const outputformat,
                                  const bam_header_t* const hin);

// last attempt to check whether an indel is supported in the given BAM file
bool is_indel_supported(knownvariant* const variant, 
                        const char* const bam_name);


#endif
