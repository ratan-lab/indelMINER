#ifndef EVIDENCE_H
#define EVIDENCE_H
    
#include "slinklist.h"
#include "readaln.h"

typedef enum evidencetype_st 
{
    SPLIT_READ,
    PAIRED_READ,
    COMPOSITE
}evidencetype;

typedef enum varianttype 
{
    INSERTION,
    DELETION
}varianttype;

typedef struct evidence_st
{
    struct evidence_st* next;
    evidencetype type;
    varianttype variantclass;
    char strand;
    uint8_t qual;
    char* qname;
    readseg* aln1;
    readseg* aln2;
    readseg* aln3;
    int32_t b1;
    int32_t b2;
    int32_t mindelsize; // only filled for PE evidence
    int32_t max; // only filled for PE evidnce
    bool isused;
}evidence;

// create a new evidence structure from this read
evidence* new_evidence(readaln* const rln, 
                       const void* const variantsegment, 
                       const varianttype vclass,
                       const evidencetype type); 

// string representation of the evidence
void print_evidence(const evidence* const evdnc);

// sort the evidence based on the chromosomal coordinates of the variation
int sort_evidence(const void* const el1, const void* const el2);

// free all resources used by evidence that is used in some variation
evidence* free_used_evidence(evidence* allevidence);
#endif
