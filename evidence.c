#include "evidence.h"

// create a new evidence structure from this read
evidence* new_evidence(readaln* const rln, 
                       const void* const variantsegment,
                       const varianttype vclass,
                       const evidencetype type) 
{
    evidence* evdnc = ckallocz(sizeof(evidence));

    evdnc->type = type;
    evdnc->variantclass = vclass;
    evdnc->qual = rln->qual;
    evdnc->strand = rln->strand;
    evdnc->isused = FALSE;    
    evdnc->qname = ckallocz(strlen(rln->qname) + 1);
    memcpy(evdnc->qname, rln->qname, strlen(rln->qname));

    if(type == SPLIT_READ){
        const readseg* varsegment = variantsegment;
        evdnc->aln1 = combine_readsegs(rln->segments, varsegment);
        evdnc->aln2 = duplicate_readseg(varsegment);
        evdnc->aln3 = combine_readsegs(varsegment->next, NULL);
        evdnc->b1 = varsegment->start;
        evdnc->b2 = varsegment->end;
    }else if(type == PAIRED_READ){
        forceassert(variantsegment == NULL);
        evdnc->aln1 = rln->segments;
    }else{
        fatalf("unknown type of evidence");
    }

    return evdnc;
}


// string representation of the evidence
void print_evidence(const evidence* const evdnc)
{
	if(evdnc->type == SPLIT_READ){
		printf("SRread evidence : %d : %d\n", evdnc->b1, evdnc->b2);
	}else if(evdnc->type == PAIRED_READ){
		printf("PEread evidence : %d : %d\n", evdnc->b1, evdnc->b2);
	}else{
		fatalf("unknown evidence type\n");
	}
}

// sort the evidence based on the chromosomal coordinates of the variation
int sort_evidence(const void* const el1, const void* const el2)
{
    evidence* e1 = *((evidence**)el1);
    evidence* e2 = *((evidence**)el2);

    if(e1->b1 == e2->b1) return e1->b2 - e2->b2;

    return e1->b1 - e2->b1;
}

static void free_evidence_block(evidence* evdnc)
{
    ckfree(evdnc->qname);
    readseg* ptr = evdnc->aln1;
    free_readsegs(&ptr);
    ptr = evdnc->aln2;
    free_readsegs(&ptr);
    ptr = evdnc->aln3;
    free_readsegs(&ptr);
}

// free all resources used by evidence that is used in some variation and return
// the list with evidence that has not been used
evidence* free_used_evidence(evidence* allevidence)
{
    evidence* remaining = NULL;

    evidence* iter = allevidence;
    evidence* tmp;
    while(iter){
        tmp = iter->next;

        if(iter->isused == TRUE){
            free_evidence_block(iter);
            ckfree(iter);
        }else{
            sladdhead(&remaining, iter);
        }
        iter = tmp;
    }

    slreverse(&remaining);
    return remaining;
}
