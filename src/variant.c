#include "variant.h"

extern uint maxdelsize;
extern uint maxpedelsize;
extern uint minbalance;
extern bool call_all_indels;

#define SUB 0
#define INS 1
#define DEL 2

// print so many bases around the deletion even
const uint neighborhood = 80;

static int sort_by_position(const void* const el1,
                            const void* const el2)
{
    variant* a = *(variant**)el1;
    variant* b = *(variant**)el2;

    forceassert(a->tid == b->tid);

    if(a->start == b->start) return a->stop - b->stop;
    return a->start - b->start;
}

static int sort_by_knownposition(const void* const el1,
                                 const void* const el2)
{
    knownvariant* a = *(knownvariant**)el1;
    knownvariant* b = *(knownvariant**)el2;

    forceassert(a->tid == b->tid);

    if(a->start == b->start) return a->bpstop - b->bpstop;
    return a->start - b->start;
}

// sort the variants based on the chromosomal coordinates
variant* sort_variants(variant* vs)
{
    slsort(&vs, sort_by_position);
    return vs;
}

knownvariant* sort_knownvariants(knownvariant* vs)
{
    slsort(&vs, sort_by_knownposition);
    return vs;
}

static char* voted_consensus(evidence** const evdnc,
                             const uint numsupport,
                             int* const maxsize)
{
    hashtable* counts = new_hashtable(4);
    uint* intcounts = ckallocz(numsupport * sizeof(uint));
    uint indx = 0;    
    int size = 0;

    uint i;
    for(i = 0; i < numsupport; i++){
        forceassert(evdnc[i]->type == SPLIT_READ);
        readseg* rsg = evdnc[i]->aln2;

        if(rsg->sequence == NULL) continue;

        int inslen = strlen(rsg->sequence);
        if(inslen > size) size = inslen;

        if(lookup_hashtable(counts,rsg->sequence, inslen) == NULL){
            add_hashtable(counts, 
                          rsg->sequence,
                          inslen, 
                          intcounts + indx);
            indx++;
        }

        int* cnt = must_find_hashtable(counts,
                                       rsg->sequence,
                                       inslen);

        *cnt = *cnt + 1;
    }
    *maxsize = size;

    // now pick the sequence that is seen the most number of times
    bin* iter;
    bin* next;
    uint maximumcount = 0;
    char* consensus = NULL;
    
    int j;
    for(j = 0; j < counts->size; j++){
        iter = counts->bins[j];
        while(iter){
            next = iter->next;
            if(*((uint*)iter->val) > maximumcount){
                maximumcount = *((uint*)iter->val); 
                consensus = iter->name;
            }
            iter = next;
        } 
    }
    
    
    char* rtconsensus = ckallocz(strlen(consensus) + 1);
    memcpy(rtconsensus, consensus, strlen(consensus));
    
    ckfree(intcounts);
    free_hashtable(&counts);
    return rtconsensus;   
}

void print_vcf_output(const variant* const variant,
                      char** const sequences,
                      const char* const bam_name UNUSED,
                      const bam_header_t* const hin)
{    
    printf("%s\t%d\t.\t", hin->target_name[variant->tid], variant->start - variant->lw);
    
    int endpos = -1;
    if(variant->type == DELETION){
        // this is the length of the reference nucs if it included the whole
        // breakpoint range
        int reflength = (variant->stop + variant->rw) 
                      - (variant->start - variant->lw - 1);
        forceassert(reflength >= 1);

        // this is the length of the alternate nucs if it included the whole
        // breakpoint range
        int altlength = (variant->start + variant->rw) 
                      - (variant->start - variant->lw - 1);
        forceassert(altlength >= 1);

        endpos = variant->start - variant->lw + reflength - altlength + 1;

        int i = 0;
        for(i = 0; i < (reflength - altlength + 1); i++){
            printf("%c",sequences[variant->tid][variant->start-variant->lw-1+i]);
        }
        printf("\t");

        printf("%c\t",sequences[variant->tid][variant->start - variant->lw - 1]);
    }else if(variant->type == INSERTION){
        int reflength = (variant->stop + variant->rw)
                      - (variant->start - variant->lw - 1);
        forceassert(reflength >= 1);

        int maxinsertsize = 0;
        char* consensus = voted_consensus(variant->evidence,
                                          variant->support,
                                          &maxinsertsize);
        int altlength = reflength + strlen(consensus) 
                      + (variant->stop + variant->rw) - (variant->start);
        forceassert(altlength >= 1);

        endpos = variant->start - variant->lw + 1;

        printf("%c\t",sequences[variant->tid][variant->start - variant->lw - 1]);
        
        printf("%c",sequences[variant->tid][variant->start - variant->lw -1]);
        printf("%s\t", consensus);
    }else{
        fatal("unhandled variant type");
    }

    printf(".\t.\t%s;",variant->type == DELETION ? "DELETION" : "INSERTION");
    if(variant->evdnctype == SPLIT_READ){
        printf("SPLIT_READ;");
    }else if(variant->evdnctype == PAIRED_READ){
        printf("PAIRED_READ;");
    }else if(variant->evdnctype == COMPOSITE){
        printf("COMPOSITE;");
    }else{
        fatalf("unknown evidence type for this variant");
    } 

    forceassert(endpos != -1);
    printf("NS=%u;END=%d;BP_END=%d", 
           variant->support, endpos, variant->stop + variant->rw + 1);

    // how many reads in forward and how many on the reverse strand support this
    // variant
    uint i;
    uint numfevidence = 0, numrevidence = 0;
    for(i = 0; i < variant->support; i++){
        if(variant->evidence[i]->strand == '+'){
            numfevidence += 1;
        }else if(variant->evidence[i]->strand == '-'){
            numrevidence += 1;
        }else{
            fatalf("unknown strand");
        }
    }
    printf(";NFS=%u;NRS=%u", numfevidence, numrevidence);
    forceassert((numfevidence + numrevidence) == variant->support);

    // look at every split read evidence and calculate the number of unique tail
    // distances in the evidence set.Also look at the evidence and find the MQ 
    //of the supporting reads.
    uint maxtaild = 100;
    char* taildistances = ckallocz(maxtaild + 1);
    uint numunique_taildistances = 0;
    uint num_pe_evidence = 0;

    int balance = INT_MAX;
    int lflank = -1, rflank = -1;

    uint mq = 0;
    uint mq30 = 0;
    uint numdiffs = 0;
    for(i = 0; i < variant->support; i++){
        mq += variant->evidence[i]->qual;
        if(variant->evidence[i]->qual >= 30) mq30 += 1;

        uint ltmp = 0;
        readseg* iter;
        for(iter = variant->evidence[i]->aln1; iter; iter = iter->next){
            switch(iter->op){
                case BAM_CMATCH: 
                    ltmp += iter->oplen;
                    uint j;
                    for(j = 0; j < iter->oplen; j++){
                        if(iter->sequence[j] != sequences[variant->tid][iter->start + j]){
                            numdiffs += 1;
                        }
                    }
                    break;
                case BAM_CPMATCH: ltmp += iter->oplen; break;
                case BAM_CMMATCH: 
                    ltmp += iter->oplen; 
                    numdiffs += iter->oplen;
                    break;
                case BAM_CINS: 
                    ltmp += iter->oplen; 
                    numdiffs += iter->oplen;
                    break;
                case BAM_CDEL: 
                    numdiffs += iter->oplen;
                    break;
                case BAM_CSOFT_CLIP: break;
                default: fatalf("unhandled BAM operation");
            }
        }

        uint rtmp = 0;
        for(iter = variant->evidence[i]->aln3; iter; iter = iter->next){
            switch(iter->op){
                case BAM_CMATCH: 
                    rtmp += iter->oplen; 
                    uint j;
                    for(j = 0; j < iter->oplen; j++){
                        if(iter->sequence[j] != sequences[variant->tid][iter->start + j]){
                            numdiffs += 1;
                        }
                    }
                    break;
                case BAM_CPMATCH: rtmp += iter->oplen; break;
                case BAM_CMMATCH: 
                    rtmp += iter->oplen; 
                    numdiffs += iter->oplen;
                    break;
                case BAM_CINS: 
                    rtmp += iter->oplen; 
                    numdiffs += iter->oplen;
                    break;
                case BAM_CDEL: 
                    numdiffs += iter->oplen;
                    break;
                case BAM_CSOFT_CLIP: break;
                default: fatalf("unhandled BAM operation");
            }
        }
        uint taild = MIN(rtmp, ltmp);
        // do I need more memory?
        if(taild > maxtaild){
            taildistances = ckreallocz(taildistances, maxtaild + 1, taild + 1);
            maxtaild = taild;
        }
        taildistances[taild] = '1';
        if(variant->evidence[i]->type == PAIRED_READ){
            num_pe_evidence += 1;
        }

        if(abs(rtmp - ltmp) < balance) {
            balance = abs(rtmp - ltmp);
            lflank = ltmp;
            rflank = rtmp;
        }
    }    

    for(uint i = 0; i < maxtaild; i++){
        if(taildistances[i] == '1') numunique_taildistances++;
    }
    numunique_taildistances += num_pe_evidence;

    printf(";UTAILS=%d;MQ=%d;MQ30=%d;DF=%d;DP=%d", 
        numunique_taildistances,
        (int)(mq * 1.0 / variant->support), 
        mq30, 
        (int)((numdiffs * 1.0 / variant->support) + 0.5),
        calculate_cov_params(bam_name, 
                             variant->tid, 
                             variant->start - variant->lw - 1, 
                             variant->stop + variant->rw + 1));
 
    printf(";BF=%d,%d", lflank, rflank);
    printf("\n");
    ckfree(taildistances);
}

static void print_insertion_output(const variant* const variant,
                                   char** const sequences,
                                   const bam_header_t* const hin)
{
    forceassert(variant->type == INSERTION);
    
    int maxinsertsize = 0;
    char* consensus UNUSED = voted_consensus(variant->evidence, 
                                             variant->support,
                                             &maxinsertsize);
    ckfree(consensus);

    char* refsequence = sequences[variant->tid];
    uint sequencelen = hin->target_len[variant->tid];

    // bases on the reference sequence
    int i, j;         
    int lpos, lmid, rmid, rpos;
    int idx1, idx2, idx3, idx4;

    idx1 = 0;
    lpos = variant->start < neighborhood 
         ? 0 : variant->start - neighborhood;

    for(i = lpos, idx2 = 0; 
        i < (int)variant->start; 
        i++, idx2++) printf("%c", toupper(refsequence[i]));

    lmid = i;
    for(idx3 = idx2, j = 0; 
        j < maxinsertsize; 
        j++, idx3++) printf("-");
    rmid = i;
    rpos = (variant->stop + neighborhood) > sequencelen 
         ? sequencelen : variant->stop + neighborhood;
    for(idx4 = idx3; i < rpos; i++, idx4++) 
        printf("%c", toupper(refsequence[i]));
    printf("\n");

    // print all the supporting reads for the indel
    uint supindx;
    for(supindx = 0; supindx < variant->support; supindx++){
        evidence* evdnc = variant->evidence[supindx];

        readseg* seg1 = evdnc->aln1;
        // proceed to the read segment which can be shown
        readseg* iter;
        for(iter = seg1; iter && (iter->end < lpos); iter = iter->next);
        if(iter == NULL) continue;

        int rstart = iter->start > lpos ? 0 : lpos - iter->start;
        if(iter->start > lpos){
            for(i = lpos; i < iter->start; i++) printf(" ");
        }

        for(; iter; iter = iter->next){
            switch(iter->op){
                case BAM_CMATCH:
                    printf("%s", iter->sequence + rstart);
                    break;
                case BAM_CPMATCH:
                    printf("%s", iter->sequence + rstart);
                    break;
                case BAM_CMMATCH:
                    printf("%s", iter->sequence + rstart);
                    break;
                case BAM_CINS:
                    break;
                case BAM_CDEL:
                    for(i = rstart; i < iter->oplen; i++) printf("-");
                    break;
                case BAM_CREF_SKIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CSOFT_CLIP:
                    break;
                case BAM_CHARD_CLIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CPAD:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                default:
                    fatalf("Unknown CIGAR operation: %d", iter->op);
            }

            rstart = 0;
        }

        // proceed to the segment of interest (which bears the indel)
        readseg* seg2 = evdnc->aln2;
        forceassert(seg2->next == NULL);

        if(seg2->op == BAM_CINS){
            for(j = 0; j < seg2->oplen; j++){
                printf("%c", tolower(seg2->sequence[j]));
            }     
            for(; j < maxinsertsize; j++){
                printf("-");
            }
        }else{
            fatal("This segment should only contain the variation");
        }
        
        // finally print the following segments
        readseg* seg3 = evdnc->aln3;
        for(iter = seg3, i = idx3; iter && (i < idx4); iter = iter->next){
            switch(iter->op){
                case BAM_CMATCH:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){
                        printf("%c", iter->sequence[j]);
                    }
                    break;
                case BAM_CPMATCH:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){
                        printf("%c", iter->sequence[j]);
                    }
                    break;
                case BAM_CMMATCH:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){
                        printf("%c", iter->sequence[j]);
                    }
                    break;
                case BAM_CINS:
                    break;
                case BAM_CDEL:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){ 
                        printf("-");
                    }
                    break;
                case BAM_CREF_SKIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CSOFT_CLIP:
                    break;
                case BAM_CHARD_CLIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CPAD:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                default:
                    fatalf("Unknown CIGAR operation: %d", iter->op);
            }
        }
        for(; i < idx4; i++) printf(" ");
        printf("%s\n", evdnc->qname);

    }
}


// the detailed output would be as following:
// ########################################################################
// 1   chr20   69849374 69849394
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAaaaaa<10>aaaaaAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
//      AAAAAAAAAAAAAAAAAAAAAAAA              AAAAAAAAAAAAAAAAAAAAAA
// AAAAAAAAAAAAAAAA                                AAAAAAAAAAAAAAAAAAAAAAA
//                     AAAAAAAAA--------------AAAAAAAAAAAAAAAAAAA
//
// |                           |              |                           |   
// lpos                        lmid           rmid                        rpos 
// idx1                        idx2           idx3                        idx4
// lpos,lmid,rmid,rpos are all coordinates on the reference sequence.
// idx? are all coordinates on the figure (0-based)
static void print_deletion_output(const variant* const variant,
                                  char** const sequences,
                                  const bam_header_t* const hin)
{
    forceassert(variant->type == DELETION);

    char* refsequence = sequences[variant->tid];
    uint sequencelen = hin->target_len[variant->tid];

    char buffer[1024];

    // bases on the reference sequence
    int i, j;         
    int numdigits; // number of digits used in the size string
    int lpos, lmid, rmid, rpos;
    int idx1, idx2, idx3, idx4;
    if((variant->stop - variant->start) < 10){
        idx1 = 0;
        lpos = variant->start < neighborhood 
             ? 0 : variant->start - neighborhood;

        for(i = lpos, idx2 = 0; 
            i < (int)variant->start; 
            i++, idx2++) printf("%c", toupper(refsequence[i]));
        lmid = i;
        for(idx3 = idx2; 
            i < (int)variant->stop; 
            i++, idx3++) printf("%c",tolower(refsequence[i]));
        rmid = i;
        rpos = (variant->stop + neighborhood) > sequencelen 
             ? sequencelen : variant->stop + neighborhood;
        for(idx4 = idx3; i < rpos; i++, idx4++) 
            printf("%c", toupper(refsequence[i]));
    }else{
        idx1 = 0;
        lpos = variant->start < neighborhood 
             ? 0 : variant->start - neighborhood;
        for(i = lpos, idx2 = 0; 
            i < (int)variant->start; 
            i++, idx2++) printf("%c", toupper(refsequence[i]));

        lmid = i;
        for(idx3 = idx2; 
            i < (int)(variant->start + 5); 
            i++, idx3++) printf("%c",tolower(refsequence[i]));
        if((variant->stop - variant->start - 10) > 0){
            printf("<%d>", (variant->stop - variant->start - 10));
            sprintf(buffer, "<%d>", (variant->stop - variant->start - 10));
            numdigits = strlen(buffer);
            idx3 += numdigits;
        }
        i = variant->stop - 5;
        for(; i < (int)variant->stop; i++, idx3++){
            printf("%c",tolower(refsequence[i]));
        }
        rmid = i;
        rpos = (variant->stop + neighborhood) > sequencelen 
             ? sequencelen : variant->stop + neighborhood;
        for(idx4 = idx3; i < rpos; i++, idx4++) 
            printf("%c", toupper(refsequence[i]));
    }
    printf("\n");

    // print all the supporting reads for the indel
    uint supindx;
    for(supindx = 0; supindx < variant->support; supindx++){
        evidence* evdnc = variant->evidence[supindx];
        if(evdnc->type == PAIRED_READ) {
            printf("%s\n", evdnc->qname);
            continue;
        }
        readseg* seg1 = evdnc->aln1;
        // proceed to the read segment which can be shown
        readseg* iter;
        for(iter = seg1; iter && (iter->end < lpos); iter = iter->next);
        forceassert(iter != NULL);

        int rstart = iter->start > lpos ? 0 : lpos - iter->start;
        if(iter->start > lpos){
            for(i = lpos; i < iter->start; i++) printf(" ");
        }

        for(; iter; iter = iter->next){
            switch(iter->op){
                case BAM_CPMATCH:
                    printf("%s", iter->sequence + rstart);
                    break;
                case BAM_CMMATCH:
                    printf("%s", iter->sequence + rstart);
                    break;
                case BAM_CMATCH:
                    printf("%s", iter->sequence + rstart);
                    break;
                case BAM_CINS:
                    break;
                case BAM_CDEL:
                    for(i = rstart; i < iter->oplen; i++) printf("-");
                    break;
                case BAM_CREF_SKIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CSOFT_CLIP:
                    break;
                case BAM_CHARD_CLIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CPAD:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                default:
                    fatalf("Unknown CIGAR operation: %d", iter->op);
            }
            rstart = 0;
        }

        // proceed to the segment of interest (which bears the indel)
        readseg* seg2 = evdnc->aln2;
        forceassert(seg2->next == NULL);

        if(seg2->op == BAM_CDEL){
            for(i = idx2; i < idx3; i++){
                printf("-");
            }
        }else{
            fatal("This segment should only contain the variation");
        }

        // finally print the following segments
        readseg* seg3 = evdnc->aln3;
        for(iter = seg3, i = idx3; iter && (i < idx4); iter = iter->next){
            switch(iter->op){
                case BAM_CPMATCH:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){
                        printf("%c", iter->sequence[j]);
                    }
                    break;
                case BAM_CMMATCH:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){
                        printf("%c", iter->sequence[j]);
                    }
                    break;
                case BAM_CMATCH:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){
                        printf("%c", iter->sequence[j]);
                    }
                    break;
                case BAM_CINS:
                    break;
                case BAM_CDEL:
                    for(j = 0; (j < iter->oplen) && (i < idx4); j++, i++){ 
                        printf("-");
                    }
                    break;
                case BAM_CREF_SKIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CSOFT_CLIP:
                    break;
                case BAM_CHARD_CLIP:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                case BAM_CPAD:
                    fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
                    break;
                default:
                    fatalf("Unknown CIGAR operation: %d", iter->op);
            }
        }

        for(; i < idx4; i++) printf(" ");
        printf("%s\n", evdnc->qname);
    } 
}

static void print_det_output(const variant* const variant,
                             char** const sequences,
                             const bam_header_t* const hin)
{
    static int indel_index = 1;

    printf("###########################################################\n");
    printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n",
        indel_index++,
        hin->target_name[variant->tid],
        variant->start,
        variant->stop,
        variant->type == DELETION ? "Deletion" : "Insertion",
        variant->start,
        variant->stop + variant->rw + 1,
//        variant->type == DELETION ? variant->stop + variant->rw + 1 : variant->stop + variant->rw,
        variant->support);

    if(variant->type == DELETION){
        print_deletion_output(variant, sequences, hin);
    }else if(variant->type == INSERTION){
        print_insertion_output(variant, sequences, hin);
    }
}


void print_variants(variant** const pvariants,
                    char** const sequences,
                    const char* const bam_name,
                    const bam_header_t* const hin,
                    const char* const outputformat,
                    const uint minsupport,
                    const uint maxdiffsallowed)
{
    // collect the variants that have the minimum support
    variant* iter = *pvariants;
    variant* selectedvariants = NULL;
    variant* filteredvariants = NULL;
    uint i;

    while(iter){
        variant* tmp = iter->next;

        // in case of insertions we  want to see evidence from both sides of the
        // breakpoint
        bool left = FALSE;
        bool right = FALSE;
        int balance = INT_MAX;
        int lflank = -1, rflank = -1;

        // check for limit on number of differences allowed
        uint numdiffs = 0;
        for(i = 0; i < iter->support; i++){
            uint ltmp = 0;
            readseg* iter2;
            for(iter2 = iter->evidence[i]->aln1; iter2; iter2 = iter2->next){
                switch(iter2->op){
                    case BAM_CMATCH: 
                        ltmp += iter2->oplen;
                        uint j;
                        for(j = 0; j < iter2->oplen; j++){
                            if(iter2->sequence[j] != sequences[iter->tid][iter2->start + j]){
                                numdiffs += 1;
                            }
                        }
                        break;
                    case BAM_CPMATCH: ltmp += iter2->oplen; break;
                    case BAM_CMMATCH: 
                        ltmp += iter2->oplen; 
                        numdiffs += iter2->oplen;
                        break;
                    case BAM_CINS: 
                        ltmp += iter2->oplen; 
                        numdiffs += iter2->oplen;
                        break;
                    case BAM_CDEL: 
                        numdiffs += iter2->oplen;
                        break;
                    case BAM_CSOFT_CLIP: 
                        numdiffs += iter2->oplen;
                        break;
                    default: fatalf("unhandled BAM operation");
                }
            }
            if(ltmp >= minbalance) left = TRUE;
    
            uint rtmp = 0;
            for(iter2 = iter->evidence[i]->aln3; iter2; iter2 = iter2->next){
                switch(iter2->op){
                    case BAM_CMATCH: 
                        rtmp += iter2->oplen; 
                        uint j;
                        for(j = 0; j < iter2->oplen; j++){
                            if(iter2->sequence[j] != sequences[iter->tid][iter2->start + j]){
                                numdiffs += 1;
                            }
                        }
                        break;
                    case BAM_CPMATCH: rtmp += iter2->oplen; break;
                    case BAM_CMMATCH: 
                        rtmp += iter2->oplen; 
                        numdiffs += iter2->oplen;
                        break;
                    case BAM_CINS: 
                        rtmp += iter2->oplen; 
                        numdiffs += iter2->oplen;
                        break;
                    case BAM_CDEL: 
                        numdiffs += iter2->oplen;
                        break;
                    case BAM_CSOFT_CLIP: 
                        numdiffs += iter2->oplen;
                        break;
                    default: fatalf("unhandled BAM operation");
                }
            }
            if(rtmp >= minbalance) right = TRUE;

            if(abs(rtmp - ltmp) < balance) {
                balance = abs(rtmp - ltmp);
                lflank = ltmp;
                rflank = rtmp;
            }
        }
        uint xnumdiffs = (int)(numdiffs * 1.0 / iter->support) + 0.5;

        if(((iter->type == DELETION && lflank >= minbalance && rflank >= minbalance) || (iter->type == INSERTION && (lflank >= minbalance || rflank >= minbalance))) && (iter->support >= minsupport) && (xnumdiffs <= maxdiffsallowed) && left && right){
            sladdhead(&selectedvariants, iter);
        }else{
            sladdhead(&filteredvariants, iter);
        }
        iter = tmp;
    }         

    slreverse(&selectedvariants);
    slreverse(&filteredvariants);

    if(call_all_indels == TRUE){
        iter = selectedvariants;
        while(iter){
            if(0 == strcmp(outputformat, "vcf")){
                print_vcf_output(iter, sequences, bam_name, hin);
            }else if(0 == strcmp(outputformat, "detailed")){
                print_det_output(iter, sequences, hin);
            }
            iter = iter->next;
        }
    }else{
        // now only print all the variants that do not overlap another variant
        iter = selectedvariants;
        variant* olapiter = NULL;
        while(iter){
            olapiter = iter->next;
    
            while(olapiter && ((olapiter->start - olapiter->lw) <= (iter->stop + iter->rw))){
                olapiter = olapiter->next;
            }
    
            const variant* tmpiter;
            uint maxsupport = 0;
            const variant* chosenone = iter;
            for(tmpiter = iter; tmpiter != olapiter; tmpiter = tmpiter->next){
                if(tmpiter->support > maxsupport){
                    maxsupport = tmpiter->support;
                    chosenone = tmpiter;
                }
            }
    
            if(0 == strcmp(outputformat, "vcf")){
                print_vcf_output(chosenone, sequences, bam_name, hin);
            }else if(0 == strcmp(outputformat, "detailed")){
                print_det_output(chosenone, sequences, hin);               
            }
    
            iter = olapiter;
        }
    }

    iter = sllast(selectedvariants);
    if(iter != NULL) {
        iter->next = filteredvariants;
    }else{
        iter = filteredvariants;
    }

    *pvariants = iter;
}

// read the variants from the VCF file
knownvariant* read_variants(const char* const vcfname, 
                       char** const sequences,
                       const int tid,
                       const char* const chromname)
{
    size_t n = 1;
    char* fptr = ckalloc(n + 1);

    int numread = 0;              // number of variants read till now
    knownvariant* variants = NULL;// the variations themselves
    char chrom[128];              // chromosome
    uint start;                   // position

    char reference[maxpedelsize];
    char alternate[maxpedelsize];

    char type[128];
    char evdnctype[128];
    char supportstring[128];
    char stopstring[128];
    char bpstopstring[128];
    uint support;
    uint stop;
    uint bpstop;
    char information[1024];

    FILE* fp = ckopen(vcfname, "r");

    while(getline(&fptr, &n, fp) != -1){
        if(fptr[0] == '#') continue;
        if(sscanf(fptr, "%s %u %*c %s %s %*c %*c %[^;];%[^;];NS=%[^;];END=%[^;];BP_END=%[^;];%s\n",  chrom, &start, reference, alternate, type, evdnctype, supportstring ,stopstring, bpstopstring, information) != 10){
            fatalf("Error in reading the variant : %s", fptr);   
        }
        support = atoi(supportstring);
        stop = atoi(stopstring);
        bpstop = atoi(bpstopstring);

        if(strcmp(chrom,chromname) != 0) continue;
        numread++;

        knownvariant* variation = ckallocz(sizeof(knownvariant));
        variation->tid = tid;
        variation->start = start;
        variation->addntlinfo = ckalloc(strlen(information)+1);
        strcpy(variation->addntlinfo, information);
        
        variation->reference = ckallocz(strlen(reference) + 1);
        sprintf(variation->reference, "%s", reference);
        variation->alternate = ckallocz(strlen(alternate) + 1);
        sprintf(variation->alternate, "%s", alternate);
        
        variation->type = strncmp(type, "DELETION", 8) == 0
                        ? DELETION : INSERTION;
        if(strcmp(evdnctype, "SPLIT_READ") == 0){
            variation->evdnctype = SPLIT_READ;
        }else if(strcmp(evdnctype, "PAIRED_READ") == 0){
            variation->evdnctype = PAIRED_READ;
        }else if(strcmp(evdnctype, "COMPOSITE") == 0){
            variation->evdnctype = COMPOSITE;
        }else{          
            fatalf("unknown evidence type for this variant");
        }

        variation->support = support;
        variation->stop = stop;
        variation->bpstop = bpstop;
        variation->sequence = sequences[tid];


        sladdhead(&variants, variation);
    }

    fclose(fp);

    variants = sort_knownvariants(variants);

    ckfree(fptr);

    fprintf(stderr, "Read %d variants for %s\n", numread, chromname);
    return variants;
}

static void move_boundaries(variant* const vs,
                            const char* const reference)
{
    int maxinsertsize;
    char* consensus = NULL;
    uint lw = 0;
    uint rw = 0;

    if(vs->type == INSERTION){
        consensus = voted_consensus(vs->evidence,
                                    vs->support,
                                    &maxinsertsize);

        while(strncmp(consensus, 
                      reference + vs->start - strlen(consensus) - lw,
                      strlen(consensus)) == 0){
            lw += strlen(consensus);            
        }
        uint shift = 0;
        while((shift < strlen(consensus)) &&
              (consensus[strlen(consensus)-shift-1]==reference[vs->start-1-lw])){
            lw++;
            shift++;
        }

//       this is the older version of this code. I will keep this around for now
//        while((lw < strlen(consensus)) &&
//              (consensus[strlen(consensus)-lw-1]==reference[vs->start-1-lw])){
//            lw++;
//        }
//        if(lw == strlen(consensus)){
//            while(strncmp(consensus,
//                  reference + vs->start - lw,
//                  strlen(consensus)) == 0){
//            lw += strlen(consensus);
//            }
//        }
 
        while(strncmp(consensus,
                      reference + vs->stop + rw,
                      strlen(consensus)) == 0){
            rw += strlen(consensus);
        }
        shift = 0;
        while(shift < strlen(consensus) && 
             (consensus[shift] == reference[vs->stop+rw])){
            rw++;
            shift++;
        }
     
//       this is the older version of this code. I will keep this around for now
//        while((rw < strlen(consensus)) &&
//              (consensus[rw] == reference[vs->stop+rw])) rw++;
//        if(rw == strlen(consensus)){
//            while(strncmp(consensus,
//                          reference + vs->stop + rw,
//                          strlen(consensus)) == 0){
//                rw += strlen(consensus);
//            }
//        }
    }else if(vs->type == DELETION){
        while(reference[vs->start-1-lw] == reference[vs->stop-1-lw]) lw++;
        while(reference[vs->start+rw] == reference[vs->stop+rw]) rw++;
    }
        
    vs->lw = lw;
    vs->rw = rw;
    ckfree(consensus);
}

static variant* mergeSRvariants(const variant* const v1,
                                  const variant* const v2)
{
    forceassert(v1->type == v2->type);
    forceassert(v1->start <= v2->start);
    forceassert(v1->tid == v2->tid) 
    forceassert(v1->evdnctype == SPLIT_READ);
    forceassert(v2->evdnctype == SPLIT_READ);

    variant* variation = ckallocz(sizeof(variant));
    variation->type = v1->type;
    variation->tid = v1->tid;
    variation->evdnctype = SPLIT_READ;
    
    variation->start = v1->start;
    variation->stop  = v1->stop;
    variation->lw = v1->lw;
    variation->rw = v1->rw;

    variation->support = v1->support + v2->support;

    variation->evidence = ckalloc(variation->support * sizeof(evidence*));
    memcpy(variation->evidence, 
           v1->evidence, 
           v1->support*sizeof(evidence*));
    memcpy(variation->evidence + v1->support, 
           v2->evidence, 
           v2->support*sizeof(evidence*));

    return variation;
}

// merge the SR variants if the exact position of the SV is unclear since bases
// at the edge of one read-half could equally well be appended to the other
// read-half. Then merge the SR and PE variants if they share the same
// breakpoints
void merge_variants(variant** pvs, 
                    const char* const reference, 
                    const bool join_sr_pe)
{
    variant* vs = *pvs;

    if(vs == NULL) return;

    // if it is an SR identified variant, then lets move the boundaries around
    // to see if  there are more than one possible breakpoints
    variant* iter1;
    for(iter1 = vs; iter1; iter1 = iter1->next){
        if(iter1->evdnctype == SPLIT_READ){
            move_boundaries(iter1, reference);
//            fprintf(stderr, 
//                    "%d:%d %d:%d\n", 
//                    iter1->start, iter1->lw, iter1->stop, iter1->rw);
        }
    }
    
    variant* pevariants = NULL;
    iter1 = vs;
    variant* tail1 = NULL;
    while(iter1){
        if(iter1->evdnctype == PAIRED_READ){
            if(tail1 == NULL){
                tail1 = iter1->next;
                sladdhead(&pevariants, iter1);
                iter1 = tail1;
                tail1 = NULL;
                vs = iter1;
                continue;
            }else{
                tail1->next = iter1->next;
                sladdhead(&pevariants, iter1);
                iter1 = tail1;
            }
        }    

        tail1 = iter1;
        iter1 = iter1->next;
    }   
    slreverse(&pevariants);     


    // lets sort the variants again
    vs = sort_variants(vs);

    // now lets combine the SR identified variants with each other if they
    // support the same breakpoints
    iter1 = vs;
    tail1 = NULL;

    while(iter1){
        variant* nvr = iter1;

        variant* iter2 = iter1->next;
        variant* tail2 = iter1;

        while(iter2 && (iter2->start <= (iter1->stop + iter1->rw))){
            forceassert(iter1->evdnctype == SPLIT_READ);
            forceassert(iter2->evdnctype == SPLIT_READ);

            if(iter2->type == iter1->type){
                if(((iter1->start - iter1->lw) == (iter2->start - iter2->lw)) &&
                   ((iter1->stop + iter1->rw) == (iter2->stop + iter2->rw))){
                    nvr = mergeSRvariants(iter1, iter2);

                    if(tail1 != NULL){
                        tail1->next = nvr;
                    }else{
                        vs = nvr;
                    }

                    if(tail2 != NULL){
                        tail2->next = iter2->next;
                    }
                    ckfree(iter2->evidence);   
                    ckfree(iter2);

                    nvr->next = iter1->next;
                    ckfree(iter1->evidence);
                    ckfree(iter1);

                    iter1 = nvr;
                    iter2 = iter1;
                }
            }
            tail2 = iter2;
            iter2 = iter2->next;
        }

        tail1 = iter1;
        iter1 = iter1->next;
    }

    if(join_sr_pe == FALSE){
        iter1 = sllast(vs);
        if(iter1 != NULL){
            iter1->next = pevariants;
        }else{
            vs = pevariants;
        }

        *pvs = sort_variants(vs);
        return;
    }
   
    // now merge the variants that could share the same breakpoints, but are
    // supported by SR and PE reads. In all these cases we trust the SR
    // breakpoints.
    variant* iter2;
    variant* candidate;

    iter1 = pevariants;
    tail1 = NULL;
    while(iter1){
        uint overlap = 0;
        candidate = NULL;        

        for(iter2 = vs; iter2; iter2 = iter2->next){
            if(iter2->start > iter1->stop) break;
            if(iter2->evdnctype == PAIRED_READ) continue;

            uint olapf = 0;
            uint olap = 0;
            if((iter1->start >= iter2->start) && 
               (iter1->start < iter2->stop)){
                olap = MIN(iter1->stop,iter2->stop) - iter1->start;
            }else if((iter2->start >= iter1->start) &&
                     (iter2->start < iter1->stop)){
                olap = MIN(iter1->stop,iter2->stop) - iter2->start;
            }

            olapf = (olap * 100.0 / (iter1->stop - iter1->start))
                  + (olap * 100.0 / (iter2->stop - iter2->start));

            if(olapf > overlap){
                overlap = olapf;
                candidate = iter2;
            }
        }

        bool tomerge = TRUE;
        if(candidate != NULL){
            int size = candidate->stop - candidate->start;
            
            for(uint i = 0; i < iter1->support; i++){
                if(size < iter1->evidence[i]->mindelsize){
                    tomerge = FALSE;
                    break;   
                }
            }
        }

        if(candidate != NULL && tomerge == TRUE){
            candidate->support += iter1->support;   
            if ((candidate->evdnctype == PAIRED_READ && 
                 iter1->evdnctype == SPLIT_READ)     || 
                (candidate->evdnctype == SPLIT_READ  && 
                 iter1->evdnctype == PAIRED_READ))  {
            candidate->evdnctype = COMPOSITE;
            }
            candidate->evidence = ckrealloc(candidate->evidence,
                                  candidate->support * sizeof(evidence*));

            uint i = candidate->support - iter1->support;
            uint j = 0;
            for(; j < iter1->support; j++){
                candidate->evidence[i++] = iter1->evidence[j];
            }
        }else{
            if(tail1 == NULL){
                tail1 = iter1->next;
                sladdhead(&vs, iter1);
                iter1 = tail1;
                tail1 = NULL;
                pevariants = iter1;
                continue;
            }else{
                tail1->next = iter1->next;
                sladdhead(&vs, iter1);
                iter1 = tail1;
            }
        }

        tail1 = iter1;
        iter1 = iter1->next;
    }

    for(iter1 = pevariants; iter1; iter1 = iter1->next){
        ckfree(iter1->evidence);
    }
    slfreelist(&pevariants);
   
    *pvs = sort_variants(vs);
}

void print_vcf_line(const knownvariant* const kiter,
                    const bam_header_t* const hin)
{   
    printf("%s\t%d\t.\t%s\t%s\t.\t.\t%s;", 
        hin->target_name[kiter->tid], 
        kiter->start, 
        kiter->reference, 
        kiter->alternate, 
        kiter->type == DELETION ? "DELETION" : "INSERTION");
    if(kiter->evdnctype == SPLIT_READ){
        printf("SPLIT_READ;");
    }else if(kiter->evdnctype == PAIRED_READ){
        printf("PAIRED_READ;");
    }else if(kiter->evdnctype == COMPOSITE){
        printf("COMPOSITE;");
    }
    printf("NS=%d;END=%d;BP_END=%d;%s", kiter->support, kiter->stop, kiter->bpstop, kiter->addntlinfo);
}

static void realign_with_indel(const char* const reference,
                               const int rstart,
                               const int rstop,
                               const readaln* const rln,
                               const int qstart,
                               const int qstop,
                               const knownvariant* const variant,
                               int* const alnsubs,
                               int* const alnindels,
                               int* const alnaligned)
{
    char* query = rln->segments->sequence;

    // lets create a fake reference sequence which has the variant in it
    char* target = copy_partial_string(reference + rstart, rstop - rstart);
    if(variant->type == DELETION){
        memmove(target + variant->start - rstart,
                target + variant->stop  - rstart - 1,
                rstop - variant->stop + 2);
    }else if(variant->type == INSERTION){
        target = ckreallocz(target, rstop - rstart + 1, rstop - rstart + 1 + strlen(variant->alternate) - 1);
        memmove(target + variant->start - rstart + strlen(variant->alternate)-1, 
                target + variant->start - rstart, 
                strlen(target + variant->start - rstart));
        memcpy(target + variant->start - rstart, variant->alternate + 1, strlen(variant->alternate) - 1);
    }else{
        fatalf("Unknown type of variant: %d\n", variant->type);
    }     

    // now align the query to the target
    char* t1 = target;
    char* t2 = query + qstart;

    uint len1 = strlen(target);
    uint len2 = qstop - qstart;
    if(strlen(t2) < len2){
        len2 = strlen(t2);
    }

    // counters in loops
    uint i, j;

    // initialize the scores
    int match    = 2;
    int mismatch = 1;
    int gopen    = 4;
    int gextend  = 1;

    // allocate the dp matrix
    int** V = ckalloc((len2 + 1)*sizeof(int*));
    for(i = 0; i <= len2; i++){
        V[i] = ckallocz((len1 + 1) * sizeof(int));
        V[i][0] = -gopen -(i*gextend);
    }
    for(j = 0; j <= len1; j++){
        V[0][j] = -gopen -(j*gextend);
    }

    // allocate E,F (two other variables to fill the dp matrix)
    int E;
    int* F = ckallocz((len1+1)*sizeof(int));

    // keep track of the indices for backtracking
    char** I = ckalloc((len2 + 1)*sizeof(char*));
    for(i = 0; i <= len2; i++){
        I[i] = ckallocz((len1 + 1) * sizeof(char));
    }

    // variables to store intermediate values if this is a sub/indel.
    int ifsub, ifdel, ifins, ifindel;

    // maximum score of the alignment
    int max_score = 0, max_i = -1, max_j = -1;

    // find the scores for the dp matrix
    for(i = 1; i <= len2; i++){
        E = 0;
        for(j = 1; j <= len1; j++){
            ifsub = V[i-1][j-1];
            ifsub = toupper(t1[j-1]) == toupper(t2[i-1]) ? \
                                        ifsub + match : ifsub - mismatch;

            ifins = MAX(F[j], V[i-1][j] - gopen) - gextend;
            F[j]  = ifins;
            
            ifdel = MAX(E, V[i][j-1] - gopen) - gextend;
            E     = ifdel;

            ifindel = MAX(ifins, ifdel);

            // priority if the scores are equal are ins, del, sub.
            V[i][j] = ifsub;
            I[i][j] = SUB;
            if(V[i][j] < ifindel){
                if(ifins >= ifdel) I[i][j] = INS; else I[i][j] = DEL;
                V[i][j] = ifindel;
            }

            // an alignment which includes more bases from the reference is
            // given preference, hence the >=
            if(V[i][j] > max_score){
                max_score = V[i][j];
                max_i     = i;
                max_j     = j;
            }
        }
    }
    ckfree(F);
  
//    // debug statements
//    printf("----------DP MATRIX--------\n");
//    for(i = 0; i <= len2; i++){
//        for(j = 0; j <= len1; j++){
//            printf("%d\t", V[i][j]);
//        }
//        printf("\n");
//    }
//    printf("----------INDEX MATRIX--------\n");
//    for(i = 0; i <= len2; i++){
//        for(j = 0; j <= len1; j++){
//            printf("%d\t", I[i][j]);
//        }
//        printf("\n");
//    }
//    // debug statements ends
 
    char* nt1 = ckallocz((len1+len2)*sizeof(char));
    char* nt2 = ckallocz((len1+len2)*sizeof(char));

    // trace back to find the best alignment 
    char direction = SUB;
    int score = max_score;
    uint a = 0, b = 0;
    i = max_i;
    j = max_j;
   
    while(score > 0){
        direction = I[i][j];

        if(SUB == direction){
            nt1[a++] = t1[j-1];
            nt2[b++] = t2[i-1];
            i--; j--;
        }else if(INS == direction){
            nt1[a++] = '-';
            nt2[b++] = t2[i-1];
            i--;
        }else{
            nt1[a++] = t1[j-1];
            nt2[b++] = '-';
            j--;
        }
        score = V[i][j];
    }


    assert(strlen(nt1) == strlen(nt2));

    int k, subs = 0, ins = 0, dels = 0, aligned = 0;
    for(k = strlen(nt1); k >= 0; k--){
        if(nt1[k] == '-'){
            ins++;
            aligned++;
        }else if(nt2[k] == '-'){
            dels++;
        }else if(nt1[k] != nt2[k]){
            subs++;
            aligned++;
        }else{
            aligned++;
        }
    }   
        
    ckfree(target);
    
    *alnsubs = subs;
    *alnindels = ins + dels;
    *alnaligned = aligned;
}


static int check_for_indel(const bam1_t* alignment, void* data)
{
    knownvariant* variant = data;
    
    // nothing to be done if this read is not aligned
    bool is_aligned = (alignment->core.flag & 0x4) == 0;
    if(is_aligned == FALSE){
        return TRUE;  
    }

    // nothing to be done if this is a secondary or supplementary alignment
    if((alignment->core.flag & 0x100) == 0x100) return TRUE;
    if((alignment->core.flag & 0x200) == 0x200) return TRUE;
    if((alignment->core.flag & 0x400) == 0x400) return TRUE;
    if((alignment->core.flag & 0x800) == 0x800) return TRUE;

    // if some read already supports this indel, then I do not need to do
    // anything else
    if(variant->diffsample_support == TRUE){
        return TRUE;
    }

    int aln1subs = 0; // how many substitutions in the current alignment
    int aln1indels = 0; // how many indels in the current alignment
    int aln1aligned = 0; // how many bases are in the alignment?

    bool overlaps = FALSE;
    readaln* rln = new_readaln(alignment);
    readseg* iter;
    int i, j, readindx;
    int qstart = -1, qstop = -1; // start,stop of the query that was aligned 
    for(iter = rln->segments, readindx = 0; iter; iter = iter->next){
        if((iter->end < (int)variant->start) || 
           (iter->start > (int)variant->stop)){
        }else{
            overlaps = TRUE;
        }
        switch(iter->op){
            case BAM_CSOFT_CLIP: 
                if(iter->next == NULL) qstop = readindx;
                readindx += iter->oplen;
                break;
            case BAM_CINS: 
                if(qstart == -1) qstart = readindx;
                if((variant->type == INSERTION) && 
                   (iter->start == (int)variant->start)){
                    // this read supports the indel which was missed when I was
                    // realigning
                    free_readaln(&rln);
                    variant->diffsample_support = TRUE;
                    return TRUE;
                }      
                readindx += iter->oplen;
                aln1indels += iter->oplen;
                aln1aligned += iter->oplen;
                break;
            case BAM_CDEL:
                if((variant->type == DELETION) && 
                   (iter->start == (int)variant->start) && 
                   (iter->end == (int)(variant->stop - 1))){
                    // this read supports the indel which was missed
                    free_readaln(&rln);
                    variant->diffsample_support = TRUE;
                    return TRUE;
                }

                aln1indels += iter->oplen;
                break;
            case BAM_CMATCH:
                if(qstart == -1) qstart = readindx;
                for(i = 0, j = iter->start; i < iter->oplen; i++, j++){
                    if(iter->sequence[i] != variant->sequence[j]){
                        aln1subs += 1;
                    }
                }
                readindx += iter->oplen;
                aln1aligned += iter->oplen;
                break;  
            default:
                fatalf("Unhandled CIGAR op: %d", iter->op); 
        }  
    }
    if(qstop == -1) qstop = aln1aligned + qstart;
//    fprintf(stderr, "%s %d %d %d %d\n", bam1_qname(alignment), aln1aligned, qstop, qstart, qstop - qstart);
    forceassert(aln1aligned == (qstop - qstart));
//    fprintf(stderr, "1. %s %d\n", bam1_qname(alignment), aln1diffs);


    // no point in doing this if the read does not overlap the variant
    if(overlaps == FALSE){
        free_readaln(&rln);
        return TRUE;
    }

    // no point in doing this if the read doesnt cover the complete variant
    int indelsize = abs(strlen(variant->alternate) - strlen(variant->reference));
    int rstart = alignment->core.pos;
    int rstop = bam_calend(&alignment->core, bam1_cigar(alignment));
    if((uint)rstop < variant->bpstop){
        free_readaln(&rln);
        return TRUE;
    }

    free_readaln(&rln);
    rln = new_unaligned_readaln(alignment, 30);
    forceassert(qstart != -1);
    forceassert(qstop != -1);

    // how many differences would I see if I changed the alignment so that the
    // read supports the indel
    int aln2subs = 0;
    int aln2indels = 0;
    int aln2aligned = 0;
    realign_with_indel(variant->sequence, rstart-indelsize, rstop+indelsize,
                       rln, qstart, qstop, variant,
                       &aln2subs, &aln2indels, &aln2aligned);

//    fprintf(stderr, "%s :: %d:%d:%d :: %d:%d:%d\n", bam1_qname(alignment), aln1subs, aln1indels, aln1aligned, aln2subs, aln2indels, aln2aligned);

    // if the number of substitutions do not increase, number of indels do not
    // increase and number of aligned bases do not decrease, then this certainly
    // could be a read which supports this indel
    if((aln2subs <= aln1subs) && 
       (aln2indels <= aln1indels) && 
       (aln2aligned >= aln1aligned)){
        variant->diffsample_support = TRUE;
    }

    return TRUE;
}



// last attempt to check whether an indel is supported in the given BAM file
bool is_indel_supported(knownvariant* const variant, 
                        const char* const bam_name)
{
    bamFile fp = bam_open(bam_name, "r");
    bam_index_t* fp_index = bam_index_load(bam_name);

    bam_fetch(fp, fp_index, variant->tid, variant->start, variant->stop, variant ,check_for_indel);

    bam_close(fp);
    bam_index_destroy(fp_index);

    return variant->diffsample_support;
}

// print the variants that we already know about along with the tag for this 
// sample (in case this variant was also found in this sample) 
knownvariant* print_knownvariants(knownvariant* const kvariants, 
                                  variant* const variants, 
                                  const char* const sample_name,
                                  const char* const bam_name,
                                  const char* const outputformat UNUSED,
                                  const bam_header_t* const hin)
{
    knownvariant* kiter;
    const variant* uiter;

    if(variants == NULL) return kvariants;
    uiter = variants;

    for(kiter = kvariants; kiter; kiter = kiter->next){
        bool is_found = FALSE;

        uint kstart = kiter->start;
        uint kstop = kiter->stop;

        const variant* first_overlap = NULL;
//        if((kstart < (uiter->start - uiter->lw)) || 
//           (kiter->evdnctype == PAIRED_READ) || 
//           (kiter->evdnctype == COMPOSITE)){
//            uiter = variants;
//        }        

        // the best way here would be to sort the variants using ->start - ->lw
        // and then using two iterators to move through the two lists. Right now
        // I do this really dumb thing of going through all the evidence for
        // every known variant. $$$
        uiter = variants;

        forceassert((uiter == variants) || 
                    (kstart >= (uiter->start - uiter->lw)));

        for(; uiter; uiter = uiter->next){
            uint ustart = uiter->start - uiter->lw;
            uint ustop = 0;

            int reflength = (uiter->stop + uiter->rw) 
                          - (uiter->start - uiter->lw - 1);
            forceassert(reflength >= 1);

            if(uiter->type == DELETION){
                int altlength = (uiter->start + uiter->rw)  
                              - (uiter->start - uiter->lw - 1);
                forceassert(altlength >= 1);
                ustop = uiter->start - uiter->lw + reflength - altlength + 1;
            }else if(uiter->type == INSERTION){
                ustop = ustart + 1;
            }
            forceassert(ustop != 0);

            if(kstart >= ustop){
            }else if(ustart >= kstop){
            }else{
                if(first_overlap == NULL) first_overlap = uiter;

                // there is some overlap
                if(((kiter->evdnctype == SPLIT_READ) || 
                    (kiter->evdnctype == COMPOSITE)) && 
                   (uiter->evdnctype == SPLIT_READ)){
                    if((kstart == ustart) && (kstop == ustop)){
                        is_found = TRUE;
                        break;
                    }
                }else if((((kiter->evdnctype == SPLIT_READ) || 
                           (kiter->evdnctype == COMPOSITE)) &&  
                          (uiter->evdnctype == PAIRED_READ)) || 
                         ((kiter->evdnctype == PAIRED_READ) && 
                          (uiter->evdnctype == SPLIT_READ)) ||
                         ((kiter->evdnctype == PAIRED_READ) &&
                          (uiter->evdnctype == PAIRED_READ))){
                    uint s = MAX(kiter->start, uiter->start);
                    uint e = MIN(kiter->bpstop, uiter->stop);
    
                    uint olap = 0;
                    if(e >= s){
                        olap = e - s;
                    }
    
                    if((olap * 100.00 / (kiter->bpstop - kiter->start)) > 50){
                        is_found = TRUE;
                        break;
                    }
                }
            }
        }

        if(uiter == NULL){
            const variant* testvar = sllast(variants);
            if(testvar->start < kstart){
                break;
            }
        }

        if(first_overlap == NULL){
            uiter = variants;
        }else{
            uiter = first_overlap;
        }

        print_vcf_line(kiter, hin); 
        if(is_found == TRUE){
            printf(";%s", sample_name);      
        }else{
            if((kiter->evdnctype == SPLIT_READ) &&
               (is_indel_supported(kiter, bam_name) == TRUE)){
                printf(";%s", sample_name);
            }
        }
        printf("\n");
    }

    return kiter;
}

