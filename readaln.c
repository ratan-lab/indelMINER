#include "readaln.h"

// convert 4-bit integer representing the base to nucleotide
char bit2char(const int encodedbase)
{
    char nuc = 'X';
    switch(encodedbase & 0xF){
        case 1: nuc = 'A'; break;
        case 2: nuc = 'C'; break;
        case 4: nuc = 'G'; break;
        case 8: nuc = 'T'; break;
        case 15: nuc = 'N'; break;
        default: fatalf("Unhandled base encoding : %d:%d\n", 
                 encodedbase, encodedbase & 0xF);
    }
    return nuc;
}

static inline readseg* allocreadseg(){
    readseg* rs = ckallocz(sizeof(readseg));
    return rs;
}

readseg* new_readseg(const char* const sequence,
                     const uint32_t cigar,
                     int* const prefindx,
                     int* const preadindx)
{

    int refindx = *prefindx;
    int readindx = *preadindx;

    readseg* rs = allocreadseg();
    rs->op    = cigar & BAM_CIGAR_MASK;
    rs->oplen = cigar >> BAM_CIGAR_SHIFT;
    rs->sequence = ckallocz(rs->oplen + 1);

    int j;
    switch(rs->op){
        case BAM_CPMATCH:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = sequence[readindx+j];
            }
            readindx += j;
            refindx += j;
            rs->end = refindx;
            break;
        case BAM_CMMATCH:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = sequence[readindx+j];
            }
            readindx += j;
            refindx += j;
            rs->end = refindx;
            break;
        case BAM_CINS:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = sequence[readindx+j];
            }
            readindx += j;
            rs->end = refindx;
            break;
        case BAM_CDEL:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = '-';
            }
            refindx += j;
            rs->end = refindx;
            break;
        case BAM_CREF_SKIP:
            fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
            break;
        case BAM_CSOFT_CLIP:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                rs->sequence[j] = sequence[readindx+j];
            }
            readindx += j;
            rs->end = refindx;
            break;
        case BAM_CHARD_CLIP:
            fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
            break;
        case BAM_CPAD:
            fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
            break;
        default:
            fatal("Unhandled cigar operation");
            break;
    }

    *prefindx = refindx;
    *preadindx = readindx;
    return rs;
}

static readseg* new_readseg_bam(const bam1_t* const alignment,
                                const uint32_t cigar,
                                int* const prefindx,
                                int* const preadindx)
{
    int refindx = *prefindx;
    int readindx = *preadindx;

    readseg* rs = allocreadseg();
    rs->op    = cigar & BAM_CIGAR_MASK;
    rs->oplen = cigar >> BAM_CIGAR_SHIFT;
    rs->sequence = ckallocz(rs->oplen + 1);

    int j;
    switch(rs->op){
        case BAM_CMATCH:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = 
                 bit2char(bam1_seqi(bam1_seq(alignment),readindx+j));
            }
            readindx += j;
            refindx += j;
            rs->end = refindx;
            break;
        case BAM_CPMATCH:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = 
                 bit2char(bam1_seqi(bam1_seq(alignment),readindx+j));
            }
            readindx += j;
            refindx += j;
            rs->end = refindx;
            break;
        case BAM_CMMATCH:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = 
                 bit2char(bam1_seqi(bam1_seq(alignment),readindx+j));
            }
            readindx += j;
            refindx += j;
            rs->end = refindx;
            break;
        case BAM_CINS:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = 
                 bit2char(bam1_seqi(bam1_seq(alignment),readindx+j));
            }
            readindx += j;
            rs->end = refindx;
            break;
        case BAM_CDEL:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                 rs->sequence[j] = '-';
            }
            refindx += rs->oplen;
            rs->end = refindx;
            break;
        case BAM_CREF_SKIP:
            fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
            break;
        case BAM_CSOFT_CLIP:
            rs->start = refindx;
            for(j = 0; j < rs->oplen; j++){
                rs->sequence[j] = 
                bit2char(bam1_seqi(bam1_seq(alignment),readindx+j));
            }
            readindx += j;
            rs->end = refindx;
            break;
        case BAM_CHARD_CLIP:
            fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
            break;
        case BAM_CPAD:
            fatalf("Implement %s:%d", __FUNCTION__, __LINE__);
            break;
        default:
            fatal("Unhandled cigar operation");
            break;
    }

    *prefindx = refindx;
    *preadindx = readindx;
    return rs;
}

// glean the information from this BAM alignment
readaln* new_readaln(const bam1_t* const alignment)
{   
    readaln* rln = ckallocz(sizeof(readaln));
    
    char* name = ckallocz(alignment->core.l_qname + 1);
    memcpy(name,bam1_qname(alignment),alignment->core.l_qname);
    rln->qname = name;
    rln->qual  = alignment->core.qual;
 
    rln->index = (alignment->core.flag & 0x40) == 0x40  
               ? '1' : '2';
    rln->strand = (alignment->core.flag & 0x10) == 0x10
                ? '-' : '+';

    bool is_aligned = (alignment->core.flag & 0x4) == 0;

    if(is_aligned){
        rln->tid = alignment->core.tid;
        rln->segments = NULL;
    
        uint32_t* cigarstring = bam1_cigar(alignment);
        int i, readindx, refindx;
    
        readindx = 0;
        refindx = alignment->core.pos; 
        for(i = 0; i < alignment->core.n_cigar; i++){
            readseg* rs = new_readseg_bam(alignment, 
                                          cigarstring[i], 
                                          &refindx,  // 0 based coordinate
                                          &readindx);
            sladdhead(&rln->segments, rs);
        }
        
        slreverse(&rln->segments);  
    }else{
        rln->tid = -1;

        uint32_t cigar = (alignment->core.l_qseq << BAM_CIGAR_SHIFT)+BAM_CSOFT_CLIP;
        int refindx = 0;
        int readindx = 0;
        rln->segments = new_readseg_bam(alignment, cigar, &refindx, &readindx);

        rln->segments->start = -1;
        rln->segments->end = -1;
    }

    return rln;
}

// create a read from this BAM alignment
readaln* new_unaligned_readaln(const bam1_t* const alignment, const uint8_t mmq)
{
    readaln* rln = ckallocz(sizeof(readaln));
    
    char* name = ckallocz(alignment->core.l_qname + 1);
    memcpy(name,bam1_qname(alignment),alignment->core.l_qname);
    rln->qname = name;
    rln->qual  = mmq;
 
    rln->index = (alignment->core.flag & 0x40) == 0x40  
               ? '1' : '2';
    rln->strand = (alignment->core.flag & 0x10) == 0x10
                ? '-' : '+';

    rln->tid = -1;

    uint32_t cigar = (alignment->core.l_qseq << BAM_CIGAR_SHIFT)+BAM_CSOFT_CLIP;
    int refindx = 0;
    int readindx = 0;
    rln->segments = new_readseg_bam(alignment, cigar, &refindx, &readindx);

    rln->segments->start = -1;
    rln->segments->end = -1;

    return rln;
}

void print_readaln(const readaln* const rln)
{
    printf("Read : %s : %d\n", rln->qname, rln->tid);
    
    readseg* iter;
    for(iter = rln->segments; iter; iter = iter->next){
        printf("%d-%d:%s\n", iter->start, iter->end, iter->sequence);
    }
    printf("\n");        
}

// combine all these readseg to one [begin, end)
readseg* combine_readsegs(const readseg* const begin, const readseg* const end)
{
    if(begin == end) return NULL;

    readseg* readsegs = NULL;

    const readseg* iter;
    for(iter = begin; iter != end; iter = iter->next){
        readseg* rs = duplicate_readseg(iter);
        sladdhead(&readsegs, rs);
    }

    if(readsegs != NULL) slreverse(&readsegs);

    return readsegs;
}

// free the resources used by this
void free_readaln(readaln** prln)
{
    readaln* rln = *prln;

    ckfree(rln->qname);
    
    readseg* iter;
    for(iter = rln->segments; iter; iter = iter->next){
        ckfree(iter->sequence);
    }    
    slfreelist(&rln->segments);
    
    ckfree(rln);
    *prln = NULL;   
}

// free the resources used by read segment
void free_readsegs(readseg** prs)
{
    readseg* rs = *prs;
    
    readseg* iter; 
    for(iter = rs; iter; iter = iter->next){
        ckfree(iter->sequence);
    }
    slfreelist(&rs); 

    *prs = NULL; 
}

// duplicate this read segment. Duplicate all members except the next pointer.
readseg* duplicate_readseg(const readseg* const rs)
{     
    readseg* dup = allocreadseg();

    dup->sequence = NULL;
    if(rs->sequence){     
        dup->sequence = ckallocz(strlen(rs->sequence) + 1);
        memcpy(dup->sequence, rs->sequence, strlen(rs->sequence));
    }
    
    dup->oplen = rs->oplen;
    dup->op = rs->op;
    dup->start = rs->start;
    dup->end = rs->end;

    return dup; 
}
// modify the segments for this readaln, as per the given cigarstring
void update_readsegs(readaln* const rln,
                     const int r1,
                     const uint32_t* cigarstring1,
                     const int numcigarops1,
                     const int index,
                     const int q2,
                     const int r2,
                     const uint32_t* cigarstring2,
                     const int numcigarops2)
{
    int i, j;
    int refindx  = r1;
    int readindx = 0;
    readseg* readsegs = NULL;
    for(i = 0, j = 0; i < numcigarops1; i++){
        int op    = (cigarstring1[i] & BAM_CIGAR_MASK);
        int oplen = (cigarstring1[i] >> BAM_CIGAR_SHIFT);
        forceassert(oplen > 0);
        if(op != BAM_CDEL) j += oplen;

        if(j <= index){
            readseg* rsg = new_readseg(rln->segments->sequence, 
                                       cigarstring1[i], 
                                       &refindx, &readindx);
            sladdhead(&readsegs, rsg);
        }

        if(j > index){
            uint32_t cstring = ((index - (j - oplen)) << BAM_CIGAR_SHIFT)+ op;
            if((index - (j - oplen)) > 0){
                readseg* rsg = new_readseg(rln->segments->sequence,
                                           cstring,
                                           &refindx, &readindx);
                sladdhead(&readsegs, rsg); 
            }
            break;
        }
    }

    int rindex = r2;
    int nextindex = index;
    if(index >= q2){
        int offset = 0;
        for(i = 0, j = 0; i < numcigarops2; i++){
            int op    = (cigarstring2[i] & BAM_CIGAR_MASK);
            int oplen = (cigarstring2[i] >> BAM_CIGAR_SHIFT);
            if(op != BAM_CDEL) j += oplen;
            if(j <= q2){
            }else if((j > q2) && (j <= index)){
                if(op != BAM_CINS){
                    offset += oplen;
                    if((j - oplen) <= q2){
                        offset -= (q2 - (j - oplen));
                    }
                }
            }else if(j > index){
                if(op != BAM_CINS){
                    if((j - oplen) <= index){
                        offset += (index - (j - oplen));
                    }
                } 
            }
        }
        
        rindex = r2 + offset;
    }else{
        uint32_t cstring = ((q2 - index) << BAM_CIGAR_SHIFT)+ BAM_CINS;
        readseg* rsg = new_readseg(rln->segments->sequence,
                                   cstring,
                                   &refindx, &readindx);
        sladdhead(&readsegs, rsg);

        nextindex += (q2 - index);
    }     
    
    // is this an deletion?
    if(refindx < rindex){
        uint32_t cstring = ((rindex - refindx) << BAM_CIGAR_SHIFT)+ BAM_CDEL;
        readseg* rsg = new_readseg(rln->segments->sequence,
                                   cstring,
                                   &refindx, &readindx);
        sladdhead(&readsegs, rsg);
    }

    for(i = 0, j = 0; i < numcigarops2; i++){
        int op    = (cigarstring2[i] & BAM_CIGAR_MASK);
        int oplen = (cigarstring2[i] >> BAM_CIGAR_SHIFT);
        if(op != BAM_CDEL) j += oplen;

        if(j > nextindex){
            uint32_t cstring = ((j - nextindex) << BAM_CIGAR_SHIFT)+ op;
            readseg* rsg = new_readseg(rln->segments->sequence,
                                       cstring,
                                       &refindx, &readindx);
            sladdhead(&readsegs, rsg);
            i++;
            break;
        }
    }

    for(; i < numcigarops2; i++){
        readseg* rsg = new_readseg(rln->segments->sequence, 
                                   cigarstring2[i], 
                                   &refindx, &readindx);
        sladdhead(&readsegs, rsg);
    }

    free_readsegs(&rln->segments);
    slreverse(&readsegs);
    rln->segments = readsegs;
}
