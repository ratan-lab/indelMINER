#include "alignment.h"

extern uint klength;
extern uint numgaps;
extern uint maxdelsize;
extern bool debug_flag;
extern FILE* debug_file;
extern uint32_t seed_mask;
extern uint ethreshold;

static const uint base2bits[] = {
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,0,0,
    0,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,1,
    0,0,0,2,0,0,0,0,0,0,
    0,0,0,0,0,0,3,0,0,0,
};

// read the seeds of length klength from the sequence. Doing a seeds[kmer] will
// return the 1 based position of the seed. Use the result to query the table
// again
static uint* read_seeds(const char* const refsequence,
                        const uint start,
                        const uint stop,
                        uint** const pptrs)
{
    uint* ptrs = *pptrs;
    const char* sequence = refsequence + start;

    // this is the array which has the seeds
    uint* seeds = ckalloc(sizeof(uint) * pow(4,klength));
    memset(seeds, EMPTY, pow(4,klength) * sizeof(uint));

    uint32_t seed;
    uint i = 0, j = 0, k = 0;
    
    seed = 0;
    for(i = 0; i < (klength - 1); i++){
        seed = (seed << 2) + base2bits[(int)sequence[i]];
    }

    for(i = klength - 1; i < (stop - start); i++){
        seed = (seed << 2) + base2bits[(int)sequence[i]];

        if(seeds[seed] == EMPTY){
            seeds[seed] = i - (klength - 2);
        }else{      
            j = seeds[seed];
            while(j != 0){
                k = j;
                j = ptrs[j];
            }
            forceassert(ptrs[k] == EMPTY);
            ptrs[k] = i - (klength - 2);
        }

        seed = (seed & seed_mask);
    }

    return seeds;
}

static void bin_diagonals(int* const diagonals,
                          const uint numdiagonals,
                          const uint* const refseeds,
                          uint* const refptrs,
                          char* const readsequence,
                          const uint zstart,
                          const uint stop,
                          const uint* const readseeds,
                          uint* const readptrs)
{
    uint readlength = stop - zstart;
    char* sequence = readsequence + zstart;
    

    uint32_t seed = 0;
    uint i = 0, j = 0;
    uint indx = 0;

    seed = 0;
    for(i = 0; i < (klength - 1); i++){
        seed = (seed << 2) + base2bits[(int)sequence[i]];
    }

    int numunique = 0;
    for(i = klength - 1; i < readlength; i++){
        seed = (seed << 2) + base2bits[(int)sequence[i]];

        if(readptrs[readseeds[seed]] == 0) numunique += 1;
        if((refseeds[seed] != EMPTY) && (readptrs[readseeds[seed]] == 0)){
            j = refseeds[seed];

            while(j != 0){
                indx  = j - 1 - (i - klength + 1); 
                indx += readlength;
                indx -= klength;
                indx += 1;
                if(indx < numdiagonals){
                    diagonals[indx] += 1;
                    j = refptrs[j];
                }else{
                    break;
                }
            }    
        }

        seed = (seed & seed_mask);
    }

    int max = 0;
    for(i = 0; i < numdiagonals; i++){
        if(diagonals[i] > max){
            max = diagonals[i];
        }
    }
    if(TRUE == debug_flag){
        fprintf(debug_file, "Unique kmers in read: %d\n", numunique);
        fprintf(debug_file, "Max kmers on a diagonal: %d\n", max);
    }
}

static void bin_bands(int* const bands,
                      const int* const diagonals,
                      const uint numdiagonals)
{
    uint i = 0, j = 0;
    for(i = 0; i < (numdiagonals - numgaps); i++){
        for(j = i; j < (i + numgaps + 1); j++){
            bands[i] += diagonals[j];
        }
    }
}

static void select_band(int* const bands,
                        const int numdiagonals UNUSED,
                        const uint left,
                        const uint right,
                        const int anchor,
                        uint* const indx1,
                        uint* const indx2,
                        const int reflength UNUSED)
{
    // find the band that has the most hits
    int max  =  0;
    int dist = INT_MAX;
    uint indx = 0;
    uint i;

    for(i = left; i < right; i++){
        if(bands[i] > max){
            max  = bands[i];
            indx = i;
            dist = abs(anchor - i);              
        }else if(bands[i] == max){
            if(abs(anchor - i) < dist){
                max  = bands[i];
                indx = i;
                dist = abs(anchor - i);              
            }
        }
    }
 
    if(TRUE == debug_flag){
        fprintf(debug_file, "Index:Max kmer counts in a band:%d:%d\n",indx,max);
    }

    if(indx == INT_MAX){
        printf("whoa. gotta investigate\n");
    }

    *indx1 = indx;
    *indx2 = indx + numgaps;
}                        


static void print_alignments(char* const reference,
                             const uint r1,
                             const uint r2,
                             char* const readsequence,
                             const uint q1,
                             const uint q2,
                             const uint readlength,
                             int* const S)
{
    uint k;
    fprintf(debug_file, "Reference:\t");
    for(k = r1 - 10; k < r2 - 1; k++) 
        fprintf(debug_file, "%c", tolower(reference[k]));
    for(k = r1 - 1; k < r2; k++) 
        fprintf(debug_file,"%c",reference[k]);
    for(k = r2; k < r2 + 10; k++) 
    fprintf(debug_file, "%c", tolower(reference[k]));
    fprintf(debug_file, "\n");

    fprintf(debug_file, "ReadSeq:\t");
    for(k=0; k < readlength;k++) 
        fprintf(debug_file,"%c", readsequence[k]);
    fprintf(debug_file, "\n");

    DISPLAY(debug_file,
       readsequence - 1 + q1 - 1,
       reference -1 + r1 - 1, 
       q2 - q1 + 1,
       r2 - r1 + 1,
       S,  
       q1,
       r1);
    fprintf(debug_file, "\n");
}

static int count_matches(const uint32_t* cigarstring1,
                         const int numcigarops1,
                         const int q1,
                         const int q2,
                         const uint32_t* cigarstring2,
                         const int numcigarops2,
                         const int q3,
                         const int q4,
                         int* const pmm)
{
    forceassert(q1 == 0);
    forceassert(q3 == q2);

    int i, j;
    int matches = 0, mm = 0;
    for(i = 0, j = q1; i < numcigarops1; i++){
        int oplen = (cigarstring1[i] >> BAM_CIGAR_SHIFT);
        int op    = (cigarstring1[i] & BAM_CIGAR_MASK);
        
        if(op != BAM_CDEL) j += oplen;

        if(j < q2){
            if(op == BAM_CPMATCH){
                matches += oplen;
            }else if(op == BAM_CMMATCH){
                mm += oplen;
            }
        }

        if(j >= q2){
            if(op == BAM_CPMATCH){
                matches += (q2 - (j - oplen));
            }else if(op == BAM_CMMATCH){
                mm += (q2 - (j - oplen));
            }
            break;
        }
    }
//    printf("Stage1: %d matches, %d mismatches\n", matches, mm);

    for(i = 0, j = 0; i < numcigarops2; i++){
        int oplen = (cigarstring2[i] >> BAM_CIGAR_SHIFT);
        int op    = (cigarstring2[i] & BAM_CIGAR_MASK);
        
        if(op != BAM_CDEL) j += oplen;

        if(j >= q3){
            if(op == BAM_CPMATCH){
                matches += (j - q3);
            }else if(op == BAM_CMMATCH){
                mm += (j - q3);
            }
            i += 1;
            break;
        }
    }
//    printf("Stage2: %d matches, %d mismatches\n", matches, mm);

    for(; i < numcigarops2; i++){
        int oplen = (cigarstring2[i] >> BAM_CIGAR_SHIFT);
        int op    = (cigarstring2[i] & BAM_CIGAR_MASK);
        
        if(op != BAM_CDEL) j += oplen;

        if(j < q4){
            if(op == BAM_CPMATCH){
                matches += oplen;
            }else if(op == BAM_CMMATCH){
                mm += oplen;
            }
        }

        if(j >= q4){
            if(op == BAM_CPMATCH){
                matches += (q4 - (j - oplen));
            }else if(op == BAM_CMMATCH){
                mm += (q4 - (j - oplen));
            }
            break;
        }
    }
//    printf("Stage3: %d matches, %d mismatches\n", matches, mm);
    *pmm = mm;
    return matches;
}


static int find_best_del_candidate(const int q1, const int q2, 
                                    const uint32_t* cigarstring1, 
                                    const int numcigarops1, 
                                    const int q3, const int q4,
                                    const uint32_t* cigarstring2, 
                                    const int numcigarops2,
                                    const int readlength)
{
    forceassert(q1 == 0);
    forceassert(q3 <= q2);

    int i, matches, bestmatches = 0, mm, bestmm = INT_MAX, index = -1;
    for(i = q3; i <= q2; i++){
        matches = count_matches(cigarstring1, numcigarops1, q1, i,
                                cigarstring2, numcigarops2, i, q4, &mm);
//        printf("Matches: %d, Mismatches: %d\n", matches, mm);
        forceassert(matches <= readlength);    

        if(matches > bestmatches){
            bestmatches = matches;
            bestmm = mm;
            index = i;
        }else if((matches == bestmatches) && (mm < bestmm)){
            bestmatches = matches;
            bestmm = mm;
            index = i;
        }
        if((matches == readlength) && (mm == 0)) break;
    }
    forceassert(index != -1);
//    printf("%d : %d : %d\n", index, bestmatches, bestmm); 

    return index;
}


// attempt diagonal alignment of refseq:zstart1-end1 to readseq:zstart2-end2
static void attempt_band_alignment(char* const refseq,
                                   const uint zstart1,
                                   const uint end1,
                                   char* const readseq,
                                   const uint zstart2,
                                   const uint end2,
                                   const int low,
                                   const int up,
                                   int* const pr1,
                                   int* const pr2,
                                   int* const pq1,
                                   int* const pq2,
                                   int* const pnumcigarops,
                                   uint32_t** const pcigarstring)
{
    if(low > up) fprintf(stderr, "low: %d up: %d\n", low, up);
    forceassert(low <= up);
    int* S  = ckallocz(((end1 - zstart1) + (end2 - zstart2))*sizeof(int)); 
    int score = local_align(readseq+zstart2, end2 - zstart2, 
                            refseq + zstart1, end1 - zstart1, 
                            low, up, pq1, pr1, pq2, pr2, S);

    if(score <= 0){
        *pr1 = 0;
        *pr2 = 0;
        *pq1 = 0;
        *pq2 = 0;
        ckfree(S);
        return;
    }
 
    fetch_cigar(readseq + zstart2 + *pq1 - 2, refseq + zstart1 + *pr1 - 2, 
                *pq2 - *pq1 + 1, *pr2 - *pr1 + 1, S,
                *pq1, *pr1, end2 - zstart2, pnumcigarops, pcigarstring);   


    if(TRUE == debug_flag){
        fprintf(debug_file, "Aligned read:%d-%d to reference:%d-%d\n", *pq1+zstart2-1, *pq2+zstart2, *pr1+zstart1-1, *pr2+zstart1);
        print_alignments(refseq+zstart1, *pr1, *pr2, 
                         readseq+zstart2, *pq1, *pq2, end2 - zstart2, S);
    }

    *pr1 = *pr1 + zstart1-1;
    *pr2 = *pr2 + zstart1;
    *pq1 = *pq1 + zstart2 - 1;
    *pq2 = *pq2 + zstart2;

    ckfree(S);
}

static void find_best_band(char* const refseq,
                           const uint zstart1,
                           const uint end1,
                           const uint anchor,
                           char* const readseq,
                           const uint zstart2,
                           const uint end2,
                           int* const plow,
                           int* const pup)
{
    uint numdiagonals = ((end1 - zstart1) - (klength - 1))
                      + ((end2 - zstart2) - (klength - 1));
    forceassert(numdiagonals > numgaps);

    forceassert(end2 >= zstart2);
    if((end2 - zstart2) < klength){
        *plow = numdiagonals - 1;
        *pup  = numdiagonals - 1;
        return;
    }

    // read the seeds from the reference sequence
    uint* refptrs  = ckallocz((end1 - zstart1) * sizeof(uint));
    uint* refseeds = read_seeds(refseq, zstart1, end1, &refptrs);

    // read the seeds from the sequence of interest
    uint* readptrs  = ckallocz((end2 - zstart2) * sizeof(uint));
    uint* readseeds = read_seeds(readseq, zstart2, end2, &readptrs);

    int* diagonals = ckallocz(sizeof(int) * numdiagonals);
    int* bands = ckallocz(sizeof(int) * numdiagonals);

    bin_diagonals(diagonals,numdiagonals,
                  refseeds,refptrs,
                  readseq, zstart2, end2, readseeds,readptrs);
    bin_bands(bands, diagonals, numdiagonals);

    uint indx1, indx2;
    select_band(bands, numdiagonals, 0, numdiagonals, anchor - zstart1, 
                &indx1, &indx2, end1 - zstart1);

    
//    printf("%u %u %u %u\n", indx1, end2,zstart2, klength);
//    printf("%u %u %u %u\n", indx2, end2,zstart2, klength);

    *plow = indx1 - ((end2 - zstart2) -klength + 1);
    *pup  = indx2 - ((end2 - zstart2) -klength + 1);
    
    ckfree(refptrs);
    ckfree(refseeds);
    ckfree(readptrs);
    ckfree(readseeds);
    ckfree(diagonals);
    ckfree(bands);
}
 
static evidence* add_evidence_from_segment(readaln* const rln, 
                                           const readseg* const ignore)
{
    readseg* iter;
    evidence* allevidence = NULL;
    for(iter = rln->segments; iter; iter = iter->next){
        if(ignore && 
           (iter->oplen == ignore->oplen) && 
           (iter->op == ignore->op) && 
           (iter->start == ignore->start) && 
           (iter->end == ignore->end)) continue;
        if(iter->op == BAM_CDEL){
            if(debug_flag){
                fprintf(debug_file, "Adding %d base deletion\n",iter->oplen);
            }
            evidence* evdnc = new_evidence(rln, iter, DELETION, SPLIT_READ);
            sladdhead(&allevidence, evdnc);
        }else if(iter->op == BAM_CINS){
            if(debug_flag){
                fprintf(debug_file, "Adding %d base insertion\n",iter->oplen);
            }
            evidence* evdnc = new_evidence(rln, iter, INSERTION, SPLIT_READ);
            sladdhead(&allevidence, evdnc);
        }
    }
    free_readsegs(&rln->segments);
    return allevidence;
}

static inline void add_prefix_soft_clip(const int cliplength, 
                                       int* const pnumcigarops,
                                       uint32_t** const pcigarstring)
{
    if(cliplength == 0) return;

    int numcigarops = *pnumcigarops;
    uint32_t* cigarstring = *pcigarstring;

    // what is the first cigar operation
    int op    = (cigarstring[0] & BAM_CIGAR_MASK);
    int oplen = (cigarstring[0] >> BAM_CIGAR_SHIFT);
    
    if(op == BAM_CSOFT_CLIP){
        oplen += cliplength;
        cigarstring[0] = (oplen << 4) + BAM_CSOFT_CLIP;
    }else{
        numcigarops += 1;
        cigarstring = ckrealloc(*pcigarstring, numcigarops * sizeof(uint32_t));
        memmove(cigarstring + 1, cigarstring, (numcigarops-1)*sizeof(uint32_t));

        cigarstring[0] = (cliplength << 4) + BAM_CSOFT_CLIP;

        *pnumcigarops = numcigarops;
        *pcigarstring = cigarstring;
    }
}

static inline void add_suffix_soft_clip(const int cliplength, 
                                        int* const pnumcigarops,
                                        uint32_t** const pcigarstring)
{
    if(cliplength == 0) return;

    int numcigarops = *pnumcigarops;
    uint32_t* cigarstring = *pcigarstring;
    forceassert(numcigarops > 0);

    // what is the first cigar operation
    int op    = (cigarstring[numcigarops - 1] & BAM_CIGAR_MASK);
    int oplen = (cigarstring[numcigarops - 1] >> BAM_CIGAR_SHIFT);
    
    if(op == BAM_CSOFT_CLIP){
        oplen += cliplength;
        cigarstring[numcigarops - 1] = (oplen << 4) + BAM_CSOFT_CLIP;
    }else{
        numcigarops += 1;
        cigarstring = ckrealloc(*pcigarstring, numcigarops * sizeof(uint32_t));

        cigarstring[numcigarops - 1] = (cliplength << 4) + BAM_CSOFT_CLIP;

        *pnumcigarops = numcigarops;
        *pcigarstring = cigarstring;
    }
}


// left2             left1          position        right1        right2
//   |                |                |              |             |
//   V                V                V              V             V
// ----------------------------------------------------------------------      
static evidence* attempt_diagonal_alignments(readaln* const rln,
                                             char* const refseq,
                                             const int32_t left1,
                                             const int32_t right1,
                                             const int32_t left2,
                                             const int32_t right2,
                                             const int32_t anchor,
                                             char* const readseq)
{
    forceassert(anchor >= left1);      
    forceassert(anchor >= left2);      
    forceassert(anchor <= right1);
    forceassert(anchor <= right2);
    forceassert(left2 >= 0);
    forceassert(right2 > 0);

    uint readlength = strlen(readseq);
    int low1, up1;
    find_best_band(refseq, left1, right1, anchor,
                   readseq, 0, readlength, &low1, &up1);

    int numcigarops1;
    int r1,r2,q1,q2;
    uint32_t* cigarstring1 = ckallocz(sizeof(uint32_t));
    attempt_band_alignment(refseq, left1, right1, 
                           readseq, 0, readlength, 
                           low1, up1, &r1, &r2, &q1, &q2, 
                           &numcigarops1, &cigarstring1);

    if(q1 == q2){
        // the read did not align
        ckfree(cigarstring1);
        return NULL;
    }

    // this is over if I was able to align the whole read here
    if((q1 == 0) && (q2 == (int)readlength)){
        update_readsegs(rln, 
                        r1, cigarstring1, numcigarops1, 
                        readlength,
                        0, -1, NULL, 0);
        ckfree(cigarstring1);
        return add_evidence_from_segment(rln, NULL);
    }
    
    // identify the first non-match base in the alignment on the read
    int i, j;
    for(i = 0, j = 0; i < numcigarops1; i++){
        int op = (cigarstring1[i] & BAM_CIGAR_MASK);
        if((i == 0) && (op == BAM_CSOFT_CLIP)) continue;
        if(op != BAM_CPMATCH) break;
        j += (cigarstring1[i] >> BAM_CIGAR_SHIFT);
    } 
    uint f_nonmatch = j;
    for(i = numcigarops1 - 1, j = 0; i >= 0; i--){
        int op = (cigarstring1[i] & BAM_CIGAR_MASK);
        if((i == (numcigarops1 - 1)) && (op == BAM_CSOFT_CLIP)) continue;
        if(op != BAM_CPMATCH) break;
        j += (cigarstring1[i] >> BAM_CIGAR_SHIFT);
    } 
    uint l_nonmatch = j;

    int low2, up2;
    int numcigarops2;
    int r3,r4,q3,q4;
    uint32_t* cigarstring2 = ckallocz(sizeof(uint32_t));
    if(r1 > anchor){
        if(q1 == 0){
            forceassert(readlength > f_nonmatch);
            if(((readlength - f_nonmatch) < ethreshold) ||
               ((right2 - r1 - f_nonmatch) < ethreshold)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            find_best_band(refseq, r1 + f_nonmatch, right2, r1, 
                           readseq, f_nonmatch, readlength, &low2, &up2);

            attempt_band_alignment(refseq, r1 + f_nonmatch, right2, 
                                   readseq, f_nonmatch, readlength, 
                                   low2, up2, 
                                   &r3, &r4, &q3, &q4, 
                                   &numcigarops2, &cigarstring2);

            if((q4 != (int)readlength) || (q3 == q4)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            add_prefix_soft_clip(f_nonmatch,&numcigarops2, &cigarstring2);
        }else if(q2 == (int)readlength){
            forceassert(readlength > l_nonmatch);
            if(((readlength - l_nonmatch) < ethreshold) || 
               ((r2 - l_nonmatch - anchor) < ethreshold)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            find_best_band(refseq, anchor, r2 - l_nonmatch, r2, 
                           readseq, 0, readlength - l_nonmatch, &low2, &up2);

            attempt_band_alignment(refseq, anchor, r2 - l_nonmatch, 
                                   readseq, 0, readlength - l_nonmatch, 
                                   low2, up2, &r3, &r4, &q3, &q4, 
                                   &numcigarops2, &cigarstring2);

            if((q3 != 0) || (q3 == q4)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            add_suffix_soft_clip(l_nonmatch,&numcigarops2, &cigarstring2);
        }else{
            ckfree(cigarstring1);
            ckfree(cigarstring2);
            return NULL;
        }
    }else if(r1 < anchor){
        if(r2 >= anchor){
            ckfree(cigarstring1);
            ckfree(cigarstring2);
            return NULL;
        }

        if(q1 == 0){
            forceassert(readlength > f_nonmatch);
            if(((readlength - f_nonmatch) < ethreshold) || 
               ((anchor - r1 - f_nonmatch) < ethreshold)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            find_best_band(refseq, r1 + f_nonmatch, anchor, r1, 
                           readseq, f_nonmatch, readlength, &low2, &up2);
            attempt_band_alignment(refseq, r1 + f_nonmatch, anchor, 
                                   readseq, f_nonmatch, readlength, 
                                   low2, up2, 
                                   &r3, &r4, &q3, &q4, 
                                   &numcigarops2, &cigarstring2);          

            if((q4 != (int)readlength) || (q3 == q4)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            add_prefix_soft_clip(f_nonmatch,&numcigarops2, &cigarstring2);
        }else if(q2 == (int)readlength){
            forceassert(readlength > l_nonmatch);
            if(((readlength - l_nonmatch) < ethreshold) || 
               ((r2 - l_nonmatch - left2) < ethreshold)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            find_best_band(refseq, left2, r2 - l_nonmatch, r2, 
                           readseq, 0 , readlength - l_nonmatch, &low2, &up2);

            attempt_band_alignment(refseq, left2, r2 - l_nonmatch,
                                   readseq, 0 , readlength - l_nonmatch, 
                                   low2, up2, &r3, &r4, &q3, &q4, 
                                   &numcigarops2, &cigarstring2);

            if((q3 != 0) || (q3 == q4)){
                ckfree(cigarstring1);
                ckfree(cigarstring2);
                return NULL;
            }
            add_suffix_soft_clip(l_nonmatch,&numcigarops2, &cigarstring2);
        }else{
            ckfree(cigarstring1);
            ckfree(cigarstring2);
            return NULL;
        }
    }else{
        // is this the place to handle insertions?
        ckfree(cigarstring1);
        ckfree(cigarstring2);
        return NULL;
    }     

    // now find the best alignment out of the possible ones
    forceassert(q1 < q2);
    forceassert(q3 < q4);

    int index = -1;
    if((q1 > q3) && (q1 <= q4)){
//        printf("%s : %d %d %d %d\n", rln->qname, q1, q2, q3, q4);
        index = find_best_del_candidate(q3,q4,cigarstring2,numcigarops2,
                                   q1,q2,cigarstring1,numcigarops1, readlength);
        update_readsegs(rln,
                        r3, cigarstring2, numcigarops2, 
                        index,
                        q1, r1, cigarstring1, numcigarops1);
    }else if((q3 > q1) && (q3 <= q2)){
//        printf("%s: %d %d %d %d\n", rln->qname, q1, q2, q3, q4);
        index = find_best_del_candidate(q1,q2,cigarstring1,numcigarops1,
                                   q3,q4,cigarstring2,numcigarops2, readlength);
        update_readsegs(rln,
                        r1, cigarstring1, numcigarops1,
                        index,
                        q3, r3, cigarstring2, numcigarops2);
    }else if((q1 > q4) && (r1 == r4)){
        update_readsegs(rln,
                        r3, cigarstring2, numcigarops2,
                        q4,
                        q1, r1, cigarstring1, numcigarops1);
    }else if((q3 > q2) && (r2 == r3)){
        update_readsegs(rln,
                        r1, cigarstring1, numcigarops1,
                        q2,
                        q3, r3, cigarstring2, numcigarops2);   
    }else{
        ckfree(cigarstring1);
        ckfree(cigarstring2);
        return NULL;
    }

    ckfree(cigarstring1);
    ckfree(cigarstring2);
    return add_evidence_from_segment(rln, NULL);
}

// align the read rln such that one segment aligns within the expected PE
// distance of the mate, and the remaining segment aligns within "maxdelsize" of
// it. If that fails then check if this is an instance of an insertion.
evidence* attempt_pe_alignment(char** const sequences,
                               const int32_t tid,
                               const int32_t position,
                               const int* const range,
                               readaln* const rln)
{
    char* refsequence = sequences[tid];
    int reflength = strlen(refsequence);

    forceassert(range[0] <= range[1]);
    int32_t left1, right1, left2, right2, distance;
    distance = range[1];
    left1    = position >= distance ? position - distance : 0;
    right1   = reflength < (position + distance) 
             ? reflength : position + distance;    

    distance = range[1] + maxdelsize;
    left2    = position >= distance ? position - distance : 0;
    right2   = reflength < (position + distance) 
             ? reflength : position + distance;    
  
    if(TRUE == debug_flag){
        fprintf(debug_file, "Attempting to align read %s (%d)\n", 
        rln->qname, position);  
    }
    evidence* evdnc = attempt_diagonal_alignments(rln,
                                          refsequence, 
                                          left1, right1, 
                                          left2, right2,
                                          position,
                                          rln->segments->sequence);

    if(TRUE == debug_flag) fprintf(debug_file, "\n\n");    

    return evdnc;
}
