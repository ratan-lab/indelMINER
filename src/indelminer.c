#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#include "bam.h"

#include "asserts.h"
#include "constants.h"
#include "errors.h"
#include "resources.h"
#include "sequences.h"
#include "hashtable.h"

#include "shared.h"
#include "bamoperations.h"
#include "evidence.h"
#include "variant.h"
#include "readaln.h"
#include "graph.h"
#include "alignment.h"

#define indelminer_VERSION 0.1

#define READCHUNK 100000

// debug file
bool debug_flag;
FILE* debug_file;

// the kmer length and number of gaps allowed in an alignment
uint maxdelsize, maxpedelsize, klength, numgaps;

// the chromosomal region to analyze
char* chromosomal_region = NULL;

// used in alignment.c
uint32_t seed_mask;

// do not identify an indel in the last "ethreshold" bases of the reads
uint ethreshold;

// if we are only trying to check if a variant is supported then it should not
// be in the last these many base pairs
uint ethreshold_vcfcheck;

// do not call indels if the mapping quality of  the read is less than this
int qthreshold;

// in case of overlapping indels only call the one with the most support
bool call_all_indels; 

// maximum differences allowed in the alignments
uint maxdiffsallowed;

typedef struct alignmentdata_st
{
    struct alignmentdata_st* next;
    hashtable* insertlengths;
    hashtable* readpairs;
    evidence* allevidence;
    int numread;
    const char* outputformat;
    char** reference_sequence;
    bam_header_t* hin;
    const char* bam_name;
    uint minsupport;
    const char* vcfname;
    knownvariant* kvariants;
    const char* sample_name;
    bool ignore_threeprime;
}alignmentdata;

// look at the nodes from [n1,n2] and estimate the breakpoint range
static uint* find_breakpoint_range(node* const n1,
                                   node* const n2)
{
    int pos;

    // find the left breakpoint
    int left = -1;
    node* iter = n1;
    while(iter){
        evidence* evdnc = iter->val;
        pos = evdnc->b1;
        if((left == -1) || (pos > left)) left = pos;

        if(iter == n2) break;
        iter = iter->next;
    }

    // find the right breakpoint
    int right = -1;
    iter = n1;
    while(iter){
        evidence* evdnc = iter->val;
        pos = evdnc->b2;
        if((right == -1) || (pos <  right)) right = pos;

        if(iter == n2) break;
        iter = iter->next;
    }

    uint* range = ckalloc(2 * sizeof(uint));
    range[0] = (uint)left;
    range[1] = (uint)right;

    return range;
}

static variant* process_evidence(evidence** pallevidence,
                                 const int32_t tid,
                                 const int marker)
{
    variant* variations = NULL;

    slsort(pallevidence, sort_evidence);
    evidence* allevidence = *pallevidence;
//    fprintf(stderr, "pairs as evidence: %d\n", slcount(allevidence));
//    for(evidence* iter = allevidence; iter; iter = iter->next){
//        print_evidence(iter);
//    }


    // names of nodes
    char name[1024];
    int index = 1;

    graph* network = new_graph();
    evidence* iter = allevidence;
    for(; iter; iter = iter->next){
    // create a node from this pair/split read. Make the connections between
    // the nodes if they intersect each other.
        if(iter->b2 >= marker){
            break;
        }     
        sprintf(name, "node_%d", index++);
        node* n = new_node(name, iter);
        add_node(network, n);
    }
    
//    print_graph(network);

    // find the connected components in the graph
    find_connected_components(network);

    // lets find the variants
    sort_nodes(network);
    node* start = NULL;
    node* end   = NULL;
    node* niter = network->node_list;
    uint support;
    while(niter){
        start = niter;
        support = 0;
        while(niter && (niter->component == start->component)){
            end   = niter;
            niter = niter->next;
            support++;
        }

        // what is the breakpoint from these nodes
        uint* range = find_breakpoint_range(start, end);

        if(range[0] <= range[1]){
            // lets collect everything about this variant and add it to the list
            variant* variation = ckallocz(sizeof(variant));
            evidence* evdnc = start->val;
            variation->type = evdnc->variantclass;  
            variation->evdnctype = evdnc->type;   
            variation->tid = tid;
            variation->start = range[0];
            variation->stop = range[1];
            variation->support = support;
            
            uint i;
            node* iter;
            variation->evidence = ckalloc(support * sizeof(evidence*));
            for(iter = start, i = 0; iter != end->next; iter = iter->next){
                evidence* evdnc = iter->val;
                variation->evidence[i] = iter->val;
                i++;
                evdnc->isused = TRUE;
            }
            forceassert(i == support);

            sladdhead(&variations, variation);   
        }else{
            uint i;
            node* iter;
            for(iter = start, i = 0; iter != end->next; iter = iter->next){
                evidence* evdnc = iter->val;
                evdnc->isused = TRUE;
            }
        }

        ckfree(range);
    }

    free_graph(&network);

    return variations;
}

static int find_marker(const hashtable* const readpairs)
{
	bin* iter;
	bin* next;
	int max = INT_MAX;

	for(int i = 0; i < readpairs->size; i++){
		iter = readpairs->bins[i];
		while(iter){
			next = iter->next;

            evidence* evdnc = iter->val;
//            print_evidence(evdnc);
            if(evdnc->aln1->start < max){
                max = evdnc->aln1->start;
            }

			iter = next;
		}
	}

    return max;
}

static int check_for_mate(const bam1_t* alignment, void* data)
{
    readaln* rln = data;

    char* qname = bam1_qname(alignment);
    if(strcmp(qname, rln->qname) == 0){
        char index = (alignment->core.flag & 0x40) == 0x40 ? '1' : '2';
        if(rln->index == index){      
            readaln* tmp = new_readaln(alignment);
            rln->qname = tmp->qname;
            rln->tid = tmp->tid;
            rln->strand = tmp->strand;
            rln->index = tmp->index;
            rln->segments = tmp->segments;

            ckfree(tmp);
        }
    }
    return TRUE;
}

static readaln* find_mate_rln(const char* const bam_name,
                             const int32_t tid,
                             const int32_t pos,
                             const char index,
                             char* const qname)
{
    bamFile fp = bam_open(bam_name, "r");
    bam_index_t* fp_index = bam_index_load(bam_name);
     
    readaln* rln = ckallocz(sizeof(readaln));
    rln->qname = qname;
    rln->index = index;

    bam_fetch(fp, fp_index, tid, pos, pos + 1, rln, check_for_mate);
   
    bam_close(fp);
    bam_index_destroy(fp_index);

    // if I did not find the mate, then I should return NULL
    if (rln->segments == NULL) {
        ckfree(rln);
        return NULL;
    }
    return rln;
}


// go through the CIGAR string of this alignment and if there is an indel,
// then create the requisite evidence structure. 
static evidence* check_variants(readaln* const rln)
{
    evidence* allevdnc = NULL;

    readseg* iter = rln->segments;
    readseg* prev = NULL;

    uint rpos = 0;
    uint tpos = 0;
    while(iter){
        if((iter->op == BAM_CPMATCH) || (iter->op == BAM_CMMATCH) || (iter->op
== BAM_CMATCH) || (iter->op == BAM_CINS)) tpos += iter->oplen;
        iter = iter->next;
    }

    iter = rln->segments;
    while(iter){
        if(iter->op == BAM_CDEL){
            if((rpos > ethreshold_vcfcheck) && 
              ((tpos - rpos) > ethreshold_vcfcheck)){
                evidence* evdnc = new_evidence(rln, iter, DELETION, SPLIT_READ);
                sladdhead(&allevdnc, evdnc);
            }
        }else if(iter->op == BAM_CINS){
            if((rpos > ethreshold_vcfcheck) && 
              ((tpos - rpos) > ethreshold_vcfcheck)){
                evidence* evdnc = new_evidence(rln, iter, INSERTION, SPLIT_READ);
                sladdhead(&allevdnc, evdnc);
            }
        }else if((iter->op == BAM_CMATCH)  || 
                 (iter->op == BAM_CPMATCH) || 
                 (iter->op == BAM_CMMATCH)){
            rpos += iter->oplen;
        }else if(iter->op == BAM_CSOFT_CLIP){
            // this is either the beginning or the end of the read
            forceassert((prev == NULL) || (iter->next == NULL));
        }else{
            fatalf("unknown cigar op");
        }

        prev = iter;
        iter = iter->next;
    }      

//    if(slcount(allevdnc) > 1){
//        evidence* iter;
//        for(iter = allevdnc; iter; iter = iter->next) iter->isused = TRUE;
//        free_used_evidence(allevdnc);
//        allevdnc = NULL;
//    } 

    return allevdnc;
}

static int fetch_func(const bam1_t* alignment, void* data)
{
    alignmentdata* alndata = data;
    hashtable* insertlengths = alndata->insertlengths;
    hashtable* readpairs = alndata->readpairs;
    char** sequences = alndata->reference_sequence;
    bam_header_t* hin = alndata->hin;

    // ignore this read if it is a duplicate or secondary alignment
    if((alignment->core.flag & 0x100) == 0x100) return 0;
    if((alignment->core.flag & 0x200) == 0x200) return 0;
    if((alignment->core.flag & 0x400) == 0x400) return 0;

    bool is_aligned = (alignment->core.flag & 0x4) == 0;
    bool is_mate_aligned = (alignment->core.flag & 0x8) == 0;
    bool is_se = (alignment->core.flag & 0x1) == 0;
    bool is_proper_pair = (alignment->core.flag & 0x2) == 0x2;
    bool is_rc = (alignment->core.flag & 0x10) == 0x10;
    bool is_mate_rc = (alignment->core.flag & 0x20) == 0x20;

    // ignore reads that are not pe
    if(is_se) return 0;

    // ignore this pair if both align on different chromosomes
    if(is_aligned && 
       is_mate_aligned && 
      (alignment->core.tid != alignment->core.mtid)) return 0;

    // what is the expected insert lengths for this pair?
    uint8_t* rg = bam_aux_get(alignment, "RG");
    char* rgname = "generic";
    if(rg != NULL){
        rgname = bam_aux2Z(rg);
    }
    int* range = must_find_hashtable(insertlengths,
                                     rgname,
                                     strlen(rgname));

    // Cases to consider
    // a) One read maps and the other read is unmapped
    // b) Both reads map in proper pair
    // c) Both reads map, but not in a proper pair
    if(is_aligned && (!is_mate_aligned)){
        // we will deal with the mate when we get there
    }else if((!is_aligned) && (is_mate_aligned)){
        // attempt to align the unaligned read such that one segment of this read
        // aligns within PE distance of the aligned read, and the other segment
        // aligns within PE distance of that.
        uint8_t* pmmq = bam_aux_get(alignment, "MQ");
        int mmq;
        if (pmmq == NULL) {
            // in this case the aligner did not tell me about the mates
            // mapping quality. Lets assume it is the same as that of 
            // this sequence      
            mmq = alignment->core.qual;
        } else {
            forceassert(pmmq[0] == 'i');
            mmq = bam_aux2i(pmmq);
        }

        if(mmq >= qthreshold){
            readaln* rln = new_unaligned_readaln(alignment, mmq);
            if(!is_mate_rc){
                rln->segments->sequence = (char*)reverse_complement_string(
                                          (uchar*)rln->segments->sequence, 
                                          strlen(rln->segments->sequence));
                rln->strand = rln->strand == '+' ? '-' : '+';
            }
    
            evidence* evdnc = attempt_pe_alignment(sequences,
                                                   alignment->core.mtid,
                                                   alignment->core.mpos,
                                                   range,
                                                   rln);
            if(evdnc != NULL){
                sladdhead(&alndata->allevidence, evdnc);
            }else{
                // $$$ this could be a candidate for longer insertions???
                free_readsegs(&rln->segments);
            }
            ckfree(rln->qname);
            ckfree(rln);
        }
    }else if(is_aligned && is_mate_aligned && is_proper_pair){
        // I only have to worry if there is something interesting going on with
        // this read. Interesting implies that there is either a soft trim
        // portion on the 5' end or a deletion.
        readaln* rln = new_readaln(alignment);

        uint numcdels = 0;
        uint numcins = 0;
        uint numcsclip = 0;
        bool is_threeprime_clip = FALSE;

        for(readseg* iter = rln->segments; iter; iter = iter->next){
            if(iter->op == BAM_CDEL) numcdels += 1;
            if(iter->op == BAM_CINS) numcins += 1;
            if(iter->op == BAM_CSOFT_CLIP) numcsclip += 1;

            if((((rln->strand == '+') && (iter->next == NULL)) || 
                ((rln->strand == '-') && (iter == rln->segments))) && 
               (iter->op == BAM_CSOFT_CLIP)){
                is_threeprime_clip = TRUE; 
            }
        }

        uint numinteresting = numcdels + numcins + numcsclip;

        if(numinteresting == 0){
            free_readaln(&rln);
        }else if(numinteresting > 0){
            // here I have a chance of finding the indel using my method. 
            evidence* evdnc = NULL;
            evidence* bwaevdnc = NULL;
    
            if(((numcsclip==0)||((numcsclip==1)&&(is_threeprime_clip==TRUE))) && 
                (numcdels == 0) && 
                (numcins == 0)){
                // do nothing
            }else{
                // lets realign the read and try to call the indels
                uint8_t* pmmq = bam_aux_get(alignment, "MQ");
                int mmq;
                if (pmmq == NULL) {
                    // in this case the aligner did not tell me about the mates
                    // mapping quality. Lets assume it is the same as that of 
                    // this sequence
                    mmq = alignment->core.qual;
                } else {
                    forceassert(pmmq[0] == 'i');
                    mmq = bam_aux2i(pmmq);
                }
                if(mmq >= qthreshold){
                    bwaevdnc = check_variants(rln);

                    free_readaln(&rln);
                    rln = new_unaligned_readaln(alignment, alignment->core.qual);
        
                    if((is_rc && is_mate_rc) || ((!is_rc) && (!is_mate_rc))){
                    rln->segments->sequence = (char*)reverse_complement_string(
                                               (uchar*)rln->segments->sequence, 
                                               strlen(rln->segments->sequence));
                        rln->strand = rln->strand == '+' ? '-' : '+';
                    }
                    
                    evdnc = attempt_pe_alignment(sequences,
                                                 alignment->core.mtid,
                                                 alignment->core.mpos,
                                                 range,
                                                 rln);
                }
            }

            if(evdnc != NULL){
                evidence* last = sllast(evdnc);
                last->next = alndata->allevidence;
                alndata->allevidence = evdnc;
                
                if(bwaevdnc != NULL){
                    bwaevdnc->isused = TRUE;
                    free_used_evidence(bwaevdnc);
                }
            }else{
                // did the aligner find any variants? If yes, then I am going to
                // include those here
                if(bwaevdnc != NULL){
                    evidence* last = sllast(bwaevdnc);
                    last->next = alndata->allevidence;
                    alndata->allevidence = bwaevdnc;
                }
                free_readsegs(&rln->segments);
            }
            ckfree(rln->qname);
            ckfree(rln);
        }
    }else if(is_aligned && is_mate_aligned && (!is_proper_pair)){
        // the only ones I am interested in this case are the ones where the
        // orientation is proper and the distance between the reads is greater
        // than that expected for this readgroup.
        if((abs(alignment->core.isize) > range[1]) && 
           (abs(alignment->core.isize) < maxpedelsize) && 
           (is_rc != is_mate_rc)){
            if(alignment->core.pos < alignment->core.mpos){
                readaln* rln = new_readaln(alignment);
                evidence* evdnc = new_evidence(rln, 
                                               NULL,
                                               DELETION, 
                                               PAIRED_READ);
                add_hashtable(readpairs, 
                              rln->qname, 
                              alignment->core.l_qname, 
                              evdnc);
                ckfree(rln->qname);
                ckfree(rln);
            }else{
                // the mate for this read should already be in the hash. If it
                // is not there, then either there is an issue with this BAM, or
                // I am running several instances of indelMINER and the
                // mate is in some other partition and will be looked at by some
                // other instance of indelMINER. However this is the instance
                // that will process this evidence, so I should get the mate
                // read here, even if this slows things down...
                evidence* evdnc = find_value_hashtable(readpairs, 
                                                       bam1_qname(alignment), 
                                                       alignment->core.l_qname);
                if(evdnc == NULL){
                    readaln* rln = find_mate_rln(alndata->bam_name,
                                                 alignment->core.mtid,
                                                 alignment->core.mpos,
                               (alignment->core.flag & 0x40) == 0x40 ? '2' : '1',
                                                 bam1_qname(alignment));
                    // if I could not find the mate, then this is an issue with
                    // the BAM file. I should just ignore this evidence in that
                    // case after I print a warning to the user.
                    if (rln == NULL) {
                        goto next; 
                    } else {
                        evdnc = new_evidence(rln,
                                             NULL,
                                             DELETION,
                                             PAIRED_READ); 
                        // the mapping quality of the evidence is the minimum of 
                        // the mapping qualities of its reads
                        if(alignment->core.qual < evdnc->qual){
                            evdnc->qual = alignment->core.qual;
                        }

                        add_hashtable(readpairs,
                                      rln->qname,
                                      alignment->core.l_qname,
                                      evdnc);
                        ckfree(rln->qname);
                        ckfree(rln);
                    }
                }   

                readaln* rln = new_readaln(alignment);
                evdnc->aln3 = rln->segments;
                ckfree(rln->qname);
                ckfree(rln);
                evdnc->b1 = ((readseg*)sllast(evdnc->aln1))->end;
                evdnc->b2 = evdnc->aln3->start;
                evdnc->mindelsize = abs(alignment->core.isize) - range[1];
                evdnc->max = range[1];

                // I probably should not include this as evidence if both 
                // reads were mapped with a low mapping quality. In most 
                // cases I should be able to look at the MQ tag and tell 
                // the mapping quality of the mate(if they have used BWA to 
                // align the reads)
                int smq = alignment->core.qual;
                uint8_t* pmmq = bam_aux_get(alignment, "MQ");
                int mmq;
                if (pmmq == NULL) {
                    // in this case the aligner did not tell me about the mates
                    // mapping quality. Lets assume it is the same as that of 
                    // this sequence
                    mmq = alignment->core.qual;
                } else {
                    forceassert(pmmq[0] == 'i');
                    mmq = bam_aux2i(pmmq);
                }
                    
                if((smq >= qthreshold) || (mmq >= qthreshold)){
                    sladdhead(&alndata->allevidence, evdnc);
                }else{
                    evdnc->isused = TRUE;
                    free_used_evidence(evdnc);
                }
next:
                remove_hashtable_entry(readpairs, 
                                       bam1_qname(alignment), 
                                       alignment->core.l_qname);
            }
        }
    }

    if((++alndata->numread % READCHUNK) == 0){
        timestamp("Read %d reads", alndata->numread);

        // If I have a read pair where I have only come across one of the reads,
        // then I should not process variants beyond that point
        int marker = find_marker(alndata->readpairs);
        if(alignment->core.pos < marker) marker = alignment->core.pos;

        // look at the evidence collected so far. Collect the evidence into
        // clusters, where each cluster supports one variant.
        variant* variations = process_evidence(&alndata->allevidence,
                                               alignment->core.tid,
                                               marker);

        // sort the variants
        variations = sort_variants(variations);
        
        if(alndata->vcfname == NULL){
            // merge the variants if the variants share the same probable 
            // breakpoints
            merge_variants(&variations, sequences[alignment->core.tid], TRUE);

            // print the variants that you have found using the evidence above
            print_variants(&variations, 
                           sequences, 
                           alndata->bam_name,
                           hin, 
                           alndata->outputformat, 
                           alndata->minsupport,
                           maxdiffsallowed);
        }else{
            // merge the variants if they share the evidence type and same
            // probable breakpoints
            merge_variants(&variations, sequences[alignment->core.tid], FALSE);

            // print the variants that we already know about along with the tag
            // for this sample (in case this variant was also found in this
            // sample) 
            alndata->kvariants = print_knownvariants(alndata->kvariants, 
                                 variations, 
                                 alndata->sample_name, 
                                 alndata->bam_name,
                                 alndata->outputformat,
                                 hin);
        }
        fflush(stdout);

        // now relinquish the resources that are no longer required
        alndata->allevidence = free_used_evidence(alndata->allevidence);
        for(variant* iter = variations; iter; iter = iter->next){
            ckfree(iter->evidence);
        }
        slfreelist(&variations);
    }
    
    return 0;
}

static void free_evidence(bin* bin)
{
    evidence* evdnc = bin->val;
    free_readsegs(&evdnc->aln1);
}

// This program is used to call indels in a single sample using a MSA cleaned
// file of alignments. 
static void indelminer(const char* const configfile,
                       const char* const reference_fasta_name,
                       const char* const vcfname,
                       const char* const sample_name,
                       const char* const bam_name,
                       const char* const ignorefile,
                       const char* const outputformat,
                       const uint minsupport,
                       const bool ignore_threeprime)
{
    // if there is a configuration file that I can use. If not, then I will have
    // to estimate all of them
    hashtable* insertlengths = new_hashtable(4);
  
    bamFile fp = bam_open(bam_name, "r");
    bam_header_t* hin = bam_header_read(fp);
    uint* averagecoverage = ckallocz(hin->n_targets * sizeof(int));
    hashtable* id2chroms = new_hashtable(5);
    for(int32_t i = 0; i < hin->n_targets; i++){
        add_hashtable_int(id2chroms, 
                          hin->target_name[i], 
                          strlen(hin->target_name[i]), 
                          i);  
    }

    int chromid = -1;
    int chromstart = -1;
    int chromstop = -1;
    if(chromosomal_region != NULL){
        bam_parse_region(hin,chromosomal_region,&chromid,&chromstart,&chromstop);
    }

    if(configfile != NULL){
        // read the values from the config file
        read_configuration(configfile,id2chroms,insertlengths,averagecoverage);
    }else{
        // estimate the insert length range for the various read groups
        estimate_insertlengths(bam_name, insertlengths, chromid);
        // estimate the average depth coverage for the various chromosomes
        estimate_average_coverage(bam_name, chromid, averagecoverage);
    }
    fprintf(stderr,"\nRead-group\tMin-value\tMax-value\n");
    fprintf(stderr,"----------\t---------\t---------\n");
    print_insertlength_ranges(insertlengths);
    fprintf(stderr,"----------\t---------\t---------\n\n");
    fprintf(stderr,"\nChromosomeID\tMean-coverage\n");
    fprintf(stderr,"-------------\t-----------\n");
    for(int32_t i = 0; i < hin->n_targets; i++){
        fprintf(stderr, "%d\t%u\n", i, averagecoverage[i]);
    }
    fprintf(stderr,"-------------\t-----------\n\n");
    timestamp("Read insertlengths and average coverage for the BAM file");

    // read the reference sequence
    char** reference_sequence = read_reference(reference_fasta_name, 
                                               hin, 
                                chromid == -1 ? NULL :hin->target_name[chromid]);
    timestamp("Read the reference sequence");
    bam_header_destroy(hin);
    bam_close(fp);

    // lets print out some header information if this is a VCF output
    if(strncmp(outputformat, "vcf", 3) == 0){
        print_vcf_preamble(indelminer_VERSION);
    }
    // more stuff if I am tagging these indels
    if(vcfname != NULL){
        printf("##INFO=<ID=%s,Number=0,Type=Flag,Description=\"The variant is also present in this sample\">\n", sample_name);
    }

    fp = bam_open(bam_name, "r");
    hin = bam_header_read(fp);
    bam_index_t* fp_index = bam_index_load(bam_name);
    if(0 == fp_index) fatal("BAM indexing file is not available.");

    alignmentdata* alndata = ckallocz(sizeof(alignmentdata));
    alndata->readpairs = new_hashtable(20);
    alndata->allevidence = NULL;
    alndata->numread = 0;
    alndata->outputformat = outputformat;
    alndata->reference_sequence = reference_sequence;
    alndata->insertlengths = insertlengths;
    alndata->hin = hin;
    alndata->bam_name = bam_name;
    alndata->minsupport = minsupport;
    alndata->vcfname = vcfname;
    alndata->kvariants = NULL;
    alndata->sample_name = sample_name;
    alndata->bam_name = bam_name;
    alndata->ignore_threeprime = ignore_threeprime;

    for(int i = 0; i < hin->n_targets; i++){
        if((chromid != -1) && (i != chromid)) continue;

        // if I am just tagging these indels, then I should read all the indels
        // on this chromosome   
        knownvariant* kvariants = NULL;
        if(vcfname != NULL){
            kvariants = read_variants(vcfname, 
                                      reference_sequence, 
                                      i, hin->target_name[i]);

            if(kvariants == NULL) continue;
            alndata->kvariants = kvariants;
        }

        // would it be wise to ignore certain regions in this chromosome?
        if(ignorefile != NULL){
        }

        if(chromid == -1){
            bam_fetch(fp,fp_index,i,0,hin->target_len[i],alndata,fetch_func);
        }else{
            forceassert(i == chromid);
            bam_fetch(fp,fp_index,i,chromstart,chromstop,alndata,fetch_func);

        }

        // look at the evidence collected so far. Collect the evidence into
        // clusters, where each cluster supports one variant.
        variant* variations = process_evidence(&alndata->allevidence,i,INT_MAX);

        // sort the variants
        variations = sort_variants(variations);

        if(vcfname == NULL){
            // merge the variants if the variants share same probable breakpoints
            merge_variants(&variations, reference_sequence[i], TRUE);

            // print the variants that you have found using the evidence above
            print_variants(&variations,
                           reference_sequence,
                           alndata->bam_name,
                           hin,
                           outputformat,
                           minsupport,
                           maxdiffsallowed);
            fflush(stdout);
        }else{
            // merge the variants if they share the evidence type and same
            // probable breakpoints
            merge_variants(&variations, reference_sequence[i], FALSE);

            // print the variants that we already know about along with the tag
            // for this sample (in case this variant was also found in this
            // sample) 
            knownvariant* leftover = print_knownvariants(alndata->kvariants, 
                                                         variations, 
                                                         sample_name, 
                                                         alndata->bam_name,
                                                         alndata->outputformat,
                                                         hin);

            knownvariant* iter;
            for(iter = leftover; iter; iter = iter->next){
                print_vcf_line(iter, hin);
                if(iter->evdnctype == SPLIT_READ && 
                   is_indel_supported(iter, bam_name) == TRUE){
                    printf(";%s", sample_name);
                }
                printf("\n");
            }

            for(iter = kvariants; iter; iter = iter->next){
                ckfree(iter->reference);
                ckfree(iter->alternate);
                ckfree(iter->addntlinfo);
            }
            slfreelist(&kvariants);
        }

        // now relinquish the resources that are no longer required
        alndata->allevidence = free_used_evidence(alndata->allevidence);
        for(variant* iter = variations; iter; iter = iter->next){
            ckfree(iter->evidence);
        }
        slfreelist(&variations);
    }

    int indx;
    for(indx = 0; indx < hin->n_targets; indx++){
        ckfree(reference_sequence[indx]);
    }
    ckfree(reference_sequence);
    bam_header_destroy(hin);
    bam_close(fp);
    bam_index_destroy(fp_index);

    free_hashtable(&id2chroms);
    free_hashtable_completely(&insertlengths, NULL);
    free_hashtable_completely(&alndata->readpairs, free_evidence);
    slfreelist(&alndata->allevidence);
    ckfree(alndata);
    ckfree(averagecoverage);
}

// print the usage and exit
static void print_help(FILE* file)
{
    fprintf(file, "\n");
    fprintf(file, 
    "Program: indelminer (Call/Tag indels from a clean BAM file)\n");
    fprintf(file, "Version: %2.2f\n\n", indelminer_VERSION);
    fprintf(file, "Usage:\n");
    fprintf(file, "\tindelminer [options] ref.fa [indels.vcf] sample=aln.bam\n");
    fprintf(file, "where the options are\n");
    fprintf(file, "\t-h   print help and return\n");
    fprintf(file, "\n");
    fprintf(file, "\t-i, read the configuration from this file\n");
    fprintf(file," \t-c, only analyze this chromosomal region [ALL]\n");
    fprintf(file," \t-t, do not call indels on 3' regions of the read\n");
    fprintf(file," \t-q, do not call indels from reads with MQ < INT [10]\n");
    fprintf(file," \t-n, disallow indel within INT bp towards the ends [10]\n");
    fprintf(file," \t    We ignore the 3' soft-clipping, since that is where\n");
    fprintf(file," \t    we expect the low quality region on the reads\n");
    fprintf(file, "\t-a, in case of overlapping indels, call all of them\n");
    fprintf(file, "\t    Default is to call the indels with most support\n");
    fprintf(file, "\t-e, minimum support for an indel [2]\n");
    fprintf(file, "\t-o, output format. vcf/detailed [vcf]\n");
    fprintf(file, "\t-s, maximum size of deletion reported using split reads [1 kbp]\n");
    fprintf(file, "\t-p, maximum size of deletion reported using PE reads [1 Mbp]\n");
    fprintf(file, "\n");
    fprintf(file, "\t-k, length of the kmers to be used in alignments[6]\n");
    fprintf(file, "\t-g, number of gaps allowed in the alignments[0]\n");
    fprintf(file, "\t-f, number of differences allowed in an alignment[6]\n");
    fprintf(file, "\n");  
    fprintf(file, "Assumptions:\n");
    fprintf(file, "\tThe BAM file is coordinate sorted\n");
    fprintf(file, "\tUnless specified in a config file, insertlengths for\n");
    fprintf(file, "\treadgroups, as well as average coverage per chromosome\n");
    fprintf(file, "\tare estimated from the BAM file (which can be slow!!!),\n");
    fprintf(file, "\tas well as lead to false negatives as some of the PE\n");
    fprintf(file, "\tevidence which is accounted for in one sample,might not\n");
    fprintf(file, "\tbe accounted for in the other\n");

}

int main(int argc, char** argv)
{
    argv0 = "indelminer";
    debug_file = stderr;

    maxdelsize = 1000;
    maxpedelsize = 1000000;
    uint minsupport = 2;
    klength = 6;
    numgaps = 0;
    char* outputformat = "vcf";
    char* ignorefile = NULL;
    char* configfile = NULL;
    bool ignore_threeprime = FALSE;
    qthreshold = 10;
    ethreshold = 10;
    ethreshold_vcfcheck = ethreshold;
    call_all_indels = FALSE;
    maxdiffsallowed = 6;

    int c;
    while(1){
        c = getopt(argc,argv,"dl:hc:e:o:k:g:x:i:s:p:tn:q:af:");

        if (c == -1) break;
        switch (c){
            case 0:
                break;
            case 'd':
                debug_flag = TRUE;
                break;
            case 'l':
                debug_file = ckopen(optarg, "w");
                break;
            case 'h':
                print_help(stdout);
                return EXIT_SUCCESS;
            case 'c':
                chromosomal_region = optarg;
                break;
            case 'e':
                if(sscanf(optarg, "%u", &minsupport) != 1)
                    fatalf("incorrect option for -e: %s\n", optarg);
                break;
            case 'o':
                outputformat = optarg;
                break;
            case 'k':
                if(sscanf(optarg, "%u", &klength) != 1)
                    fatalf("incorrect option for -k: %s\n", optarg);
                break;
            case 'f':
                if(sscanf(optarg, "%u", &maxdiffsallowed) != 1)
                    fatalf("incorrect option for -f: %s\n", optarg);
                break;
            case 'g':
                if(sscanf(optarg, "%u", &numgaps) != 1)
                    fatalf("incorrect option for -g: %s\n", optarg);
                break;
            case 'x':
                ignorefile = optarg;
                break;
            case 'i':
                configfile = optarg;
                break;
            case 's':
                if(sscanf(optarg, "%d", &maxdelsize) != 1)
                    fatalf("incorrect option for -s: %s\n", optarg);
                break;
            case 'p':
                if(sscanf(optarg, "%d", &maxpedelsize) != 1)
                    fatalf("incorrect option for -p: %s\n", optarg);
                break;
            case 't':
                ignore_threeprime = TRUE;
                break;
            case 'n':
                if(sscanf(optarg, "%d", &ethreshold) != 1)
                    fatalf("incorrect option for -n: %s\n", optarg);
                if(ethreshold < klength) ethreshold = klength;
                ethreshold_vcfcheck = ethreshold;
                break;
            case 'q':
                if(sscanf(optarg, "%d", &qthreshold) != 1)
                    fatalf("incorrect option for -q: %s\n", optarg);
                break;
            case 'a':
                call_all_indels = TRUE;
                break;
            case '?':
                break;
            default:
                print_help(stderr);
                return EXIT_FAILURE;
        }
    }

    forceassert(maxdelsize > 0);
    forceassert((klength > 1) && (klength < 16));
    forceassert((strcmp(outputformat, "vcf") == 0) ||
                (strcmp(outputformat, "detailed") == 0)); 

    // do I have the correct arguments?
    if(argc == optind){  
        print_help(stderr);
        return EXIT_FAILURE;
    }
    forceassert(argc - optind > 1);

    // clock in 
    t0 = time(0);
 
    // print the arguments being used
    char* fasta_reference = argv[optind++];
    
    // if we have a VCF file, then we tag the indels only
    char* vcfname = NULL;

    char* ptr = argv[optind++];
    if(strchr(ptr, '=') == NULL){
        vcfname = ptr;
        minsupport = 1;
        outputformat = "vcf";
    }      

    // lets find the indels here
    if(vcfname != NULL) ptr = argv[optind++];
    char* samplename = ptr;
    while(*ptr != '=') ptr++;
    forceassert(*ptr == '=');
    *ptr = 0;
    char* samplealn = ++ptr;  

    fprintf(stderr, "Reference fasta file: %s\n", fasta_reference);
    fprintf(stderr, "Chromosomal region  : %s\n", 
        chromosomal_region == NULL ? "ALL" : chromosomal_region);
    fprintf(stderr, "BAM file            : %s\n", samplealn);
    if(vcfname != NULL){
        fprintf(stderr, "VCF file            : %s\n", vcfname);
    }
    
    seed_mask = pow(2,2*(klength-1)) - 1;
 
    // if we are checking for support of indels, then I should be lenient
    if(vcfname != NULL) ethreshold_vcfcheck = 0;

    indelminer(configfile, 
               fasta_reference, 
               vcfname, 
               samplename,
               samplealn, 
               ignorefile,
               outputformat, 
               minsupport,
               ignore_threeprime);

    // print the stats used by me
    print_usage();

    return EXIT_SUCCESS;
}
