#include "bamoperations.h"

static void print_name_range(bin* bin)
{
    int* val = bin->val;
    fprintf(stderr, "%s\t%d\t%d\n", bin->name, *val, *(val+1));
}

// print the names of the readgroups and the insert ranges
void print_insertlength_ranges(hashtable* const hash)
{
	func_hashtable(hash, print_name_range);
}

static int fetch_ins_func(const bam1_t* alignment, void* data)
{
    hashtable* insertlengths = data;

    // is this alignment for a read that was part of a proper pair?
    if((alignment->core.flag & 0x1) == 0x0) return 0;
    if((alignment->core.flag & 0x4) == 0x4) return 0;
    if((alignment->core.flag & 0x2) == 0x0) return 0;    

    // ignore duplicates and secondary alignments
    if((alignment->core.flag & 0x100) == 0x100) return 0; 
    if((alignment->core.flag & 0x200) == 0x200) return 0; 
    if((alignment->core.flag & 0x400) == 0x400) return 0; 

//    // only consider that ones which map uniquely
//    if(alignment->core.qual <= 0) return 0;

    // we can safely ignore the sides where isize < 0. 
    if(alignment->core.isize < 0) return 0;

    // get the read group information
    uint8_t* rg =  bam_aux_get(alignment, "RG");

    char* rgname = "generic";
    if(rg != NULL){
        forceassert(rg[0] == 'Z');
        rgname = bam_aux2Z(rg);
    }
    
    int32_t isize  = alignment->core.isize;
    if((alignment->core.mpos - alignment->core.pos) < 0) return 0;
    if(isize < (alignment->core.mpos - alignment->core.pos)) return 0;

    if(lookup_hashtable(insertlengths, rgname, strlen(rgname)) == NULL){
        int32_t* range = ckalloc(2*sizeof(int32_t));
        range[0] = isize;
        range[1] = isize;
        add_hashtable(insertlengths,rgname,strlen(rgname),range);
    }else{
        int32_t* range=must_find_hashtable(insertlengths,rgname,strlen(rgname));
        if(range[0] > isize) range[0] = isize;
        if(range[1] < isize) range[1] = isize;
    }

    return 0;
}

void estimate_insertlengths(const char* const bam_name,
                            hashtable* const insertlengths,
                            const int chromid)
{
    // open the BAM file
    bamFile fp = bam_open(bam_name, "r");

    // read the BAM header
    bam_header_t* hin = bam_header_read(fp);

    // load the index to the BAM file
    bam_index_t* idx = bam_index_load(bam_name);
    if(0 == idx) fatal("BAM indexing file is not available.");

    int32_t i; 
    for(i = 0; i < hin->n_targets; i++){
        if((chromid != -1) && (i != chromid)) continue;

        bam_fetch(fp,idx,i,0,hin->target_len[i],insertlengths,fetch_ins_func);
    }

    bam_index_destroy(idx);
    bam_header_destroy(hin);
    bam_close(fp);       
}

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
    bam_plbuf_t *buf = (bam_plbuf_t*)data;
    bam_plbuf_push(b, buf);
    return 0;
}

// callback function for bam_plbuf_init()
static int cov_pileup_func(uint32_t tid UNUSED, 
                           uint32_t pos UNUSED, 
                           int n, 
                           const bam_pileup1_t *pl UNUSED, 
                           void *data)
{
    covsum* cvsum = (covsum*) data;
    cvsum->covsum += n;
    if(n > 0) cvsum->numcov += 1;

    return 0;
}

// observed coverage on the chromosome(s)
void estimate_average_coverage(const char* const bam_name,
                               const int chromid,
                               uint* const meancovs)
{
    // open the BAM file
    bamFile fp = bam_open(bam_name, "r");

    // read the header
    bam_header_t* hin = bam_header_read(fp);

    bam_index_t *idx;
    bam_plbuf_t *buf;

    // load the index
    idx = bam_index_load(bam_name);
    if(0 == idx) fatal("BAM indexing file is not available.");

    // iterate through the chromosomes and find the average coverage on each
    int32_t i;
    for(i = 0; i < hin->n_targets; i++){
        if((chromid != -1) && (i != chromid)) continue;

        // fill in the coverage buffer
        covsum* cvsum = ckallocz(sizeof(covsum));
        buf = bam_plbuf_init(cov_pileup_func, cvsum);
        bam_fetch(fp, idx, i, 0, hin->target_len[i], buf, fetch_func);
        bam_plbuf_push(0, buf);
        bam_plbuf_destroy(buf);  

        meancovs[i] = floor(cvsum->covsum * 1.0 / cvsum->numcov);
        ckfree(cvsum);
    }

    bam_index_destroy(idx);
    bam_header_destroy(hin);
    bam_close(fp);
}
