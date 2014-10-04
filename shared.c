#include "shared.h"

extern char* argv0;

void read_configuration(const char* const filename, 
                        hashtable* const id2chroms,
                        hashtable* insertlengths,
                        uint* const averagecoverage)
{
    size_t n = 1;
    char* fptr = ckalloc(n + 1);

    char rgname[128]; // name of the library
    uint minsize; // minimum size in a proper pair
    uint maxsize; // maximum size in a proper pair

    char chrom[128]; // chromosome
    uint meancov;   // mean coverage on that chromosome

    FILE* fp = ckopen(filename,"r");
    while(getline(&fptr, &n, fp) != -1){
        if(strncmp(fptr, "IL", 2) == 0){
            if(sscanf(fptr, "IL %s %u %u\n", rgname, &minsize, &maxsize) != 3){
                fatalf("error in reading the insert length range: %s", fptr);
            }
            int32_t* range = ckalloc(2*sizeof(int32_t));
            range[0] = minsize;
            range[1] = maxsize;
            add_hashtable(insertlengths,rgname,strlen(rgname),range);
        } else if(averagecoverage && (strncmp(fptr, "RC", 2) == 0)){
            if(sscanf(fptr, "RC %s %u\n", chrom, &meancov) != 2){
                fatalf("error in reading the mean coverage: %s", fptr);
            }
            averagecoverage[must_find_hashtable_int(id2chroms, 
                                                    chrom, 
                                                    strlen(chrom))] = meancov;
        } else {
            fatalf("unknown tag in configuration: %s", fptr);
        }
    }
    
    ckfree(fptr);
    fclose(fp);
}

char** read_reference(const char* const reference_fasta_name,
                      const bam_header_t* const hin,
                      const char* const chromosomes)
{
    char** sequences = ckallocz(hin->n_targets * sizeof(char*));

    // go through the reference fasta sequence and read all the sequences
    sequence* sp = read_fasta_sequence(reference_fasta_name);
    
    int indx = 0;
    while(sp){
        if((chromosomes != NULL) &&
           (strcmp(hin->target_name[indx],chromosomes) != 0)){
            sp = get_next_sequence(sp);
            indx += 1;
            continue;
        }

        char* seq = ckallocz(sp->slen + 1);

        uint i;
        for(i = 0; i < sp->slen; i++){
            sp->sequence[i] = toupper(sp->sequence[i]);
        }

        memcpy(seq, sp->sequence, sp->slen);
        sequences[indx] = seq;

        sp = get_next_sequence(sp);
        indx += 1;
    }
    forceassert(indx == hin->n_targets);

    close_fasta_sequence(sp);

    return sequences;
}

void print_vcf_preamble(const float majorversion)
{
    printf("##fileformat=VCFv4.1\n");
    printf("##%sVersion=%2.2f\n", argv0, majorversion);

    // format of the info fields
    // ##INFO=<ID=ID,Number=number,Type=type,Description=”description”>
    // Possible Types for INFO fields are: Integer, Float, Flag, Character, and
    // String. The Number entry is an Integer that describes the number of
    // values that can be included with the INFO field.
    printf("##INFO=<ID=INSERTION,Number=0,Type=Flag,Description=\"Indicates that the variant is an insertion.\">\n");
    printf("##INFO=<ID=DELETION,Number=0,Type=Flag,Description=\"Indicates that the variant is a deletion.\">\n");
    printf("##INFO=<ID=SPLIT_READ,Number=0,Type=Flag,Description=\"Indicates that at least one split read supports this variant.\">\n");
    printf("##INFO=<ID=PAIRED_READ,Number=0,Type=Flag,Description=\"Indicates that at least one PE read supports this variant.\">\n");
    printf("##INFO=<ID=COMPOSITE,Number=0,Type=Flag,Description=\"Indicates that at least one split read and at least one PE read supports this variant.\">\n");
    printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of reads supporting the variant\">\n");
    printf("##INFO=<ID=END,Number=1,Type=Integer,Description=\"end position of the variant described in this record\">\n");
    printf("##INFO=<ID=BP_END,Number=1,Type=Integer,Description=\"possible 3' end of the breakpoint described in this record\">\n");
    printf("##INFO=<ID=NFS,Number=1,Type=Integer,Description=\"Number of reads supporting the variant on the forward strand\">\n");
    printf("##INFO=<ID=NRS,Number=1,Type=Integer,Description=\"Number of reads supporting the variant on the forward strand\">\n");
    printf("##INFO=<ID=UTAILS,Number=1,Type=Integer,Description=\"The number of unique tail distances in supporting reads for this variant\">\n");
    printf("##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS mapping quality of the reads covering the breakpoints\">\n");
    printf("##INFO=<ID=MQ30,Number=1,Type=Integer,Description=\"Number of reads with mapping quality greater than or equal to 30, covering the breakpoints\">\n");
    printf("##INFO=<ID=DF,Number=1,Type=Integer,Description=\"Average number of other differences on reads supporting the reported variant\">\n");
    printf("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Average read depth across the breakpoints\">\n");
    printf("##INFO=<ID=BF,Number=2,Type=Integer,Description=\"Flanks from the split read or pair best sorrounding the variant\">\n");


    // FILTERs that have been applied to the data should be described as:
    // ##FILTER=<ID=ID,Description=”description”>

    // Genotype fields specified in the FORMAT field should be described as
    // follows:
    // ##FORMAT=<ID=ID,Number=number,Type=type,Description=”description”>
    // Possible Types for FORMAT fields are: Integer, Float, Character, and 
    // String (this field is otherwise defined precisely as the INFO field)

    // Symbolic alternate alleles for imprecise structural variants:
    // ##ALT=<ID=type,Description=description>
    // 
    // The ID field indicates the type of structural variant, and can be a
    // colon-separated list of types and subtypes. ID values are case sensitive
    // strings and may not contain whitespace or angle brackets. 
    // The first level type must be  one of the following:
    // 
    // DEL Deletion relative to the reference
    // INS Insertion of novel sequence relative to the reference
    // DUP Region of elevated copy number relative to the reference
    // INV Inversion of reference sequence
    // CNV Copy number variable region (may be both deletion and duplication)
    //
    // The CNV category should not be used when a more specific category can be
    // applied. Reserved subtypes include:
    // 
    // DUP:TANDEM Tandem duplication
    // DEL:ME Deletion of mobile element relative to the reference
    // INS:ME Insertion of a mobile element relative to the reference
    // In addition, it is highly recommended (but not required) that the header
    // include tags describing the reference and contigs backing the data 
    // contained in the file.  These tags are based on the SQ field from the 
    // SAM spec; all tags are optional (see the VCF example above).
    // 
    // Breakpoint assemblies for structural variations may use an external file:
    // 
    // ##assembly=url
}

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
    bam_plbuf_t *buf = (bam_plbuf_t*)data;
    bam_plbuf_push(b, buf);
    return 0;
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid UNUSED, 
                       uint32_t pos, 
                       int n, 
                       const bam_pileup1_t *pl,
                       void *data)
{
    covdata* cvdt = (covdata*)data;
    int i;
    if((pos >= cvdt->begin) && (pos < cvdt->end)){
        for(i = 0; i < n; i++){
            if(pl[i].is_del == 0){
                cvdt->coverage[pos-cvdt->begin] += 1;
            }
        }
    }
    return 0;
}

uint calculate_cov_params(const char* const bam_name,
                          const int32_t tid,
                          const int32_t start,
                          const int32_t stop)
{
    bamFile fp = bam_open(bam_name, "r");
    bam_index_t* fp_index = bam_index_load(bam_name);

    bam_plbuf_t *buf;

    covdata* cvdt = ckallocz(sizeof(covdata));
    cvdt->tid = tid;
    cvdt->begin = start;
    cvdt->end   = stop;
    cvdt->coverage = ckallocz((cvdt->end - cvdt->begin) * sizeof(uint32_t));
    
    buf = bam_plbuf_init(pileup_func, cvdt);
    bam_fetch(fp, fp_index, tid, start, stop, buf, fetch_func);
    bam_plbuf_push(0, buf);
    bam_plbuf_destroy(buf);  

    // calculate the mean coverage in the region of the putative deletion
    uint i, covsum;
    for(i = 0, covsum = 0; i < (cvdt->end - cvdt->begin); i++){
        covsum += cvdt->coverage[i];
    }
  
    uint avgcov = floor(covsum * 1.0/(cvdt->end - cvdt->begin));
    ckfree(cvdt->coverage);
    ckfree(cvdt);

    bam_close(fp);   
    bam_index_destroy(fp_index);
    return avgcov;
}

static int calcrms(const bam1_t* const alignment, void* data)
{
    mapqcov* mqds = *(mapqcov**)data;

    // ignore this read if it is a duplicate or secondary alignment
    if((alignment->core.flag & 0x100) == 0x100) return TRUE;
    if((alignment->core.flag & 0x200) == 0x200) return TRUE;
    if((alignment->core.flag & 0x400) == 0x400) return TRUE;

    mqds->mqsqsum += (alignment->core.qual * alignment->core.qual);
    mqds->numcov  += 1;
    if(alignment->core.qual == 0) mqds->numcov0 += 1;

    return TRUE;
}

void calculate_mq_params(bamFile* const pfp,
                         bam_index_t* const fp_index,
                         const int32_t tid, 
                         const int32_t start,
                         const int32_t stop, 
                         uint* const mqrms, 
                         uint* const mq0)
{
    bamFile fp = *pfp;  
 
    mapqcov* mqds = ckallocz(sizeof(mapqcov));

    bam_fetch(fp,
              fp_index,
              tid, start, stop,
              &mqds,
              calcrms);

    *mqrms = sqrt(mqds->mqsqsum*1.0/mqds->numcov);
    *mq0 = mqds->numcov0;

    ckfree(mqds);
   
    return;
}
