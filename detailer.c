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

static void detailer(const char* const vcfname,
                     const char* const alnname)
{
    size_t n = 1;
    ssize_t num_read;
    char* fptr = ckalloc(n+1);


    FILE* fp = ckopen(vcfname, "r");
    bamFile bfp = bam_open(alnname, "r");
    bam_index_t* bfp_index = bam_index_load(alnname);
    bam_header_t* bfp_hin = bam_header_read(bfp);

    // we need a chrom -> tid mapping
    hashtable* id2chroms = new_hashtable(5);
    for(int32_t i = 0; i < bfp_hin->n_targets; i++){
        add_hashtable_int(id2chroms, 
                          bfp_hin->target_name[i], 
                          strlen(bfp_hin->target_name[i]), 
                          i);  
    }

    char chrom[128];
    uint position;
    char information[1024];

    while((num_read = getline(&fptr, &n, fp)) != -1){
        if(fptr[0] == '#'){
            printf("%s", fptr);
            continue;
        }

        if(sscanf(fptr, "%s %u %*c %*s %*s %*c %*c %s\n", 
                        chrom, &position, information) != 3){
            fatalf("error in reading the VCF line %s", fptr);
        }

        char* ptr = strstr(information, "END=");
        if(ptr == NULL) fatalf("no end coordinates for the variant");
        uint stop = atoi(ptr+4);

        uint meancov = calculate_cov_params(&bfp, bfp_index, 
                        must_find_hashtable_int(id2chroms,chrom,strlen(chrom)),
                       position - 1, 
                       stop);
        
        uint mqrms;
        uint mq0;
        calculate_mq_params(&bfp, bfp_index,
            must_find_hashtable_int(id2chroms, chrom, strlen(chrom)),
            position - 1,
            stop,
            &mqrms,
            &mq0);

        fptr[num_read - 1] = 0;
        printf("%s;MQ=%u;MQ0=%u;DP=%u\n", fptr, mqrms, mq0, meancov);

    }       

    bam_index_destroy(bfp_index);
    bam_header_destroy(bfp_hin);
    bam_close(bfp);

    fclose(fp);
}

int main(int argc, char** argv)
{
    argv0 = "detailer";

    if(argc != 3){
        fatalf("detailer indels.vcf alignments.bam");
    }

    detailer(argv[1], argv[2]);

    return EXIT_SUCCESS;
}
