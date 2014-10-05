#ifndef GLOBALALIGN_H
#define GLOBALALIGN_H

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <inttypes.h>

#include "bam.h"

#include "errors.h"
#include "memalloc.h"
#include "files.h"
#include "readaln.h"

#define MININT -9999999
#define DIGIT 10.0

int ALIGN(char* A, 
          char* B, 
          int M, 
          int N, 
          int low, 
          int up, 
          int W[][128], 
          int G,
          int H, 
          int* S);

int DISPLAY(FILE* F,
            char* A,
            char* B,
            int M,
            int N, 
            int* S,
            int AP,
            int BP);

int fetch_cigar(char* A,
                char* B,
                int M,
                int N, 
                int* S,
                int AP,
                int BP,
                const int readlength,
                int* const pnumops,
                uint32_t** pcigar);

#endif
