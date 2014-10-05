#ifndef LOCALALIGN_H
#define LOCALALIGN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "asserts.h"
#include "constants.h"
#include "errors.h"
#include "memalloc.h"
#include "resources.h"
#include "globalalign.h"

int  local_align(char* seq1,
                 const int seq1len,
                 char* seq2,
                 const int seq2len,
                 const int indx1,
                 const int indx2,
                 int* const psi,
                 int* const psj,
                 int* const pei,
                 int* const pej,
                 int* const S);

#endif
