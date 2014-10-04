#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <inttypes.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>

#include "constants.h"
#include "memalloc.h"
#include "evidence.h"
#include "variant.h"
#include "localalign.h"

#define EMPTY 0

#define TORIGHT 0
#define TOLEFT  1
#define DONTCARE 2

evidence* attempt_pe_alignment(char** const sequences,
                               const int32_t tid,
                               const int32_t position,
                               const int* const range,
                               readaln* const rln);
#endif
