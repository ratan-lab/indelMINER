#ifndef MEMALLOC_H
#define MEMALLOC_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "errors.h"
#include "constants.h"

// Routines for memory allocation
void* ckalloc(const size_t size);
void* ckallocz(const size_t size);
void* ckrealloc(void * p, const size_t size);
void* ckreallocz(void * p, const size_t oldsize, const size_t newsize);

void ckfree(void* ptr);
#endif
