#ifndef RESOURCES_H
#define RESOURCES_H

// for calculating and printing resources
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include<stdlib.h>
#include<stdio.h>
#include<stdarg.h>
#include<string.h>
#include<time.h>

#include "asserts.h"

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))

// clock in time
extern time_t t0;

// timestamp
void timestamp(const char* const fmt, ...);

// print the resources used
void print_usage();

#endif
