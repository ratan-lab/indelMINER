#ifndef CONSTANTS_H
#define CONSTANTS_H

// exit codes
#define MEM_ALLOC_FAILURE 2
#define MEM_SET_FAILURE 3

// bool constants
typedef enum bool_em
{
    FALSE = 0,
    TRUE  = 1
}bool;

// some systems do not have uint
typedef unsigned int uint;

// ignore warnings for these local variables
# define UNUSED __attribute__((unused))

#endif
