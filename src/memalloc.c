#include "memalloc.h"

// make a call to malloc, and panic if the allocation request fails. If enabled
// at compile time a special  memory  debugging  capabilities replaces the
// normal version of ckalloc, which aids in detecting memory overwrites and
// leaks. If the compiler flag MEM_DEBUG is set, then the debugging capabilities
// are invoked
void* ckalloc(const size_t size)
{
    void* ptr;

    if((ptr = malloc(size)) == NULL){
        fprintf(stderr, "Error in allocating %zd bytes.\n", size);
        perror(NULL);
        exit(MEM_ALLOC_FAILURE);
    }
    return ptr;
}

// fill the allocated memory with '0' before returning a pointer to it
void* ckallocz(const size_t size)
{
    void* ptr = ckalloc(size);

    if(memset(ptr, 0, size) != ptr){    
        fprintf(stderr, "Error in initializing the area\n");
        exit(MEM_SET_FAILURE);
    }

    return ptr;
}

/*reallocate the memory to the given size*/
void *ckrealloc(void * p, const size_t size)
{
    p = p ? realloc(p, size) : malloc(size);
    if (!p){
        fatal("ckrealloc failed");
    }
    return p;
}

// reallocate the memory to given size. If the new memory size is more than the
// old memory size, then initialize the new memory to 0
void* ckreallocz(void * p, const size_t oldsize, const size_t newsize)
{
    p = ckrealloc(p, newsize);

    // if the oldsize was more or the same, this is just equivalent to ckrealloc
    if(newsize <= oldsize) return p;

    memset(p + oldsize, 0, newsize - oldsize);

    return p;
}

// free the memory pointed to by this pointer. If MEM_DEBUG is set, then check
// the guards.
void ckfree(void* ptr)
{
    if(ptr) free(ptr);
}
