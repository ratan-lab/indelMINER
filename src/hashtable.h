#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "asserts.h"
#include "memalloc.h"
#include "slinklist.h"
#include "superfasthash.h"

typedef struct bin_st
{
	struct bin_st* next;
	char* name;
	void* val;
}bin;

typedef struct hashtable_st 
{	
	int po2size;	/*size of the hash table in power of 2*/
	int size;		/*size of the hash table*/
	int mask;		/*the mask to be applied*/

	int elcount;	/*number of elements in the hash table*/
	int collisions;	/*number of collisions in this hash table*/

	bin** bins;		/*individual bins*/
}hashtable;

/*allocate a new hash table*/
hashtable* new_hashtable(const int po2size);

/*add the following substring of length len to the hashtable. Return the 
 * bin corresponding to it*/
bin* add_hashtable(hashtable* const hash,  /*the hashtable*/
				   const char* const name, /*the string*/
				   const int length, 	   /*use 'length' characters of name*/
				   void* val);				/*the value*/

bin* add_hashtable_int(hashtable* const hash, /*the hashtable*/
				      const char* const name, /*the string*/
				      const int length, 	  /*use 'length' characters name*/
				      const int val);	      /*the value*/

/*look up the hash table to see an entry corresponding to the string of length
 * len exists. Return the bin if it does exist*/
void* lookup_hashtable(hashtable* const hash, 
				       const char* const name, 
				       const int len);

/*look in the hashtable and return the value associated with the string of
 * length len. Returns NULL is the bin doesnt exist*/	
void* find_value_hashtable(hashtable* const hash,
						   const char* const name,
						   const int len);

/* find a value in the hashtable or die if the key does not exist in the hash*/
void* must_find_hashtable(hashtable* const hash,
						  const char* const name,
						  const int len);

/* find a value in the hashtable or die if the key does not exist in the hash*/
int must_find_hashtable_int(hashtable* const hash,
				            const char* const name,
						    const int len);

/*remove this entry from the hashtable and return the value associated with it*/
void* remove_hashtable_entry(hashtable* const hash,
							 const char* const name,
							 const int len);

/*free all the resources allocated to the hash table*/
void free_hashtable(hashtable** const phash);

/*free all the resources allocated to the hash table. Also for each bin, free
 * the value in the bin. Before doing that use the provided function to free the
 * values */
void free_hashtable_completely(hashtable** const phash, void (*func)(bin* b));

/* apply this function to every non-null bin in the hashtable */
void func_hashtable(hashtable* const hash, void (*func)(bin* b));


/* for debugging only. print all the names in the hashtable */
void print_hashtable(hashtable* const hash);
#endif
