#include "hashtable.h"

/*allocate a new hash table*/
hashtable* new_hashtable(const int po2size)
{
	int max_size = po2size;
	if(po2size > 24){
		fprintf(stderr,"hash po2size should not exceed 24, using 24\n");
		max_size = 24;
	}

	hashtable* hashtable = ckallocz(sizeof(struct hashtable_st));
	hashtable->po2size = max_size;
	hashtable->size = (1 << hashtable->po2size);
	hashtable->mask = hashtable->size - 1;

	hashtable->bins  = ckallocz(hashtable->size*sizeof(struct bin_st*));

	return hashtable;
}

/*add the following substring of length len to the hashtable. Return the 
 * bin corresponding to it*/
bin* add_hashtable(hashtable* const hash,  /*the hashtable*/
				   const char* const name, /*the string*/
				   const int length, 	   /*use 'length' characters of name*/
				   void* val)				/*the value*/
{
    pre(name != NULL);

	uint32_t index =  superfasthash(name, length);
	index = index & hash->mask;
	
	char* hname = ckalloc(length+1);
	memcpy(hname, name, length);
	hname[length] = '\0';

	bin* bin = ckallocz(sizeof(struct bin_st));
	bin->name = hname;
	bin->val = val;
	if(hash->bins[index]){
		hash->collisions++;
	}
	bin->next = hash->bins[index];
	hash->bins[index] = bin;
	hash->elcount++;
	
    post(bin != NULL);
	return bin;
}

bin* add_hashtable_int(hashtable* const hash,
                       const char* const name,
                       const int length,
                       const int val)
{
    return add_hashtable(hash, name, length, val + NULL);
} 

/*look up the hash table to see an entry corresponding to the string of length
 * len exists. Return the bin if it does exist*/
void* lookup_hashtable(hashtable* const hash, 
				       const char* const name, 
				       const int len)
{
    pre(name != NULL);

	uint32_t index = superfasthash(name, len);
	index = index & hash->mask;
	
	bin* iter = NULL;
	bin* bin = NULL;
	if(hash->bins[index] != NULL){
		for(iter = hash->bins[index];iter; iter = iter->next){
			if(strncmp(iter->name, name, len) == 0){
				bin = iter;
			}
		}
	}
	return bin;
}

/*look in the hashtable and return the value associated with the string of
 * length len. Returns NULL is the bin doesnt exist*/	
void* find_value_hashtable(hashtable* const hash,
						   const char* const name,
						   const int len)
{
    pre(name != NULL);

	bin* bin = lookup_hashtable(hash, name, len);
	if(bin == NULL){
		return NULL;
	}
	return bin->val;
}

/* find a value in the hashtable or die if the key does not exist in the hash*/
void* must_find_hashtable(hashtable* const hash,
						  const char* const name,
						  const int len)
{
    pre(name != NULL);

	bin* bin = lookup_hashtable(hash, name, len);
	if(bin == NULL){
		fatalf("did not find %s in the hash", name);
	}
	return bin->val;
}

int must_find_hashtable_int(hashtable* const hash,
                            const char* const name,
                            const int len)
{
    return must_find_hashtable(hash,name,len) - NULL;
}

/*remove this entry from the hashtable and return the value associated with it*/
void* remove_hashtable_entry(hashtable* const hash,
							 const char* const name,
							 const int len)
{
    pre(name != NULL);

	uint32_t index = superfasthash(name, len);
	index = index & hash->mask;
	
	bin* iter = NULL;
	bin* bin = NULL;
	if(hash->bins[index] != NULL){
		for(iter = hash->bins[index];iter; iter = iter->next){
			if(strncmp(iter->name, name, len) == 0){
				bin = iter;
				break;
			}
		}
	}

	void* val = NULL;

	if(bin){
		val = bin->val;
		slremove(&hash->bins[index], bin);
		ckfree(bin->name);
		ckfree(bin);
		hash->elcount--;
	}
	return val;
}

/* apply this function to every non-null bin in the hashtable */
void func_hashtable(hashtable* const hash, 
					void (*func)(bin* bin))
{
    pre(func != NULL);

	int i;
	bin* iter;
	bin* next;
	
	for(i = 0; i < hash->size; i++){
		iter = hash->bins[i];
		while(iter){
			next = iter->next;
			func(iter);
			iter = next;
		}
	}
}

static void free_hash(bin* bin)
{
	ckfree(bin->name);
	ckfree((char*)bin);
}

static void comp_free_hash(bin* bin)
{
	ckfree(bin->name);
	ckfree(bin->val);
	ckfree((char*)bin);
}

/*free all the resources allocated to the hash table*/
void free_hashtable(hashtable** const phash)
{
	if(*phash != NULL){
		func_hashtable(*phash, free_hash);

		hashtable* hash = *phash;
		ckfree(hash->bins);
		ckfree(hash);
		*phash = NULL;
	}
}

/*free all the resources allocated to the hash table. Also for each bin, free
 * the value in the bin.*/
void free_hashtable_completely(hashtable** const phash, void (*func)(bin* bin))
{	
	if(*phash != NULL){
        if(func != NULL) func_hashtable(*phash, func);
		func_hashtable(*phash, comp_free_hash);

		hashtable* hash = *phash;
		ckfree(hash->bins);
		ckfree(hash);
		*phash = NULL;
	}
}

static void print_names(bin* bin)
{
    pre(bin != NULL);
	printf("%s\n", bin->name);
}

/* for debugging only. print all the names in the hashtable */
void print_hashtable(hashtable* const hash)
{
	func_hashtable(hash, print_names);
}

