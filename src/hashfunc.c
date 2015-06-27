#include "hashfunc.h"

/*
**	FNV-1a
**	https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function
*/
uint32_t hashfunc(const char* data, int len){
	static const uint32_t PRIME = 16777619;
	static const uint32_t OFFSET = 2166136261;
	int i;
	uint32_t result = OFFSET;
	for(i=len-1; i>=0; i--){
		result ^= data[i];
		result *= PRIME;
	}
	return result;
}

/*
**	DJB2
**	http://www.cse.yorku.ca/~oz/hash.html
*/
uint32_t djb2(const char* data, int len){
	int i;
	uint32_t result = 5381;
	for(i=len-1; i>=0; i--){
		result += (result << 5) + data[i];
	}
	return result;
}
