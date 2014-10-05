#include "strings.h"

char* copy_partial_string(const char* const string, size_t numchars)
{
    char* copy = ckallocz(numchars + 1);
    memcpy(copy, string, numchars);
    return copy;
}


char* copy_string(const char* const string)
{
    return copy_partial_string(string, strlen(string));
}
