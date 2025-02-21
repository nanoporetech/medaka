#ifndef _MEDAKA_COMMON_H
#define _MEDAKA_COMMON_H

#include <stdint.h>

#include "khash.h"
#include "kvec.h"


// bam tag used for datatypes
static const char datatype_tag[] = "DT";


/** Simple integer min/max
 * @param a
 * @param b
 *
 * @returns the min/max of a and b
 *
 */
static inline int max ( int a, int b ) { return a > b ? a : b; }
static inline int min ( int a, int b ) { return a < b ? a : b; }


/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg);


/** Reallocates memory with a message on failure.
 *
 *  @param ptr pointer to realloc.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xrealloc(void *ptr, size_t size, char* msg);


/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, int position, int length);


/** Format a uint32_t to a string
 *
 * @param value to format.
 * @param dst destination char.
 * @returns length of string.
 *
 */
size_t uint8_to_str(uint8_t value, char *dst);

/** Swap two strings
 *
 *  @param a first string
 *  @param b second string
 *  @returns a plp_data pointer.
 *
 */
void swap_strings(char** a, char** b);

/** Format an array values as a comma seperate string
 *
 * @param values integer input array
 * @param length size of input array
 * @param result output char buffer of size 4 * length * sizeof char
 * @returns void
 *
 * The output buffer size comes from:
 *    a single value is max 3 chars
 *    + 1 for comma (or \0 at end)
 */
void format_uint8_array(uint8_t* values, size_t length, char* result);

// Simple container for strings
typedef struct string_set {
    size_t n;
    char **strings;
} string_set;


/** Destroys a string set
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_string_set(string_set strings);


/** Retrieves contents of key-value tab delimited file.
 *
 *  @param fname input file path.
 *  @returns a string_set
 *
 *  The return value can be free'd with destroy_string_set.
 *  key-value pairs are stored sequentially in the string set
 *
 */
string_set read_key_value(char * fname);



#endif
