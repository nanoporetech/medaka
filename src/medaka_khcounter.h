#ifndef _MEDAKA_KHCOUNTER_H
#define _MEDAKA_KHCOUNTER_H

#include "khash.h"

typedef struct kh_counter_stats_t {
    size_t sum;
    size_t max;
} kh_counter_stats_t;

KHASH_MAP_INIT_STR(KH_COUNTER, size_t)

// create a counter
static inline khash_t(KH_COUNTER) *kh_counter_init() {
    khash_t(KH_COUNTER) *h = kh_init(KH_COUNTER);
    return h;
}

// Get a value from a counter 
size_t kh_counter_val(khash_t(KH_COUNTER) *hash, char *key);

// Clean up a counter
void kh_counter_destroy(khash_t(KH_COUNTER) *hash);

// Increment a counter by one
size_t kh_counter_increment(khash_t(KH_COUNTER) *hash, char *key);

// Increment a counter by a given amount
size_t kh_counter_add(khash_t(KH_COUNTER) *hash, char *key, size_t val);

// Retrieve statistics on counter
kh_counter_stats_t kh_counter_stats(khash_t(KH_COUNTER) *hash);

// Print contents of a counter
void kh_counter_print(khash_t(KH_COUNTER) *hash);

#endif
