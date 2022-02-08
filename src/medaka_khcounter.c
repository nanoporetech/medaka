// Wrap khash to make it more consise to use

#include <stdio.h>
#include "khash.h"
#include "medaka_khcounter.h"
#include "medaka_common.h"

/* Implementation of a counter of strings (increasing only)
 *
 * khash_t(KH_COUNTER) *h = kh_init(KH_COUNTER);
 * kh_counter_increment(h, "one");
 * kh_counter_increment(h, "two");
 * kh_counter_increment(h, "two");
 * kh_counter_add(h, "three", 2);
 * kh_counter_increment(h, "three", 1);
 * kh_counter_print(h);
 * kh_counter_destroy(h);
 *
 */

const int KH_COUNTER = 99;
KHASH_MAP_INIT_STR(KH_COUNTER, size_t)

size_t kh_counter_val(khash_t(KH_COUNTER) *hash, char *key) {
    khiter_t k = kh_get(KH_COUNTER, hash, key);
    size_t val = k != kh_end(hash) ? kh_val(hash, k) : 0;
    return val;
}

size_t kh_counter_add(khash_t(KH_COUNTER) *hash, char *key, size_t val) {
    // note: key is copied so no need for caller to hold on to it
    int ret;
    khiter_t k = kh_put(KH_COUNTER, hash, key, &ret);
    if (ret == 1) { // new key
        kh_key(hash, k) = strdup(key);
        kh_value(hash, k) = val;
    } else if (ret == 0) {  // exists
        // get value and add
        size_t cur = kh_counter_val(hash, key);
        kh_value(hash, k) = cur + val;
    } else {
        // shouldnt get here - previously deleted key
    }
    return ret;
}

size_t kh_counter_increment(khash_t(KH_COUNTER) *hash, char *key) {
    kh_counter_add(hash, key, 1);
}

kh_counter_stats_t kh_counter_stats(khash_t(KH_COUNTER) *hash) {
    kh_counter_stats_t stats = { .sum=0, .max=0 };
    for (khiter_t k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            const char *key = kh_key(hash, k);
            size_t val = kh_value(hash, k);
            stats.sum += val;
            stats.max = max(stats.max, val);
        }
    }
    return stats;
}

void kh_counter_destroy(khash_t(KH_COUNTER) *hash) {
    for (khiter_t k = 0; k < kh_end(hash); ++k){
        if (kh_exist(hash, k)) {
            free((char*) kh_key(hash, k));
        }
    }
    kh_destroy(KH_COUNTER, hash);
}

void kh_counter_print(khash_t(KH_COUNTER) *hash) {
    for (khiter_t k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            const char *key = kh_key(hash, k);
            size_t val = kh_value(hash, k);
            printf("%s -> %lu\n", key, val);
        }
    }
    kh_counter_stats_t stats = kh_counter_stats(hash);
    printf("max: %lu, sum: %lu\n", stats.max, stats.sum);
}


// Demonstrates usage
int main(int argc, char *argv[]) {
    khash_t(KH_COUNTER) *h = kh_init(KH_COUNTER);
    kh_counter_increment(h, "one");
    kh_counter_increment(h, "two");
    kh_counter_increment(h, "two");
    kh_counter_add(h, "three", 2);
    kh_counter_increment(h, "three");
    kh_counter_print(h);
    kh_counter_destroy(h);
}
