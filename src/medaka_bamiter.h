#ifndef MBAMITER_H
#define MBAMITER_H

#include <stdbool.h>
#include "htslib/sam.h"

// parameters for bam iteration
typedef struct {
    htsFile *fp;
    bam_hdr_t *hdr;
    hts_itr_t *iter;
    int min_mapQ;
    char tag_name[2];
    int tag_value;
    bool keep_missing;
} mplp_data;

// iterator for reading bam
int read_bam(void *data, bam1_t *b);

#endif
