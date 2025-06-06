#define _GNU_SOURCE
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "htslib/sam.h"

#include "khash.h"
#include "kvec.h"
#include "medaka_bamiter.h"
#include "medaka_common.h"
#include "medaka_counts.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table

// For recording bad reads and skipping processing
KHASH_SET_INIT_STR(BADREADS)


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param buffer_cols number of pileup columns.
 *  @param num_dtypes number of datatypes in pileup.
 *  @param num_homop maximum homopolymer length to consider.
 *  @param featlen number of features per column.
 *  @param fixed_size if not zero data matrix is allocated as fixed_size * n_cols, ignoring other arguments
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t buffer_cols, size_t num_dtypes, size_t num_homop, size_t featlen, size_t fixed_size) {
    assert(buffer_cols >= n_cols);
    plp_data data = xalloc(1, sizeof(_plp_data), "plp_data");
    data->buffer_cols = buffer_cols;
    data->num_dtypes = num_dtypes;
    data->num_homop = num_homop;
    data->featlen = featlen;
    data->n_cols = n_cols;
    if (fixed_size != 0) {
        assert(buffer_cols == n_cols);
        data->matrix = xalloc(fixed_size * n_cols, sizeof(size_t), "matrix");
    } else {
        data->matrix = xalloc(featlen * num_dtypes * buffer_cols * num_homop, sizeof(size_t), "matrix");
    }
    data->major = xalloc(buffer_cols, sizeof(size_t), "major");
    data->minor = xalloc(buffer_cols, sizeof(size_t), "minor");
    return data;
}


/** Enlarge the internal buffers of a pileup data structure.
 *
 *  @param pileup a plp_data pointer.
 *  @param buffer_cols number of pileup columns for which to allocate memory
 *
 */
void enlarge_plp_data(plp_data pileup, size_t buffer_cols) {
    assert(buffer_cols > pileup->buffer_cols);
    size_t old_size = pileup->featlen * pileup->num_dtypes * pileup->num_homop * pileup->buffer_cols;
    size_t new_size = pileup->featlen * pileup->num_dtypes * pileup->num_homop * buffer_cols;

    pileup->matrix = xrealloc(pileup->matrix, new_size * sizeof(size_t), "matrix");
    pileup->major = xrealloc(pileup->major, buffer_cols * sizeof(size_t), "major");
    pileup->minor = xrealloc(pileup->minor, buffer_cols * sizeof(size_t), "minor");
    // zero out new part of matrix
    for (size_t i = old_size; i < new_size; ++i) {
        pileup->matrix[i] = 0;
    }
    pileup->buffer_cols = buffer_cols;
}


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data) {
    free(data->matrix);
    free(data->major);
    free(data->minor);
    free(data);
}


/** Prints a pileup data structure.
 *
 *  @param pileup a pileup structure.
 *  @param dtypes datatype prefix strings.
 *  @returns void
 *
 */
void print_pileup_data(plp_data pileup, char *dtypes[]){
    fprintf(stdout, "pos\tins\t");
    if (pileup->num_dtypes > 1) {  //TODO header for multiple dtypes and num_homop > 1
        for (size_t i = 0; i < pileup->num_dtypes; ++i) {
            for (size_t j = 0; j < pileup->featlen; ++j){
                fprintf(stdout, "%s.%c\t", dtypes[i], plp_bases[j]);
            }
        }
    } else {
        for (size_t k = 0; k < pileup->num_homop; ++k) {
            for (size_t j = 0; j < pileup->featlen; ++j){
                fprintf(stdout, "%c.%lu\t", plp_bases[j], k+1);
            }
        }
    }
    fprintf(stdout, "depth\n");
    for (size_t j = 0; j < pileup->n_cols; ++j) {
        int s = 0;
        fprintf(stdout, "%zu\t%zu\t", pileup->major[j], pileup->minor[j]);
        for (size_t i = 0; i < pileup->num_dtypes * pileup->featlen * pileup->num_homop; ++i){
            const size_t c = pileup->matrix[j * pileup->num_dtypes * pileup->featlen * pileup->num_homop + i];
            s += c;
            fprintf(stdout, "%zu\t", c);
        }
        fprintf(stdout, "%d\n", s);
    }
}


float* _get_weibull_scores(const bam_pileup1_t *p, const size_t indel, const size_t num_homop, khash_t(BADREADS) *bad_reads) {
    // Create homopolymer scores using Weibull shape and scale parameters.
    // If prerequisite sam tags are not present an array of zero counts is returned.
    float* fraction_counts = xalloc(num_homop, sizeof(float), "weibull_counts");
    static const char* wtags[] = {"WL", "WK"};  // scale, shape
    double wtag_vals[2] = {0.0, 0.0};
    for (size_t i=0; i < 2; ++i) {
        uint8_t *tag = bam_aux_get(p->b, wtags[i]);
        if (tag == NULL) {
            char* read_id = bam_get_qname(p->b);
            khiter_t k;
            k = kh_get(BADREADS, bad_reads, read_id);
            if (k == kh_end(bad_reads)) { // a new bad read
                int putret;
                k = kh_put(BADREADS, bad_reads, read_id, &putret);
                fprintf(stderr, "Failed to retrieve Weibull parameter tag '%s' for read %s.\n",
                        wtags[i], read_id);
            }
            return fraction_counts;
        }
        uint32_t taglen = bam_auxB_len(tag);
        if (p->qpos + indel >= taglen) {
            fprintf(stderr, "%s tag was out of range for %s position %lu. taglen: %i\n",
                    wtags[i], bam_get_qname(p->b), p->qpos + indel, taglen);
            return fraction_counts;
        }
        wtag_vals[i] = bam_auxB2f(tag, p->qpos + indel);
    }

    // found tags, fill in values
    float scale = wtag_vals[0];  //wl
    float shape = wtag_vals[1];  //wk
    for (size_t x = 1; x < num_homop + 1; ++x) {
        float a = pow((x-1) / scale, shape);
        float b = pow(x / scale, shape);
        fraction_counts[x-1] = fmax(0.0, -exp(-a) * expm1(a - b));
    }
    return fraction_counts;
}


/** Generates medaka-style feature data in a region of a bam.
 *
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
 *  @param num_dtypes number of datatypes in bam.
 *  @param dtypes prefixes on query names indicating datatype.
 *  @param num_homop maximum homopolymer length to consider.
 *  @param tag_name by which to filter alignments.
 *  @param tag_value by which to filter data.
 *  @param keep_missing alignments which do not have tag.
 *  @param weibull_summation use predefined bam tags to perform homopolymer partial counts.
 *  @returns a pileup data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 *  If num_dtypes is 1, dtypes should be NULL; all reads in the bam will be
 *  treated equally. If num_dtypes is not 1, dtypes should be an array of
 *  strings, these strings being prefixes of query names of reads within the
 *  bam file. Any read not matching the prefixes will cause exit(1).
 *
 *  If tag_name is not NULL alignments are filtered by the (integer) tag value.
 *  When tag_name is given the behaviour for alignments without the tag is
 *  determined by keep_missing.
 *
 */
plp_data calculate_pileup(
        const char *region, const bam_fset* bam_set, size_t num_dtypes, char *dtypes[],
        size_t num_homop, const char tag_name[2], const int tag_value, const bool keep_missing,
        bool weibull_summation, const char *read_group, const int min_mapQ) {
    if (num_dtypes == 1 && dtypes != NULL) {
        fprintf(stderr, "Recieved invalid num_dtypes and dtypes args.\n");
        exit(1);
    }
    const size_t dtype_featlen = featlen * num_dtypes * num_homop;

    // extract `chr`:`start`-`end` from `region`
    //   (start is one-based and end-inclusive),
    //   hts_parse_reg below sets return value to point
    //   at ":", copy the input then set ":" to null terminator
    //   to get `chr`.
    int start, end;
    char *chr = xalloc(strlen(region) + 1, sizeof(char), "chr");
    strcpy(chr, region);
    char *reg_chr = (char *) hts_parse_reg(chr, &start, &end);
    // start and end now zero-based end exclusive
    if (reg_chr) {
        *reg_chr = '\0';
    } else {
        fprintf(stderr, "Failed to parse region: '%s'.\n", region);
    }

    // open bam etc.
    // this is all now deferred to the caller
    htsFile *fp = bam_set->fp;
    hts_idx_t *idx = bam_set->idx;
    sam_hdr_t *hdr = bam_set->hdr;

    // setup bam interator
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr; data->iter = bam_itr_querys(idx, hdr, region);
    data->min_mapQ = min_mapQ; memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    data->keep_missing = keep_missing; data->read_group = read_group;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);
    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    // allocate output assuming one insertion per ref position
    int n_cols = 0;
    size_t buffer_cols = 2 * (end - start);
    plp_data pileup = create_plp_data(n_cols, buffer_cols, num_dtypes, num_homop, featlen, 0);

    // get counts
    size_t major_col = 0;  // index into `pileup` corresponding to pos
    n_cols = 0;            // number of processed columns (including insertions)
    khash_t(BADREADS) *no_rle_tags = kh_init(BADREADS);  // maintain set of reads without rle tags

    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {
        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if (pos < start) continue;
        if (pos >= end) break;
        n_cols++;

        // find maximum insert
        int max_ins = 0;
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->indel > 0 && max_ins < p->indel) max_ins = p->indel;
        }

        // reallocate output if necessary
        if (n_cols + max_ins > pileup->buffer_cols) {
            float cols_per_pos = (float) (n_cols + max_ins) / (1 + pos - start);
            // max_ins can dominate so add at least that
            buffer_cols = max_ins + max(2 * pileup->buffer_cols, (int) cols_per_pos * (end - start));
            enlarge_plp_data(pileup, buffer_cols);
        }

        // set major/minor position indexes, minors hold ins
        for (int i = 0; i <= max_ins; ++i ) {
            pileup->major[major_col / dtype_featlen + i] = pos;
            pileup->minor[major_col / dtype_featlen + i] = i;
        }

        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) continue;

            // find to which datatype the read belongs
            int dtype = 0;
            if (num_dtypes > 1){
                bool failed = false;
                char *tag_val;
                uint8_t *tag = bam_aux_get(p->b, datatype_tag);
                if (tag == NULL) { // tag isn't present
                    failed = true;
                } else {
                    tag_val = bam_aux2Z(tag);
                    failed = errno == EINVAL;
                }
                if (!failed) {
                    bool found = false;
                    for (dtype = 0; dtype < num_dtypes; ++dtype){
                        if(strcmp(dtypes[dtype], tag_val) == 0) {
                            found = true;
                            break;
                        }
                    }
                    failed = !found;
                }
                if (failed) {
                    fprintf(stderr, "Datatype not found for %s.\n", bam_get_qname(p->b));
                    exit(1);
                }
            }

            int base_i, min_minor = 0;
            int max_minor = p->indel > 0 ? p->indel : 0;
            if (p->is_del) {
                // deletions are kept in the first layer of qscore stratification, if any
                int qstrat = 0;
                base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
                //base = plp_bases[base_i];
                pileup->matrix[major_col + featlen * dtype * num_homop + featlen * qstrat + base_i] += 1;
                min_minor = 1; // in case there is also an indel, skip the major position
            }
            // loop over any query bases at or inserted after pos
            int query_pos_offset = 0;
            for (int minor = min_minor; minor <= max_minor; ++minor, ++query_pos_offset){
                int base_j = bam1_seqi(bam1_seq(p->b), p->qpos + query_pos_offset);
                //base = seq_nt16_str[base_j];
                if bam_is_rev(p->b) {
                    base_j += 16;
                }

                base_i = num2countbase[base_j];
                if (base_i != -1) {  // not an ambiguity code
                    size_t partial_index =
                        major_col + dtype_featlen * minor              // skip to column
                        + featlen * dtype * num_homop                  // skip to datatype
                        //+ featlen * qstrat                           // skip to qstrat/homop
                        + base_i;                                      // the base

                    if (weibull_summation) {
                        float* fraction_counts = _get_weibull_scores(p, query_pos_offset, num_homop, no_rle_tags);
                        for (size_t qstrat = 0; qstrat < num_homop; ++qstrat) {
                            static const int scale = 10000;
                            pileup->matrix[partial_index + featlen * qstrat] +=
                                scale * fraction_counts[qstrat];
                        }
                        free(fraction_counts);
                    } else {
                        int qstrat = 0;
                        if (num_homop > 1) {
                            // want something in [0, num_homop-1]
                            qstrat = min(bam_get_qual(p->b)[p->qpos + query_pos_offset], num_homop);
                            qstrat = max(0, qstrat - 1);
                        }
                        pileup->matrix[partial_index + featlen * qstrat] += 1;
                    }
                }
            }
        }
        major_col += (dtype_featlen * (max_ins+1));
        n_cols += max_ins;
    }
    kh_destroy(BADREADS, no_rle_tags);
    pileup->n_cols = n_cols;

    bam_itr_destroy(data->iter);
    bam_mplp_destroy(mplp);
    free(data);
    free(plp);
    free(chr);

    return pileup;
}
