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
#include "htslib/faidx.h"

#include "khash.h"
#include "kvec.h"
#include "medaka_bamiter.h"
#include "medaka_common.h"
#include "medaka_counts.h"
#include "medaka_khcounter.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str
#define bam_nt16_table seq_nt16_table

// For recording bad reads and skipping processing
KHASH_SET_INIT_STR(BADREADS)

/** Swap two strings
 *
 *  @param a first string
 *  @param b second string
 *  @returns a plp_data pointer.
 *
 */
void swap_strings(char** a, char** b) {
    char *temp = *a;
    *a = *b;
    *b = temp;
}

size_t base2_index(char c) {
    if (c == 'A') return 0;
    else if (c == 'C') return 1;
    else if (c == 'G') return 2;
    else if (c == 'T') return 3;
    else if (c == 'a') return 9;
    else if (c == 'c') return 10;
    else if (c == 'g') return 11;
    else if (c == 't') return 12;
    else return 0;
}

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
void format_uint8_array(uint8_t* values, size_t length, char* result) {
    size_t len = 0;
    for (size_t i = 0; i < length; ++i) {
        len += uint8_to_str(values[i], result + len);
        strcpy(result + len, ",");
        len += 1;
    }
    result[len-1] = '\0';
}


/** Destroys a string set
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_string_set(string_set strings) {
    for(size_t i = 0; i < strings.n; ++i) {
        free(strings.strings[i]);
    }
    free(strings.strings);
}


/** Retrieves contents of key-value tab delimited file.
 *
 *  @param fname input file path.
 *  @returns a string_set
 *
 *  The return value can be free'd with destroy_string_set.
 *  key-value pairs are stored sequentially in the string set
 *
 */
string_set read_key_value(char * fname) {
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    kvec_t(char*) strings;
    kv_init(strings);
    
    fp = fopen(fname, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while ((read = getdelim(&line, &len, '\t', fp)) != -1) {
        line[read - 1] = '\0';
        char *key = NULL; swap_strings(&key, &line);
        kv_push(char*, strings, key);
        read = getline(&line, &len, fp);
        line[read - 1] = '\0';
        char *value = NULL; swap_strings(&value, &line);
        kv_push(char*, strings, value);
    }
    free(line);
    // move strings into a simpler container (awkward to pass kvec_t through cffi)
    string_set my_strings;
    my_strings.n = strings.n;
    my_strings.strings = strings.a;
    return(my_strings);
}


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param buffer_cols number of pileup columns.
 *  @param feature_length length of feature vector.
 *  @param num_dtypes number of datatypes in pileup.
 *  @param num_homop maximum homopolymer length to consider.
 *  @param fixed_size if not zero data matrix is allocated as fixed_size * n_cols, ignoring other arguments
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t buffer_cols, size_t feature_length, size_t num_dtypes, size_t num_homop, size_t fixed_size) {
    assert(buffer_cols >= n_cols);
    plp_data data = xalloc(1, sizeof(_plp_data), "plp_data");
    data->buffer_cols = buffer_cols;
    data->num_dtypes = num_dtypes;
    data->num_homop = num_homop;
    data->n_cols = n_cols;
    if (fixed_size != 0) {
        assert(buffer_cols == n_cols);
        data->matrix = xalloc(fixed_size * n_cols, sizeof(int), "matrix");
    } else {
        data->matrix = xalloc(feature_length * num_dtypes * buffer_cols * num_homop, sizeof(size_t), "matrix");
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
void enlarge_plp_data(plp_data pileup, size_t buffer_cols, size_t feature_length) {
    assert(buffer_cols > pileup->buffer_cols);
    size_t old_size = feature_length * pileup->num_dtypes * pileup->num_homop * pileup->buffer_cols;
    size_t new_size = feature_length * pileup->num_dtypes * pileup->num_homop * buffer_cols;

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
 *  @param num_dtypes number of datatypes in the pileup.
 *  @param dtypes datatype prefix strings.
 *  @param num_homop maximum homopolymer length to consider.
 *  @returns void
 *
 */
void print_pileup_data(plp_data pileup, size_t num_dtypes, char *dtypes[], size_t num_homop){
    fprintf(stdout, "pos\tins\t");
    if (num_dtypes > 1) {  //TODO header for multiple dtypes and num_homop > 1
        for (size_t i = 0; i < num_dtypes; ++i) {
            for (size_t j = 0; j < featlen; ++j){
                fprintf(stdout, "%s.%c\t", dtypes[i], plp_bases[j]);
            }
        }
    } else {
        for (size_t k = 0; k < num_homop; ++k) {
            for (size_t j = 0; j < featlen; ++j){
                fprintf(stdout, "%c.%lu\t", plp_bases[j], k+1);
            }
        }
    }
    fprintf(stdout, "depth\n");
    for (size_t j = 0; j < pileup->n_cols; ++j) {
        int s = 0;
        fprintf(stdout, "%zu\t%zu\t", pileup->major[j], pileup->minor[j]);
        for (size_t i = 0; i < num_dtypes * featlen * num_homop; ++i){
            size_t c = pileup->matrix[j * num_dtypes * featlen * num_homop + i];
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
        bool weibull_summation, const char *read_group) {
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
    data->min_mapQ = 1; memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    data->keep_missing = keep_missing; data->read_group = read_group;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);
    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    // allocate output assuming one insertion per ref position
    int n_cols = 0;
    size_t buffer_cols = 2 * (end - start);
    plp_data pileup = create_plp_data(n_cols, buffer_cols, featlen, num_dtypes, num_homop, 0);

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
            float cols_per_pos = (float) (n_cols + max_ins) / (pos - start);
            // max_ins can dominate so add at least that
            buffer_cols = max_ins + max(2 * pileup->buffer_cols, (int) cols_per_pos * (end - start));
            enlarge_plp_data(pileup, buffer_cols, featlen);
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



/** Generates clair3-style pileup feature data in a region of a bam.
 *
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
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

/**
 * The pileup input is 594 integers – 33 genome positions wide with 18 features at each position –
 * 
 * A+, C+, G+, T+, I_S+, I^1 S+, D_S+, D^1_S+, D_R+, A-, C-, G-, T-, I_S-, I^1_S-, D_S-, D^1_S-, and D_R-
 *
 * A, C, G, T, I, D, +, - means the count of read support of the four nucleotides: insertion,
 * deletion, positive strand, and negative strand. Superscript “1” means only the indel with the
 * highest read support is counted (i.e., all indels are counted if without “1“). Subscript “S”/“R” means
 * the starting/non-starting position of an indel. For example, a 3bp deletion with the most reads support
 * will have the first deleted base counted in either D1_S+ or D1_S-, and the second and third deleted bases
 * counted in either D_R+ or D_R-. The design was determined experimentally, but the rationale is that for
 * 1bp indels that are easy to call, look into the differences between the “S” counts, but reduce the
 * quality if the “R” counts and discrepancy between positions increase.
 *
 */

plp_data calculate_clair3_pileup(
        const char *region, const bam_fset* bam_set, const char * fasta_path) {
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
    // TODO: the code in bam iterator might not be correct for clair3
    //       e.g. mapQ filters, primary/secondary/supplementary, ...
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr; data->iter = bam_itr_querys(idx, hdr, region);
    data->min_mapQ = 1;
    //memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    //data->keep_missing = keep_missing; data->read_group = read_group;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);
    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    // allocate output, clair3 doesn't have insertion columns, so buffer remains fixed
    int n_cols = 0;
    size_t buffer_cols = end - start;
    plp_data pileup = create_plp_data(n_cols, buffer_cols, featlenclair3, 1, 1, 0);

    // get counts
    size_t major_col = 0;  // index into `pileup` corresponding to pos
    n_cols = 0;            // number of processed columns (including insertions, which clair3 doesn't have ;))

    // for recording deletion lengths
    size_t del_buf_size = 32;
    size_t* dels_f = xalloc(del_buf_size, sizeof(size_t), "dels_f");
    size_t* dels_r = xalloc(del_buf_size, sizeof(size_t), "dels_r");

    faidx_t* fai = fai_load(fasta_path);
    int len = 0;
    char *ref_seq = NULL;
    ref_seq = faidx_fetch_seq(fai, chr, start, end, &len);

    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {
        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if (pos < start) continue;
        if (pos >= end) break;
        n_cols++;

        memset(dels_f, 0, del_buf_size * sizeof(size_t));
        memset(dels_r, 0, del_buf_size * sizeof(size_t));

        // we still need this as positions might not be contiguous
        pileup->major[major_col / featlenclair3] = pos;
        pileup->minor[major_col / featlenclair3] = 0;

        // counters for insertion strings
        khash_t(KH_COUNTER) *ins_counts_f = kh_init(KH_COUNTER);
        khash_t(KH_COUNTER) *ins_counts_r = kh_init(KH_COUNTER);

        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) continue;

            if (p->indel < 0) {
                // there's a deletion starting on next genomic position,
                // record the length here and finalise after the read loop
                //  - actually deleted bases get recorded in next block
                size_t d = (size_t) -1 * p->indel;
                if (d >= del_buf_size) {
                    size_t new_size = max(d, 2 * del_buf_size);
                    dels_f = xrealloc(dels_f, new_size*sizeof(size_t), "dels_f");
                    memset(dels_f, 0, del_buf_size * sizeof(size_t));
                    dels_r = xrealloc(dels_r, new_size*sizeof(size_t), "dels_r");
                    memset(dels_r, 0, del_buf_size * sizeof(size_t));
                    del_buf_size = new_size;
                } 
                if (bam_is_rev(p->b)) {
                    dels_r[d - 1] += 1;
                } else {
                    dels_f[d - 1] += 1;
                }
            }

            // handle ref_base/sub/del
            int base_i;
            if (p->is_del) {
                // we've been deleted, +1 to DR
                base_i = bam_is_rev(p->b) ? c3_rev_del : c3_fwd_del;
            } else {
                // just a base
                int base_j = bam1_seqi(bam1_seq(p->b), p->qpos);
                if bam_is_rev(p->b) { base_j += 16; }
                base_i = num2countbaseclair3[base_j];
            }
            pileup->matrix[major_col + base_i] += 1;

            // handle insertion
            //  - build insert string then hash
            if (p->indel > 0) {
                // clair3 puts the reference base in the allele (a la VCF), we don't know
                // that so just take the inserted bases. Take care to handle case of
                // reference base being deleted when building insertion string
                size_t first = p->is_del ? 0 : 1;
                char* indel = (char*) xalloc(p->indel + 1, sizeof(char), "indel");
                for (size_t i = 0, j = first; j < p->indel + first; ++i, ++j) {
                    indel[i] = seq_nt16_str[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
                }
                indel[p->indel] = '\0'; // actually should be 0 already
                //fprintf(stderr, "pos: %i\tp->indel: %i\tp->is_del: %i\tindel: '%s'\n", pos, p->indel, p->is_del, indel);
                if (bam_is_rev(p->b)) {
                    kh_counter_increment(ins_counts_r, indel);
                } else {
                    kh_counter_increment(ins_counts_f, indel);
                }
                free(indel);
            }
        }

        // finalise deletions: DS (all) and D1S (best)
        //
        // forward
        size_t best_count = 0;
        size_t all_count = 0;
        for (size_t i = 0; i < del_buf_size; ++i) {
            size_t d = dels_f[i];
            all_count += d;
            best_count = max(best_count, d);
        }
        pileup->matrix[major_col + c3_fwd_del_all] = all_count;
        pileup->matrix[major_col + c3_fwd_del_best] = best_count;
        // reverse
        best_count = 0;
        all_count = 0;
        for (size_t i = 0; i < del_buf_size; ++i) {
            size_t d = dels_r[i];
            all_count += d;
            best_count = max(best_count, d);
        }
        pileup->matrix[major_col + c3_rev_del_all] = all_count;
        pileup->matrix[major_col + c3_rev_del_best] = best_count;

        
        // finalise IS and I1S
        // TODO: eventually we need to return all these
        //       - is that sensible, could be a lot
        //
        // forward
        kh_counter_stats_t stats = kh_counter_stats(ins_counts_f);
        pileup->matrix[major_col + c3_fwd_ins_all] = stats.sum;
        pileup->matrix[major_col + c3_fwd_ins_best] = stats.max;
        kh_counter_destroy(ins_counts_f);
        // reverse
        stats = kh_counter_stats(ins_counts_r);
        pileup->matrix[major_col + c3_rev_ins_all] = stats.sum;
        pileup->matrix[major_col + c3_rev_ins_best] = stats.max;
        kh_counter_destroy(ins_counts_r);

        //https://github.com/HKU-BAL/Clair3/blob/main/preprocess/CreateTensorPileup.py#L163
        int offset = pos - start;
        char ref_base = ref_seq[offset];
        int ref_offset_forward = base2_index(ref_base);
        int ref_offset_reverse = ref_offset_forward + reverse_pos_start;
        size_t forward_sum = 0;
        size_t reverse_sum = 0;
        for (size_t i = 0; i < 4; i++) {
            forward_sum += pileup->matrix[major_col + i];
            reverse_sum += pileup->matrix[major_col + i + reverse_pos_start];
        }
        pileup->matrix[major_col + ref_offset_forward] = -1 * forward_sum;
        pileup->matrix[major_col + ref_offset_reverse] = -1 * reverse_sum;
        // move to next position
        major_col += featlenclair3;
    }
    pileup->n_cols = n_cols;

    bam_itr_destroy(data->iter);
    bam_mplp_destroy(mplp);
    fai_destroy(fai);
    free(data);
    free(plp);
    free(chr);

    return pileup;
}


// Demonstrates usage
int main(int argc, char *argv[]) {
    if(argc < 3) {
        fprintf(stderr, "Usage %s <bam> <region>.\n", argv[0]);
        exit(1);
    }
    const char *bam_file = argv[1];
    const char *reg = argv[2];

    size_t num_dtypes = 1;
    char **dtypes = NULL;
    if (argc > 3) {
        num_dtypes = argc - 3;
        dtypes = &argv[3];
    }
    char tag_name[2] = "";
    int tag_value = 0;
    bool keep_missing = false;
    size_t num_homop = 5;
    bool weibull_summation = false;
    const char* read_group = NULL;

    bam_fset* bam_set = create_bam_fset(bam_file);

    //plp_data pileup = calculate_pileup(
    //    reg, bam_set, num_dtypes, dtypes,
    //    num_homop, tag_name, tag_value, keep_missing,
    //    weibull_summation, read_group);
    //print_pileup_data(pileup, num_dtypes, dtypes, num_homop);
    plp_data pileup = calculate_clair3_pileup(
        reg, bam_set, read_group);
    fprintf(stdout, "pileup is length %zu, with buffer of %zu columns\n", pileup->n_cols, pileup->buffer_cols);
    destroy_plp_data(pileup);
    destroy_bam_fset(bam_set);
    exit(0);
}
