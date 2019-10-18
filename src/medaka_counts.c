#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"

#include "medaka_bamiter.h"
#include "medaka_common.h"
#include "medaka_counts.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param num_dtypes number of datatypes in pileup.
 *  @param num_homop maximum homopolymer length to consider.
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t buffer_cols, size_t num_dtypes, size_t num_homop) {
    assert(buffer_cols >= n_cols);
    plp_data data = xalloc(1, sizeof(*data), "plp_data");
    data->buffer_cols = buffer_cols;
    data->num_dtypes = num_dtypes;
    data->num_homop = num_homop;
    data->n_cols = n_cols;
    data->matrix = xalloc(featlen * num_dtypes * buffer_cols * num_homop, sizeof(size_t), "matrix");
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
    size_t old_size = featlen * pileup->num_dtypes * pileup->num_homop * pileup->buffer_cols;
    size_t new_size = featlen * pileup->num_dtypes * pileup->num_homop * buffer_cols;
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


float* _get_weibull_scores(const bam_pileup1_t *p, const size_t indel, const size_t num_homop) {
    static const char* wtags[] = {"WL", "WK"};  // scale, shape
    float* fraction_counts = xalloc(num_homop, sizeof(float), "weibull_counts");
    double wtag_vals[2] = {0.0, 0.0};
    for (size_t i=0; i < 2; ++i) {
        uint8_t *tag = bam_aux_get(p->b, wtags[i]);
        if (tag == NULL) {
            fprintf(stderr, "Failed to retrieve tag for Weibull parameter %s.\n",
                    wtags[i]);
            exit(1);
        }
        uint32_t taglen = bam_auxB_len(tag);
        if (p->qpos + indel >= taglen) {
            fprintf(stderr, "%s tag was out of range for %s position %lu. taglen: %i\n",
                    wtags[i], bam_get_qname(p->b), p->qpos + indel, taglen);
            exit(1);
        }
        wtag_vals[i] = bam_auxB2f(tag, p->qpos + indel);

        /* TODO: understand why errno comes back set even when value found ok,
         *       could then skip our own test above
        if (errno == ERANGE) {
            uint32_t taglen = bam_auxB_len(tag);
            if (errno == EINVAL) {
                fprintf(stderr, "%s tag was not a binary tag\n", wtags[i]);
            } else {
                fprintf(stderr, "%s tag was out of range for %s position %lu. taglen: %i\n",
                        wtags[i], bam_get_qname(p->b), p->qpos + indel, taglen);
                exit(1);
            }
        }
        */
    }

    float wl = wtag_vals[0];
    float wk = wtag_vals[1];
    for (size_t x = 1; x < num_homop + 1; ++x) {
        fraction_counts[x-1] = exp(-pow(x / wl, wk)) - exp(-pow((x + 1) / wl, wk));
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
        const char *region, const char *bam_file, size_t num_dtypes, char *dtypes[],
        size_t num_homop, const char tag_name[2], const int tag_value, const bool keep_missing,
        bool weibull_summation) {
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
    htsFile *fp = hts_open(bam_file, "rb");
    hts_idx_t *idx = sam_index_load(fp, bam_file);
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || idx == 0 || fp == 0) {
        hts_close(fp); hts_idx_destroy(idx); bam_hdr_destroy(hdr);
        free(chr);
        fprintf(stderr, "Failed to read .bam file '%s'.", bam_file);
        exit(1);
    }

    // setup bam interator
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr; data->iter = bam_itr_querys(idx, hdr, region);
    data->min_mapQ = 1; memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    data->keep_missing = keep_missing;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);
    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    // allocate output assuming one insertion per ref position
    int n_cols = 0;
    size_t buffer_cols = 2 * (end - start);
    plp_data pileup = create_plp_data(n_cols, buffer_cols, num_dtypes, num_homop);

    // get counts
    size_t major_col = 0;  // index into `pileup` corresponding to pos
    n_cols = 0;            // number of processed columns (including insertions)
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
            float cols_per_pos = (float) n_cols / (pos - start);
            buffer_cols = max(2 * buffer_cols, cols_per_pos * (end - start));
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

            int base_i;
            if (p->is_del) {
                // deletions are kept in the first layer of qscore stratification, if any
                int qstrat = 0;
                base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
                //base = plp_bases[base_i];
                pileup->matrix[major_col + featlen * dtype * num_homop + featlen * qstrat + base_i] += 1;
            } else { // handle pos and any following ins
                int max_j = p->indel > 0 ? p->indel : 0;
                for (int j = 0; j <= max_j; ++j){
                    int base_j = bam1_seqi(bam1_seq(p->b), p->qpos + j);
                    //base = seq_nt16_str[base_j];
                    if bam_is_rev(p->b) {
                        base_j += 16;
                    }

                    base_i = num2countbase[base_j];
                    if (base_i != -1) {  // not an ambiguity code
                        size_t partial_index = 
                            major_col + dtype_featlen * j  // skip to column
                            + featlen * dtype * num_homop  // skip to datatype
                            //+ featlen * qstrat           // skip to qstrat/homop
                            + base_i;                      // the base

                        if (weibull_summation) {
                            float* fraction_counts = _get_weibull_scores(p, j, num_homop);
                            for (size_t qstrat = 0; qstrat < num_homop; ++qstrat) {
                                static const int scale = 10000;
                                pileup->matrix[partial_index + featlen * qstrat] += 
                                    scale * fraction_counts[qstrat];
                            }
                            free(fraction_counts);
                        } else {
                            int qstrat = min(bam_get_qual(p->b)[p->qpos + j], num_homop) - 1;
                            pileup->matrix[partial_index + featlen * qstrat] += 1;
                        }
                    }
                }
            }
        }
        major_col += (dtype_featlen * (max_ins+1));
        n_cols += max_ins;
    }
    pileup->n_cols = n_cols;

    bam_itr_destroy(data->iter);
    bam_mplp_destroy(mplp);
    free(data);
    free(plp);
    free(chr);

    hts_close(fp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);

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
    size_t num_homop = 2;
    bool weibull_summation = true;

    plp_data pileup = calculate_pileup(
        reg, bam_file, num_dtypes, dtypes,
        num_homop, tag_name, tag_value, keep_missing,
        weibull_summation);
    print_pileup_data(pileup, num_dtypes, dtypes, num_homop);
    fprintf(stdout, "pileup is length %zu, with buffer of %zu columns\n", pileup->n_cols, pileup->buffer_cols);
    destroy_plp_data(pileup);
    exit(0); 
}

