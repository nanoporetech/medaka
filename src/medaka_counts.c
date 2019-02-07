#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "htslib/sam.h"

#include "medaka_counts.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str


/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg){
    void *res = calloc(num, size);
    if (res == NULL){
        fprintf(stderr, "Failed to allocate mem for %s\n", msg);
        exit(1);
    }
    return res;
}


/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, int position, int length) {
   char *ptr;
   size_t i;

   ptr = malloc(length + 1);

   for (i = 0 ; i < length ; i++) {
      *(ptr + i) = *(string + position);
      string++;
   }

   *(ptr + i) = '\0';
   return ptr;
}


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param num_dtypes number of datatypes in pileup.
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t num_dtypes) {
    plp_data data = xalloc(1, sizeof(*data), "plp_data");
    data->n_cols = n_cols;
    data->counts = xalloc(featlen * num_dtypes * n_cols, sizeof(size_t), "count");
    data->major = xalloc(n_cols, sizeof(size_t), "major");
    data->minor = xalloc(n_cols, sizeof(size_t), "minor");
    return data;
}


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data) {
    free(data->counts);
    free(data->major);
    free(data->minor);
    free(data);
}


// parameters for bam iteration
typedef struct {
    htsFile *fp;
    bam_hdr_t *hdr;
    hts_itr_t *iter;
    int min_mapQ;
    const char tag_name[2];
    int tag_value;
    bool keep_missing;
} mplp_data;


// iterator for reading bam
static int read_bam(void *data, bam1_t *b) {
    mplp_data *aux = (mplp_data*) data;
    uint8_t *tag;
    bool check_tag = (strcmp(aux->tag_name, "") != 0);
    int ret;
    while (1) {
        ret = aux->iter ? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if (ret<0) break;
        // only take primary alignments
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        // filter by mapping quality
        if ((int)b->core.qual < aux->min_mapQ) continue;
        // filter by tag
        if (check_tag) {
            tag = bam_aux_get((const bam1_t*) b, aux->tag_name);
            if (tag == NULL){ // tag isn't present or is currupt
                if (aux->keep_missing) {
                    break;
                } else {
                    continue;
                }
            }
            int tag_value = bam_aux2i(tag);
            if (errno == EINVAL) continue; // tag was not integer
            if (tag_value != aux->tag_value) continue;
        }
        break;
    }
    return ret;
}


/** Generates medaka-style feature data in a region of a bam.
 *  
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
 *  @param num_dtypes number of datatypes in bam.
 *  @param dtypes prefixes on query names indicating datatype.
 *  @returns a pileup counts data pointer.
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
plp_data calculate_pileup(const char *region, const char *bam_file, size_t num_dtypes, char *dtypes[], const char tag_name[2], const int tag_value, const bool keep_missing) { 

    if (num_dtypes == 1 && dtypes != NULL) {
        fprintf(stderr, "Recieved invalid num_dtypes and dtypes args.\n");
        exit(1);
    }
    //fprintf(stderr, "%u\n", num_dtypes);
    //for (size_t i=0; i<num_dtypes; ++i) {
    //    fprintf(stderr, "%s\n", dtypes[i]);
    //}
    const size_t dtype_featlen = featlen * num_dtypes;

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
        fprintf(stderr, "Failed to read .bam file '%s'.", bam_file);
        exit(1);
    }

    // setup bam interator
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr; data->iter = bam_itr_querys(idx, hdr, region);
    data->min_mapQ = 1; memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;

    bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)& data);
    const bam_pileup1_t **plp = xalloc(1, sizeof(bam_pileup1_t *), "pileup");
    int ret, pos, tid, n_plp;

    // first iterate to find the total number of cols (including inserts)
    int n_cols = 0;
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {
        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if(pos < start) continue;
        if(pos >= end) break;
        n_cols++;

        int max_ins = 0;
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->indel > 0 && max_ins < p->indel) max_ins = p->indel;
        }
        n_cols += max_ins;
    }
    plp_data pileup = create_plp_data(n_cols, num_dtypes);

    // reset and iterate to get counts
    bam_itr_destroy(data->iter);
    data->iter = bam_itr_querys(idx, hdr, region);
    bam_mplp_destroy(mplp);
    mplp = bam_mplp_init(1, read_bam, (void **)& data);
    size_t major_col = 0; // col of `pileup` corresponding to pos
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, &n_plp, plp) > 0)) {
        const char *c_name = data->hdr->target_name[tid];
        if (strcmp(c_name, chr) != 0) continue;
        if(pos < start) continue;
        if(pos >= end) break;
           
        // find maximum insert, again...
        int max_ins = 0;
		for (int i = 0; i < n_plp; ++i) {
		    const bam_pileup1_t *p = plp[0] + i;
		    if (p->indel > 0 && max_ins < p->indel) max_ins = p->indel;
		}

        // set major/minor position indexes, minors hold ins
        for (int i = 0; i <= max_ins; ++i ) {
            pileup->major[major_col / dtype_featlen + i] = pos;
            pileup->minor[major_col / dtype_featlen + i] = i;
        }

        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) {
                continue;
            }

            // find to which datatype the read belongs
            int dtype = 0;
            if (num_dtypes > 1){
                char *qname = bam_get_qname(p->b);
                bool found = false;
                for (dtype = 0; dtype < num_dtypes; ++dtype){
                    char *aux = substring(qname, 0, strlen(dtypes[dtype]));
                    if(strcmp(dtypes[dtype], aux) == 0) {
                        found = true;
                        free(aux);
                        break;
                    }
                    free(aux);
                }
                if (!found) {
                    fprintf(stderr, "Datatype not found for %s.\n", qname);
                    exit(1);
                }
            }

            int base_i;
            if (p->is_del) {
                base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
                //base = plp_bases[base_i]; 
                pileup->counts[major_col + featlen * dtype + base_i] += 1;
            } else { // handle pos and any following ins
                int max_j = p->indel > 0 ? p->indel : 0;
                for (int j = 0; j <= max_j; ++j){
                    int base_j = bam1_seqi(bam1_seq(p->b), p->qpos + j);
                    //base = seq_nt16_str[base_j];
                    if bam_is_rev(p->b) {
                        base_j += 16;
                    }
                    base_i = num2countbase[base_j];
                    pileup->counts[major_col + dtype_featlen * j + featlen * dtype + base_i] += 1;
                }
            }
        }
        major_col += (dtype_featlen * (max_ins+1));
    }

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


/** Prints a pileup data structure.
 *
 *  @param pileup a pileup counts structure.
 *  @param num_dtypes number of datatypes in the pileup.
 *  @param dtypes datatype prefix strings.
 *  @returns void
 *
 */
void print_pileup_data(plp_data pileup, size_t num_dtypes, char *dtypes[]){
    fprintf(stdout, "pos\tins\t");
    if (num_dtypes > 1) {
        for (size_t i = 0; i < num_dtypes; ++i) {
            for (size_t j = 0; j < featlen; ++j){
                fprintf(stdout, "%s.%c\t", dtypes[i], plp_bases[j]);
            }
        }
    } else {
        for (size_t j = 0; j < featlen; ++j){
            fprintf(stdout, "%c\t", plp_bases[j]);
        }
    }
    fprintf(stdout, "depth\n");
    for (size_t j = 0; j < pileup->n_cols; ++j) {
        int s = 0;
        fprintf(stdout, "%zu\t%zu\t", pileup->major[j], pileup->minor[j]);
        for (size_t i = 0; i < num_dtypes * featlen; ++i){
            size_t c = pileup->counts[j * num_dtypes * featlen + i];
            s += c;
            fprintf(stdout, "%zu\t", c);
        }
        fprintf(stdout, "%d\n", s);
    }
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

    plp_data pileup = calculate_pileup(
        reg, bam_file, num_dtypes, dtypes,
        tag_name, tag_value, keep_missing);
    print_pileup_data(pileup, num_dtypes, dtypes);
    destroy_plp_data(pileup);
    exit(0); 
}

