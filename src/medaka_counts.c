#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str


// lazyness
void *xalloc(size_t num, size_t size, char* msg){
    void *res = calloc(num, size);
    if (res == NULL){
        fprintf(stderr, "Failed to allocate mem for %s\n", msg);
        exit(1);
    }
    return res;
}


// medaka-style base encoding
const char plp_bases[] = "XacgtACGTdD";
const size_t featlen = 11;
const size_t fwd_del = 10;
const size_t rev_del = 9;

// convert 16bit IUPAC (+16 for strand) to plp_base index
size_t num2countbase[32] = {
 0, 5, 6, 0, 7, 0, 0, 0,
 8, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 2, 0, 3, 0, 0, 0,
 4, 0, 0, 0, 0, 0, 0, 0,
}; 


// medaka-style feature data
typedef struct {
    size_t n_cols;
    size_t *counts;
    size_t *major;
    size_t *minor;
} _plp_data;

typedef _plp_data *plp_data; 


// constructor for feature data
plp_data create_plp_data(size_t n_cols) {
    plp_data data = xalloc(1, sizeof(*data), "plp_data");
    data->n_cols = n_cols;
    data->counts = xalloc(featlen * n_cols, sizeof(size_t), "count");
    data->major = xalloc(n_cols, sizeof(size_t), "major");
    data->minor = xalloc(n_cols, sizeof(size_t), "minor");
    return data;
}


// destructor for feature data
void destroy_plp_data(plp_data data) {
    free(data->counts);
    free(data->major);
    free(data->minor);
    free(data);
}


// some necessities
typedef struct {
    htsFile *fp;
    bam_hdr_t *hdr;
    hts_itr_t *iter;
    int min_mapQ;
} mplp_data;


static int read_bam(void *data, bam1_t *b) {
    mplp_data *aux = (mplp_data*) data;
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        break;
    }
    return ret;
}


// Generate medaka-style feature data in a region of a bam
plp_data calculate_pileup(const char *region, const char *bam_file) { 

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
    //fprintf(stderr, "%s  %d  %d\n", region, start, end); 
   
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
    data->fp = fp; data->hdr = hdr; data->min_mapQ = 1;
    data->iter = bam_itr_querys(idx, hdr, region);
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
    plp_data pileup = create_plp_data(n_cols);

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
            pileup->major[major_col / featlen + i] = pos;
            pileup->minor[major_col / featlen + i] = i;
        }

        // loop through all reads at this position
        for (int i = 0; i < n_plp; ++i) {
            const bam_pileup1_t *p = plp[0] + i;
            if (p->is_refskip) {
                continue;
            }
            int base_i;
            if (p->is_del) {
                base_i = bam_is_rev(p->b) ? rev_del : fwd_del;
                //base = plp_bases[base_i]; 
                pileup->counts[major_col + base_i] += 1;
            } else { // handle pos and any following ins
                int max_j = p->indel > 0 ? p->indel : 0;
                for (int j = 0; j <= max_j; ++j){
                    int base_j = bam1_seqi(bam1_seq(p->b), p->qpos + j);
                    //base = seq_nt16_str[base_j];
                    if bam_is_rev(p->b) {
                        base_j += 16;
                    }
                    base_i = num2countbase[base_j];
                    pileup->counts[major_col + featlen * j + base_i] += 1;
                }
            }
        }
        major_col += (featlen * (max_ins+1));
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


void print_pileup_data(plp_data pileup){
    fprintf(stdout, "pos\tins\t");
    for (size_t j = 0; j < featlen; ++j){
        fprintf(stdout, "%c\t", plp_bases[j]);
    }
    fprintf(stdout, "depth\n");
    for (size_t j = 0; j < pileup->n_cols; ++j) {
        int s = 0;
        fprintf(stdout, "%zu\t%zu\t", pileup->major[j], pileup->minor[j]);
        for (size_t i = 0; i < featlen; ++i){
            size_t c = pileup->counts[j * featlen + i];
            s += c;
            fprintf(stdout, "%zu\t", c);
        }
        fprintf(stdout, "%d\n", s);
    }
}


int main(int argc, char *argv[]) {
    if(argc != 3) {
        fprintf(stderr, "Usage %s <bam> <region>.\n", argv[0]);
        exit(1);
    }
    const char *bam_file = argv[1];
    const char *reg = argv[2];

    plp_data pileup = calculate_pileup(reg, bam_file);
    print_pileup_data(pileup);
    destroy_plp_data(pileup);
    exit(0); 
}

