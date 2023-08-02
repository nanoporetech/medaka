#include <assert.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"

#include "kvec.h"
#include "kstring.h"

#include "medaka_trimbam.h"
#include "medaka_bamiter.h"
#include "medaka_common.h"

#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define bam1_seqi(s, i) (bam_seqi((s), (i)))
#define bam_nt16_rev_table seq_nt16_str


void upper_string(char s[]) {
   int c = 0;
   while (s[c] != '\0') {
      if (s[c] >= 'a' && s[c] <= 'z') {
         s[c] = s[c] - 32;
      }
      c++;
   }
}


/** Contructs a set of sequences.
 *
 * @returns a trimmed_reads pointer.
 *
*/
trimmed_reads create_trimmed_reads(void) {
    trimmed_reads reads = xalloc(1, sizeof(*reads), "trimmed_reads");
    kv_init(reads->sequences);
    kv_init(reads->is_reverse);
    return reads;
}


/** Destroys set of sequences.
 *
 * @returns a trimmed_reads pointer.
 *
*/
void destroy_trimmed_reads(trimmed_reads reads) {
    for (size_t i=0; i<reads->sequences.n; ++i){
        free(reads->sequences.a[i]);
    }
    kv_destroy(reads->sequences);
    kv_destroy(reads->is_reverse);
    free(reads);
}


/** Adds a read to a set of sequences.
 *
 * @param reads the set of reads.
 * @param read new sequence.
 * @param is_rev is sequence reversed.
 *
*/
void add_read(trimmed_reads reads, char * read, bool is_rev){
    kv_push(char *, reads->sequences, read);
    kv_push(bool, reads->is_reverse, is_rev);
}


/** Find start and end positions in read which cause it
 *  to be trimmed to be trimmed to specified reference
 *  coordinates.
 *
 *  @param record htslib bam1_t alignment.
 *  @param rstart start reference position.
 *  @param rend end reference position.
 *  @param partial allow partial overlaps.
 *  @param qstart (out) query start position.
 *  @param qend (out) query end position.
 *  @returns string with corresponding reference sequence.
 *
 *  If the read does not completely span the reference region one
 *  or both of qstart and qend will be set to -1. Additionally on
 *  any other error the return value will be NULL.
 */
char * trim_read(bam1_t *record, int rstart, int rend, bool partial, int *qstart, int *qend){
    uint32_t *cigar = bam_get_cigar(record);
    const bam1_core_t *c = &record->core;
    *qstart = -1;
    *qend = -1;
    char *ref_chunk = NULL;
 
    // check we span start
    if (record->core.pos > rstart) {
        if (partial) {
            *qstart = 0;
        } else {
            return ref_chunk;
        }
    }

    // read pos is an index into the original sequence that is present in the FASTQ
    // on the strand matching the reference
    int read_pos = 0;
    int ref_pos = c->pos;
    char *qname = bam_get_qname(record);
    kstring_t ref = {0, 0, NULL};

    // Disable calculation of reference for now, since we don't need it.
    // Though to make the rest of the code happy we need to push at least
    // one character. see also <-> mark below
    kputc('N', &ref);
    
    //uint8_t *tag = bam_aux_get((const bam1_t*) record, "cs");
    //if (tag == NULL) { // tag isn't present
    //    fprintf(stderr, "Could not read 'cs' tag from read %s\n", qname);
    //    return ref_chunk;
    //}
    //char* cs = bam_aux2Z(tag); // should check this

    //// 1:=(match), 2:-(del), 3:+(ins), 4:*(subs), 5:skip next(sub)
    //size_t flag = 0;
    //for (size_t i = 0; i < strlen(cs); ++i){
    //    // not handling intron '~'
    //    char c = cs[i];
    //    if (c == '=') {
    //        flag = 1;
    //    } else if (c == '-') {
    //        flag = 2;
    //    } else if (c == '+') {
    //        flag = 3;
    //    } else if (c == '*') {
    //        flag = 4;
    //    } else if (flag == 5){
    //        continue;
    //    } else {
    //        if (flag == 1 || flag == 2) {
    //            kputc(c, &ref);
    //        } else if (flag == 4) {
    //            kputc(c, &ref);
    //            flag = 5;
    //        }
    //    }
    //}

    bool found_start = false;
    bool found_end = false;
    int cigar_len = 0;
    int cigar_op = 0;

    for (int ci = 0; ci < c->n_cigar; ++ci) {
        cigar_len = cigar[ci] >> 4;
        cigar_op = cigar[ci] & 0xf;

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;
 
        // Process match between the read and the reference
        bool is_aligned = false;
        if (cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            is_aligned = true;
            read_inc = 1;
            ref_inc = 1;
        } else if (cigar_op == BAM_CDEL) {
            ref_inc = 1;   
        } else if (cigar_op == BAM_CREF_SKIP) {
            // erm don't handle this one
            ref_inc = 1;
            fprintf(stderr, "Unhandled cigar op, %d (REF_SKIP), in read %s\n", cigar_op, qname);
            return ref_chunk;
        } else if (cigar_op == BAM_CINS) {
            read_inc = 1;
        } else if (cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1;
        } else if (cigar_op == BAM_CHARD_CLIP) {
            read_inc = 0;
        } else {
            fprintf(stderr, "Unhandled cigar op, %d, in read %s\n", cigar_op, qname);
            return ref_chunk;
        }

        // Iterate over the pairs of aligned bases
        for (int j = 0; j < cigar_len; ++j) {
            if (is_aligned) {
                if (!found_start) {
                    if (ref_pos == rstart){
                        *qstart = read_pos;
                        found_start = true;
                    } else if (ref_pos > rstart) {
                        // we've just moved past, take the previous.
                        // initial check ensures we covered start pos
                        *qstart = read_pos - 1;
                        found_start = true;
                    }
                }
                if (!found_end) {
                    if (ref_pos == rend) {
                        *qend = read_pos;
                        found_end = true;
                    } else if (ref_pos > rend) {
                        // we've just moved past, take the previous
                        *qend = read_pos - 1;
                        found_end = true;
                    }
                }
            }
            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
    if (*qend == -1 && partial) {
        *qend = read_pos;
        // correct for soft-clipping
        if (cigar_op == BAM_CSOFT_CLIP) {
            *qend -= cigar_len;
        }
    }


    // correct
    size_t st = record->core.pos > rstart ? 0 : rstart - record->core.pos;
    size_t length = min(strlen(ref.s) - st, rend - rstart);
    st = 0; length = 1; //<-> as part of ref disablement
    ref_chunk = substring(ref.s, st, length);
    upper_string(ref_chunk);
    free(ref.s);
    return ref_chunk;
}



/** Get reads overlapping a region, trimmed to the region.
 *
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
 *  @param num_dtypes number of datatypes in bam.
 *  @param dtypes prefixes on query names indicating datatype.
 *  @param tag name used for read filtering.
 *  @param tag value used for read filtering.
 *  @param keep_missing, if true keep reads without the tag.
 *  @param partial, whether to keep partially spanning reads. 
 *  @param read_group used for read filtering.
 *
 *  @returns void
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
trimmed_reads retrieve_trimmed_reads(
    const char *region, const char *bam_file, size_t num_dtypes, char *dtypes[],
    const char tag_name[2], const int tag_value, const bool keep_missing, 
    const bool partial, const char *read_group, const int min_mapq){

    if (num_dtypes == 1 && dtypes != NULL) {
        fprintf(stderr, "Recieved invalid num_dtypes and dtypes args.\n");
        exit(1);
    }

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
        exit(1);
    }
   
    // open bam etc. 
    htsFile *fp = hts_open(bam_file, "rb");
    hts_idx_t *idx = sam_index_load(fp, bam_file);
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == 0 || idx == 0 || fp == 0) {
        hts_close(fp); hts_idx_destroy(idx); sam_hdr_destroy(hdr);
        free(chr);
        fprintf(stderr, "Failed to read .bam file '%s'.", bam_file);
        exit(1);
    }

    // setup bam interator
    mplp_data *data = xalloc(1, sizeof(mplp_data), "pileup init data");
    data->fp = fp; data->hdr = hdr; data->iter = bam_itr_querys(idx, hdr, region);
    data->min_mapQ = min_mapq; memcpy(data->tag_name, tag_name, 2); data->tag_value = tag_value;
    data->keep_missing = keep_missing; data->read_group = read_group;

    // get trimmed reads
    bam1_t *record = bam_init1();
    //kvec_t(char*) sequences; kv_init(sequences);
    //kvec_t(bool) is_reverse; kv_init(is_reverse);
    trimmed_reads reads = create_trimmed_reads();
    //TODO: support multiple data types
    char * ref = xalloc(1, sizeof(char), "chr"); 
    while(read_bam((void *) data, record) > 0){
        int qstart, qend;
        char * ref_chunk = trim_read(record, start, end, partial, &qstart, &qend);
        //if (qstart < 0) qstart = 0;
        //if (qend < 0 ) qend = record->core.l_qseq;
 
        if (qstart >=0 && qend >=0 && ref_chunk != NULL) {
            if (strlen(ref_chunk) > strlen(ref)) {
                free(ref); ref = ref_chunk;
            }
            uint8_t *q = bam_get_seq(record);
            //uint32_t len = record->core.l_qseq;
            if (qend - qstart > 1) {
                char *qseq = xalloc(qend - qstart + 1, sizeof(char), "seq");
                int j = 0;
                for(int i=qstart; i<qend ; i++){
                    qseq[j++] = seq_nt16_str[bam_seqi(q, i)];
                }
                //kv_push(char *, sequences, qseq);
                //kv_push(bool, is_reverse, bam_is_rev(record));
                add_read(reads, qseq, bam_is_rev(record));
            }
        } else {
            if (ref_chunk != NULL) {
               free(ref_chunk);
            }
            //fprintf(stderr, "dont like read %s: %i %i, %s\n", bam_get_qname(record), qstart, qend, ref_chunk);
        }
    }
    bam_destroy1(record);
    add_read(reads, ref, 0);

    bam_itr_destroy(data->iter);
    free(data);
    free(chr);

    hts_close(fp);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    return reads;
}


// Demonstrates usage
int _main(int argc, char *argv[]) {
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
    bool partial = true;
    const char* read_group = NULL;
    const int min_mapq = 1;

    trimmed_reads reads = retrieve_trimmed_reads(
        reg, bam_file, num_dtypes, dtypes,
        tag_name, tag_value, keep_missing, partial, read_group, min_mapq);
    for (size_t i=0; i<reads->sequences.n; ++i){
        fprintf(stderr, "%i  %s\n", reads->is_reverse.a[i], reads->sequences.a[i]);
    }
    destroy_trimmed_reads(reads);
    exit(0); 
}

