#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "medaka_common.h"
#include "fastrle.h"


KSEQ_INIT(gzFile, gzread)

void rle(char *in, int inlen, size_t block_size, char *out, char *outruns) {
    if (block_size > 94) {
        fprintf(stderr, "Invalid block length %i (> 94)", block_size); exit(1);
    }
    char c = in[0];
    size_t l = 1;
    size_t oi = 0;
    for (size_t i=1; i<inlen; ++i) {
        if ((in[i] != c) || (l == block_size)) {
            out[oi] = c;
            // note, we are skipping quality 0 (!) for sanity
            outruns[oi] = l + 33;
            c = in[i];
            l = 1;
            oi++;
        } else {
            l++;
        }
    }
    out[oi] = c;
    outruns[oi] = l + 33;
    oi++;
    out[oi] = '\0';
    outruns[oi] = '\0';
}


size_t fastrle(char *fname, size_t block_size) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fname, "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        char* out = xalloc(seq->seq.l + 1, sizeof(char), "out");
        char* outruns = xalloc(seq->seq.l + 1, sizeof(char), "outlen");
        rle(seq->seq.s, seq->seq.l, block_size, out, outruns);
        printf("@%s\n%s\n+\n%s\n", seq->name.s, out, outruns);
        free(out);
        free(outruns);
    }
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}
