#define _GNU_SOURCE
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "medaka_bamiter.h"
#include "medaka_counts.h"
#include "medaka_read_matrix.h"

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
    const char tag_name[2] = "";
    const int tag_value = 0;
    const bool keep_missing = false;
    const size_t num_homop = 5;
    const bool weibull_summation = false;
    const char* read_group = NULL;
    const int min_mapQ = 1;
    const bool row_per_read = false;
    const bool include_dwells = false;
    const bool include_haplotype = false;
    const int max_reads = 100;

    bam_fset* bam_set = create_bam_fset(bam_file);

    plp_data pileup = calculate_pileup(
        reg, bam_set, num_dtypes, dtypes,
        num_homop, tag_name, tag_value, keep_missing,
        weibull_summation, read_group, min_mapQ);
    print_pileup_data(pileup, dtypes);
    fprintf(stdout, "pileup is length %zu, with buffer of %zu columns\n", pileup->n_cols, pileup->buffer_cols);
    destroy_plp_data(pileup);

    read_aln_data read_aln = calculate_read_alignment(
        reg, bam_set, num_dtypes, dtypes,
        tag_name, tag_value, keep_missing,
        read_group, min_mapQ, row_per_read,
        include_dwells, include_haplotype, max_reads);
    print_read_aln_data(read_aln);
    fprintf(stdout, "pileup is length %zu with %zu reads, with buffer of %zu positions and %zu reads\n", read_aln->n_pos, read_aln->n_reads, read_aln->buffer_pos, read_aln->buffer_reads);
    destroy_read_aln_data(read_aln);

    destroy_bam_fset(bam_set);
    return 0;
}
