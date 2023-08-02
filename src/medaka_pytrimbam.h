
// a set of strings (with which cffi is happy, avoid exposing kvec_t)
typedef struct {
    size_t n_seqs;
    char **seqs;
    _Bool *is_rev;
    void *t_reads; 
} _seq_set;

typedef _seq_set *seq_set;


/** Python wrapper to retrieve_trimmed_reads 
 *
 * @param seqs as seq_set structure.
 *
 * @returns a seq_set with pointers to underlying kvecs
 *
 */
seq_set PY_retrieve_trimmed_reads(
    const char *region, const char *bam_file, size_t num_dtypes, char *dtypes[],
    const char tag_name[2], const int tag_value, const bool keep_missing, 
    const bool partial, const char *read_group, const int min_mapq);


/** Cleans up the result of PY_retrieve_trimmed_reads
 *
 * @param seqs as seq_set structure.
 *
 * @returns void
 *
 */
void PY_destroy_reads(seq_set reads);

