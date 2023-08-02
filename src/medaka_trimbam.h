#ifndef _MEDAKA_TRIMBAM_H
#define _MEDAKA_TRIMBAM_H

// an dynamic set of sequences
typedef struct _trimmed_reads {
    kvec_t(char*) sequences;
    kvec_t(bool) is_reverse;
} _trimmed_reads;

typedef _trimmed_reads *trimmed_reads;


/** Contructs a set of sequences.
 *
 * @returns a trimmed_reads pointer.
 *
*/
trimmed_reads create_trimmed_reads(void);


/** Destroys set of sequences.
 *
 * @returns a trimmed_reads pointer.
 *
*/
void destroy_trimmed_reads(trimmed_reads reads);


/** Adds a read to a set of sequences.
 *
 * @param reads the set of reads.
 * @param read new sequence.
 * @param is_rev is sequence reversed.
 *
*/
void add_read(trimmed_reads reads, char * read, bool is_rev);


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
    const bool partial, const char *read_group, const int min_mapq);
#endif
