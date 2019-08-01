// medaka-style feature data
typedef struct {
    size_t buffer_cols;
    size_t num_dtypes;
    size_t n_cols;
    size_t *counts;
    size_t *major;
    size_t *minor;
} _plp_data;

typedef _plp_data *plp_data; 


// medaka-style base encoding
static const char plp_bases[] = "XacgtACGTdD";
static const size_t featlen = 11; // len of the above
static const size_t fwd_del = 10; // position of D
static const size_t rev_del = 9;  // position of d


// convert 16bit IUPAC (+16 for strand) to plp_bases index
static size_t num2countbase[32] = {
 0, 5, 6, 0, 7, 0, 0, 0,
 8, 0, 0, 0, 0, 0, 0, 0,
 0, 1, 2, 0, 3, 0, 0, 0,
 4, 0, 0, 0, 0, 0, 0, 0,
}; 


/** Constructs a pileup data structure.
 *
 *  @param n_cols number of pileup columns.
 *  @param buffer_cols number of pileup columns for which to allocate memory
 *  @param num_dtypes number of datatypes in pileup.
 *  @see destroy_plp_data
 *  @returns a plp_data pointer.
 *
 *  The return value can be freed with destroy_plp_data.
 *
 */
plp_data create_plp_data(size_t n_cols, size_t buffer_cols, size_t num_dtypes);


/** Enlarge the internal buffers of a pileup data structure.
 *
 *  @param pileup a plp_data pointer.
 *  @param buffer_cols number of pileup columns for which to allocate memory
 *
 */
void enlarge_plp_data(plp_data pileup, size_t buffer_cols);


/** Destroys a pileup data structure.
 *
 *  @param data the object to cleanup.
 *  @returns void.
 *
 */
void destroy_plp_data(plp_data data);


/** Prints a pileup data structure.
 *
 *  @param pileup a pileup counts structure.
 *  @param num_dtypes number of datatypes in the pileup.
 *  @param dtypes datatype prefix strings.
 *  @returns void
 *
 */
void print_pileup_data(plp_data pileup, size_t num_dtypes, char *dtypes[]);


/** Generates medaka-style feature data in a region of a bam.
 *  
 *  @param region 1-based region string.
 *  @param bam_file input aligment file.
 *  @param num_dtypes number of datatypes in bam.
 *  @param dtypes prefixes on query names indicating datatype.
 *  @param tag_name by which to filter alignments
 *  @param tag_value by which to filter data
 *  @param keep_missing alignments which do not have tag
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
plp_data calculate_pileup(const char *region, const char *bam_file, size_t num_dtypes, char *dtypes[], const char tag_name[2], const int tag_value, const _Bool keep_missing); 
