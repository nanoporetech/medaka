// medaka-style feature data
typedef struct {
    size_t n_cols;
    size_t *counts;
    size_t *major;
    size_t *minor;
} _plp_data;

typedef _plp_data *plp_data; 

plp_data create_plp_data(size_t n_cols);
void destroy_plp_data(plp_data data);
void print_pileup_data(plp_data pileup);
plp_data calculate_pileup(const char *region, const char *bam_file); 
