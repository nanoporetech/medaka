import os

from cffi import FFI

#TODO: configure this better
htslib_dir=os.path.join('submodules', 'samtools-1.3.1', 'htslib-1.3.1')

libraries=['z', 'lzma', 'bz2', 'pthread']
library_dirs=[htslib_dir]
src_dir='src'

ffibuilder = FFI()
ffibuilder.set_source("libmedaka",
    r"""
    #include "medaka_counts.h"
    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=[src_dir, htslib_dir],
    sources=[os.path.join(src_dir, 'medaka_counts.c')],
    extra_compile_args=['-std=c99', '-msse3', '-O3'],
    extra_objects=['libhts.a']
)


ffibuilder.cdef("""

typedef struct {
    size_t n_cols;
    size_t *counts;
    size_t *major;
    size_t *minor;
} _plp_data;

typedef _plp_data *plp_data; 

plp_data calculate_pileup(const char *region, const char *bam_file);
void print_pileup_data(plp_data pileup);
plp_data create_plp_data(size_t n_cols);
void destroy_plp_data(plp_data data);

"""
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
