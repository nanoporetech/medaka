import itertools
import os
import platform

from cffi import FFI

dir_path = os.path.dirname(os.path.realpath(__file__))

#samver is pulled from this file in the Makefile
samver = "1.14"
htslib_dir = os.path.join(dir_path, 'submodules', 'samtools-{}'.format(samver), 'htslib-{}'.format(samver))
deflatever = "1.10"
deflate_dir = os.path.join(dir_path, 'submodules', 'libdeflate-{}'.format(deflatever))

libraries=['m', 'z', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
library_dirs=[htslib_dir]
if os.getenv('WITHDEFLATE') == "1":
    print("Using deflate")
    libraries.append('deflate')
    library_dirs.append(deflate_dir)
src_dir='src'

extra_compile_args = ['-std=c99', '-O3']
extra_link_args = []
if platform.machine() in {"aarch64", "arm64"}:
    if platform.system() == "Darwin":
        pass
    else:
        extra_compile_args.append("-march=armv8-a+simd")
else:
    extra_compile_args.append("-mtune=haswell")

if os.getenv('MEDAKA_COMPILE_ARGS') is not None:
    extra_compile_args.append(os.getenv('MEDAKA_COMPILE_ARGS'))
if os.getenv('MEDAKA_LINK_ARGS') is not None:
    extra_link_args.append(os.getenv('MEDAKA_LINK_ARGS'))

print("Extra compile args: ", extra_compile_args)
print("Extra link args: ", extra_link_args)


ffibuilder = FFI()
ffibuilder.set_source("libmedaka",
    r"""
    #include "kvec.h"
    #include "medaka_bamiter.h"
    #include "medaka_common.h"
    #include "medaka_counts.h"
    #include "medaka_trimbam.h"
    #include "medaka_pytrimbam.h"
    #include "medaka_read_matrix.h"
    #include "medaka_rnn_variants.h"
    #include "fastrle.h"
    #include "kseq.h"

    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=[src_dir, htslib_dir],
    sources=[
        os.path.join(src_dir, x) for x in (
            'medaka_bamiter.c', 'medaka_common.c', 'medaka_counts.c',
            'fastrle.c', 'medaka_trimbam.c', 'medaka_pytrimbam.c',
            'medaka_read_matrix.c',
            'medaka_rnn_variants.c')],
    extra_compile_args=extra_compile_args,
    extra_objects=['libhts.a'],
    extra_link_args=extra_link_args,
)

cdef = [
    "typedef struct { ...; } bam_fset;"
    "bam_fset* create_bam_fset(char* fname);"
    "void destroy_bam_fset(bam_fset* fset);"
]
for header in ('medaka_counts.h','fastrle.h', 'medaka_pytrimbam.h', 'medaka_read_matrix.h', 'medaka_rnn_variants.h'):
    with open(os.path.join(src_dir, header), 'r') as fh:
        # remove directives
        lines = ''.join(x for x in fh.readlines() if not x.startswith('#'))
        cdef.append(lines)

ffibuilder.cdef('\n\n'.join(cdef))


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
