import itertools
import os

from cffi import FFI

#TODO: configure this better
samver="1.10"
htslib_dir=os.path.join('submodules', 'samtools-{}'.format(samver), 'htslib-{}'.format(samver))

libraries=['m', 'z', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
library_dirs=[htslib_dir]
src_dir='src'

ffibuilder = FFI()
ffibuilder.set_source("libmedaka",
    r"""
    #include "kvec.h"
    #include "medaka_bamiter.h"
    #include "medaka_common.h"
    #include "medaka_counts.h"
    #include "medaka_trimbam.h"
    #include "medaka_pytrimbam.h"
    #include "fastrle.h"
    #include "kseq.h"

    """,
    libraries=libraries,
    library_dirs=library_dirs,
    include_dirs=[src_dir, htslib_dir],
    sources=[os.path.join(src_dir, x) for x in ('medaka_bamiter.c', 'medaka_common.c', 'medaka_counts.c', 'fastrle.c', 'medaka_trimbam.c', 'medaka_pytrimbam.c')],
    extra_compile_args=['-std=c99', '-msse3', '-O3'],
    extra_objects=['libhts.a']
)

cdef = []
for header in ('medaka_counts.h','fastrle.h', 'medaka_pytrimbam.h'):
    with open(os.path.join(src_dir, header), 'r') as fh:
        # remove directives
        lines = ''.join(x for x in fh.readlines() if not x.startswith('#'))
        cdef.append(lines)

ffibuilder.cdef('\n\n'.join(cdef))


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
