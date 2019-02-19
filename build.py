import os

from cffi import FFI

#TODO: configure this better
samver="1.9"
htslib_dir=os.path.join('submodules', 'samtools-{}'.format(samver), 'htslib-{}'.format(samver))

libraries=['z', 'lzma', 'bz2', 'pthread', 'curl', 'crypto']
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

with open(os.path.join(src_dir, 'medaka_counts.h'), 'r') as fh:
    ffibuilder.cdef(fh.read())


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
