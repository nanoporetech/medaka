#!/bin/bash
set -e -x
export MANYLINUX=1

PACKAGE_NAME='medaka'

cd /io

yum install -y zlib-devel bzip2 bzip2-devel xz-devel
rm -rf libhts.a bincache/*
make scripts/mini_align clean_htslib libhts.a


# Compile wheels
for minor in 4 5 6; do
    PYBIN="/opt/python/cp3${minor}-cp3${minor}m/bin"
    # auditwheel/issues/102
    "${PYBIN}/pip" install --upgrade cffi setuptools pip wheel==0.31.1
    "${PYBIN}/pip" wheel . -w wheelhouse/
done


# Bundle external shared libraries into the wheels
for whl in "wheelhouse/${PACKAGE_NAME}"*.whl; do
    auditwheel repair "${whl}" -w wheelhouse/
done


# Install packages
for minor in 4 5 6; do
    PYBIN="/opt/python/cp3${minor}-cp3${minor}m/bin"
    "${PYBIN}/pip" install "${PACKAGE_NAME}" --no-index -f /io/wheelhouse
    "${PYBIN}/medaka_counts" --print medaka/test/data/test_reads.bam Consensus_Consensus_Consensus_Consensus_utg000001l:10000-10010
done

cd wheelhouse && ls | grep -v "${PACKAGE_NAME}.*manylinux" | xargs rm
