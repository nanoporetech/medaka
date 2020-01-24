#!/bin/bash
# Usage: ./build-wheels.sh <workdir> <pyminorversion1> <pyminorversion2> ...
set -e -x
export MANYLINUX=1

PACKAGE_NAME='medaka'
parasail=parasail-1.1.17-py2.py3-none-linux_x86_64.whl

workdir=$1
shift

echo "Changing cwd to ${workdir}"
cd ${workdir}

yum install -y zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel
rm -rf libhts.a bincache/*
make scripts/mini_align clean libhts.a
mkdir -p wheelhouse
curl -L https://ont-research.s3-eu-west-1.amazonaws.com/${parasail} \
    -o ${parasail}

# Compile wheels
for minor in $@; do
    PYBIN="/opt/python/cp3${minor}-cp3${minor}m/bin"
    # auditwheel/issues/102
    "${PYBIN}/pip" install --upgrade cffi setuptools pip wheel==0.31.1
    "${PYBIN}/pip" install ${parasail}
    "${PYBIN}/pip" wheel . -w ./wheelhouse/
done


# Bundle external shared libraries into the wheels
for whl in "wheelhouse/${PACKAGE_NAME}"*.whl; do
    auditwheel repair "${whl}" -w ./wheelhouse/
done


# Install packages
for minor in $@; do
    PYBIN="/opt/python/cp3${minor}-cp3${minor}m/bin"
    "${PYBIN}/pip" install ${parasail}
    "${PYBIN}/pip" install "${PACKAGE_NAME}" --no-index -f ./wheelhouse
    "${PYBIN}/medaka_counts" --print medaka/test/data/test_reads.bam utg000001l:10000-10010
done

cd wheelhouse && ls | grep -v "${PACKAGE_NAME}.*manylinux" | xargs rm
