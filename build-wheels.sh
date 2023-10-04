#!/bin/bash
# Usage: ./build-wheels.sh <workdir> <pyminorversion1> <pyminorversion2> ...
set -e -x
export MANYLINUX=1
export MEDAKA_DIST=1
export WITHDEFLATE=1

arch=$(uname -m)

workdir=$1
shift

echo "Changing cwd to ${workdir}"
cd ${workdir}

if [[ -z "${PACKAGE_NAME}" ]]; then
    PACKAGE_NAME='medaka'
fi
PACKAGE_FILE_NAME=${PACKAGE_NAME//-/_}
sed -i "s/__dist_name__ = 'medaka'/__dist_name__ = '${PACKAGE_NAME}'/" setup.py


# some many linux containers are centos-based, others are debian!
if [ -f /etc/centos-release ]; then
    yum install -y zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel ncurses-devel
else
    apt update
    apt install -y zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-gnutls-dev libssl-dev libffi-dev
fi

rm -rf libhts.a bincache/*
make scripts/mini_align clean libhts.a
mkdir -p wheelhouse
# this only exists after building libhts
LIBDEFLATE=$(ls -d ${PWD}/submodules/libdeflate-*/)
echo "LIBDEFLATE = '${LIBDEFLATE}'"

echo "PYTHON VERSIONS AVAILABLE"
ls /opt/python/

# Compile wheels
for minor in $@; do
    if [[ "${minor}" == "8" ]]  || [[ "${minor}" == "9" ]] || [[ "${minor}" == "10" ]]; then
        PYBIN="/opt/python/cp3${minor}-cp3${minor}/bin"
    else
        PYBIN="/opt/python/cp3${minor}-cp3${minor}m/bin"
    fi
    # auditwheel/issues/102
    "${PYBIN}"/pip install --upgrade setuptools pip wheel==0.31.1 cffi==1.15.0
    "${PYBIN}"/pip wheel --prefer-binary --no-dependencies . -w ./wheelhouse/
done


# Bundle external shared libraries into the wheels
export LD_LIBRARY_PATH=$PWD/libdeflate
for whl in "wheelhouse/${PACKAGE_FILE_NAME}"*.whl; do
    LD_LIBRARY_PATH=${LIBDEFLATE} auditwheel repair "${whl}" -w ./wheelhouse/
done
unset LD_LIBRARY_PATH


## Install packages
if [[ "${DO_COUNT_TEST}" == "1" ]]; then
    for minor in $@; do
        if [[ "${minor}" == "8" || "${minor}" == "9" ]] || [[ "${minor}" == "10" ]]; then
            PYBIN="/opt/python/cp3${minor}-cp3${minor}/bin"
        else
            PYBIN="/opt/python/cp3${minor}-cp3${minor}m/bin"
        fi
        "${PYBIN}"/pip install --use-pep517 --prefer-binary -r requirements.txt
        "${PYBIN}"/pip install "${PACKAGE_NAME}" --no-index -f ./wheelhouse
        "${PYBIN}"/medaka_counts --print medaka/test/data/test_reads.bam utg000001l:10000-10010
    done
fi

mkdir wheelhouse-final
cp wheelhouse/${PACKAGE_FILE_NAME}*manylinux* wheelhouse-final
