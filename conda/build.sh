export INCLUDE_PATH="${PREFIX}/include"
export LIBRARY_PATH="${PREFIX}/lib"
export LDFLAGS="${LDFLAGS} -L${PREFIX}/lib"

export CFLAGS="-I${PREFIX}/include ${LDFLAGS}"

# disable Makefile driven build of htslib.a
sed -i.bak "s/'build_ext': HTSBuild//" setup.py

# just link to htslib
sed -i.bak 's/extra_objects.*//' build.py
sed -i.bak 's/^libraries=\[/libraries=\["hts",/' build.py

export MEDAKA_DIST=1
export SKIP_LFS_CHECK=1
#export CONDA_OVERRIDE_CUDA="11.8"
make scripts/mini_align
${PYTHON} -m pip install . --no-deps --ignore-installed -vvv
