# Builds a cache of binaries which can just be copied for CI
BINARIES=samtools minimap2 tabix bgzip bcftools
BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)
OS := $(shell uname)
ARCH := $(shell arch)

ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
LDFLAGS="-L/opt/homebrew/opt/openssl@3/lib"
else
SEDI=sed -i
endif

ifeq ($(ARCH), $(filter $(ARCH), aarch64 arm64))
TUNE=
MM2ARGS=arm_neon=1 aarch64=1
ifeq ($(OS), Darwin) # macos arm
export HDF5_DIR=/opt/homebrew/Cellar/hdf5/1.12.1
.package-reqs: pyprep-m1
endif

else # x64
TUNE=-mtune=haswell
MM2ARGS2=
LDFLAGS=
endif

PYTHON ?= python3
VENV ?= venv
DOCKERTAG ?= "medaka/medaka:latest"
COVFAIL = 80

VERSION := $(shell grep "__version__" medaka/__init__.py | awk '{gsub("\"","",$$3); print $$3}')

CFLAGS ?= -fpic -O3 -std=c99
LIBS ?= -lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto
STATIC_HTSLIB ?= htslib/libhts.a
HTS_CONF_ARGS ?=
HTS_CONF_ENV ?= CFLAGS="$(CFLAGS) $(TUNE)"

# set WITHDFLATE=1 to build with deflate. Enabling this will
# require setting LD_LIBRARY_PATH=submodules/libdeflate-<ver>/
# at runtime. This doesn't apply to wheel which will bundle
# the library (as curl, crypto, lzma, etc)
WITHDEFLATE ?= 
DEFLATEVER=$(shell sed -n 's/deflatever = "\(.*\)"/\1/p' build.py)
DEFLATE = $(PWD)/submodules/libdeflate-${DEFLATEVER}
DEFLATEREQ =
ifeq ($(WITHDEFLATE), 1)
CFLAGS += -I$(DEFLATE) -L$(DEFLATE)
LIBS += -ldeflate
HTS_CONF_ARGS += --with-libdeflate
HTS_CONF_ENV += LDFLAGS="-L$(DEFLATE)"
DEFLATEREQ = submodules/libdeflate-${DEFLATEVER}/libdeflate.so.0
endif

submodules/libdeflate-$(DEFLATEVER)/libdeflate.so.0:
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd submodules \
		&& curl -L -o libdeflate-v${DEFLATEVER}.tar.gz https://github.com/ebiggers/libdeflate/archive/refs/tags/v${DEFLATEVER}.tar.gz \
		&& tar -xzf libdeflate-v${DEFLATEVER}.tar.gz \
		&& cd libdeflate-${DEFLATEVER} \
		&& make


SAMVER=$(shell sed -n 's/samver = "\(.*\)"/\1/p' build.py)
submodules/samtools-$(SAMVER)/Makefile:
	cd submodules; \
		curl -L -o samtools-${SAMVER}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
		tar -xjf samtools-${SAMVER}.tar.bz2; \
		rm samtools-${SAMVER}.tar.bz2


libhts.a: submodules/samtools-$(SAMVER)/Makefile $(DEFLATEREQ)
	# this is required only to add in -fpic so we can build python module
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd submodules/samtools-${SAMVER}/htslib-${SAMVER}/ \
		&& $(HTS_CONF_ENV) ./configure $(HTS_CONF_ARGS) \
		&& make -j 4
	cp submodules/samtools-${SAMVER}/htslib-${SAMVER}/$@ $@


$(BINCACHEDIR)/samtools: | libhts.a $(BINCACHEDIR)
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	# copy our hack up version of tview
	${SEDI} 's/tv->is_dot = 1;/tv->is_dot = 0;/' submodules/samtools-${SAMVER}/bam_tview.c
	cd submodules/samtools-${SAMVER} && make -j 4
	cp submodules/samtools-${SAMVER}/$(@F) $@


$(BINCACHEDIR)/tabix: | libhts.a $(BINCACHEDIR)
	cp submodules/samtools-${SAMVER}/htslib-${SAMVER}/$(@F) $@


$(BINCACHEDIR)/bgzip: | libhts.a $(BINCACHEDIR)
	cp submodules/samtools-${SAMVER}/htslib-${SAMVER}/$(@F) $@


.PHONY: clean_htslib
clean_htslib:
	cd submodules/samtools-${SAMVER} && make clean || exit 0
	cd submodules/samtools-${SAMVER}/htslib-${SAMVER} && make clean || exit 0


MINIMAPVER=2.17
$(BINCACHEDIR)/minimap2: | $(BINCACHEDIR)
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	cd submodules; \
		curl -L -o minimap2-${MINIMAPVER}.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAPVER}/minimap2-${MINIMAPVER}.tar.bz2; \
		tar -xjf minimap2-${MINIMAPVER}.tar.bz2; \
	    rm -rf minimap2-${MINIMAPVER}.tar.bz2
	cd submodules/minimap2-${MINIMAPVER} && make ${MM2ARGS}
	cp submodules/minimap2-${MINIMAPVER}/minimap2 $@


$(BINCACHEDIR)/bcftools: | $(BINCACHEDIR)
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	if [ ! -d submodules/bcftools-v${SAMVER} ]; then \
		cd submodules; \
		curl -L -o bcftools-v${SAMVER}.tar.bz2 https://github.com/samtools/bcftools/releases/download/${SAMVER}/bcftools-${SAMVER}.tar.bz2; \
		tar -xjf bcftools-v${SAMVER}.tar.bz2; \
		cd bcftools-${SAMVER}; \
		make; \
	fi
	cp submodules/bcftools-${SAMVER}/bcftools $@


$(BINCACHEDIR)/vcf2fasta: | $(BINCACHEDIR)
	cd src/vcf2fasta && g++ -std=c++11 \
		-I./../../submodules/samtools-${SAMVER}/htslib-${SAMVER}/ vcf2fasta.cpp \
		./../../submodules/samtools-${SAMVER}/htslib-${SAMVER}/libhts.a \
		$(CFLAGS) $(LDFLAGS) $(LIBS) \
		-o $(@F)
	cp src/vcf2fasta/$(@F) $@


scripts/mini_align:
	@echo "\x1b[1;33mMaking $(@F)\x1b[0m"
	curl https://raw.githubusercontent.com/nanoporetech/pomoxis/master/scripts/mini_align -o $@
	chmod +x $@


venv: ${VENV}/bin/activate
IN_VENV=. ./${VENV}/bin/activate

$(VENV)/bin/activate:
	test -d $(VENV) || $(PYTHON) -m venv $(VENV) --prompt "medaka"
	${IN_VENV} && pip install pip wheel --upgrade

$(VENV)/bin/%: $(BINCACHEDIR)/%
	cp $< $@


.PHONY: check_lfs
check_lfs: venv
	${IN_VENV} && python -c "from setup import check_model_lfs; check_model_lfs()"

.PHONY: pyprep
pyprep-m1: venv
	@echo "\x1b[1;33mInstalling prerequisites with homebrew\x1b[0m"
	brew install pkg-config hdf5@1.12 openssl@3
	@echo "\x1b[1;33mCompiling and installing mappy from custom repository\x1b[0m"
	${IN_VENV} && pip install git+https://github.com/cjw85/minimap2.git@1951fe908ae83efd4f459e060d29b21254d54993


.PHONY: .package-reqs
.package-reqs: venv check_lfs scripts/mini_align libhts.a | $(addprefix $(VENV)/bin/, $(BINARIES))


.PHONY: install
install: .package-reqs
	@echo "\x1b[1;33mInstalling medaka\x1b[0m"
	${IN_VENV} && LDFLAGS=$(LDFLAGS) pip install .


.PHONY: develop
develop: .package-reqs 
	@echo "\x1b[1;33mInstalling medaka in development mode\x1b[0m"
	${IN_VENV} && WITHDEFLATE=$(WITHDEFLATE) LDFLAGS=$(LDFLAGS) pip install -e .


.PHONY: test
test: install
	@echo "\x1b[1;33mRunning tests on Python package\x1b[0m"
	${IN_VENV} && pip install pytest pytest-cov flake8 flake8-rst-docstrings flake8-docstrings flake8-import-order flake8-forbid-visual-indent
	# TODO: add these exclusions back in after outstanding PRs
	${IN_VENV} && flake8 medaka --import-order-style google --application-import-names medaka,libmedaka --exclude \
		medaka/test/,medaka/medaka.py \
		--statistics
	${IN_VENV} && pytest medaka --doctest-modules \
		--cov=medaka --cov-report html --cov-report term \
		--cov-fail-under=${COVFAIL} --cov-report term-missing


# mainly here for the Dockerfile
.PHONY: install_root
install_root: check_lfs scripts/mini_align libhts.a | $(addprefix $(VENV)/bin/, $(BINARIES)) 
	pip3 install pip --upgrade
	pip3 install -r requirements.txt
	pip3 install .
	# copy binaries to PATH
	cp $| $$(dirname $$(which pip3))


.PHONY: docker
docker: clean
	@echo "\x1b[1;33mBuilding docker tag: '$(DOCKERTAG)'\x1b[0m"
	docker build -t $(DOCKERTAG) .


.PHONY: clean
clean: clean_htslib
	(${IN_VENV} && python setup.py clean) || echo "Failed to run setup.py clean"
	rm -rf libhts.a libmedaka.abi3.so venv* build dist/ medaka.egg-info/ __pycache__ medaka.egg-info
	find . -name '*.pyc' -delete


.PHONY: mem_check
mem_check: pileup
	valgrind --error-exitcode=1 --tool=memcheck ./pileup medaka/test/data/test_reads.bam utg000001l:5000-5500 || (ret=$$?; rm mem_test.bam* && exit $$ret)
	rm -rf mem_test.bam*


pileup: libhts.a
	gcc -pthread  -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Isubmodules/samtools-${SAMVER}/htslib-${SAMVER} \
		src/medaka_common.c src/medaka_counts.c src/medaka_bamiter.c libhts.a \
		$(CFLAGS) $(LDFLAGS) $(LIBS) \
		-o $(@) -std=c99 -msse3 -O3


trim_reads: libhts.a
	gcc -pthread -pg -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Isubmodules/samtools-${SAMVER}/htslib-${SAMVER} \
		src/medaka_common.c src/medaka_trimbam.c src/medaka_bamiter.c libhts.a \
		$(CFLAGS) $(LDFLAGS) $(LIBS) \
		-o $(@) -std=c99 -msse3 -O3


.PHONY: wheels
wheels:
	docker run -v `pwd`:/io quay.io/pypa/manylinux2010_x86_64 /io/build-wheels.sh /io 7 8 9

.PHONY: build_env
build_env: pypi_build/bin/activate
IN_BUILD=. ./pypi_build/bin/activate
pypi_build/bin/activate:
	test -d pypi_build || $(PYTHON) -m venv pypi_build --prompt "pypi"
	${IN_BUILD} && pip install pip --upgrade
	${IN_BUILD} && pip install --upgrade pip setuptools twine wheel readme_renderer[md]

.PHONY: sdist
sdist: pypi_build/bin/activate scripts/mini_align submodules/samtools-$(SAMVER)/Makefile
	${IN_BUILD} && MEDAKA_DIST=1 python setup.py sdist


# Documentation
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
DOCSRC = docs

.PHONY: docs
docs: venv
	${IN_VENV} && pip install sphinx==3.5.4 sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	rm -rf docs/modules.rst docs/medaka.rst  
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll


