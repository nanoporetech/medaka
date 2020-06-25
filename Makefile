
# Builds a cache of binaries which can just be copied for CI
BINARIES=samtools minimap2 tabix bgzip racon bcftools
BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif

PYTHON ?= python3

binaries: $(addprefix $(BINCACHEDIR)/, $(BINARIES))

SAMVER=1.10
submodules/samtools-$(SAMVER)/Makefile:
	cd submodules; \
		curl -L -o samtools-${SAMVER}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
		tar -xjf samtools-${SAMVER}.tar.bz2; \
		rm samtools-${SAMVER}.tar.bz2


libhts.a: submodules/samtools-$(SAMVER)/Makefile
	# this is required only to add in -fpic so we can build python module
	@echo Compiling $(@F)
	cd submodules/samtools-${SAMVER}/htslib-${SAMVER}/ && CFLAGS=-fpic ./configure && make
	cp submodules/samtools-${SAMVER}/htslib-${SAMVER}/$@ $@


$(BINCACHEDIR)/samtools: | libhts.a $(BINCACHEDIR)
	@echo Making $(@F)
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
	@echo Compiling $(@F)
	cd submodules; \
		curl -L -o minimap2-${MINIMAPVER}.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAPVER}/minimap2-${MINIMAPVER}.tar.bz2; \
		tar -xjf minimap2-${MINIMAPVER}.tar.bz2; \
	    rm -rf minimap2-${MINIMAPVER}.tar.bz2
	cd submodules/minimap2-${MINIMAPVER} && make
	cp submodules/minimap2-${MINIMAPVER}/minimap2 $@


$(BINCACHEDIR)/bcftools: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -d submodules/bcftools-v${SAMVER} ]; then \
		cd submodules; \
		curl -L -o bcftools-v${SAMVER}.tar.bz2 https://github.com/samtools/bcftools/releases/download/${SAMVER}/bcftools-${SAMVER}.tar.bz2; \
		tar -xjf bcftools-v${SAMVER}.tar.bz2; \
		cd bcftools-${SAMVER}; \
		make; \
	fi
	cp submodules/bcftools-${SAMVER}/bcftools $@


RACONVER=1.3.1
$(BINCACHEDIR)/racon: | $(BINCACHEDIR)
	@echo Making $(@F)
	@echo GCC is $(GCC)
	if [ ! -e submodules/racon-v${RACONVER}.tar.gz ]; then \
	  cd submodules; \
	  curl -L -o racon-v${RACONVER}.tar.gz https://github.com/isovic/racon/releases/download/${RACONVER}/racon-v${RACONVER}.tar.gz; \
	  tar -xzf racon-v${RACONVER}.tar.gz; \
	fi
	cd submodules/racon-v${RACONVER}; \
		rm -rf build; \
		mkdir build; \
		cd build; \
		cmake -DCMAKE_BUILD_TYPE=Release ..; \
		make;
	cp submodules/racon-v${RACONVER}/build/bin/racon $@


$(BINCACHEDIR)/vcf2fasta: | $(BINCACHEDIR)
	cd src/vcf2fasta && g++ -std=c++11 \
		-I./../../submodules/samtools-${SAMVER}/htslib-${SAMVER}/ vcf2fasta.cpp \
		./../../submodules/samtools-${SAMVER}/htslib-${SAMVER}/libhts.a \
		-lz -llzma -lbz2 -lpthread \
		-o $(@F)
	cp src/vcf2fasta/$(@F) $@


scripts/mini_align:
	@echo Making $(@F)
	curl https://raw.githubusercontent.com/nanoporetech/pomoxis/master/scripts/mini_align -o $@
	chmod +x $@


venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=$(PYTHON) --prompt "(medaka) "
	${IN_VENV} && pip install pip --upgrade


.PHONY: check_lfs
check_lfs: venv
	${IN_VENV} && python -c "from setup import check_model_lfs; check_model_lfs()"


.PHONY: install
install: venv check_lfs scripts/mini_align libhts.a | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && pip install -r requirements.txt
	${IN_VENV} && MEDAKA_BINARIES=1 python setup.py install


.PHONY: test
test: install
	${IN_VENV} && pip install pytest pytest-cov flake8 flake8-rst-docstrings flake8-docstrings flake8-import-order
	# TODO: add these exclusions back in after outstanding PRs
	${IN_VENV} && flake8 medaka --import-order-style google --application-import-names medaka,libmedaka --exclude \
		medaka/test/,medaka/medaka.py \
		--statistics
	${IN_VENV} && pytest medaka --doctest-modules \
		--cov=medaka --cov-report html --cov-report term \
		--cov-fail-under=80 --cov-report term-missing


.PHONY: clean
clean: clean_htslib
	(${IN_VENV} && python setup.py clean) || echo "Failed to run setup.py clean"
	rm -rf libhts.a libmedaka.abi3.so venv build dist/ medaka.egg-info/ __pycache__ medaka.egg-info
	find . -name '*.pyc' -delete


.PHONY: mem_check
mem_check: install pileup
	${IN_VENV} && python -c "import medaka.test.test_counts as tc; tc.create_simple_bam('mem_test.bam', tc.simple_data['calls'])"
	valgrind --error-exitcode=1 --tool=memcheck ./pileup mem_test.bam ref:1-8 || (ret=$$?; rm mem_test.bam* && exit $$ret)
	rm -rf mem_test.bam*


pileup: libhts.a
	gcc -pthread  -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Isubmodules/samtools-${SAMVER}/htslib-${SAMVER} \
		src/medaka_common.c src/medaka_counts.c src/medaka_bamiter.c libhts.a \
		-lm -lz -llzma -lbz2 -lpthread -lcurl -lcrypto \
		-o $(@) -std=c99 -msse3 -O3


trim_reads: libhts.a
	gcc -pthread -pg -g -Wall -fstack-protector-strong -D_FORTIFY_SOURCE=2 -fPIC -std=c99 -msse3 -O3 \
		-Isrc -Isubmodules/samtools-${SAMVER}/htslib-${SAMVER} \
		src/medaka_common.c src/medaka_trimbam.c src/medaka_bamiter.c libhts.a \
		-lz -llzma -lbz2 -lpthread -lcurl -lcrypto \
		-o $(@) -std=c99 -msse3 -O3


.PHONY: wheels
wheels:
	docker run -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/build-wheels.sh /io 5 6

.PHONY: build
build: pypi_build/bin/activate
IN_BUILD=. ./pypi_build/bin/activate
pypi_build/bin/activate:
	test -d pypi_build || virtualenv pypi_build --python=python3 --prompt "(pypi) "
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
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	rm -rf docs/modules.rst docs/medaka.rst  
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll


