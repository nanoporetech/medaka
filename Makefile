.PHONY: docs

# Builds a cache of binaries which can just be copied for CI
BINARIES=samtools minimap2 mini_align vcf2fasta
BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif

binaries: $(addprefix $(BINCACHEDIR)/, $(BINARIES))

SAMVER=1.3.1
$(BINCACHEDIR)/samtools: | $(BINCACHEDIR)
	# TODO: make this a bit nicer, we're only doing this for tview
	@echo Making $(@F)
	# tar.bz is not a dependency, since that would cause it to be fetched
	#   even when installing from $(BINCACHEDIR)
	if [ ! -e submodules/samtools-${SAMVER}.tar.bz2 ]; then \
	  cd submodules; \
	  wget https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
	fi
	cd submodules && tar -xjf samtools-${SAMVER}.tar.bz2
	# copy our hack up version of tview
	${SEDI} 's/tv->is_dot = 1;/tv->is_dot = 0;/' submodules/samtools-${SAMVER}/bam_tview.c
	cd submodules/samtools-${SAMVER} && make
	cp submodules/samtools-${SAMVER}/samtools $@

$(BINCACHEDIR)/minimap2: | $(BINCACHEDIR)
	@echo Making $(@F)
	wget https://github.com/lh3/minimap2/releases/download/v2.11/minimap2-2.11_x64-linux.tar.bz2 
	tar -xvf minimap2-2.11_x64-linux.tar.bz2
	cp minimap2-2.11_x64-linux/minimap2 $@
	rm -rf minimap2-2.11_x64-linux.tar.bz2 minimap2-2.11_x64-linux

$(BINCACHEDIR)/mini_align: | $(BINCACHEDIR)
	@echo Making $(@F)
	curl https://raw.githubusercontent.com/nanoporetech/pomoxis/master/scripts/mini_align -o $@
	chmod +x $@


$(BINCACHEDIR)/vcf2fasta: | $(BINCACHEDIR)
	cd src/vcf2fasta && g++ -std=c++11 \
		-I./../../submodules/samtools-${SAMVER}/htslib-${SAMVER}/ vcf2fasta.cpp \
		./../../submodules/samtools-${SAMVER}/htslib-${SAMVER}/libhts.a \
		-lz -llzma -lbz2 -lpthread \
		-o $(@F)
	cp src/vcf2fasta/$(@F) $@


libhts.a: $(BINCACHEDIR)/samtools
	# this is required only to add in -fpic so we can build python module
	@echo Compiling $(@F)
	cd submodules/samtools-${SAMVER}/htslib-${SAMVER}/ && CFLAGS=-fpic ./configure && make
	cp submodules/samtools-${SAMVER}/htslib-${SAMVER}/$@ $@


venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=python3 --prompt "(medaka) "
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

install: venv libhts.a | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && python setup.py install

test: install
	${IN_VENV} && pip install nose
	${IN_VENV} && python setup.py nosetests

clean:
	. ${VENV} && python setup.py clean || echo "Failed to run setup.py clean"
	rm -rf libhts.a venv build dist/ medaka.egg-info/ __pycache__ medaka.egg-info
	find . -name '*.pyc' -delete

docker: binaries libhts.a
	mkdir for_docker && cp -r medaka scripts bincache setup.py build.py requirements.txt for_docker 
	docker build -t medaka .
	rm -rf for_docker

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

DOCSRC = docs

docs: venv
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	rm -rf docs/modules.rst docs/medaka.rst  
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll


