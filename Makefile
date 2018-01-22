.PHONY: docs

# Builds a cache of binaries which can just be copied for CI
BINARIES=samtools
BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)
OS := $(shell uname)
ifeq ($(OS), Darwin)
SEDI=sed -i '.bak'
else
SEDI=sed -i
endif


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


venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install numpy # needs to get done before other things
	${IN_VENV} && pip install -r requirements.txt

install: venv | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && python setup.py install

test: install
	${IN_VENV} && python -m unittest discover

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
