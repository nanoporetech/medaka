
![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)


Medaka
======

[![Build Status](https://travis-ci.org/nanoporetech/medaka.svg?branch=master)](https://travis-ci.org/nanoporetech/medaka)

[![](https://img.shields.io/pypi/v/medaka.svg)](https://pypi.org/project/medaka/)
[![](https://img.shields.io/pypi/wheel/medaka.svg)](https://pypi.org/project/medaka/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/medaka)
[![](https://img.shields.io/conda/pn/bioconda/medaka.svg)](https://anaconda.org/bioconda/medaka)



`medaka` is a tool to create a consensus sequence from nanopore sequencing data.
This task is performed using neural networks applied from a pileup of individual
sequencing reads against a draft assembly. It outperforms graph-based methods
operating on basecalled data, and can be competitive with state-of-the-art
signal-based methods, whilst being much faster.

© 2018 Oxford Nanopore Technologies Ltd.

Features
--------

  * Requires only basecalled data. (`.fasta` or `.fastq`)
  * Improved accurary over graph-based methods (e.g. Racon).
  * 50X faster than Nanopolish (and can run on GPUs).
  * Benchmarks are provided [here](https://nanoporetech.github.io/medaka/benchmarks.html).
  * Includes extras for implementing and training bespoke correction
    networks.
  * Works on Linux and MacOS.
  * Open source (Mozilla Public License 2.0).

Tools to enable the creation of draft assemblies can be found in a sister
project [pomoxis](https://github.com/nanoporetech/pomoxis).

Documentation can be found at https://nanoporetech.github.io/medaka/.


Installation
------------

Medaka can be installed in one of several ways.

**Installation with conda**

Perhaps the simplest way to start using medaka on both Linux and MacOS is
through conda; medaka is available via the
[bioconda](https://anaconda.org/bioconda/medaka) channel:

    conda install -c bioconda medaka

**Installation with pip**
  
For those who prefer python's native pacakage manager, medaka is also available
on pypi and can be installed using pip:

    pip install medaka

On Linux platforms this will install a precompiled binary, on MacOS (and other)
platforms this will fetch and compile a source distribution.

We recommend using medaka within a virtual environment, viz.:

    virtualenv medaka --python=python3 --prompt "(medaka) "
    . medaka/bin/activate
    pip install medaka

Using this method requires the user to provide several binaries:

 * [samtools](https://github.com/samtools/samtools),
 * [minimap2](https://github.com/lh3/minimap2),
 * [tabix](https://github.com/samtools/htslib), and
 * [bgzip](https://github.com/samtools/htslib)

and place these within the `PATH`. `samtools/bgzip/tabix` version 1.9 and
`minimap2` version 2.11 are recommended as these are those used in development
of medaka.

**Installation from source**

Medaka can be installed from its source quite easily on most systems.

 > Before installing medaka it may be required to install some
 > prerequisite libraries, best installed by a package manager. On Ubuntu
 > theses are:
 > * gcc
 > * zlib1g-dev
 > * libbz2-dev
 > * liblzma-dev
 > * libffi-dev
 > * libncurses5-dev
 > * libcurl4-gnutls-dev
 > * libssl-dev
 > * make
 > * wget
 > * python3-all-dev
 > * python-virtualenv

A Makefile is provided to fetch, compile and install all direct dependencies
into a python virtual environment. To set-up the environment run:

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

Using this method both `samtools` and `minimap2` are built from source and need
not be provided by the user.

**Using a GPU**

All installation methods will allow medaka to be used with CPU resource only.
To enable the use of GPU resource it is necessary to install the
`tensorflow-gpu` package. To outline, this can be achieved with:

    pip uninstall tensorflow
    pip install tensorflow-gpu

However, note that The `tensorflow-gpu` GPU package is compiled against a
specific version of the NVIDIA CUDA library; users are directed to the 
[tensorflow installation](https://www.tensorflow.org/install/gpu) pages
for further information.


Usage
-----

`medaka` can be run using its default settings through the `medaka_consensus`
program. An assembly in `.fasta` format and basecalls in `.fasta` or `.fastq`
formats are required. The program uses both `samtools` and `minimap2`. If
medaka has been installed using the from-source method these will be present
within the medaka environment, otherwise they will need to be provided by the user.

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m r94

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. When `medaka_consensus` has finished running, the consensus
will be saved to `${OUTDIR}/consensus.fasta`.

   **It is crucially important to specify the correct model, `-m` in the
   above, according to the basecaller used. Allowed values can be found by
   running `medaka consensus --help`. The default model is appropriate for
   basecallers using the transducer algorithm (Albacore or Guppy<2.1.3). For
   Guppy versions >=2.1.3 where the flip-flop algorithm has been used, users
   should select the highest numbered model equal to or less than the Guppy
   version used for basecalling.**

### Origin of the draft sequence

Medaka has been trained to correct draft sequences processed through
[racon](https://github.com/isovic/racon), specifically `racon` run four times
iteratively with:

    racon -m 8 -x -6 -g -8 -w 500 ...

Processing a draft sequence from alternative sources (e.g. the output of
[canu](https://github.com/marbl/canu) or
[wtdbg2](https://github.com/ruanjue/wtdbg2)) may lead to poorer results
even when the draft is of a superior quality than that obtained from `racon`.

The [walkthrough](https://nanoporetech.github.io/medaka/walkthrough.html#walkthrough)
outlines one recommended workflow rapid construction of a draft for input into
`medaka`. A second approach would be to run `canu` followed by `racon` applied
twice iteratively before entry into `medaka`.


Acknowledgements
----------------

We thank [Joanna Pineda](https://github.com/jopineda) and
[Jared Simpson](https://github.com/jts) for providing htslib code samples which aided
greatly development of the optimised feature generation code, and for testing the
version 0.4 release candidates.

Help
----

**Licence and Copyright**

© 2018 Oxford Nanopore Technologies Ltd.

`medaka` is distributed under the terms of the Mozilla Public License 2.0.

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the 
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.
