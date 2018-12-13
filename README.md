
![Oxford Nanopore Technologies logo](images/ONT_logo_590x106.png)


Medaka
======

[![Build Status](https://travis-ci.org/nanoporetech/medaka.svg?branch=master)](https://travis-ci.org/nanoporetech/medaka)

`medaka` is a tool to create a consensus sequence of nanopore sequencing data.
This task is performed using neural networks applied a pileup of individual
sequencing reads against a draft assembly. It outperforms graph-based methods
operating on basecalled data, and can be competitive with state-of-the-art
signal-based methods whilst being much faster.


Features
--------

  * Requires only basecalled data. (`.fasta` or `.fastq`)
  * Improved accurary over graph-based methods (e.g. Racon).
  * 50X faster than Nanopolish.
  * Benchmarks are provided [here](https://nanoporetech.github.io/medaka/benchmarks.html).
  * Includes extras for implementing and training bespoke correction
    networks.
  * Works on Linux (MacOS and Windows support is untested).
  * Open source (Mozilla Public License 2.0).

Documentation can be found at https://nanoporetech.github.io/medaka/.


Installation
------------

There are currently two installation methods for medaka, detailed below.

**Installation with pip**
  
Medaka can be installed on Linux using the python package manager, pip:

    pip install medaka

We recommend using medaka within a virtual environment, viz.:

    virtualenv medaka --python=python3 --prompt "(medaka) "
    . medaka/bin/activate
    pip install medaka

Using this method requires the user to provide a
[samtools](https://github.com/samtools/samtools) and
[minimap2](https://github.com/lh3/minimap2) binary and place these
within the `PATH`.

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
 > * make
 > * wget
 > * python3-all-dev
 > * python-virtualenv

A Makefile is provided to fetch, compile and install all direct dependencies
into a python virtual environment. To setup the environment run:

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

Using this method both `samtools` and `minimap2` are built from source and need
not be provided by the user.

Usage
-----

`medaka` can be run using its default settings through the `medaka_consensus`
program. An assembly in `.fasta` format and basecalls in `.fasta` or `.fastq`
format are required. The program uses both `samtools` and `minimap2`. If
medaka has been installed using the from-source method these will be present
within the medaka environment, else they will need to be provided by the user.

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC}

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. When `medaka_consensus` has finished running, the consensus
will be saved to `${OUTDIR}/consensus.fasta`.

Acknowledgements
----------------

We thank [Joanna Pineda](https://github.com/jopineda) and
[Jared Simpson](https://github.com/jts) for providing htslib code samples which aided
greatly development of the optimised feature generation code, and for testing the
version 0.4 release candidates.
