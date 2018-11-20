
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
  * [Benchmarks](https://nanoporetech.github.io/medaka/benchmarks.html)
  * Includes extras for implementing and training bespoke correction
    networks.
  * Works on Linux (MacOS and Windows support is untested).
  * Open source (Mozilla Public License 2.0).

Documentation can be found at https://nanoporetech.github.io/medaka/.


Installation
------------
  
Medaka should be installed inside a virtual environment. A Makefile is
provided to fetch, compile and install all direct dependencies into an
environment.

To setup the environment run:

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate


Usage
-----

`medaka` can be run using its default settings through the `medaka_consensus`
program. An assembly in `.fasta` format and basecalls in `.fasta` or `.fastq`
format are required. 

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC}

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. When `medaka_consensus` has finished running, the consensus
will be saved to `${OUTDIR}/consensus.fasta`.
