
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
  * Methylation aggregation from Guppy `.fast5` files.
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

    conda create -n medaka -c conda-forge -c bioconda medaka

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
`minimap2` version 2.17 are recommended as these are those used in development
of medaka.

**Installation from source**

Medaka can be installed from its source quite easily on most systems.

 Before installing medaka it may be required to install some
 prerequisite libraries, best installed by a package manager. On Ubuntu
 theses are:
 >     bzip2 g++ zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev
 >     libcurl4-gnutls-dev libssl-dev curl make cmake wget python3-all-dev
 >     python-virtualenv
 In addition it is required to install and set up git LFS before cloning
 the repository.

A Makefile is provided to fetch, compile and install all direct dependencies
into a python virtual environment. To set-up the environment run:

    # Note: certain files are stored in git-lfs, https://git-lfs.github.com/,
    #       which must therefore be installed first.
    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

Using this method both `samtools` and `minimap2` are built from source and need
not be provided by the user.


**Using a GPU**

All installation methods will allow medaka to be used with CPU resource only.
To enable the use of GPU resource it is necessary to install the
`tensorflow-gpu` package. Unfortunately depending on your python version it
may be necessary to modify the requirements of the `medaka` package for it
to run without complaining. Using the source code from github a working
GPU-powered `medaka` can be configured with:

    # Note: certain files are stored in git-lfs, https://git-lfs.github.com/,
    #       which must therefore be installed first.
    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    sed -i 's/tensorflow/tensorflow-gpu/' requirements.txt
    make install

However, note that The `tensorflow-gpu` GPU package is compiled against
specific versions of the NVIDIA CUDA and cuDNN libraries; users are directed to the
[tensorflow installation](https://www.tensorflow.org/install/gpu) pages
for further information. cuDNN can be obtained from the
[cuDNN Archive](https://developer.nvidia.com/rdp/cudnn-archive), whilst CUDA
from the [CUDA Toolkit Archive](https://developer.nvidia.com/cuda-toolkit-archive).

Depending on your GPU, `medaka` may show out of memory errors when running.
To avoid these the inference batch size can be reduced from the default
value by setting the `-b` option when running `medaka_consensus`. A value
`-b 100` is suitable for 11Gb GPUs.

For users with RTX series GPUs it may be required to additionally set an
environment variable to have `medaka` run without failure:

    export TF_FORCE_GPU_ALLOW_GROWTH=true

In this situation a further reduction in batch size may be required.


Usage
-----

`medaka` can be run using its default settings through the `medaka_consensus`
program. An assembly in `.fasta` format and basecalls in `.fasta` or `.fastq`
formats are required. The program uses both `samtools` and `minimap2`. If
medaka has been installed using the from-source method these will be present
within the medaka environment, otherwise they will need to be provided by
the user.

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m r941_min_high_g303

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. For the selection of the model (`-m r941_min_high_g303` in the
example above) see the Model section following.

When `medaka_consensus` has finished running, the consensus will be saved to
`${OUTDIR}/consensus.fasta`.

Models
------

For best results it is important to specify the correct model, `-m` in the
above, according to the basecaller used. Allowed values can be found by
running `medaka tools list\_models`.


Medaka models are named to indicate i) the pore type, ii) the sequencing
device (MinION or PromethION), iii) the basecaller variant, and iv) the
basecaller version, with the format:

    {pore}_{device}_{caller variant}_{caller version}

For example the model named `r941_min_fast_g303` should be used with data from
MinION (or GridION) R9.4.1 flowcells using the fast Guppy basecaller version
3.0.3. By contrast the model `r941_prom_hac_g303` should be used with PromethION
data and the high accuracy basecaller (termed "hac" in Guppy configuration
files). Where a version of Guppy has been used without an exactly corresponding
medaka model, the medaka model with the highest version equal to or less than
the guppy version should be selected.

Methylation Calling
-------------------

`medaka` includes a basic workflow for aggregating Guppy basecalling results
for Dcm, Dam, and CpG methylation. The workflow is currently very preliminary
and subject to change and improvement.

Aggregating the information from Guppy outputs is a two stage process, first
the basecalling results are extracted `.fast5` files and placed in a `.bam`
file:

    FAST5PATH=guppy/workspace
    REFERENCE=grch38.fasta
    OUTBAM=meth.bam
    medaka methylation guppy2sam ${FAST5PATH} ${REFERENCE} \
        --workers 74 --recursive \
        | samtools sort -@ 8 | samtools view -b -@ 8 > ${OUTBAM}
    samtools sort ${OUTBAM}

This program will extract both the basecall sequence and methylation scores,
align the basecall to the reference, and store results in a standard format.
In this preliminary workflow the methylation scores are stored in two SAM
tags, 'MC' and 'MA', one each for 5mC and 6mA respectively. The tags are
8bit integer array-values, one value per basecall position. This is a
different form to that proposed in the current
[hts-specs proposition](https://github.com/samtools/hts-specs/pull/418/files),
but allows for more trivial parsing.

The second step is to aggregate the reference-aligned information to produce
a simple tabular summary of read methylation counts:

    BAM=meth.bam
    REFERENCE=grch38.fasta
    REGION=chr20:500000-1000000
    OUTPUT=meth.tsv
    medaka methylation call --meth cpg ${BAM} ${REFERENCE} ${REGION} ${OUTPUT}

Here the option `--meth cpg` indicates that loci containing the sequence
motif `CG` should be examined for 5mC presence. Other choices are
`dcm` for which the motifs `CCAGG` and `CCTGG` are examined for 5mC and `dam`
(`GATC`) for 6mA.

The output file is a simple tab-delimited text file with columns:
'ref.name', 'position', 'motif', 'fwd.meth.count', 'rev.meth.count',
'fwd.canon.count', and 'rev.canon.count'. Here fwd./ref. indicate counts on the
two DNA strands and meth./canon. indicate counts for methylated and
canonical bases. Note that the position field records the position of the
first base in the motif recorded.

### Origin of the draft sequence

Medaka has been trained to correct draft sequences processed through
[racon](https://github.com/isovic/racon), specifically `racon` run four times
iteratively with:

    racon -m 8 -x -6 -g -8 -w 500 ...

Processing a draft sequence from alternative sources (e.g. the output of
[canu](https://github.com/marbl/canu) or
[wtdbg2](https://github.com/ruanjue/wtdbg2)) may lead to different results.

The [documentation](https://nanoporetech.github.io/medaka/draft_origin.html)
provides a discussion and some guidance on how to obtain a draft sequence.


Acknowledgements
----------------

We thank [Joanna Pineda](https://github.com/jopineda) and
[Jared Simpson](https://github.com/jts) for providing htslib code samples which aided
greatly development of the optimised feature generation code, and for testing the
version 0.4 release candidates.

We thank [Devin Drown](https://github.com/devindrown) for
[working through](https://github.com/nanoporetech/medaka/issues/70)
use of `medaka` with his RTX 2080 GPU.

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
