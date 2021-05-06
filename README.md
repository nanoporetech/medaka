
![Oxford Nanopore Technologies logo](https://github.com/nanoporetech/medaka/raw/master/images/ONT_logo_590x106.png)


Medaka
======

[![Build Status](https://travis-ci.org/nanoporetech/medaka.svg?branch=master)](https://travis-ci.org/nanoporetech/medaka)

[![](https://img.shields.io/pypi/v/medaka.svg)](https://pypi.org/project/medaka/)
[![](https://img.shields.io/pypi/wheel/medaka.svg)](https://pypi.org/project/medaka/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/medaka)
[![](https://img.shields.io/conda/pn/bioconda/medaka.svg)](https://anaconda.org/bioconda/medaka)


`medaka` is a tool to create consensus sequences and variant calls from
nanopore sequencing data.  This task is performed using neural networks applied
a pileup of individual sequencing reads against a draft assembly. It
outperforms graph-based methods operating on basecalled data, and can be
competitive with state-of-the-art signal-based methods whilst being much
faster.

© 2018- Oxford Nanopore Technologies Ltd.

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

For creating draft assemblies we recommend [Flye](https://github.com/fenderglass/Flye).

Documentation can be found at https://nanoporetech.github.io/medaka/.


Installation
------------

Medaka can be installed in one of several ways.

**Installation with conda**

Perhaps the simplest way to start using medaka is
through conda; medaka is available via the
[bioconda](https://anaconda.org/bioconda/medaka) channel:

    conda create -n medaka -c conda-forge -c bioconda medaka

Occasionally the conda releases lag behind the source code and
[PyPI](https://pypi.org/project/medaka/) releases.

**Installation with pip**

For those who prefer Python's native pacakage manager, medaka is also available
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

and place these within the `PATH`. `samtools/bgzip/tabix` version 1.11 and
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

Since version 1.1.0 `medaka` uses Tensorflow 2.2, prior versions used Tensorflow 1.4.
For `medaka` 1.1.0 and higher installation from source or using `pip` can make
immediate use of GPUs. However, note that the `tensorflow` package is compiled against
specific versions of the NVIDIA CUDA and cuDNN libraries; users are directed to the
[tensorflow installation](https://www.tensorflow.org/install/source#gpu) pages
for further information. cuDNN can be obtained from the
[cuDNN Archive](https://developer.nvidia.com/rdp/cudnn-archive), whilst CUDA
from the [CUDA Toolkit Archive](https://developer.nvidia.com/cuda-toolkit-archive).

For `medaka` prior to version 1.1.0, to enable the use of GPU resource it is
necessary to install the `tensorflow-gpu` package. Using the source code from github
a working GPU-powered `medaka` can be configured with:

    # Note: certain files are stored in git-lfs, https://git-lfs.github.com/,
    #       which must therefore be installed first.
    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    sed -i 's/tensorflow/tensorflow-gpu/' requirements.txt
    make install
    
*GPU Usage notes*

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

**Human variant calling**

Variant calling on human samples can be performed with the `medaka_variant` workflow:

    medaka_variant -f <REFERENCE.fasta> -i <reads.bam>

which requires simply a path to the human reference sequence, and a `.bam` file
of reads aligned to the reference.

**Bacterial (ploidy-1) variant calling**

Variant calling for monoploid samples is enabled through the `medaka_haploid_variant`
workflow:

    medaka_haploid_variant <reads.fastq> <ref.fasta>
    
which requires the reads as a `.fasta` or `.fastq` and a reference sequence as a
`.fasta` file.


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

### Origin of the draft sequence

Medaka has been trained to correct draft sequences processed through
[racon](https://github.com/isovic/racon) with:

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

© 2018- Oxford Nanopore Technologies Ltd.

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
