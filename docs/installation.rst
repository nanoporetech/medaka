
.. _installation:

Getting Started
===============

Medaka can be installed in one of several ways.

**Installation with conda**

Perhaps the simplest way to start using medaka on both Linux and MacOS is
through conda; medaka is available via the
`bioconda <https://anaconda.org/bioconda/medaka>`_ channel:

.. code-block:: bash

    conda install -c bioconda medaka


**Installation with pip**
  
Medaka can be installed using the python package manager, pip:

.. code-block:: bash

    pip install medaka

On Linux platforms this will install a precompiled binary, on MacOS (and other)
platforms this will fetch and compile a source distribution.

We recommend using medaka within a virtual environment, viz.:

.. code-block:: bash

    virtualenv medaka --python=python3 --prompt "(medaka) "
    . medaka/bin/activate
    pip install medaka

.. note::

    Using this method requires the user to provide several binaries:

    `samtools <https://github.com/samtools/samtools>`_,
    `minimap2 <https://github.com/lh3/minimap2>`_,
    `tabix <https://github.com/samtools/htslib>`_, and
    `bgzip <https://github.com/samtools/htslib>`_,

    and place these within the `PATH`. `samtools` version 1.9 and `minimap2`
    version 2.11 are recommended as these are those used in development of
    medaka.


**Installation from source**

Medaka can be installed from its source quite easily on most systems.

.. note::

    Before installing medaka it may be required to install some
    prerequisite libraries, best installed by a package manager. On Ubuntu
    theses are:
    
    gcc zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev
    libcurl4-gnutls-dev libssl-dev curl make wget python3-all-dev python-virtualenv

A Makefile is provided to fetch, compile and install all direct dependencies
into a python virtual environment. To setup the environment run:

.. code-block:: bash

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

Using this method both ``samtools`` and ``minimap2`` are built from source and need
not be provided by the user.


**Using a GPU**

All installation methods will allow medaka to be used with CPU resource only.
To enable the use of GPU resource it is necessary to install the
``tensorflow-gpu`` package. Unfortunately depending on your python version it
may be necessary to modify the requirements of the ``medaka`` package for it
to run without complaining. Using the source code from github a working
GPU-powered ``medaka`` can be configured with:

.. code-block:: bash

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    sed -i 's/tensorflow/tensorflow-gpu/' requirements.txt
    make install

However, note that The ``tensorflow-gpu`` GPU package is compiled against a
specific version of the NVIDIA CUDA library; users are directed to the 
`tensorflow installation <https://www.tensorflow.org/install/gpu>`_ pages
for further information.

.. _sequence_correction:

Sequence correction
-------------------
 
After installing the software (see :ref:`installation`), `medaka` can be run
using its default settings through the `medaka_consensus` program. An
assembly in `.fasta` format and basecalls in `.fasta` or `.fastq` format are
required (see :ref:`basecalling_and_draft_assembly` for an detailed example
of one method of obtaining these). The program uses both `samtools` and `minimap2`.
If medaka has been installed using the from-source method these will be present
within the medaka environment, else they will need to be provided by the user.

.. code-block:: bash

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m r94

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. When `medaka_consensus` has finished running, the consensus
will be saved to `${OUTDIR}/consensus.fasta`.

.. warning::

    It is crucially important to specify the correct model, ``-m`` in the
    above, according to the basecaller used. Allowed values can be found by
    running ``medaka tools list\_models``.
    
    For guppy v3.0.3 models are named similarly to their basecalling counterparts
    with a "fast" and "high accuracy" model, for example ``r941_min_fast`` and
    ``r941_min_high``. The medaka models are equal in speed regardless of basecaller
    speed/accuracy.
    
    For guppy versions >=2.1.3 where the flip-flop algorithm has been used, users
    should select the highest numbered model equal to or less than the guppy
    version used for basecalling. There are two models here: ``r941_flip213`` and
    ``r941_flip235``
    
    A final model ``r941_trans`` is available where a basecaller with the transducer
    algorithm has been used (Albacore or Guppy<2.1.3).


Improving parallelism
~~~~~~~~~~~~~~~~~~~~~

The ``medaka_consensus`` program is good for simple datasets but perhaps not
optimal for running large datasets at scale. examples. A higher level of
parallelism can be achieved by running independently the component steps
of ``medaka_consensus``. The program performs three tasks:

1. alignment or reads to input assembly (via ``mini_align`` which is a thin
   veil over ``minimap2``)
2. running of consensus algorithm across assembly regions
   (``medaka consensus``, note no underscore!)
3. aggregation of the results of 2. to create consensus sequences
   (``medaka stitch``)

The three steps are discrete, and can be split apart an run independently. In
most cases, Step 2. is the bottleneck and can be trivially parallelized. The
``medaka consensus program`` can be supplied a ``--regions``
argument which will restrict its action to particular assembly sequences from
the ``.bam`` file output in Step 1. Therefore individual jobs can be run for batches
of assembly sequences simultaneously. In the final step, ``medaka stitch``
can take as input one or more of the ``.hdf`` files output by Step 2.

So in summary something like this is possible:

.. code-block:: bash

    # align reads to assembly
    mini_align -i basecalls.fasta -r assembly.fasta -P -m \
        -p calls_to_draft.bam -t <threads>
    # run lots of jobs like this, change model as appropriate
    mkdir results
    medaka consensus calls_to_draft.bam results/contigs1-4.hdf \
        --model r941_flip235 --batch 200 --threads 8 \
        --region contig1 contig2 contig3 contig4
    ...
    # wait for jobs, then collate results
    medaka stitch results/*.hdf polished.assembly.fasta

It is not recommended to specify a value of ``--threads`` greater than 8 for
``medaka consensus`` since the compute scaling efficiency is poor beyond this.
Note also than ``medaka consensus`` may been seen to use resource equivalent to
``<threads> + 4`` as an additional 4 threads are used for reading and preparing
input data.

Improving parallelism
~~~~~~~~~~~~~~~~~~~~~

The ``medaka_consensus`` program is good for simple datasets but perhaps not
optimal for running large datasets at scale. examples. A higher level of
parallelism can be achieved by running independently the component steps
of ``medaka_consensus``. The program performs three tasks:

1. alignment or reads to input assembly (via ``mini_align`` which is a thin
   veil over ``minimap2``)
2. running of consensus algorithm across assembly regions
   (``medaka consensus``, note no underscore!)
3. aggregation of the results of 2. to create consensus sequences
   (``medaka stitch``)

The three steps are discrete, and can be split apart an run independently. In
most cases, Step 2. is the bottleneck and can be trivially parallelized. The
``medaka consensus program`` can be supplied a ``--regions``
argument which will restrict its action to particular assembly sequences from
the ``.bam`` file output in Step 1. Therefore individual jobs can be run for batches
of assembly sequences simultaneously. In the final step, ``medaka stitch``
can take as input one or more of the ``.hdf`` files output by Step 2.

So in summary something like this is possible:

.. code-block:: bash

    # align reads to assembly
    mini_align -i basecalls.fasta -r assembly.fasta -P -m \
        -p calls_to_draft.bam -t <threads>
    # run lots of jobs like this, change model as appropriate
    mkdir results
    medaka consensus calls_to_draft.bam results/contigs1-4.hdf \
        --model r941_flip235 --batch 200 --threads 8 \
        --region contig1 contig2 contig3 contig4
    ...
    # wait for jobs, then collate results
    medaka stitch results/*.hdf polished.assembly.fasta

It is not recommended to specify a value of ``--threads`` greater than 8 for
``medaka consensus`` since the compute scaling efficiency is poor beyond this.
Note also than ``medaka consensus`` may been seen to use resource equivalent to
``<threads> + 4`` as an additional 4 threads are used for reading and preparing
input data.

Origin of the draft sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Medaka has been trained to correct draft sequences processed through
`racon <https://github.com/isovic/racon>`_), specifically `racon` run four times
iteratively with:

    racon -m 8 -x -6 -g -8 -w 500 ...

Processing a draft sequence from alternative sources (e.g. the output of
`canu <https://github.com/marbl/canu>`_ or
`wtdbg2 <https://github.com/ruanjue/wtdbg2>`_) may lead to poorer results
even when the draft is of a superior quality than that obtained from `racon`.

The [walkthrough](https://nanoporetech.github.io/medaka/walkthrough.html#walkthrough)
outlines one recommended workflow rapid construction of a draft for input into
`medaka`. A second approach would be to run `canu` followed by `racon` applied
twice iteratively before entry into `medaka`.
