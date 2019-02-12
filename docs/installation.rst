
.. _installation:

Installation
============

There are currently two installation methods for medaka, detailed below.

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

    and place these within the `PATH`. `samtools` version 1.3.1 and `minimap2`
    version 2.11 are recommended as these are those used in development of
    medaka.


**Installation from source**

Medaka can be installed from its source quite easily on most systems.

.. note::

    Before installing medaka it may be required to install some
    prerequisite libraries, best installed by a package manager. On Ubuntu
    theses are:
    
    gcc zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev make wget
    python3-all-dev python-virtualenv

A Makefile is provided to fetch, compile and install all direct dependencies
into a python virtual environment. To setup the environment run:

.. code-block:: bash

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

Using this method both `samtools` and `minimap2` are built from source and need
not be provided by the user.


**Using a GPU**

All installation methods will allow medaka to be used with CPU resource only.
To enable the use of GPU resource it is necessary to install the
`tensorflow-gpu` package. In outline this can be achieved with:

.. code-block:: bash

    pip uninstall tensorflow
    pip install tensorflow-gpu

However, note that The `tensorflow-gpu` GPU package is compiled against a
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

   It is crucially important to specify the correct model, :code:`-m` in the
   above, according to the basecaller used. Allowed values can be found by
   running :code:`medaka consensus --help`. For example to run medaka with a
   model suitable for the flip-flop basecaller in Guppy use :code:`-m r941_flip`.


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
