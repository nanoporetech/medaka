
.. _installation:

Installation
============

Medaka has been tested on Linux, specifically Ubuntu 14 and Ubuntu 16.

Medaka should be installed inside a `virtual environment
<https://docs.python.org/3/tutorial/venv.html>`_. A Makefile is provided to
fetch, compile and install all direct dependencies into an environment.

To setup the environment run:

.. code-block:: bash

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    source ./venv/bin/activate

A simple `DOCKERFILE` is provided also should users wish to deploy the code
in a container, but this is by no means a requirement.

.. _sequence_correction:


Sequence correction
-------------------
 
After installing the software (see :ref:`installation`), `medaka` can be run
using its default settings through the `medaka_consensus` program. An
assembly in `.fasta` format and basecalls in `.fasta` or `.fastq` format are
required (see :ref:`basecalling_and_draft_assembly` for an detailed example
of one method of obtaining these):

.. code-block:: bash

    source ${MEDAKA}  # i.e. medaka/venv/bin/activate
    NPROC=$(nproc)
    BASECALLS=basecalls.fa
    DRAFT=draft_assm/assm_final.fa
    OUTDIR=medaka_consensus
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC}

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. When `medaka_consensus` has finished running, the consensus
will be saved to `${OUTDIR}/consensus.fasta`.
