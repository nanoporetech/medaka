
.. _installation:

Installation
============


There are currently two installation methods for medaka, detailed below.

**Installation with pip**
  
Medaka can be installed on Linux using the python package manager, pip:

.. code-block:: bash

    pip install medaka

We recommend using medaka within a virtual environment, viz.:

.. code-block:: bash

    virtualenv pomoxis --python=python3 --prompt "(medaka) "
    . pomoxis/bin/activate
    pip install medaka

Using this method requires the user to provide a
`samtools <https://github.com/samtools/samtools>`_ and
`minimap2 <https://github.com/lh3/minimap2>`_ binary and place these
within the `PATH`.


**Installation from source**

Medaka can be installed from its source quite easily on most systems.

.. note::

    Before installing pomoxis it may be required to install some
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
    medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC}

The variables `BASECALLS`, `DRAFT`, and `OUTDIR` in the above should be set
appropriately. When `medaka_consensus` has finished running, the consensus
will be saved to `${OUTDIR}/consensus.fasta`.
