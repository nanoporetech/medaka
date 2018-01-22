Getting Started
===============

`medaka` demonstrates a framework for error correcting sequencing data,
particularly aimed at nanopore sequencing. Tools are provided for both training
and inference. The code exploits the `keras <https://keras.io>`_ deep learning
library.

Tools to enable the creation of draft assembies can be found in a sister
project `pomoxis <https://github.com/nanoporetech/pomoxis>`_.

A simple method presented to inspire further ideas is to align input sequences
and a truth sequence to a common baseline (e.g. reads and a reference to a
draft assembly), and extract 'feature vectors' for input into a neural network
classifier.

The currently implemented neural networks are limited in their ability to
correct errors: to predict an output base from a window of oriented base
features. They make no use of primary signal data from nanopore experiments,
as for instance performed by `nanopolish <https://github.com/jts/nanopolish>`_.
Future projects will implement correction schemes working directly from signal
data of multiple reads. For more details see :doc:`future`.

For a complete examples starting from signal data, calculating basecalls,
through forming a draft assembly, to training and using a consensus network
see the :doc:`walkthrough`.


.. _creating_software_env:

Creating a Software Environment
-------------------------------

While not strictly necessary to use `medaka`, the :doc:`walkthrough` makes use of
several other open-source tools from Oxford Nanopore Technologies. To obtain
these run the following:

.. code-block:: bash

    MEDAKAHOME=${PWD}
    git clone https://github.com/nanoporetech/pomoxis --recursive
    git clone https://github.com/nanoporetech/medaka
    git clone https://github.com/nanoporetech/scrappie
    
    # While it is possible to install pomoxis and medaka into the same
    #   virtual environment, we will install each package into its own
    #   environment for simplicity. For more details see the readme for
    #   each of the packages.

    cd pomoxis && make install && cd ..
    cd medaka && make install && cd ..
    cd scrappie && mkdir build && cd build && cmake .. && make && cd ../../

    POMOXIS=${MEDAKAHOME}/pomoxis/venv/bin/activate
    MEDAKA=${MEDAKAHOME}/medaka/venv/bin/activate

If any of the above steps fail consult the install instructions of the
individual packages. 
 
In order to train effectively consensus models it is desirable to use a system
with a GPU. In order to use `medaka` training with a GPU, one should further
install the `tensorflow-gpu` package from pip into the medaka virtual
environment with:

.. code-block:: bash

    source ${MEDAKA}
    pip install tensorflow-gpu

Depending on your environment, specifically the versions of `CUDA` and `cuDNN`
that one has installed, it may be necessary to use versions of this package other
than the latest. `medaka` has been used with all release versions of `tensorflow`
after and including version 1.0.0.


.. _sequence_correction:

Sequence correction
-------------------
 
After installing all the necessary software (see :ref:`creating_software_env`),
`medaka` can be run using its default settings through the `medaka_consensus`
program. An assembly in `.fasta` format and basecalls in `.fasta` or `.fastq`
format are required, see :ref:`basecalling_and_draft_assembly`.

.. code-block:: bash

    source ${MEDAKA}  # i.e. medaka/bin/activate as above
    NPROC=$(nproc)
    medaka_consensus -i basecalls.fa -d draft_assembly.fa -o medaka_consensus -t ${NPROC} -p ${POMOXIS}

When `medaka_consensus` has finished running, the consensus will be saved to
`medaka_consensus/consensus.fasta`.
