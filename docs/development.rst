Development
===========

As well as providing the useful functionality of creating consensus calls from
basecall data, `medaka` demonstrates a framework for both training and
inference. The code exploits the `keras <https://keras.io>`_ deep learning
library.

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

For a complete example starting basecalls, through forming a draft assembly,
to training and using a consensus network see the :doc:`walkthrough`.

For more details on the internals of `medaka` see :doc:`medaka`.


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
    
    # While it is possible to install pomoxis and medaka into the same
    #   virtual environment, we will install each package into its own
    #   environment for simplicity. For more details see the readme for
    #   each of the packages.

    cd pomoxis && make install && cd ..
    cd medaka && make install && cd ..

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
    pip install tensorflow-gpu=={version}

where `{version}` should be replaced by the version found in the
`requirements.txt` file found within the software repository.

Depending on your environment, specifically the versions of `CUDA` and `cuDNN`
that one has installed, it may be necessary to compile from source tensorflow;
the precompiled versions are naturally linked against specific versions of
these libraries.


API reference
-------------

The medaka library comprising feature generation, data labelling, and batching
code is detailed below.

.. toctree::
   :maxdepth: 3
      
   medaka
