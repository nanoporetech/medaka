Welcome to Medaka's documentation!
==================================

Medaka demonstrates a framework for error correcting sequencing data,
particularly aimed at nanopore sequencing. Tools are provided for both training
and inference. The code exploits the `keras <https://keras.io>`_ deep learning
library.

The framework provided has had limited testing and benchmarking, but it is
being released in an effort to allow researchers to experiment without having
to write much of the tedious data preparation code. Researchers can extend the
tools provided to achieve better results that those obtained currently.

See :doc:`examples` for further herlp in getting started and
:ref:`FutureDirections` for inspirational ideas and some background into
nanopore sequence correction.

Tools to enable the creation of draft assembies can be found in a sister
project `pomoxis <https://github.com/nanoporetech/pomoxis>`_.


Installation
------------

Medaka should be installed inside a virtual environment. A Makefile is
provided to fetch, compile and install all direct dependencies into an
environment.

To setup the environment run:

.. code-block:: bash

    git clone https://github.com/nanoporetech/medaka.git
    cd medaka
    make install
    . ./venv/bin/activate

See :doc:`examples` for common workflows.

Contents
--------

.. toctree::
   :maxdepth: 2

   examples
   future

Full API reference
------------------

.. toctree::
   :maxdepth: 3
      
   medaka


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

