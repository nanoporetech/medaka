Welcome to Medaka's documentation!
==================================

Medaka demonstrates a framework for rapid error correction of sequencing data,
particularly aimed at nanopore sequencing. Tools are provided for both training
and inference. The code exploits the `keras <https://keras.io>`_ deep learning
library.

The framework provided has had limited testing and benchmarking, but it is
being released in an effort to allow researchers to experiment without having
to write much of the tedious data preparation code. Researchers can extend the
tools provided to achieve better results that those obtained currently.

See :doc:`examples` for standard usage in correcting a draft assembly,
:doc:`walkthrough` for a complete overview of training a model and its use,
and :doc:`future` for inspirational ideas and some background into nanopore
sequence correction.

Development of medaka as a "base-space" consensus tool has been paused to focus
on methods which exploit the nanopore signal data. Nevertheless researchers may
find medaka useful as a method to generate quickly sufficiently accurate
consensus sequences in many use cases; see the :doc:`benchmarks` page for further
details.

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
   walkthrough
   benchmarks
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

