Medaka
======

`medaka` is a tool to create a consensus sequence of nanopore sequencing data.
This task is performed using neural networks applied a pileup of individual
sequencing reads against a draft assembly. It outperforms graph-based methods
operating on basecalled data, and can be competitive with state-of-the-art
signal-based methods whilst being much faster.

As input `medaka` accepts reads in either a `.fastq` or a `.fastq` file. It
requires a draft assembly as a `.fasta`.

Features
--------

  * Requires only basecalled data.
  * Improved accurary over graph-based methods (e.g. Racon).
  * 50X faster than Nanopolish.
  * Includes extras for implementing and training bespoke correction
    networks.
  * Works on Linux and MacOS (Windows support is untested).
  * Open source (Mozilla Public License 2.0).


Tools to enable the creation of draft assembies can be found in a sister
project `pomoxis <https://github.com/nanoporetech/pomoxis>`_.



Table of contents
-----------------

.. toctree::
   :maxdepth: 2

   installation
   benchmarks
   walkthrough
   development
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

