Medaka
======

© 2018 Oxford Nanopore Technologies Ltd.

`medaka` is a tool to create a consensus sequence of nanopore sequencing data.
This task is performed using neural networks applied a pileup of individual
sequencing reads against a draft assembly. It outperforms graph-based methods
operating on basecalled data, and can be competitive with state-of-the-art
signal-based methods whilst being much faster.

As input `medaka` accepts reads in either a `.fasta` or a `.fastq` file. It
requires a draft assembly as a `.fasta`.

`medaka` is distributed under the terms of the Mozilla Public License 2.0.


Features
--------

  * Requires only basecalled data (`.fasta` or `.fastq`).
  * Improved accurary over graph-based methods (e.g. Racon).
  * 50X faster than Nanopolish (and can run on GPUs)..
  * Benchmarks are provided (see :doc:`benchmarks`).
  * Includes extras for implementing and training bespoke correction
    networks.
  * Works on Linux (MacOS and Windows support is untested).
  * Open source (Mozilla Public License 2.0).

Tools to enable the creation of draft assembies can be found in a sister
project `pomoxis <https://github.com/nanoporetech/pomoxis>`_.

.. admonition:: Research Release

    Research releases are provided as technology demonstrators to provide early
    access to features or stimulate Community development of tools. Support for
    this software will be minimal and is only provided directly by the developers.
    Feature requests, improvements, and discussions are welcome and can be
    implemented by forking and pull requests. However much as we would
    like to rectify every issue and piece of feedback users may have, the 
    developers may have limited resource for support of this software. Research
    releases may be unstable and subject to rapid iteration by Oxford Nanopore
    Technologies.


Acknowledgements
................

We thank `Joanna Pineda <https://github.com/jopineda>`_ and
`Jared Simpson <https://github.com/jts>`_ for providing htslib code samples
which aided greatly development of the optimised feature generation code,
and for testing the version 0.4 release candidates.

We thank `Devin Drown <https://github.com/devindrown>`_ for
`working through <https://github.com/nanoporetech/medaka/issues/70>`_
use of ``medaka`` with his RTX 2080 GPU.


Table of contents
-----------------

.. toctree::
   :maxdepth: 2

   installation
   benchmarks
   snp
   walkthrough
   draft_origin
   development
   future
   

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

