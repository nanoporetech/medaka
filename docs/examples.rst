Medaka Examples
===============

Medaka demonstrates a framework for error correcting sequencing data,
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
Future versions will implement correction schemes working directly from signal
data of multiple reads. For more details see :ref:`FutureDirections`.

.. _SequenceCorrection:

Sequence correction
-------------------



Neural network training
-----------------------

