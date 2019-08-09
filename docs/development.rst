Development
===========

As well as providing the useful functionality of creating consensus calls from
basecall data, ``medaka`` demonstrates a framework for both training and
inference. The code exploits the `keras <https://keras.io>`_ deep learning
library.

A simple method presented to inspire further ideas is to align input sequences
and a truth sequence to a common baseline (e.g. reads and a reference to a
draft assembly), and extract 'feature vectors' for input into a neural network
classifier.


API reference
-------------

The medaka library comprising feature generation, data labelling, and batching
code is detailed below.

.. toctree::
   :maxdepth: 3
      
   medaka
