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

`medaka consensus` uses a neural network error model to fix errors in input sequences.

.. code-block:: bash

	usage: medaka consensus [-h] [--encoding ENCODING]
							[--output_fasta OUTPUT_FASTA] [--start START]
							[--end END]
							(--pileupdata PILEUPDATA | --feature_file FEATURE_FILE | --alignments reads.bam ref.fasta ref_name)
							model

	positional arguments:
	  model                 Model .hdf file from training.

	optional arguments:
	  -h, --help            show this help message and exit
	  --encoding ENCODING   Model label encoding .json, used only if encoding not
							in .hdf model
	  --output_fasta OUTPUT_FASTA
							Polished consensus output file.
	  --start START         Reference position at which to start, only used with
							--alignments.
	  --end END             Reference position at which to end, only used with
							--alignments.
	  --pileupdata PILEUPDATA
							Pileup input data.
	  --feature_file FEATURE_FILE
							Pregenerated features as stored during training.
	  --alignments reads.bam ref.fasta ref_name
							Input alignments, reference fasta and reference name
							(within fasta).


Neural network training
-----------------------

Training an error model is a two-stage process.

Generate a training data HDF5 file using `medaka prepare`.

.. code-block:: bash

	usage: medaka prepare [-h] [--start START] [--end END] [--chunk_len CHUNK_LEN]
						  [--truth TRUTH]
						  ref_fasta ref_name output bam

	positional arguments:
	  ref_fasta             .fasta reference corresponding to bams.
	  ref_name              Name of reference within ref_fasta.
	  output                Output .hdf.
	  bam                   Input alignments (aligned to ref).

	optional arguments:
	  -h, --help            show this help message and exit
	  --start START         Start reference coordinate.
	  --end END             End reference coordinate.
	  --chunk_len CHUNK_LEN
							Width of pileup chunks (in ref coords) to produce.
	  --truth TRUTH         Input bam of truth aligned to ref to label data.

Then supply the training data file to `medaka train`.

.. code-block:: bash
    
    Train model using preprocessed training data.

	usage: medaka train [-h] [--train_name TRAIN_NAME] [--model MODEL]
						[--features]
						pileupdata

	positional arguments:
	  pileupdata            Path for training data.

	optional arguments:
	  -h, --help            show this help message and exit
	  --train_name TRAIN_NAME
							Name for training run.
	  --model MODEL         Model definition and initial weights .hdf.
	  --features            Stop after generating features.
