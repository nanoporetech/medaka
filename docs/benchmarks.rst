Benchmarks
==========

The following demonstrates the utility of a recurrent neural network trained to
perform a consensus call from a pileup in which basecalls and the draft
assembly have been compressed using run-length
encoding (as  implemented with :ref:`polishing_with_compressed_hp`).  The network
receives counts of run-length encodings within each column of a pileup obtained
by aligning the compressed basecalls to the compressed draft assembly. 

Results were obtained using the model released with `Medaka`.  This model was
trained on e.coli, yeast and human data. Training data was basecalled with Guppy
0.3.0.  Draft assemblies here were created using the
`mini_assemble <https://nanoporetech.github.io/pomoxis/examples.html#fast-de-novo-assembly>`_
pipeline in `pomoxis <https://github.com/nanoporetech/pomoxis>`_. 

Error statistics were calculated using the 
`Pomoxis <https://github.com/nanoporetech/pomoxis>`_ script `stats_from_bam` after
aligning 10kb chunks of the consensus to the reference. 

Comparison of `Medaka` and `Nanopolish` 
-------------------------------------------------

Evaluation of the model was performed using the `Medaka` e.coli
:doc:`walkthrough` dataset (not part of the model-training data), basecalled with 
`Scrappie <https://github.com/nanoporetech/scrappie>`_ using the rgrgr_r94 model. The
pileup had a median depth of ~80. 
`Nanopolish <https://github.com/jts/nanopolish>`_ was run with homopolymer corrections but
without methylation corrections. `Medaka` and `Nanopolish` were run on the same hardware. 
`Medaka` was also run on a Macbook Pro and took 45 minutes running on two cores.

+-----------------+--------+------------+
|                 | medaka | nanopolish |
+=================+========+============+
| Q(Accuracy)     |  30.53 |  30.80     |
+-----------------+--------+------------+
| Q(Identity)     |  45.35 |  42.27     |
+-----------------+--------+------------+
| Q(Deletion)     |  31.82 |  31.64     |
+-----------------+--------+------------+
| Q(Insertion)    |  37.03 |  40.60     |
+-----------------+--------+------------+
| Runtime (hours) |  0.17  |  3.0 2     |
+-----------------+--------+------------+
| Cores           |  4     |  32        |
+-----------------+--------+------------+
| Core hours      |  0.67  |  96.53     |
+-----------------+--------+------------+

`Medaka` delivers similar results to `Nanopolish` in a fraction of the time. 


Evaluation Across Samples and Depths
------------------------------------

Evaluation of the model was performed using e.coli (using a set of reads not
included in model training), human chromosome 21 (not part of the model
training data) and `klebsiella
<https://github.com/rrwick/Basecalling-comparison>`_. 
The e.coli and human reads were basecalled with `Guppy` version 0.3.0,
while klebsiella reads were basecalled with `Scrappie
<https://github.com/nanoporetech/scrappie>`_ using the rgrgr_r94 model.  The
draft assemblies here were created at multiple depths using the `mini_assemble
<https://nanoporetech.github.io/pomoxis/examples.html#fast-de-novo-assembly>`_
pipeline in `Pomoxis <https://github.com/nanoporetech/pomoxis>`_, which
assembles a draft using `Miniasm <https://github.com/lh3/miniasm>`_ and
performs four rounds of `Racon <https://github.com/isovic/racon>`_ polishing. 

+------------------+-----------------+------------------+
| Data set         | Racon Error (%) | Medaka Error (%) |
+==================+=================+==================+
| E.coli 25 X      |       0.47      |       0.19       |
+------------------+-----------------+------------------+
| E.coli 57 X      |       0.37      |       0.10       |
+------------------+-----------------+------------------+
| Klebsiella 25 X  |       0.86      |       0.63       |
+------------------+-----------------+------------------+
| Klebsiella 57 X  |       0.72      |       0.45       |
+------------------+-----------------+------------------+
| Human chr21 31 X |       1.01      |       0.48       |
+------------------+-----------------+------------------+


Medaka reduces the total error in the `Racon` consensus by roughly a factor of two. 
