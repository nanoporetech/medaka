Benchmarks
==========

The following demonstrates the utility of a recurrent neural network trained to
perform a consensus call from a pileup in which basecalls and the draft
assembly have been reduced using run-length encoding (as demonstrated in
:ref:`sequence_correction`). The network receives counts of base
run-lengths within each column of a pileup obtained by aligning the encoded
basecalls to the encoded draft assembly. 

Results were obtained using the default model provided with `medaka`. This model
was trained using data obtained from E.coli, S.cerevisaie and H.sapiens samples.
Training data were basecalled using Guppy 0.3.0. Draft assemblies were created
using the `mini_assemble <https://nanoporetech.github.io/pomoxis/examples.html#fast-de-novo-assembly>`_
pipeline in `pomoxis <https://github.com/nanoporetech/pomoxis>`_. 

Error statistics were calculated using the 
`pomoxis <https://github.com/nanoporetech/pomoxis>`_ program `stats_from_bam` after
aligning 10kb chunks of the consensus to the reference. 


Comparison of `medaka` and `nanopolish` 
---------------------------------------

Evaluation of the model was performed using the `medaka` E.coli
:doc:`walkthrough` dataset. These data we not used to train the model.
Basecalling was performed with 
`scrappie <https://github.com/nanoporetech/scrappie>`_ using the `rgrgr_r94`
model. The pileup had a median depth of ~80-fold.
`nanopolish <https://github.com/jts/nanopolish>`_ was run with homopolymer
correction but without methylation correction. `medaka` and `nanopolish` were
run on the same hardware. 

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
| runtime (hours) |  0.17  |  3.0       |
+-----------------+--------+------------+
| cpu threads [1]_|  4     |  32        |
+-----------------+--------+------------+
| cpu hours       |  0.67  |  96.53     |
+-----------------+--------+------------+

.. [1] A CPU using SMT was used.

For this dataset `medaka` delivers similar results to `nanopolish` in a
fraction of the time. 


Evaluation across samples and depths
------------------------------------

Evaluation of the model was performed using E.coli, H.sapiens chromosome 21,
and `K.pneumoniae <https://github.com/rrwick/Basecalling-comparison>`_. 
The E.coli and human reads were basecalled with `Guppy` version 0.3.0,
while the Klebsiella reads were basecalled with `scrappie
<https://github.com/nanoporetech/scrappie>`_ using the `rgrgr_r94` model. The
draft assemblies here were created at multiple depths using the `mini_assemble
<https://nanoporetech.github.io/pomoxis/examples.html#fast-de-novo-assembly>`_
pipeline in `pomoxis <https://github.com/nanoporetech/pomoxis>`_.

+---------------------+-----------------+------------------+
| Data set            | Racon Error (%) | Medaka Error (%) |
+=====================+=================+==================+
| E.coli 25X          |       0.47      |       0.19       |
+---------------------+-----------------+------------------+
| E.coli 57X          |       0.37      |       0.10       |
+---------------------+-----------------+------------------+
| K.pneumoniae 25X    |       0.86      |       0.63       |
+---------------------+-----------------+------------------+
| K.pneumoniae 57X    |       0.72      |       0.45       |
+---------------------+-----------------+------------------+
| H.sapiens chr21 31X |       1.01      |       0.48       |
+---------------------+-----------------+------------------+

`medaka` reduces the total error in the `racon` consensus by roughly a factor of two. 
