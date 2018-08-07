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
aligning 100kb chunks of the consensus to the reference. Reported metrics are median values over all chunks. 


Comparison of `medaka` and `nanopolish` 
---------------------------------------

Evaluation of the model was performed using the `medaka` E.coli
:doc:`walkthrough` dataset. These data we not used to train the model.
Basecalling was performed with 
`scrappie <https://github.com/nanoporetech/scrappie>`_ using the `rgrgr_r94`
model. The pileup had a median depth of ~80-fold.
`nanopolish v0.10.1 <https://github.com/jts/nanopolish>`_ was run with homopolymer
correction but without methylation correction. `medaka` and `nanopolish` were
run on the same hardware.  

+-----------------+--------+------------+
|                 | medaka | nanopolish |
+=================+========+============+
| Q(Accuracy)     |  31.55 |  31.22     |
+-----------------+--------+------------+
| Q(Identity)     |  46.02 |  45.23     |
+-----------------+--------+------------+
| Q(Deletion)     |  32.60 |  31.81     |
+-----------------+--------+------------+
| Q(Insertion)    |  39.21 |  43.01     |
+-----------------+--------+------------+
| runtime (hours) |   0.27 |  1.05      |
+-----------------+--------+------------+
| CPU cores       |   4    |  48        |
+-----------------+--------+------------+
| CPU hours       |   1.06 |  50.4      |
+-----------------+--------+------------+

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

+---------------------------+-----------------+------------------+----------------------+
| Data set                  | Racon Error (%) | Medaka Error (%) | Nanopolish Error (%) |
+===========================+=================+==================+======================+
| E.coli 25X                |       0.291     |       0.125      |       0.127          |
+---------------------------+-----------------+------------------+----------------------+
| E.coli 50X                |       0.247     |       0.079      |       0.080          |
+---------------------------+-----------------+------------------+----------------------+
| E.coli  75X               |       0.238     |       0.065      |       0.069          |
+---------------------------+-----------------+------------------+----------------------+
| E.coli 100X               |       0.228     |       0.063      |       0.060          |
+---------------------------+-----------------+------------------+----------------------+
| E.coli 150X               |       0.231     |       0.065      |       0.057          |
+---------------------------+-----------------+------------------+----------------------+
| E.coli 200X               |       0.231     |       0.061      |       0.052          |
+---------------------------+-----------------+------------------+----------------------+
| L.brevis 250X             |       0.293     |       0.055      |       0.047          |
+---------------------------+-----------------+------------------+----------------------+
| K.pneumoniae* 200X        |       0.576     |       0.318      |       0.086          |
+---------------------------+-----------------+------------------+----------------------+
* native (non-PCR'd) data. Nanopolish was run with the --methylation-aware=dcm,dam option.

`medaka` produces similar results to Nanopolish (on PCR'd data) in a fraction of the time. 
