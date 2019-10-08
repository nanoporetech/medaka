.. _Benchmarks:

Benchmarks
==========

The following demonstrates the utility of ``medaka``'s neural network in forming an
improved consensus from a pileup of reads.

Results were obtained using the default models provided with ``medaka``. These models
were trained using data obtained from a set of prokaroyotic samples.

Error statistics were calculated using the `pomoxis
<https://github.com/nanoporetech/pomoxis>`_ program ``assess_assembly`` after
aligning 100kb chunks of the consensus to the reference. Reported metrics are
median values over all chunks. 


Comparison of `medaka` and `nanopolish` 
---------------------------------------

In this comparison *E. coli* data from the :doc:`walkthrough` was used.
These data were not used to train the model. Basecalling was performed using
``Guppy v3.2.3``. Basecalled reads were trimmed using `porechop
<https://github.com/rrwick/Porechop>`_ to remove adapters, and assembly was
performed using `canu v1.8 <https://github.com/marbl/canu>`_. The assembly was
optionally corrected using `racon v1.2.1 <https://github.com/isovic/racon>`_ before being passed
to ``medaka`` or ``nanopolish``. A gpu-enabled version (commit ``896b8066``) of
`nanopolish <https://github.com/jts/nanopolish>`_ was run from
`PR661 <https://github.com/jts/nanopolish/pull/661>`_.

The workflow used here optionally includes a **single** iteration of ``racon``, see
:ref:`draftorigin` for further details.

.. table::
    Consensus calling of a 90-fold coverage *E. coli* data set using various methods. ``Medaka`` is seen
    to work effectively, with or without prior application of racon.

    +--------------------+--------+-------------------+--------------+----------+--------------+
    |                    | *canu* | *medaka no racon* | *racon (x1)* | *medaka* | *nanopolish* |
    +--------------------+--------+-------------------+--------------+----------+--------------+
    | Q(accuracy)        |   25.4 |              36.9 |         28.8 |     37.2 |         31.4 |
    +--------------------+--------+-------------------+--------------+----------+--------------+
    | Q(substitution)    |   51.0 |              53.0 |         47.9 |     53.0 |         47.0 |
    +--------------------+--------+-------------------+--------------+----------+--------------+
    | Q(deletion)        |   25.4 |              37.3 |         29.2 |     37.7 |         31.7 |
    +--------------------+--------+-------------------+--------------+----------+--------------+
    | Q(insertion)       |   57.2 |              48.2 |         41.6 |     47.0 |         47.0 |
    +--------------------+--------+-------------------+--------------+----------+--------------+
    | CPU time / hr:min  |        |              0:05 |         0:26 |     0:05 |         1:54 |
    +--------------------+--------+-------------------+--------------+----------+--------------+

For ``racon`` and ``medaka`` the reported CPU times include read to draft
overlapping performed by ``minimap2``. ``medaka`` was run using an 
NVIDIA 1080Ti GPU, the total execution time for ``medaka`` itself was
35 seconds. The wallclock time for ``nanopolish`` execution was 17 minutes.

A particular advantage of ``medaka`` over other methods is its improved
accuracy in recovering homopolymer lengths.

.. image:: images/hp_acc.png
    :align: center

Above the main plot we show homopolymer frequencies from H.sapiens Chrom. 1,
adapted from `Statistical analysis of simple repeats in the human genome <http://dirac.cnrs-orleans.fr/~piazza/PB/files/DNA.pdf>`_.

Evaluation across samples and depths
------------------------------------

The comparison below illustrates results at various coverage depths for a
collection of further organisms. Assemblies were performed as above with
canu and racon, using the ``Guppy v3.0.3`` high accuracy basecaller and
``medaka v0.6.5``.

.. image:: images/cov_acc.png
    :align: center

