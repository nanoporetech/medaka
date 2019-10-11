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

In this comparison R9.4.1 *E. coli* data from the :doc:`walkthrough` were used.
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
collection of further organisms using both the R9.4.1 and R10 pores. Assemblies were
performed as above with canu and racon, using the ``Guppy v3.0.3`` high accuracy
basecaller and ``medaka v0.6.5``.

.. image:: images/cov_acc.png
    :align: center

Complex organism genome assembly
--------------------------------

The 100Mbase genome of the nematode worm *C. elegans* was assembled from a
`100-fold coverage R10 dataset <https://ont-research.s3-eu-west-1.amazonaws.com/r10_celegans.fq.gz>`_,
with median read length of 34kb (all reads >20kb) using the
`canu <https://canu.readthedocs.io/en/latest/index.html>`_,
`flye <https://github.com/fenderglass/Flye>`_ and
`shasta <https://github.com/chanzuckerberg/shasta>`_ assemblers.
In this experiment, training sets for flip-flop basecallers and medaka models
included a subset of the *C. elegans* genome. For brevity ``medaka`` was
applied to drafts direct from the assemblers *with no intermediate racon consensus polishing*.


Assembly metrics
****************

The table below presents basic assembly metrics for the three assemblers. Run times are
for illustrative purposes only and are not direct comparisons. The assemblers were run
on different hardware appropriate to their predicted run time: *canu*: a 6 node cluster of
commodity servers, *flye*: a single server using 20 CPU threads, *shasta*: a single
server using 46 CPU threads. Assembly metrics were calculated via
`quast <https://www.ncbi.nlm.nih.gov/pubmed/23422339>`_,
removing contigs less than 5kb in length `(Kolmogorov et. al., 2019) <https://www.nature.com/articles/s41587-019-0072-8>`_.
Accuracy metrics were calculated as described above for the bacterial datasets.

+-----------------+--------+--------+----------+
|                 | *canu* | *flye* | *shasta* |
+-----------------+--------+--------+----------+
| run time        | ~24 hr | ~13 hr | ~34 min  |
+-----------------+--------+--------+----------+
| #contigs        | 37     | 35     | 66       |
+-----------------+--------+--------+----------+
| length / Mb     | 105    | 106    | 102      |
+-----------------+--------+--------+----------+
| NG50 / Mb       | 5.8    | 7      | 3.6      |
+-----------------+--------+--------+----------+
| NGA50 / Mb      | 2.4    | 1.8    | 1.8      |
+-----------------+--------+--------+----------+
| misassemblies   | 221    | 163    | 93       |
+-----------------+--------+--------+----------+
| genome fraction | 99.7%  | 94.5%  | 99.2%    |
+-----------------+--------+--------+----------+
+-----------------+--------+--------+----------+
| Q(accuracy)     |   34.6 |   29.2 |     32.4 |
+-----------------+--------+--------+----------+
| Q(substitution) |   47.5 |   41.2 |     44.0 |
+-----------------+--------+--------+----------+
| Q(deletion)     |   39.2 |   36.4 |     37.2 |
+-----------------+--------+--------+----------+
| Q(insertion)    |   41.1 |   33.2 |     38.5 |
+-----------------+--------+--------+----------+

The results show ``canu`` to give the overall highest accuracy assembly, though
``shasta`` gives similar quality with a fraction of the compute resource. However,
the ``shasta`` assembly is somewhat more fragmented than the other two assemblies
albeit with a lower misassembly rate.

The accuracy metrics can be calculated also outside of repetitive regions of the
genome by use of `RepeatMasker <http://www.repeatmasker.org/species/ce.html>`_:

+-----------------+--------+--------+----------+
|                 | *canu* | *flye* | *shasta* |
+-----------------+--------+--------+----------+
| Q(accuracy)     |   37.3 |   32.7 |     35.2 |
+-----------------+--------+--------+----------+
| Q(substitution) |   50.0 |   43.0 |     47.0 |
+-----------------+--------+--------+----------+
| Q(deletion)     |   40.9 |   38.0 |     39.2 |
+-----------------+--------+--------+----------+
| Q(insertion)    |   44.0 |   34.9 |     40.9 |
+-----------------+--------+--------+----------+

Finally `busco <https://busco.ezlab.org>`_ can be used
to measure the gene content of the three assemblies. Results are reported in the usual
"busco notation", C: complete, S: complete and single copy, D: complete and duplicated,
F: fragmented, M: missing. To aid comparison we also give scores for the reference sequence
used for the assessments above.

+-----------+-------------------------------------------------+
|           | BUSCO scores                                    |
+-----------+-------------------------------------------------+
| reference | ``C:98.6%[S:98.0%,D:0.6%],F:0.8%,M:0.6%,n:982`` |
+-----------+-------------------------------------------------+
| canu      | ``C:99.3%[S:98.6%,D:0.7%],F:0.6%,M:0.1%,n:982`` |
+-----------+-------------------------------------------------+
| flye      | ``C:98.5%[S:98.1%,D:0.4%],F:1.1%,M:0.4%,n:982`` |
+-----------+-------------------------------------------------+
| shasta    | ``C:99.0%[S:98.5%,D:0.5%],F:0.9%,M:0.1%,n:982`` |
+-----------+-------------------------------------------------+

Curiously it is found that the canu and shasta assemblies show a higher complete count than
the reference sequence (`WBcel235 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6/>`_).


Future Work
***********

We have illustrated ``medaka``'s role in efficiently creating quality consensus
sequences, and the role of long nanopore reads in creating highly contiguous assemblies
with modest compute requirements. The assemblies are in agreement with recent studies
from `Tyson et. al. <https://genome.cshlp.org/content/12/5/669.full.html>`_ and
`Yoshimura et. al. <https://genome.cshlp.org/content/29/6/1009.full>`_ in suggesting
that the *C. elegans* reference has missing sequence. The assemblies above are not the
final word on assembly with nanopore sequencing data; in contrast to latter preceding
reference we have expended relatively little effort in the assembly process, relying
entirely on off-the-shelf automated methods.

Aside from improving assembly methodology we are exploring improved chemistries, 
basecallers, and consensus algorithms. Our current algorithmic focus is on so-called
run-length encoded (RLE) methods for both basecalling and consensus. Early work
in this area has shown promise in further reducing homopolymer error in particular.
