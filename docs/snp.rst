Variant calling
===============

Medaka implements a small variant calling pipeline to call single nucleotide
polymorphisms (SNPs) together with insertions and deletions (Indels) from Nanopore
basecalls. The pipeline leverages a diploid-aware neural network, read haplotype
tagging via `WhatsHap <https://whatshap.readthedocs.io>`_, and haplotype consensus
using ``medaka``'s usual consensus network.

The following benchmarking is performed on chromosome 21 of the HG001 (NA12878) sample,
according to `GA4GH best practices <https://www.nature.com/articles/s41587-019-0054-x>`_
using `hap.py <https://github.com/Illumina/hap.py>`_ with the
`GIAB truth set <http://jimb.stanford.edu/giab-resources/>`_ stratified to exclude
repeat regions as defined by the 
`GA4GH process <https://github.com/jzook/genome-data-integration/tree/master/NISTv3.3.2/filtbeds/GRCh38>`_.
Variants were filtered by their ``QUAL`` values (separately for SNP and Indel classes)
to maximise the F1 score.

SNP and Indel calling
---------------------

*Last updated December 2019*

``Medaka``'s variant calling pipeline first aligns all reads to a reference sequence,
creates a read pileup and uses a recurrent neural network to predict a pair of bases for
every reference locus. The predictions are combined with the reference
sequence to create candidate variants under an independence assumption between
loci; no attempt is made at this point to combine predicted bases across reference
loci. This approach alone gives a reasonable estimation of variants:

.. table::
    Single nucleotide polymorphism calling from the first stage of the
    ``medaka variant`` pipeline for chromosome 21 of the NA12878 sample.
    Comparison is made with the `GIAB <http://jimb.stanford.edu/giab-resources/>`_
    high confidence callset using `hap.py <https://github.com/Illumina/hap.py>`_.

   +--------------------+-----------+---------+----------+
   |                    | Precision | Recall  | F1 score |
   +--------------------+-----------+---------+----------+
   | medaka first stage |    0.9929 |  0.9392 |    0.965 |
   +--------------------+-----------+---------+----------+

Single nucleotide polymorphisms account for 85% of small variants in the human genome,
leaving a good proportion of insertion and deletion variants. It is therefore
desirable to be able to detect indel variants. To improve on the above results,
``medaka``'s pipeline phases the recovered variants
using `WhatsHap <https://whatshap.readthedocs.io>`_. This process allows
recovery of maternal and paternal haplotypes, and the assignment of each read to one of them. Having
partitioned reads into their haplotype, ``medaka consensus`` is run
independently for each haplotype to calculate haplotype specific consensus sequences.

.. image:: images/phased_snps.png

This second pass is a task in which we know ``medaka`` excels, see :ref:`Benchmarks`.
From the multiple haplotype consensus sequences it is a simple task to reconstruct
variant calls for both homo- and hetero-zygous sites.

.. table::
    Small variant calling at 100-fold coverage using `medaka variant` pipeline.
    Comparison is made with the `GIAB <http://jimb.stanford.edu/giab-resources/>`_
    high confidence callset using `hap.py <https://github.com/Illumina/hap.py>`_.

    +------------------+-------+-----------+---------+----------+
    |                  | Class | Precision | Recall  | F1 score |
    +------------------+-------+-----------+---------+----------+
    | medaka variant   | SNP   |    0.9967 |  0.9954 |    0.996 |
    +                  +-------+-----------+---------+----------+
    |                  | Indel |    0.9558 |  0.9174 |    0.936 |
    +------------------+-------+-----------+---------+----------+
    | clair            | SNP   |    0.9931 |  0.9950 |    0.994 |
    +                  +-------+-----------+---------+----------+
    |                  | Indel |    0.9094 |  0.8092 |    0.856 |
    +------------------+-------+-----------+---------+----------+
    | nanopolish       | SNP   |    0.9938 |  0.9662 |    0.980 |
    +------------------+-------+-----------+---------+----------+

Shown also are results from `clair <https://github.com/HKU-BAL/Clair>`_,
and older results for SNP calling from ``nanopolish``. We note that Clair
and medaka are competitive in terms of SNP calling but that medaka
excels Clair for indel calling. In contrast to medaka's haplotype-consensus
approach, clair attempts to call indels directly from the read pileup;
medaka's factorisation of variant calling problem would appear to be more
successful. This seems to be particularly true when comparing heterozygous
indels to homozygous cases:


.. table::
    Indel calling comparison of medaka and clair. Within the test
    region there are approximately equal numbers of homozygous and
    heterozygous cases.

    +------------------+-------+-----------+---------+----------+
    |                  | Indel | Precision | Recall  | F1 score |
    +------------------+-------+-----------+---------+----------+
    | medaka variant   | Hom   |    0.9790 |  0.9553 |    0.967 |
    +                  +-------+-----------+---------+----------+
    |                  | Het   |    0.8995 |  0.8738 |    0.886 |
    +------------------+-------+-----------+---------+----------+
    | clair            | Hom   |    0.9566 |  0.9175 |    0.937 |
    +                  +-------+-----------+---------+----------+
    |                  | Het   |    0.8707 |  0.7332 |    0.796 |
    +------------------+-------+-----------+---------+----------+


R10.3 SNP and Indel calling
-----------------------------

*Last updated February 2020*

The benchmarks were performed as described above, but on chromosome 20 of
the HG002 (NA24385) sample at 60-fold coverage.

.. table::
    Small variant calling at 60-fold coverage using `medaka variant` pipeline.
    Comparison is made with the `GIAB <http://jimb.stanford.edu/giab-resources/>`_
    high confidence callset using `hap.py <https://github.com/Illumina/hap.py>`_.

    +------------------+-------+-----------+---------+----------+
    |                  | Class | Precision | Recall  | F1 score |
    +------------------+-------+-----------+---------+----------+
    | R9.4.1           | SNP   |    0.9898 |  0.9931 |    0.992 |
    +                  +-------+-----------+---------+----------+
    |                  | Indel |    0.9199 |  0.8634 |    0.891 |
    +------------------+-------+-----------+---------+----------+
    | R10.3            | SNP   |    0.9892 |  0.9946 |    0.992 |
    +                  +-------+-----------+---------+----------+
    |                  | Indel |    0.9508 |  0.9385 |    0.945 |
    +------------------+-------+-----------+---------+----------+



Performing Variant Calling
--------------------------

.. note::

    The ``medaka_variant`` (and ``medaka_consensus``) pipeline only operate on
    a ``.bam`` alignment file containing a single sample (value of the `RG`
    alignment tag). It will refuse to run in the case of two read groups
    being present.

The pipeline described above is implemented in the ``medaka_variant`` program:

.. code-block:: bash

    source ${MEDAKA}
    medaka_variant -f <REFERENCE.fasta> -b <reads.bam>

This will run all steps of the process described above, finally outputting a
phased ``.vcf`` variant file.

.. note::

    Variants output by ``medaka_variant`` are filtered using default quality thresholds. 
    Users may wish to apply alternative filtering for their datasets and applications 
    using the -U, -P and -N options to ``medaka_variant``.


Further Improvements
--------------------

Medaka SNP and indel calling is an extremely active area of research.
We anticipate that the final step of variant calling, which amounts to a
consensus call on a relatively pure haplotyped set of reads, will benefit
(along with consensus calling) from further innovations in feature encoding,
network architecture and training strategy. Future releases may also exploit
direct snp and variant calling from mixed read populations.

