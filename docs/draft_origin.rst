
.. _draftorigin:

Origin of the draft sequence
============================

*August 8th 2019*

Much has been discussed in the
`popular press <https://twitter.com/ErichMSchwarz/status/1158808236557914112>`_
on the nature and origin of the expected draft sequence for ``medaka``. Here
we give some background of previous workflow recommendations by exploring how medaka
works and what this may imply for how it should be used.

.. note::

    Readers are encouraged to read all of the following to appreciate the concerns
    involved in choosing a particular workflow. For those looking for a single 
    recommendation it is advised to run ``racon`` at least once before running
    ``medaka``. Further application of ``racon`` can be fruitless.


Medaka's Algorithm
------------------

All ``medaka`` releases to date have heavily leveraged the idea of a read pileup
against a draft sequence to create input data for a neural network. ``medaka`` is
trained by creating such pileups and annotating each column with the truth base
(or gap) via alignment of reference sequences to the same draft sequences.

A rather general idea in the field of machine learning is the generalization
of models to data outside the set used for training, either to notionally
similar data or to characteristically new forms. The following discussion
highlights ways in which the input data to ``medaka`` can become similar or
markedly different form that to which the models were trained.

Aligners
........

The creation of a read pileup to a draft sequence requires the use of an
alignment procedure. The first concern for the performance of ``medaka`` is then
that results may (at least conceptually if not practically) be sensitive to the
aligner used to align reads. Furthermore there could be a
sensitivity to the exact alignment parameters used with any one aligner.

Within the convenience program ``medaka_consensus`` the aligner used is the
same as that during the training of ``medaka`` models, namely ``minimap2``. The
core consensus program ``medaka consensus`` will however take as input an
arbitrarily created ``.bam`` alignment file.

This sensitivity to the aligner results largely stems from how gaps are
calculated and recorded in the alignment: a different aligner with a different
parameterization can lead to a different gap structure in the pileup. This will
lead to a different input into the neural network, possibly one which the
network is not expecting. By way of example imagine a machine learning
algorithm trained to recognise different apple cultivars from photographs. When
presented with a photo of a pear, what is this algorithm to do? It certainly is
not going to give a meaningful answer to a human using the system. The
situation with ``medaka`` is likely not as severe as comparing apples and
pears, although notions of what constitutes a valid input remain.

Basecallers
...........

It is not simply the choice of aligner which can lead to characteristically
different neural network inputs. Naturally the choice of basecaller is a
consideration: different basecallers will produce calls with different
characteristic error profiles (with respect to a truth sequence). Again this
will effect the nature of the input into a neural network and potentially lead
to a lack of generalization across basecallers. It is for this reason that
``medaka`` conservatively provides models for a variety of basecallers
and basecaller models.

The draft sequence
..................

Perhaps more subtle is the effect of the draft input sequence. Variability in
draft sequences errors, for example that caused by the use of different
assembly methodologies or programs, can lead to different neural network
inputs. The crucial point is that if the relative rates of substitution,
insertion, and deletion *between the draft and reads* is different to those
presented to the algorithm during its training, performance may be degraded.

Sequencing depth
................

The structure of a read pileup is also affected by the sequencing depth. This
is not altogether apparent in alignment viewers such as IGV which tend to
`squeeze-out` insertions in their view. Consider the case where an insertion
occurs in one read with respect to the draft sequence `or` to other reads.
The current implementation in medaka creates an additional pileup column with
this single base. At higher sequencing depths more of these insertion columns
derived from random insertion events are present: the pileup becomes progressively
stretched with respect to the truth at higher sequencing depths. This effect is
explicitely mitigated in the techniques used to train ``medaka`` models.

Discussion
..........

It is for these and similar considerations that previous recommendations
regarding how to produce a draft sequence have been somewhat draconian. The
recommendations were based on the specific way in which the ``medaka`` models
are trained, with the anticipation that users would produce neural network
inputs that looked like apples, and not pears! These recommendations were by
their nature pessimistic, the inference error in presenting ``medaka`` with a
slightly bruised apple appears minimal.

Devin Drown usefully provided an example of this in the twitter conversation
where he tested the effect of running racon a varying number of times before
using medaka. The result being that the resultant accuracy was unchanged. In a
single tweet there is not space to pose or answer questions such as: is the
iterated racon procedure actually producing a different result? Naturally if
the draft does not evolve during iteration no variation in the output of
``medaka`` is to be expected.

.. raw:: html

    <blockquote class="twitter-tweet" data-conversation="none"><p lang="en" dir="ltr">Using the scripts for assembly benchmark from <a href="https://twitter.com/nanopore?ref_src=twsrc%5Etfw">@nanopore</a> (<a href="https://t.co/ezulj7RFN2">https://t.co/ezulj7RFN2</a>) I also calculated consensus Q-scores. For this example, I achieved &gt; Q35 <a href="https://t.co/XLPdvUFZYr">pic.twitter.com/XLPdvUFZYr</a></p>&mdash; Devin Drown (@devindrown) <a href="https://twitter.com/devindrown/status/1159223021246095360?ref_src=twsrc%5Etfw">August 7, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

One point to note is that with older basecallers than Devin has used, where
assembler output was of lower quality, the draft quality would improve with
subsequent applications of ``racon``. With the latest basecallers improved
accuracy of drafts can make mute the concerns noted above; ``medaka``
may now be less sensitive to variations in precise properties of the draft


How should I create my draft sequence?
--------------------------------------

Given the above discussion the question remains, what is the best way to create
the draft sequence in order to achieve optimal post-``medaka`` results?  As
noted the previous recommendations were based on how ``medaka`` models
are trained. Growing empirical evidence has shown that ``medaka`` is not overly
sensitive to the nature of the draft input. Devin's tweet does show there is
still some dependence on the general quality level. More debugging of this
example is warranted but a general effect could be the reason the single
application of ``racon`` is better than the no-``racon`` case: an
improvement in draft quality leads to ``medaka`` viewing a higher effective
coverage as more reads align the the improved draft. It should be stressed that
there is a difference between general quality and the notions of draft nature
and error profiles discussed previously. Two drafts could have very similar
overall quality but have very different error profiles.

There is a known case where is it possible to `over-polish` a draft sequence,
as noted by Sergey Koren:

.. raw:: html

    <blockquote class="twitter-tweet" data-conversation="none"><p lang="en" dir="ltr">Why 4x, have seen good evidence &gt;1 racon corrupts repeats in assembly, especially of diploids. Examples for canu of correctly resolved BACs: 0=74%, 1 = 77%, 2=70%, 3=63%. For miniasm: 0= 21%, 1=60%, 2=52%. Identity does go up but very slightly. I&#39;d vote for 1 round.</p>&mdash; Sergey Koren (@sergekoren) <a href="https://twitter.com/sergekoren/status/1158850929338310657?ref_src=twsrc%5Etfw">August 6, 2019</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

One reason for this is that reads can be calculated as overlappin the wrong
copy of a repeat; the iterated application of ``racon`` (or ``medaka``) can
lead to an averaging of these repeats rather than maintaining their distinct
sequence.

Recommendations
...............

It is hard to recommend a single best-practice workflow to obtain optimal results
on all datasets. However it is still the case that running ``racon`` at least
once before giving a dataset to ``medaka`` is recommended. The primary effect
of doing this is to normalize the error profile of the draft regardless of the
assembler used. As the basecall quality and so assembler output quality
improves this normalization step could become redundant when the distinction
between quality differences and error profiles becomes blurred.

During training of ``medaka`` non-default parameter values are used with ``racon``.
The arguments presented above regarding the relative indel frequencies between
reads and draft likely still have relevance. It is therefore recommend to use these
with racon before using ``medaka``:

.. code:: bash

    racon -m 8 -x -6 -g -8 -w 500 ...

The effect of deviating from this prescription has not been explored with recent
basecallers, it may well be the case that ``medaka`` is not overly sensitive to
changes to these values.

Future developments
...................

An open question is how to remove entirely the dependence of ``medaka`` results
on the nature of the draft sequence, at least to the point that all that is
required is a reasonably accurate draft. Methods to remove the draft from the
calculation entirely and examination of the read data at face value are being
explored. Additionally methods to more fully utilise the rich output of
assemblers are in development, to avoid some of the more subtle issues. As a
separate, simpler possiblity, users may wish to train ``medaka`` models against
the output of an assembler of their choice.

Medaka is actively developed both for consensus and variant calling, with continuous
improvements being released under and open source model via
`github <https://github.com/nanoporetech/medaka>`_,
`pypi <https://pypi.org/project/medaka/>`_, and
`bioconda <https://bioconda.github.io/recipes/medaka/README.html>`_.

