History and Future Directions
=============================

Medaka has grown out of some early work on the study of errors in basecalls
from recurrent neural network (RNN) basecallers of Oxford Nanopore
Technologies' (ONT) data. It is relatively straightforward to enumerate
dominant error classes in reads and correct for some of these. For example it
is well known that deletion of bases forms a large proportion of all errors;
occurring both within and outside homopolymers. Perhaps not so well known is
that these deletions can show a 'strand' bias, being present in one orientation
of reads but not the other. This fact can be exploited. Results from ad-hoc,
aggressive repairing of deletions was presented by ONT at London Calling 2017.

Having corrected dominant error classes it becomes progressively harder to
achieve further improvement by manual inspection. At this point persuing
automated methods to learn and correct for errors becomes prudent. Medaka
attempts this by implementing an interface to prepare training data from alignment
of basecalls and a truth set to a common baseline. This baseline may be
anything, perhaps a draft assembly or chosen single molecule basecalls. Once
these data have been collated they may be used directly or prepared further to
exploit known errors. For example one may reduce or augment the with data counts
of orientated bases and indels as this is known to be a relavant consideration.

Medaka allows researchers to experiment with their own consensus ideas without
having to write much of the tedious data preparation code. Researchers can
extend the tools provided to achieve better results that those obtained
currently. There are certainly obvious extension which can be made, including:

* use of quality scores from the basecaller,
* relax contraints on how features are selected: there are currently some
  shortcuts taken which could better handle alignment subtleties,
* modifying the feature vectors to take wider windows of data,
* altering the size, structure, and connections between the RNN layers
* use of convolutional layers.

Using lower level data directly (not basecalls alone) can provide a more
powerful method; this is something which is being actively researched. Early
examples include `nanopolish <https://github.com/jts/nanopolish>`_ and
`poreseq <https://github.com/tszalay/poreseq>`_ which take the approach of
iteratively refining a candidate sequence by alignment of nanopore event data
to the candidate sequence. Where their data models predict discrepancies
between the observed data and the candidate sequence, these methods mutate the
candidate to achieve better alignment scores. In this way the methods maximise
the probability of having observed all input event data under their model.

Neural network classifiers have been shown (in the context of basecalling) to
model better the primary raw, and secondary event, data from nanopore devices
than the generative Hidden Markov Models used in both of the above methods.
It is natural therefore to desire to encorporate neural networks into a
procedure to error correct candidate sequences from primary and secondary
data; this is the focus of current research. Further it may be possible that
such an approach need not be iterative: a single pass on the inputs could be
performed to achieve results in quicker time.

Nevertheless researchers should find medaka useful as a method to generate
quickly accurate consensus sequences in many use cases; see the
:ref:`Benchmarks` page for more details.
