Benchmarks
==========

The following demonstrates the utility of a recurrent neural network trained
to make corrections of a draft assembly by considering counts of bases and
indels in pileup columns obtained from reads aligned to the draft, as
implemented with :ref:`SequenceCorrection`.

All models here were trained from reads aligning to a region of E.coli and
tested on a distinct, non-overlapping region. In all cases
`scrappie <https://github.com/nanoporetech/scrappie>`_ was used to perform
the basecalling; the 6mer transducer model was trained specifically for this
experiment. The draft assemblies here were created using the
`mini_assemble <https://nanoporetech.github.io/pomoxis/examples.html#fast-de-novo-assembly>`_
pipeline in `pomoxis <https://github.com/nanoporetech/pomoxis>`_. Statistics
were calculated using `alignqc <https://www.healthcare.uiowa.edu/labs/au/AlignQC/>`_.

+----------------------------+-----------------+
|                            | Raw transducer  |
|                            +--------+--------+
|                            |  draft | medaka |
+===================+========+========+========+
| **Total Error %** |        |  0.269 |  0.081 |
+-------------------+--------+--------+--------+
| **Mismatch**      |        |  0.007 |  0.011 |
+-------------------+--------+--------+--------+
| **Deletion**      | Total  |  0.222 |  0.045 |
+                   +--------+--------+--------+
|                   | Non-HP |  0.011 |  0.005 |
+                   +--------+--------+--------+
|                   | HP     |  0.211 |  0.040 |
+-------------------+--------+--------+--------+
| **Insertion**     | Total  |  0.040 |  0.026 |
+                   +--------+--------+--------+
|                   | Non-HP |  0.019 |  0.004 |
+                   +--------+--------+--------+
|                   | HP     |  0.021 |  0.022 |
+-------------------+--------+--------+--------+

Medaka reduces the total error by more than a factor of three. This is achieved
mainly through reducing the homopolymer (three or more bases) deletion error. 

This rather simple correction model does not reach the level of
nanopolish, but does reduce the runtime of nanopolish considerably.

Future versions of medaka will aim to improve on the above results with the
aim to surpass the nanopolish results whilst also improving runtime. See
:ref:`FutureDirections` for more information.

