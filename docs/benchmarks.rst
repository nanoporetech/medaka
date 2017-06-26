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

+----------------------------+---------------------+-----------------+-----------------------------+
|                            | 5mer non-transducer | 5mer transducer | 6mer transducer             |
|                            +------------+--------+--------+--------+-------+--------+------------+
|                            |      draft | medaka |  draft | medaka | draft | medaka | nanopolish |
+===================+========+============+========+========+========+=======+========+============+
| **Total Error %** |        |      0.648 |  0.320 |  0.578 |  0.255 | 0.203 |  0.169 |      0.074 |
+-------------------+--------+------------+--------+--------+--------+-------+--------+------------+
| **Mismatch**      |        |      0.008 |  0.009 |  0.019 |  0.015 | 0.017 |  0.005 |      0.007 |
+-------------------+--------+------------+--------+--------+--------+-------+--------+------------+
| **Deletion**      | Total  |      0.638 |  0.277 |  0.532 |  0.184 | 0.157 |  0.127 |      0.042 |
+                   +--------+------------+--------+--------+--------+-------+--------+------------+
|                   | Non-HP |      0.092 |  0.021 |  0.109 |  0.046 | 0.008 |  0.004 |      0.004 |
+                   +--------+------------+--------+--------+--------+-------+--------+------------+
|                   | HP     |      0.546 |  0.256 |  0.423 |  0.138 | 0.148 |  0.123 |      0.038 |
+-------------------+--------+------------+--------+--------+--------+-------+--------+------------+
| **Insertion**     | Total  |      0.002 |  0.034 |  0.027 |  0.056 | 0.029 |  0.037 |      0.025 |
+                   +--------+------------+--------+--------+--------+-------+--------+------------+
|                   | Non-HP |      0.002 |  0.019 |  0.025 |  0.034 | 0.029 |  0.016 |      0.014 |
+                   +--------+------------+--------+--------+--------+-------+--------+------------+
|                   | HP     |      0.000 |  0.015 |  0.002 |  0.022 | 0.001 |  0.020 |      0.011 |
+-------------------+--------+------------+--------+--------+--------+-------+--------+------------+

In the best case (5mer non-transducer) medaka reduces the total error by half.
This is achieved mainly through reducing the homopolymer (three or more bases)
deletion error. A similar level of improvement is seen in the 5mer
transducer case, whilst the relative gain is smaller in the 6mer case where
the quality of the draft assembly is much improved compare to the previous
cases. This rather simple correction model does not reach the level of
nanopolish, but does reduce the runtime of nanopolish considerably.

Future versions of medaka will aim to improve on the above results with the
aim to surpass the nanopolish results whilst also improving runtime. See
:ref:`FutureDirections` for more information.

