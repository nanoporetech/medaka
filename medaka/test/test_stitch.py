import operator
import os
import unittest
import tempfile

import numpy as np

import medaka.stitch

class TestStitch(unittest.TestCase):

    def test_fill_gaps(self):

        # create gapped consensus contigs
        bases = 'ATGCN'
        def rand_seq(n):
            return  ''.join(np.random.choice(list(bases), n, replace=True))
        cases = [
            ('chr1',  # name
             'chr1',  # info (after gap-filling stitching)
             rand_seq(20),  # draft seq
             ((3, 7), (11, 15)), # polished contig boundaries (end exclusive)
             ((0, 3), (7, 11), (15, 20))),  # gaps
            ('chr2',
             'chr2',
             rand_seq(25),
             ((0, 5), (6, 7), (10, 12), (16, 25)),
             ((5, 6), (7, 10), (12, 16))),
            ('chr3',  # no gap-filling required
             'chr3',
             rand_seq(50),
             ((0, 50),),
             ())
        ]

        _, draft = tempfile.mkstemp()
        contigs = []
        with open(draft, 'w') as fh:
            for ref_name, _, full_seq, bounds, gaps, in cases:
                fh.write('>{}\n{}\n'.format(ref_name, full_seq))
                contigs.extend([
                    ('{}_segment{}'.format(ref_name, i),
                    # -1 as medaka sample names are end inclusive
                    '{}:{}.0-{}.0'.format(ref_name, start,  end - 1),
                    full_seq[start:end])
                    for i, (start, end) in enumerate(bounds)
                ])

        results, gap_trees = medaka.stitch.fill_gaps(contigs, draft)

        self.assertEqual(len(results), len(cases))
        for i, (case, result) in enumerate(zip(cases, results)):
            for item, descr in (0, 'chrom'), (1, 'info'), (2, 'sequence'):
                msg = 'Failed case {} {}'.format(i, descr)
                self.assertEqual(case[item], result[item], msg=msg)
            # check gaps are correct
            sorted_gaps = tuple([(i.begin, i.end) for i in sorted(
                gap_trees[case[0]], key=operator.attrgetter('begin'))])
            self.assertEqual(case[-1], sorted_gaps, msg=msg.format(i, 'gaps'))

        os.remove(draft)
