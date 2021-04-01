import operator
import os
import unittest
from unittest.mock import patch
import tempfile

import numpy as np

from medaka.common import Region
import medaka.stitch

class TestStitch(unittest.TestCase):

    def test_fill_gaps(self):
        bases = 'ATGCN'
        def rand_seq(n):
            return  ''.join(np.random.choice(list(bases), n, replace=True))
        cases = [
            (
                'no_gaps',  # no gap-filling required
                rand_seq(50),
                ((0, 50),),
                ()),
            (
                'inside_gaps',
                rand_seq(25),
                ((0, 5), (6, 7), (10, 12), (16, 25)),
                ((5, 6), (7, 10), (12, 16))),
            (
                'outside_gaps',  # name
                rand_seq(20),  # draft seq
                ((3, 7), (11, 15)),  # polished contig boundaries (end exclusive)
                ((0, 3), (7, 11), (15, 20))),  # implied gaps
        ]

        _, draft = tempfile.mkstemp()
        contigs = []
        with open(draft, 'w') as fh:
            for ref_name, full_seq, bounds, gaps, in cases:
                fh.write('>{}\n{}\n'.format(ref_name, full_seq))
                contigs.extend([
                    ((ref_name, start,  end - 1), [full_seq[start:end]])
                    for i, (start, end) in enumerate(bounds)])

        results, gap_trees = medaka.stitch.fill_gaps(contigs, draft)

        self.assertEqual(len(results), len(cases))
        for i, (case, result) in enumerate(zip(cases, results)):
            # check filled sequence
            self.assertEqual(
                case[1], ''.join(result[1]), msg="Sequence for case '{}'.".format(case[0]))
            # check gaps are correct
            msg = 'Gaps incorrect for case {}.'.format(i,)
            sorted_gaps = tuple([(i.begin, i.end) for i in sorted(
                gap_trees[case[0]], key=operator.attrgetter('begin'))])
            self.assertEqual(case[-1], sorted_gaps, msg=msg)
        os.remove(draft)

    def test_010_fill_neighbours(self):
        seq = 'ATCG'
        contigs = iter([
            (('chr1', 0, 1000), [seq] * 2),
            (('chr1', 1001, 2000), [seq] * 2),  # join with above
            (('chr2', 0, 1000), [seq] * 2),
            (('chr2', 2000, 3000), [seq] * 2),  # not joined
            (('chr3', 3001, 1000), [seq] * 2)])  # not joined
        expected = [
            (('chr1', 0, 2000), [seq] * 4),
            (('chr2', 0, 1000), [seq] * 2),
            (('chr2', 2000, 3000), [seq] * 2),  # not joined
            (('chr3', 3001, 1000), [seq] * 2)]  # not joined
        contigs = list(medaka.stitch.collapse_neighbours(contigs))
        self.assertEqual(contigs, expected, 'Fill neighbours')


class RegressionStitch(unittest.TestCase):

    @patch('pysam.FastaFile')
    def test_001_cases(self, FastaFile):
        class Args:
            def __init__(self, **kwargs):
                for k, v in kwargs.items():
                    setattr(self, k, v)
        files = {1:"utg1190", 2:"scaffold_117"}
        FastaFile.return_value = Args(references=['scaffold_117', 'utg1190'], lengths=[45079626, 16772114])
        temp = tempfile.NamedTemporaryFile()
        args = Args(draft="", threads=1, output=temp.name, fillgaps=False)

        for fid, region in files.items():
            fname = os.path.join(os.path.dirname(__file__), "data", "test_stitch_{}.hdf".format(fid))
            args.inputs = fname
            args.regions=[Region.from_string(region)]
            try:
                medaka.stitch.stitch(args)
            except Exception as e:
                self.fail("Stitching raise and Exception:\n {}".format(e))
