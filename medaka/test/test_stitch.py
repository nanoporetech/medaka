import operator
import os
import random
import unittest
from unittest.mock import patch
import tempfile

import numpy as np
import pysam

from medaka.common import Region
import medaka.stitch

class TestStitch(unittest.TestCase):

    def test_010_fill_gaps(self):
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

    def test_020_fill_neighbours(self):
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


class Args:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

class MyFastaFile:

    def __init__(*args, **kwargs):
        pass

    @property
    def references(self):
        return ['utg1190', 'scaffold_117']

    @property
    def lengths(self):
        return [16772114, 45079626]

    def fetch(self, reference, start=None, end=None, region=None):
        if start is None:
            start = 0
        if end is None:
            end = self.lengths[self.references.index(reference)]
        return "".join(random.choices('ACGT', k=end - start))



class RegressionStitch(unittest.TestCase):

    @patch('pysam.FastaFile', MyFastaFile)
    def _run_one(self, expected, fillgaps=False, regions=None):
        args = Args(draft="", threads=1, regions=regions, fillgaps=fillgaps)

        outputs = list()
        for fid, (region, exp), in enumerate(zip(MyFastaFile().references, expected), 1):
            with tempfile.NamedTemporaryFile(delete=False) as temp:
                args.output = temp.name
                args.inputs = os.path.join(os.path.dirname(__file__), "data", "test_stitch_{}.hdf".format(fid))
                try:
                    medaka.stitch.stitch(args)
                except Exception as e:
                    self.fail("Stitching raise and Exception:\n {}".format(e))
                outputs.append(temp)
        return outputs


    def _check(self, files, expected):
        for file, exp in zip(files, expected):
            with pysam.FastaFile(file.name) as fh:
                self.assertEqual(exp, fh.references)


    def test_001_cases(self):
        # each file should have just one entry
        expected = [
            ['utg1190_0'],
            ['scaffold_117_0']]
        files = self._run_one(expected, regions=None, fillgaps=False)
        self._check(files, expected)


    def test_002_copy_missing(self):
        # each file should have both entries, note how the missing
        # contig comes second
        expected = [
            ['utg1190', 'scaffold_117'],
            ['scaffold_117', 'utg1190']]
        files = self._run_one(expected, regions=None, fillgaps=True)
        self._check(files, expected)
