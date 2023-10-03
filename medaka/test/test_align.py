import array
import tempfile
import unittest

import parasail

import medaka.align


@unittest.skipIf(parasail.dnafull is None, "Using fake parasail.")
class ParasailAlignment(unittest.TestCase):
    """Check medaka.align.parasail_alignment function."""

    ref = 'AGCATGTTAGATAAGATA'

    def test_alignment_001(self):
        """Test alignment of read starting mid-reference and clip at end.

        AGCATGTTAGATAA**GATA
        TTAGATAAAGGATActg
        """
        seq = 'TTAGATAAAGGATACTG'
        rstart, cigar = medaka.align.parasail_alignment(seq, self.ref)
        expected = (6, '8=2I4=3S')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')

    def test_alignment_002(self):
        """Test alignment of read starting mid-reference with clip at begining.

        AGCATGTTAGATAA*GATA
        aaaAGATAAGGATA
        """
        seq = 'AAAAGATAAGGATA'
        rstart, cigar = medaka.align.parasail_alignment(seq, self.ref)
        expected = (8, '3S6=1I4=')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')

    def test_alignment_003(self):
        """Test alignment with clipping in both ends.

        AGCATGTTAGATAAGATA
        ccCATGTTAGATAAGcc
        """
        seq = 'CCCATGTTAGATAAGCC'
        rstart, cigar = medaka.align.parasail_alignment(seq, self.ref)
        expected = (2, '2S13=2S')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')

    def test_alignment_004(self):
        """Test alignment of read starting before the reference.

        AGCATGTTAGATAAGATA
        ctgAGCATGT
        """
        seq = 'CTGAGCATGT'
        rstart, cigar = medaka.align.parasail_alignment(seq, self.ref)
        expected = (0, '3S7=')
        self.assertEqual(cigar, expected[1], 'cigars differ!')
        self.assertEqual(rstart, expected[0], 'rstart wrong!')


class CigarProcessing(unittest.TestCase):
    """Test processing cigar strings."""

    cigar1 = '4S2=1D40=1X29=2D12=1I10=3S'
    cigar2 = '4S2=1D40=1X29=2D12=1I10='
    cigar3 = '2=1D40=1X29=2D12=1I10=3S'
    cigar4 = '2=1D40=1X29=2D12=1I10='
    cigars = (cigar1, cigar2, cigar3, cigar4)

    def test_cigar_ops_from_start_end(self):
        split1 = (('4', 'S'), ('2', '='), ('1', 'D'), ('40', '='), ('1', 'X'),
                 ('29', '='), ('2', 'D'), ('12', '='), ('1', 'I'), ('10', '='),
                 ('3', 'S'))
        self.assertEqual(split1,
                         tuple(medaka.align.cigar_ops_from_start(self.cigar1)))
        self.assertEqual(split1[::-1],
                         tuple(medaka.align.cigar_ops_from_end(self.cigar1)))

    def test_trim_cigar_start(self):
        base_cigar = '1=1D40=1X29=2D12=1I10='
        cases = [('4I', 4, 0), # cigar prefix, trim_q, rstart_offset
                 ('2D', 0, 2),
                 ('3X', 3, 3),
                 ('1I2X', 3, 2)]
        for prefix, expt_q_trim, expt_rstart_offset in cases:
            cigar = prefix + base_cigar
            trimmed, q_trim, r_offset = medaka.align.trim_cigar(cigar,
                                                                start=True)
            self.assertEqual(base_cigar, trimmed)
            self.assertEqual(expt_q_trim, q_trim)
            self.assertEqual(expt_rstart_offset, r_offset)

    def test_no_added_clipping(self):
        """Do not add anything, cigar should return untouched."""
        for cigar in self.cigars:
            got = medaka.align.add_extra_clipping(cigar, 0, 0)
            self.assertEqual(cigar, got)

    def test_add_start(self):
        """Add clipping only at the beginning of the cigar."""
        start_clipped = 4
        expected = (
            '8S2=1D40=1X29=2D12=1I10=3S',
            '8S2=1D40=1X29=2D12=1I10=',
            '4S2=1D40=1X29=2D12=1I10=3S',
            '4S2=1D40=1X29=2D12=1I10=')
        got = tuple(
            medaka.align.add_extra_clipping(cigar, start_clipped, 0)
            for cigar in self.cigars)
        self.assertEqual(expected, got)

    def test_add_end(self):
        """Add clipping at the end of the cigar."""
        end_clipped = 8
        expected = (
            '4S2=1D40=1X29=2D12=1I10=11S',
            '4S2=1D40=1X29=2D12=1I10=8S',
            '2=1D40=1X29=2D12=1I10=11S',
            '2=1D40=1X29=2D12=1I10=8S')
        got = tuple(
            medaka.align.add_extra_clipping(cigar, 0, end_clipped)
            for cigar in self.cigars)
        self.assertEqual(expected, got)

    def test_add_both(self):
        """Add clipping both ends of the cigar."""
        start_clipped = 6
        end_clipped = 9
        expected = (
            '10S2=1D40=1X29=2D12=1I10=12S',
            '10S2=1D40=1X29=2D12=1I10=9S',
            '6S2=1D40=1X29=2D12=1I10=12S',
            '6S2=1D40=1X29=2D12=1I10=9S')
        got = tuple(
            medaka.align.add_extra_clipping(cigar, start_clipped, end_clipped)
            for cigar in self.cigars)
        self.assertEqual(expected, got)


class InitialiseAlignment(unittest.TestCase):
    """Test `medaka.align.initialise_alignment`."""

    input_kwargs = {
        'query_name': 'test',
        'reference_id': 0,
        'reference_start': 2,
        'query_sequence': 'GCCCTGTTGATCTT',
        'cigarstring': "1S3=1D8=2S",
        'flag': 0,
        'mapping_quality': 60}
    len_query = len(input_kwargs['query_sequence'])
    input_kwargs['query_qualities'] = array.array('B', [20] * len_query)

    def test_inputs(self):
        """Test inputs are correctly passed to alignment."""
        alignment = medaka.align.initialise_alignment(**self.input_kwargs)
        for key, expected in self.input_kwargs.items():
            got = getattr(alignment, key)
            self.assertEqual(expected, got)

    def test_derived(self):
        """Test arguments derived from inputs."""
        alignment = medaka.align.initialise_alignment(**self.input_kwargs)
        expected_kwargs = {
            'query_alignment_start': 1,
            'query_alignment_end': 12,
            'query_alignment_sequence': 'CCCTGTTGATC'}

        for key, expected in expected_kwargs.items():
            got = getattr(alignment, key)
            self.assertEqual(got, expected)


class ChunkedEdlibAlignment(unittest.TestCase):
    """Test `medaka.align.chunked_edlib_align."""

    rseq = 'ATCGGGATAG'
    qseq = 'ATCGGGCTAG'

    def test_no_edits_no_chunking(self):
        """Test a short set aligns with no edits to itself in one chunk."""

        seq = 'ATCGGGATAG'
        alns = list(medaka.align.chunked_edlib_align(seq, seq, 'contig1',
                                                     chunk_size=2*len(seq)))
        self.assertEqual(len(alns), 1)
        aln = alns[0]
        self.assertEqual(aln.cigarstring, '{}='.format(len(seq)))


    def test_no_edits_chunking(self):
        """Test arguments derived from inputs."""
        seq = 'ATCGGGATAGATTG'
        # alignments have min and max sizes of 2 and chunk_size + 1 due to
        # chunks overlapping by 1 base.
        cases = [
            #seq_len, [(rstart, rend, chunk length)]
            (5, [(0, 5, 5)]),
            (6, [(0, 5, 5), (4, 6, 2)]),
            (7, [(0, 5, 5), (4, 7, 3)]),
            (8, [(0, 5, 5), (4, 8, 4)]),
            (9, [(0, 5, 5), (4, 9, 5)]),
            (10, [(0, 5, 5), (4, 10, 6)]),
            (11, [(0, 5, 5), (4, 9, 5), (8, 11, 3)]),
            (12, [(0, 5, 5), (4, 9, 5), (8, 12, 4)]),
            (13, [(0, 5, 5), (4, 9, 5), (8, 13, 5)]),
            (14, [(0, 5, 5), (4, 9, 5), (8, 14, 6)]),
        ]
        for l, exp_alns in cases:
            s = seq[:l]
            alns = list(medaka.align.chunked_edlib_align(s, s, 'contig1',
                                                         chunk_size=5))
            self.assertEqual(len(alns), len(exp_alns))
            aln_tups = [(a.reference_start, a.reference_end,
                     a.query_alignment_length) for a in alns]
            self.assertEqual(aln_tups, exp_alns)


    def test_indels_at_ends(self):
        """Test how various alignment modes handle indels at seq ends."""

        r = '*ATCGGGATAG*'.replace('*', '')
        q = 'GATCGGGATAGC'.replace('*', '')
        lr = len(r)
        lq = len(q)
        cases = [
            # query, ref, mode, cigar, qstart, qend, rstart, rend, chnk cigars
            (q, r, 'HWT', '10=', 1, lq - 1,  0, lr, ['3=', '4=', '4=', '2=']),
            (q, r, 'HW', '1I10=1I', 0, lq,  0, lr,
             ['1I3=', '4=', '4=', '2=1I']),
            (q, r, 'NW', '1I10=1I', 0, lq,  0, lr,
             ['1I3=', '4=', '4=', '2=1I']),
            # now swap around query and ref, as edlib is global in query
            # the NW case differs.
            (r, q, 'HWT', '10=', 0, lr, 1, lq - 1, ['4=', '4=', '4=']),
            (r, q, 'HW', '10=', 0, lr, 1, lq - 1, ['4=', '4=', '4=']),
            (r, q, 'NW', '1D10=1D', 0, lr,  0, lq, ['1D4=', '4=', '4=1D']),
        ]
        for i, (q, r, mode, cigar, qstart, qend, rstart, rend, chnk_cigs) in \
                enumerate(cases):
            msg = 'Failed case {} unchunked.'.format(i)
            alns = list(medaka.align.chunked_edlib_align(q, r, 'contig1',
                                                         chunk_size=1000,
                                                         mode=mode))
            self.assertEqual(len(alns), 1, msg=msg)
            a = alns[0]
            qtrim = q[qstart:qend]
            self.assertEqual(a.reference_start, rstart, msg=msg)
            self.assertEqual(a.reference_end, rend, msg=msg)
            self.assertEqual(a.query_alignment_sequence, qtrim, msg=msg)
            self.assertEqual(a.query_alignment_start, 0, msg=msg)
            self.assertEqual(a.query_alignment_end, len(qtrim), msg=msg)
            self.assertEqual(a.cigarstring, cigar, msg=msg)

            msg = 'Failed case {} chunked.'.format(i)
            alns = list(medaka.align.chunked_edlib_align(q, r, 'contig1',
                                                         chunk_size=4,
                                                         mode=mode))
            self.assertEqual(len(alns), len(chnk_cigs), msg=msg)
            self.assertEqual(alns[0].reference_start, rstart, msg=msg)
            self.assertEqual(alns[-1].reference_end, rend, msg=msg)
            self.assertListEqual([a.cigarstring for a in alns], chnk_cigs)


    def test_snps_at_ends(self):
        """Test how various alignment modes handle snps at seq ends."""

        r = 'CATCGGGATAGT'
        q = 'GATCGGGATAGC'
        lr = len(r)
        lq = len(q)
        cases = [
            # query, ref, mode, cigar, qstart, qend, rstart, rend, chnk cigars
            (q, r, 'HWT', '10=', 1, lq - 1,  1, lr - 1,
             ['3=', '4=', '4=', '2=']),
            # with edlib, the cost of aligning an indel is same as snp..
            # hence HW probably not best for robust alignment at ends.
            (q, r, 'HW', '1X10=1I', 0, lq,  0, lr - 1,
             ['1X3=', '4=', '4=', '2=1I']),
            (q, r, 'NW', '1X10=1X', 0, lq,  0, lr,
             ['1X3=', '4=', '4=', '2=1X']),
        ]
        for i, (q, r, mode, cigar, qstart, qend, rstart, rend, chnk_cigs) in \
                enumerate(cases):
            msg = 'Failed case {} unchunked.'.format(i)
            alns = list(medaka.align.chunked_edlib_align(q, r, 'contig1',
                                                         chunk_size=1000,
                                                         mode=mode))
            self.assertEqual(len(alns), 1, msg=msg)
            a = alns[0]
            qtrim = q[qstart:qend]
            self.assertEqual(a.reference_start, rstart, msg=msg)
            self.assertEqual(a.reference_end, rend, msg=msg)
            self.assertEqual(a.query_alignment_sequence, qtrim, msg=msg)
            self.assertEqual(a.query_alignment_start, 0, msg=msg)
            self.assertEqual(a.query_alignment_end, len(qtrim), msg=msg)
            self.assertEqual(a.cigarstring, cigar, msg=msg)

            msg = 'Failed case {} chunked.'.format(i)
            alns = list(medaka.align.chunked_edlib_align(q, r, 'contig1',
                                                         chunk_size=4,
                                                         mode=mode))
            self.assertEqual(len(alns), len(chnk_cigs), msg=msg)
            self.assertEqual(alns[0].reference_start, rstart, msg=msg)
            self.assertEqual(alns[-1].reference_end, rend, msg=msg)
            self.assertListEqual([a.cigarstring for a in alns], chnk_cigs)
