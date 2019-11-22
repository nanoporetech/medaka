"""Module to mock data needed for several tests to run"""
import array
import os
import tempfile

import h5py
import numpy as np
import pysam

from medaka import common

#  Ref          A    C    A    T    *    G    A    T    G
#
# Basecall1:   2A   1C   4A   5T        1G   1A   2T   1G
# Basecall2:   0A   1C   4A    *        1G   1A   1T   2G
# Basecall3:   2a   1c   4a   5t   1a   1g   1a   2t   1g
# Basecall4:   2a   1c   4a   1c        1g   1a   2t   1g
#
# The zero entry is deliberate to trigger edge case, it should be
# treated as a 1, previously caused invalid memory access.

simple_data = {
    'ref': 'ACATGATG',
    'truth': {
        'query_name': 'truth',
        'seq': 'ACATAGATCTG',  # the A from third below and another CT
        'quality':  array.array('B', [2, 1, 4, 5, 1, 1, 1, 2, 1, 1, 1]),
        'cigarstring': '4=1I3=2I1=',
        'flag': 0,
        'tags': {'MD': '8'}
    },
    'calls': [
        {
            'query_name': 'basecall_1',
            'seq': 'ACATGATG',  # exactly ref
            'quality': array.array('B', [2, 1, 4, 5, 1, 1, 2, 1]),
            'cigarstring': '8=',
            'flag': 0,
            'tags': {
                'cs': '=ACATGATG',
                'AA': 1,
                'WL': [1.5, 0.5, 3.5, 4.5, 0.5, 0.5, 1.5, 0.5],
                'WK': [1e3] * 8,  # sharply peaked, on wl + 0.5
                'DT': 'r9'}
            },
        {
            'query_name': 'basecall_2',
            'seq': 'ACAGATG',  # deletion of T in the middle
            'quality': array.array(
                'B', [0, 1, 4, 1, 1, 1, 2]),  # zero, see above
            'cigarstring': '3=1D4=',
            'flag': 0,
            'tags': {
                'cs': '=ACA-t=GATG',
                'AA': 1,
                'WL': [1.0] * 7,
                'WK': [1.0] * 7,
                'DT': 'r9'}
            },
        {
            'query_name': 'basecall_3',
            'seq': 'ACATAGATG',  # insertion of A in the middle
            'quality':  array.array('B', [2, 1, 4, 5, 1, 1, 1, 2, 1]),
            'cigarstring': '4=1I4=',
            'flag': 16,
            'tags': {
                'cs': '=ACAT+a=GATG',
                'AA': 2,
                'WL': [1.0] * 9,
                'WK': [1.0] * 9,
                'DT': 'r9'}
            },
        {
            'query_name': 'basecall_4',
            'seq': 'ACACGATG',  # substitution T->C
            'quality': array.array('B', [2, 1, 4, 1, 1, 1, 2, 1]),
            'cigarstring': '3=1X4=',
            'flag': 16,
            'tags': {
                # no 'AA' tag to test "keep_missing" behaviour of medaka.features.PileupCounts
                'cs': '=ACA*tc=GATG',
                'WL': [1.0] * 8,
                'WK': [1.0] * 8,
                'DT': 'r10'}
            }
    ]
}


def create_simple_bam(fname, calls):
    """Create a small bam file with RLE encoding coded in the qscores."""
    ref_len = len(simple_data['ref'])

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': ref_len, 'SN': 'ref'}]}

    tmp_file = '{}.tmp'.format(fname)
    with pysam.AlignmentFile(
            tmp_file, 'wb', reference_names=['ref', ],
            reference_lengths=[ref_len, ], header=header) as output:
        for index, basecall in enumerate(calls):
            a =common.initialise_alignment(
                query_name=basecall['query_name'],
                reference_id=0,
                reference_start=0,
                query_sequence=basecall['seq'],
                cigarstring=basecall['cigarstring'],
                flag=basecall['flag'],
                query_qualities=basecall['quality'],
                tags=basecall['tags'])

            output.write(a)

    pysam.sort("-o", fname, tmp_file)
    os.remove(tmp_file)
    pysam.index(fname)


def mock_fast5_file():
    """Create a fast5 file with fake shape/scale rle data"""
    fast5_file = tempfile.NamedTemporaryFile(suffix='.fast5').name
    data_path = (
        'read_{}/Analyses/Basecall_1D_000/'
        'BaseCalled_template/RunlengthBasecall')

    with h5py.File(fast5_file, 'w') as h:
        for call in simple_data['calls']:
            read_id = call['query_name']

            bases = call['seq']
            shapes = call['tags']['WL']
            scales = call['tags']['WK']

            # (basecall, WL, WK) defined above wrt ref.
            if call['flag'] == 16:
                bases = common.reverse_complement(bases)
                shapes = shapes[::-1]
                scales = shapes[::-1]

            arr = np.fromiter(
                [(b, s, sc) for (b, s, sc) in zip(bases, shapes, scales)],
                dtype=[('base', 'S1'), ('shape', '>f4'), ('scale', '>f4')])
            h.create_dataset(data_path.format(read_id), data=arr)

    return fast5_file


def mock_summary_file(fast5_fname):
    """Create a summary file to link read_id and fast5 filename"""
    read_ids = [x['query_name'] for x in simple_data['calls']]

    summary_file = tempfile.NamedTemporaryFile(suffix='.txt').name
    with open(summary_file, 'w') as output:
        output.write('read_id\tfilename\n')
        for read_id in read_ids:
            output.write('{}\t{}\n'.format(read_id, fast5_fname))

    return summary_file
