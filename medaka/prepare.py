"""
Prepare data for error model training.
"""

import argparse
from functools import partial
from multiprocessing import Pool
import h5py
import numpy as np

from medaka.util.pileup import (prepare_training_data, get_reference_names,
                                get_reference_length)
from medaka.util.sequences import segment_limits

def get_parser():
    parser = argparse.ArgumentParser(
        description="""Generate training data from bam files.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--threads', type=int, default=1,
                        help='number of processes')
    parser.add_argument('truth_to_ref',
                        help='bam file with true sequence aligned to reference.')
    parser.add_argument('reads_to_ref',
                        help='bam file with reads aligned to reference.')
    parser.add_argument('output_name', type=str, default='',
                        help='name for output dataset.')
    return parser


def multi_prepare_reference(reads_bam, truth_bam, reference, threads=1):
    """prepare reference segments in parallel"""
    ref_length = get_reference_length(reads_bam, reference)
    limits = segment_limits(ref_length, overlap_len=0)
    func = partial(prepare_training_data, reads_bam, truth_bam, reference)
    pool = Pool(threads)
    for d, l in pool.imap(func, limits):
        yield (d, l)
    pool.close()
    pool.join()


def main():
    args = get_parser().parse_args()
    references = get_reference_names(args.truth_to_ref)

    data = []
    label = []

    for reference in references:
        for X, y in multi_prepare_reference(args.reads_to_ref,
                                            args.truth_to_ref,
                                            reference, threads=args.threads):

            data.append(X)
            label.append(y)

    data = np.vstack(data)
    label = np.vstack(label)
    names = ('data', 'label')
    data_items = (data, label)

    output_filename = '_'.join([args.output_name, 'training_data.h5'])

    with h5py.File(output_filename, 'w') as h:
        for name, item in zip(names, data_items):
            h.create_dataset(name, data=item)


if __name__ == '__main__':
    main()
