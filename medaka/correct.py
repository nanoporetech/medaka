"""
Correct a reference sequence using a reads-to-reference sequence
alignment and a pre-trained error model.
"""

import sys
import argparse
from functools import partial
from multiprocessing import Pool

from Bio import SeqIO

import keras.backend as K
from keras.models import load_model

from medaka.util.pileup import (bam_to_feature_array, get_reference_length,
                                get_reference_names)
from medaka.util.reshape import trim_to_step_multiple
from medaka.util.generators import serve_data_batch
from medaka.util.sequences import (model_output_to_fasta, assemble_sequences,
                                   segment_limits)


def correct_reference_segment(bam, reference, model, limits=(None, None)):
    """Correct a segment of a reference sequence

    :param bam: (sorted indexed) bam with read alignment to assembly
    :param reference: name of reference to process
    :param model: keras model object
    :param limits: int start and end points of segment
    :returns: SeqRecord corrected sequence
    """
    _, data = bam_to_feature_array(bam, reference, start=limits[0], end=limits[1])
    batch, window, feature = model.input_shape
    data, steps = trim_to_step_multiple(data, batch, window)
    data_generator = serve_data_batch(data, batch, window)
    predictions = model.predict_generator(data_generator, steps)
    segment_name = "%s:%d-%d" % (reference, limits[0], limits[1])
    corrected_seq = model_output_to_fasta(predictions, segment_name)
    return corrected_seq


def multi_correct_reference(bam, reference, model, threads=4):
    """Correct reference segments in parallel
    """
    ref_length = get_reference_length(bam, reference)
    limits = segment_limits(ref_length)
    func = partial(correct_reference_segment, bam, reference, model)
    pool = Pool(threads)
    for p in pool.imap(func, limits):
        yield p
    pool.close()
    pool.join()


def get_parser():
    parser = argparse.ArgumentParser(
        description="""Take an alignment of reads to a reference sequence
                     and produce a corrected reference sequence.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--threads', type=int, default=4,
                        help='number of processes')
    parser.add_argument('bam', help='alignment of reads to reference')
    parser.add_argument('model', help='model .h5 filepath')
    return parser


def main():
    K.set_learning_phase(0)
    sys.setrecursionlimit(2000)

    args = get_parser().parse_args()
    model = load_model(args.model)
    references = get_reference_names(args.bam)

    corrected_sequences = []

    for reference in references:
        for corrected in multi_correct_reference(args.bam, reference, model,
                                            threads=args.threads):
            corrected_sequences.append(corrected)

    corrected_references = assemble_sequences(corrected_sequences)
    SeqIO.write(corrected_references,
        'corrected_reference.fa', 'fasta')


if __name__ =='__main__':
    main()
