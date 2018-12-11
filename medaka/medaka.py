import argparse
import logging
import os
from pkg_resources import resource_filename

import  numpy as np

from medaka.inference import train, predict
from medaka.stitch import stitch
from medaka.features import create_labelled_samples, create_samples

default_model = os.path.join(resource_filename(__package__, 'data'), 'medaka_model.hdf5')

def _log_level():
    """Parser to set logging level and acquire software version/commit"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #parser.add_argument('--version', action='version', version=get_version())

    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument('--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument('--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only).')

    return parser


def _chunking_feature_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('bam', help='Input alignments.')
    parser.add_argument('--model', default=default_model, help='Model definition.')
    parser.add_argument('--batch_size', type=int, default=5, help='Inference batch size.')
    parser.add_argument('--regions', default=None, nargs='+', help='Genomic regions to analyse.')
    parser.add_argument('--chunk_len', type=int, default=10000, help='Chunk length of samples.')
    parser.add_argument('--chunk_ovlp', type=int, default=1000, help='Overlap of chunks.')
    parser.add_argument('--read_fraction', type=float, help='Fraction of reads to keep',
        nargs=2, metavar=('lower', 'upper'))
    parser.add_argument('--rle_ref', default=None, help='Encoded reference file (required only for some model types.')
    parser.add_argument('--max_label_len', type=int, default=1, help='Maximum label length.')
    return parser


def feature_gen_dispatch(args):
    if hasattr(args, 'truth'):
        create_labelled_samples(args)
    else:
        create_samples(args)


def main():
    parser = argparse.ArgumentParser('medaka',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True


    # Transformation of sequence data
    #TODO: is this needed?
    #pparser = subparsers.add_parser('compress',
    #    help='Transform basecalls and draft assemblies.',
    #    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #pparser.set_defaults(func=compress)
    #pparser.add_argument('input', help='.fasta/fastq file.')
    #pparser.add_argument('--output', default=None, help='Output file, default it stdout.')
    #pparser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel execution.')

    # Creation of feature files
    fparser = subparsers.add_parser('features',
        help='Create features for inference.',
        parents=[_log_level(), _chunking_feature_args()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    fparser.set_defaults(func=feature_gen_dispatch)
    fparser.add_argument('output', help='Output features file.')
    fparser.add_argument('--truth', help='Bam of truth aligned to ref to create features for training.')
    fparser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel execution.')

    # Training program
    tparser = subparsers.add_parser('train',
        help='Train a model from features.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    tparser.set_defaults(func=train)
    tparser.add_argument('features', nargs='+', help='Path for training data.')
    tparser.add_argument('--train_name', type=str, default='keras_train', help='Name for training run.')
    tparser.add_argument('--model', help='Model definition and initial weights .hdf, or .yml with kwargs to build model.')
    tparser.add_argument('--max_label_len', type=int, default=1, help='Maximum label length.')
    tparser.add_argument('--epochs', type=int, default=5000, help='Maximum number of trainig epochs.')
    tparser.add_argument('--validation_split', type=float, default=0.2, help='Fraction of data to validate on.')
    tparser.add_argument('--batch_size', type=int, default=200, help='Training batch size.')
    tparser.add_argument('--max_samples', type=int, default=np.inf, help='Only train on max_samples.')
    tparser.add_argument('--mini_epochs', type=int, default=1, help='Reduce fraction of data per epoch by this factor')
    tparser.add_argument('--balanced_weights', action='store_true', help='Balance label weights.')
    tparser.add_argument('--seed', type=int, help='Seed for random batch shuffling.')

    # Consensus from bam input
    cparser = subparsers.add_parser('consensus',
        help='Run inference from a trained model and alignments.',
        parents=[_log_level(), _chunking_feature_args()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cparser.set_defaults(func=predict)
    cparser.add_argument('output', help='Output file.')
    cparser.add_argument('--threads', type=int, default=1, help='Number of threads used by inference.')

    # Consensus from features input
    cfparser = subparsers.add_parser('consensus_from_features',
        help='Run inference from a trained model on existing features.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cfparser.add_argument('features', nargs='+', help='Pregenerated features (from medaka features).')
    cfparser.add_argument('--model', default=default_model, help='Model definition.')
    cfparser.add_argument('--ref_rle', default=None, help='Encoded reference file (required only for some model types.')

    # Post-processing of consensus outputs
    sparser = subparsers.add_parser('stitch',
        help='Stitch together output from medaka consensus into final output.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sparser.set_defaults(func=stitch)
    sparser.add_argument('inputs', nargs='+', help='Consensus .hdf files.')
    sparser.add_argument('output', help='Output .fasta.', default='consensus.fasta')
    sparser.add_argument('--regions', default=None, nargs='+', help='Limit stitching to these reference names')

    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger(__package__)
    logger.setLevel(args.log_level)

    #TODO: do common argument validation here: e.g. rle_ref being present if
    #      required by model
    args.func(args)

    #TODO: subcommand to print extract model / feature yaml data and print to screen / dump to txt yaml file.


if __name__ == '__main__':
    main()
