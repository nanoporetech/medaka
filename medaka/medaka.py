import argparse
import logging
import numpy as np
import yaml

from medaka.inference import train, predict
from medaka.stitch import stitch
from medaka.common import write_yaml_data
from medaka.features import choose_feature_func, compress


def _consensus_feat():
    """Parser with arguments common to generating consensus features from alignments."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Consensus features arguments.',
        add_help=False)

    parser.set_defaults(func=choose_feature_func)
    parser.add_argument('bam', help='Input alignments.')
    parser.add_argument('reference', help='Input reference .fasta corresponding to `bam`.')
    parser.add_argument('-r', '--regions', nargs='+', default=None, type=str, help='Regions in samtools format.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for parallel execution.')
    parser.add_argument('--chunk_len', type=int, default=10000, help='Chunk length of samples.')
    parser.add_argument('--chunk_ovlp', type=int, default=1000, help='Overlap of chunks.')
    parser.add_argument('--read_fraction', type=float, help='Fraction of reads to keep',
        nargs=2, metavar=('lower', 'upper'))
    return parser


def _train_feat():
    """Parser with arguments for generating training features from alignments."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Consensus features arguments.',
        add_help=False)

    parser.add_argument('output', help='Output features file.')
    parser.add_argument('-T', '--truth', help='Bam of truth aligned to ref to create features for training.')
    parser.add_argument('--batch_size', type=int, default=None, help='Write training batches rather than samples')
    parser.add_argument('--max_label_len', type=int, default=10, help='Max label length (only used if writing batches).')
    parser.add_argument('-m', '--model', help='Create features expected by this model (.yml or .hdf)')
    return parser


def _predict():
    """Parser with arguments for running consensus."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Consensus.',
        add_help=False)

    parser.set_defaults(func=predict)
    parser.add_argument('model', help='Model .hdf file from training.')
    parser.add_argument('--model_yml', help='Model yml containing label encoding and model options, required only if training ended prematurely.')
    parser.add_argument('--output_prefix', default='consensus_chunks', help='Consensus sequence output fasta file.')
    parser.add_argument('--output_probs', default=None, help='Consensus probabilities output hdf file.')
    parser.add_argument('--batch_size', type=int, default=5, help='Prediction batch size.')

    return parser


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('medaka', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True


    pparser = subparsers.add_parser('compress', help='Compress basecalls / draft assemblies.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pparser.set_defaults(func=compress)
    pparser.add_argument('input', help='.fasta/fastq file.')
    pparser.add_argument('-o', '--output', default=None, help='Output file, default it stdout.')
    pparser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for parallel execution.')

    fparser = subparsers.add_parser('features', help='Create features for inference.',
                                    parents=[_consensus_feat(), _train_feat()],
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    tparser = subparsers.add_parser('train', help='Train a model from features.',
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
    tparser.add_argument('--balanced_weights', action='store_true', help='Balance label weights.')

    cparser = subparsers.add_parser('consensus', help='Run inference from a trained model and alignments.',
                                    parents=[_consensus_feat(), _predict()],
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    cfparser = subparsers.add_parser('consensus_from_features', help='Run inference from a trained model on existing features.',
                                    parents=[_predict()],
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cfparser.add_argument('--features', nargs='+', help='Pregenerated features (saved with --output_features option).')
    cfparser.add_argument('--ref_names', nargs='+', help='Only load these references from pregenerated features.', default=None)


    sparser = subparsers.add_parser('stitch', help='Stitch basecalled chunks together.',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sparser.set_defaults(func=stitch)
    sparser.add_argument('inputs', nargs='+', help='Consensus sequence (.fasta) or probabilities (.hdf) files.')
    sparser.add_argument('output', help='Output .fasta.', default='consensus.fasta')
    sparser.add_argument('--ref_names', nargs='+', help='Limit stitching to these reference names', default=None)
    sparser.add_argument('--model_yml', help='Model yml containing label encoding, required only if consensus ended prematurely.')

    fparser = subparsers.add_parser('fix', help='Add data to a hdf file.')
    fparser.set_defaults(func=fix)
    fparser.add_argument('hdfs', nargs='+', help='hdf files to add data to.')
    fparser.add_argument('yml', help='yml file containing data to add.')

    args = parser.parse_args()
    args.func(args)

def fix(args):
    with open(args.yml) as fh:
        yml_str = fh.read()
    d = yaml.load(yml_str)
    for hdf in args.hdfs:
        write_yaml_data(hdf, d)


if __name__ == '__main__':
    main()
