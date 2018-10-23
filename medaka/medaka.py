import argparse
import logging
import numpy as np
import yaml

from medaka.inference import train, predict, prepare
from medaka.stitch import stitch
from medaka.common import write_yaml_data


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('medaka')
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    pparser = subparsers.add_parser('prepare', help='Create pileup feature datasets and write to hdf.')
    pparser.set_defaults(func=prepare)
    pparser.add_argument('ref_fasta', help='.fasta reference corresponding to bams.')
    # pparser.add_argument('bams', nargs='+', help='Input alignments (all aligned to same ref).')
    pparser.add_argument('bam', help='Input alignments (aligned to ref).')
    pparser.add_argument('output', help='Output .hdf.')
    pparser.add_argument('--ref_name', default=None, type=str, help='Name of reference within ref_fasta.')
    pparser.add_argument('--start', default=0, type=int, help='Start reference coordinate.')
    pparser.add_argument('--end', type=int, help='End reference coordinate.')
    pparser.add_argument('--sample_length', type=int, default=1000, help='Number of pileup columns in each feature sample.')
    pparser.add_argument('--overlap', default=100, type=int, help='Overlap of neighbouring samples in pileup columns.')
    pparser.add_argument('--truth', default=None, help='Input bam of truth aligned to ref to label data.')


    tparser = subparsers.add_parser('train', help='Train a model from features.')
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

    cparser = subparsers.add_parser('consensus', help='Run inference from a trained model.')
    cparser.set_defaults(func=predict)
    cparser.add_argument('model', help='Model .hdf file from training.')
    cparser.add_argument('--model_yml', help='Model yml containing label encoding and model options, required only if training ended prematurely.')
    cparser.add_argument('--output_fasta', default='consensus_chunks.fa', help='Consensus sequence output fasta file.')
    cparser.add_argument('--output_probs', default=None, help='Consensus probabilities output hdf file.')
    cparser.add_argument('--start', default=0, type=int, help='Reference position at which to start, only used with --alignments.')
    cparser.add_argument('--end', default=None, type=int, help='Reference position at which to end, only used with --alignments.')
    cparser.add_argument('--overlap', default=10000, type=int, help='Overlap of neighbouring samples in pileup columns.')
    cparser.add_argument('--sample_length', type=int, default=50000, help='Number of pileup columns in each consensus sample.')
    ingroup = cparser.add_mutually_exclusive_group(required=True)
    ingroup.add_argument('--pileup', help='Pileup input data (saved with prepare)')
    ingroup.add_argument('--features', nargs='+', help='Pregenerated features (saved with --output_features option).')
    ingroup.add_argument('--alignments', help='Input alignments, reference fasta and reference name (within fasta).', nargs=3,
                         metavar=('reads.bam', 'ref.fasta', 'ref_name'))
    cparser.add_argument('--batch_size', type=int, default=5, help='Prediction batch size.')
    cparser.add_argument('--ref_names', nargs='+', help='Only load these references from pregenerated features.', default=None)

    sparser = subparsers.add_parser('stitch', help='Stitch basecalled chunks together.')
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
