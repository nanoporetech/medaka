import argparse
import logging

from medaka.inference import train, predict
from medaka.tview import prepare


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('medaka')
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    pparser = subparsers.add_parser('prepare', help='Create (and merge) pileup datasets.')
    pparser.set_defaults(func=prepare)
    pparser.add_argument('ref_fasta', help='.fasta reference corresponding to bams.')
    pparser.add_argument('ref_name', help='Name of reference within ref_fasta.')
    pparser.add_argument('output', help='Output .hdf.')
    # pparser.add_argument('bams', nargs='+', help='Input alignments (all aligned to same ref).')
    pparser.add_argument('bam', help='Input alignments (aligned to ref).')
    pparser.add_argument('--start', default=0, type=int, help='Start reference coordinate.')
    pparser.add_argument('--end', type=int, help='End reference coordinate.')
    pparser.add_argument('--chunk_len', type=int, default=20000, help='Width of pileup chunks (in ref coords) to produce.')
    pparser.add_argument('--truth', default=None, help='Input bam of truth aligned to ref to label data.')


    tparser = subparsers.add_parser('train', help='Train a model from pileup data.')
    tparser.set_defaults(func=train)
    tparser.add_argument('pileupdata', help='Path for training data.')
    tparser.add_argument('--train_name', type=str, default='keras_train', help='Name for training run.')
    tparser.add_argument('--model', help='Model definition and initial weights .hdf.')
    tparser.add_argument('--features', action='store_true', help='Stop after generating features.')

    cparser = subparsers.add_parser('consensus', help='Run inference from a trained model.')
    cparser.set_defaults(func=predict)
    cparser.add_argument('model', help='Model .hdf file from training.')
    cparser.add_argument('--encoding', help='Model label encoding .json, used only if encoding not in .hdf model')
    cparser.add_argument('--output_fasta', default='basecalls.fasta', help='Polished consensus output file.')
    cparser.add_argument('--start', default=0, type=int, help='Reference position at which to start, only used with --alignments.')
    cparser.add_argument('--end', default=None, type=int, help='Reference position at which to end, only used with --alignments.')
    cparser.add_argument('--overlap', default=250, type=int, help='Overlap of neighbouring polished chunks in pileup columns.')
    ingroup = cparser.add_mutually_exclusive_group(required=True)
    ingroup.add_argument('--pileupdata', help='Pileup input data.')
    ingroup.add_argument('--feature_file', help='Pregenerated features as stored during training.')
    ingroup.add_argument('--alignments', help='Input alignments, reference fasta and reference name (within fasta).', nargs=3,
                         metavar=('reads.bam', 'ref.fasta', 'ref_name'))

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
