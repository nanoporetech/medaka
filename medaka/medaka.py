import argparse
import logging
import numpy as np

from medaka.inference import train, predict
from medaka.tview import prepare
from medaka.stitch import stitch


def main():
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    parser = argparse.ArgumentParser('medaka')
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    pparser = subparsers.add_parser('prepare', help='Create (and merge) pileup datasets.')
    pparser.set_defaults(func=prepare)
    pparser.add_argument('ref_fasta', help='.fasta reference corresponding to bams.')
    # pparser.add_argument('bams', nargs='+', help='Input alignments (all aligned to same ref).')
    pparser.add_argument('bam', help='Input alignments (aligned to ref).')
    pparser.add_argument('output', help='Output .hdf.')
    pparser.add_argument('--ref_name', default=None, type=str, help='Name of reference within ref_fasta.')
    pparser.add_argument('--start', default=0, type=int, help='Start reference coordinate.')
    pparser.add_argument('--end', type=int, help='End reference coordinate.')
    pparser.add_argument('--chunk_len', type=int, default=200000, help='Width of pileup chunks (in ref coords) to produce.')
    pparser.add_argument('--truth', default=None, help='Input bam of truth aligned to ref to label data.')


    tparser = subparsers.add_parser('train', help='Train a model from pileup data.')
    tparser.set_defaults(func=train)
    tparser.add_argument('pileupdata', help='Path for training data.')
    tparser.add_argument('--train_name', type=str, default='keras_train', help='Name for training run.')
    tparser.add_argument('--model', help='Model definition and initial weights .hdf.')
    tparser.add_argument('--features', action='store_true', help='Stop after generating features.')
    tparser.add_argument('--max_label_len', type=int, default=1, help='Maximum label length.')
    tparser.add_argument('--epochs', type=int, default=5000, help='Maximum number of trainig epochs.')
    tparser.add_argument('--validation_split', type=float, default=0.2, help='Fraction of data to validate on.')
    tparser.add_argument('--batch_size', type=int, default=200, help='Training batch size.')
    #TODO: plumb these in
    #tparser.add_argument('--min_depth', type=int, default=15, help='Exclude samples in which any pileup column has insufficient depth.')
    #tparser.add_argument('--max_depth', type=int, default=np.inf, help='Exclude samples in which any pileup column has excessive depth.')
    #TODO: add options to include/exclude regions of pileupdata by ref_name/start/end ?

    cparser = subparsers.add_parser('consensus', help='Run inference from a trained model.')
    cparser.set_defaults(func=predict)
    cparser.add_argument('model', help='Model .hdf file from training.')
    cparser.add_argument('--encoding', help='Model label encoding .json, used only if encoding not in .hdf model')
    cparser.add_argument('--output_fasta', default='consensus_chunks.fa', help='Consensus output file.')
    cparser.add_argument('--output_features', default=None, help='Save features for later reuse.')
    cparser.add_argument('--output_probs', default=None, help='Label probabilities output file.')
    cparser.add_argument('--start', default=0, type=int, help='Reference position at which to start, only used with --alignments.')
    cparser.add_argument('--end', default=None, type=int, help='Reference position at which to end, only used with --alignments.')
    cparser.add_argument('--overlap', default=10000, type=int, help='Overlap of neighbouring samples in pileup columns.')
    cparser.add_argument('--sample_length', type=int, default=50000, help='Number of pileup columns in each consensus sample.')
    ingroup = cparser.add_mutually_exclusive_group(required=True)
    ingroup.add_argument('--pileup', help='Pileup input data (saved with prepare)')
    ingroup.add_argument('--features', help='Pregenerated features (saved with --output_features option).')
    ingroup.add_argument('--alignments', help='Input alignments, reference fasta and reference name (within fasta).', nargs=3,
                         metavar=('reads.bam', 'ref.fasta', 'ref_name'))
    cparser.add_argument('--batch_size', type=int, default=100, help='Prediction batch size.')

    sparser = subparsers.add_parser('stitch', help='Stitch basecalled chunks together.')
    sparser.set_defaults(func=stitch)
    sparser.add_argument('fastas', nargs='+', help='.fasta files containing chunks to stitch.')
    sparser.add_argument('output', help='output .fasta.', default='stitched_basecalls.fasta')
    sparser.add_argument('--min_overlap', help='minimum overlap required by stitching.', default=50, type=int)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
