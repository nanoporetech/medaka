import argparse
import logging
import os
from pkg_resources import resource_filename
import yaml

import  numpy as np

from medaka.datastore import DataStore
from medaka.inference import train, predict
from medaka.stitch import stitch, snps, merge_vcfs
from medaka.features import create_labelled_samples, create_samples

model_store = resource_filename(__package__, 'data')

model_dict = {
  'r94': 'medaka_model.hdf5',
  'r941_flip': 'r941_flip_model.hdf5'
}
model_dict = {k:os.path.join(model_store, v) for k,v in model_dict.items()}
default_model = 'r94'


class ResolveModel(argparse.Action):
    """Resolve model filename or ID into filename"""
    def __init__(self, option_strings, dest, default=None, required=False, help='Model file.'):
        super().__init__(
            option_strings, dest, nargs=1, default=default, required=required,
            help='{} {{{}}}'.format(help, ', '.join(model_dict.keys()))
        )

    def __call__(self, parser, namespace, values, option_string=None):
        val = values[0]
        if not os.path.exists(val):
            # try lookup
            try:
                val = model_dict[val]
            except:
                raise RuntimeError(
                    "Filepath for '--{}' argument does not exist and is not a known model ID ({})".format(
                        self.dest, val)
                )
            #TODO: verify the file is a model?
        setattr(namespace, self.dest, val)



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
    parser.add_argument('--model', action=ResolveModel, default=model_dict[default_model],
                        help='Model definition, default is equivalent to {}.'.format(default_model))
    parser.add_argument('--batch_size', type=int, default=200, help='Inference batch size.')
    parser.add_argument('--regions', default=None, nargs='+', help='Genomic regions to analyse.')
    parser.add_argument('--chunk_len', type=int, default=10000, help='Chunk length of samples.')
    parser.add_argument('--chunk_ovlp', type=int, default=1000, help='Overlap of chunks.')
    parser.add_argument('--read_fraction', type=float, help='Fraction of reads to keep',
        nargs=2, metavar=('lower', 'upper'))
    parser.add_argument('--rle_ref', default=None, help='Encoded reference file (required only for some model types.')
    return parser


def feature_gen_dispatch(args):
    if hasattr(args, 'truth'):
        create_labelled_samples(args)
    else:
        create_samples(args)


def hdf2yaml(args):
    with DataStore(args.input) as ds, open(args.output, 'w') as fh:
        yaml.dump(ds.meta, fh)


def yaml2hdf(args):
    with DataStore(args.output, 'a') as ds, open(args.input) as fh:
        ds.update_meta(yaml.load(fh))


def main():
    from medaka import __version__
    parser = argparse.ArgumentParser('medaka',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    parser.add_argument('--version', action='version',
        version='%(prog)s {}'.format(__version__))

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
    tparser.add_argument('features', nargs='+', help='Paths to training data.')
    tparser.add_argument('--train_name', type=str, default='keras_train', help='Name for training run.')
    tparser.add_argument('--model', action=ResolveModel, help='Model definition and initial weights .hdf, or .yml with kwargs to build model.')
    tparser.add_argument('--max_label_len', type=int, default=1, help='Maximum label length.')
    tparser.add_argument('--epochs', type=int, default=5000, help='Maximum number of trainig epochs.')
    tparser.add_argument('--batch_size', type=int, default=200, help='Training batch size.')
    tparser.add_argument('--max_samples', type=int, default=np.inf, help='Only train on max_samples.')
    tparser.add_argument('--mini_epochs', type=int, default=1, help='Reduce fraction of data per epoch by this factor')
    tparser.add_argument('--balanced_weights', action='store_true', help='Balance label weights.')
    tparser.add_argument('--seed', type=int, help='Seed for random batch shuffling.')
    tparser.add_argument('--threads_io', type=int, default=1, help='Number of threads for parallel IO.')
    tparser.add_argument('--device', type=int, default=0, help='GPU device to use.')

    vgrp = tparser.add_mutually_exclusive_group()
    vgrp.add_argument('--validation_split', type=float, default=0.2, help='Fraction of data to validate on.')
    vgrp.add_argument('--validation_features', nargs='+', default=None, help='Paths to validation data')

    # Consensus from bam input
    cparser = subparsers.add_parser('consensus',
        help='Run inference from a trained model and alignments.',
        parents=[_log_level(), _chunking_feature_args()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cparser.set_defaults(func=predict)
    cparser.add_argument('output', help='Output file.')
    cparser.add_argument('--threads', type=int, default=1, help='Number of threads used by inference.')
    cparser.add_argument('--save_features', action='store_true', default=False,
                         help='Save features with consensus probabilities.')
    tag_group = cparser.add_argument_group('filter tag', 'Filtering alignments by an integer valued tag.')
    tag_group.add_argument('--tag_name', type=str, help='Two-letter tag name.')
    tag_group.add_argument('--tag_value', type=int, help='Value of tag.')
    tag_group.add_argument('--tag_keep_missing', action='store_true', help='Keep alignments when tag is missing.')

    # Consensus from features input
    cfparser = subparsers.add_parser('consensus_from_features',
        help='Run inference from a trained model on existing features.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cfparser.add_argument('features', nargs='+', help='Pregenerated features (from medaka features).')
    cfparser.add_argument('--model', action=ResolveModel, default=default_model, help='Model definition.')
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

    pparser = subparsers.add_parser('snp',
        help='Decode probabilities as dipoloid SNPs.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pparser.set_defaults(func=snps)
    pparser.add_argument('ref_fasta', help='Reference sequence .fasta file.')
    pparser.add_argument('inputs', nargs='+', help='Consensus .hdf files.')
    pparser.add_argument('output', help='Output .vcf.', default='snps.vcf')
    pparser.add_argument('--regions', default=None, nargs='+', help='Limit stitching to these reference names')
    pparser.add_argument('--threshold', default=0.04, type=float, help='Threshold for considering secondary calls.')
    pparser.add_argument('--ref_vcf', default=None, help='Reference vcf to compare to.')

    # Tools
    toolparser = subparsers.add_parser('tools',
        help='tools sub-command.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    toolsubparsers = toolparser.add_subparsers(title='tools', description='valid tool commands', help='additional help', dest='tool_command')

    # Dump model/feature meta to yaml
    hparser = toolsubparsers.add_parser('hdf2yaml',
        help='Dump medaka meta in a hdf to yaml.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    hparser.set_defaults(func=hdf2yaml)
    hparser.add_argument('input', help='Input .hdf file.')
    hparser.add_argument('output', help='Output .yaml file.', default='meta.yaml')

    # Create model .hdf containing model/feature meta from yaml
    yparser = toolsubparsers.add_parser('yaml2hdf',
        help='Dump medaka meta in a yaml to hdf.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    yparser.set_defaults(func=yaml2hdf)
    yparser.add_argument('input', help='Input .yaml file.')
    yparser.add_argument('output', help='Output .hdf, will be appended to if it exists.', default='meta.hdf')

    # merge two haploid VCFs into a diploid VCF.
    yparser = toolsubparsers.add_parser('merge_vcfs',
        help='Merge two haploid VCFs into a phased diploid VCF.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    yparser.set_defaults(func=merge_vcfs)
    yparser.add_argument('vcf1', help='Input .vcf file.')
    yparser.add_argument('vcf2', help='Input .vcf file.')
    yparser.add_argument('vcfout', help='Output .vcf.')


    # Create model .hdf containing model/feature meta from yaml
    yparser = toolsubparsers.add_parser('yaml2hdf',
        help='Dump medaka meta in a yaml to hdf.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    yparser.set_defaults(func=yaml2hdf)
    yparser.add_argument('input', help='Input .yaml file.')
    yparser.add_argument('output', help='Output .hdf, will be appended to if it exists.', default='meta.hdf')

    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger(__package__)
    logger.setLevel(args.log_level)

    #TODO: do common argument validation here: e.g. rle_ref being present if
    #      required by model
    args.func(args)


if __name__ == '__main__':
    main()
