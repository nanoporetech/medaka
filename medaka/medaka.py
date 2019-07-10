import argparse
import logging
import os
from pkg_resources import resource_filename
import yaml

import numpy as np
import pysam

import medaka.datastore
import medaka.features
import medaka.inference
import medaka.labels
import medaka.stitch
import medaka.variant
import medaka.vcf

model_store = resource_filename(__package__, 'data')
allowed_models = [
    'r941_trans', 'r941_flip213', 'r941_flip235',
    'r941_min_fast', 'r941_min_high', 'r941_prom_fast', 'r941_prom_high',
]
default_model = 'r941_min_high'
model_dict = {
    k:os.path.join(model_store, '{}_model.hdf5'.format(k))
    for k in allowed_models
}


class ResolveModel(argparse.Action):
    """Resolve model filename or ID into filename"""
    def __init__(self, option_strings, dest, default=None, required=False, help='Model file.'):
        super().__init__(
            option_strings, dest, nargs=1, default=default, required=required,
            help='{} {{{}}}'.format(help, ', '.join(allowed_models))
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


class CheckBam(argparse.Action):
    """Check a bam has < 2 samples (RG tag)"""

    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.exists(values):
            raise RuntimeError(
                "Filepath for '--{}' argument does not exist ({})".format(
                    self.dest, values)
            )
        with pysam.AlignmentFile(values) as bam:
            header_dict = bam.header.as_dict()
            if 'RG' in header_dict and len(header_dict['RG']) > 1:
                raise RuntimeError('The bam {} contains more than one read group.'.format(values))
        setattr(namespace, self.dest, values)


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


def _model_arg():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--model', action=ResolveModel, default=model_dict[default_model],
                        help='Model definition, default is equivalent to {}.'.format(default_model))
    return parser


def _chunking_feature_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False,
        parents=[_model_arg()],
    )
    parser.add_argument('bam', help='Input alignments.', action=CheckBam)
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
        medaka.features.create_labelled_samples(args)
    else:
        medaka.features.create_samples(args)


def hdf2yaml(args):
    with medaka.datastore.DataStore(args.input) as ds, open(args.output, 'w') as fh:
        yaml.dump(ds.meta, fh)


def yaml2hdf(args):
    with medaka.datastore.DataStore(args.output, 'a') as ds, open(args.input) as fh:
        ds.update_meta(yaml.unsafe_load(fh))


def print_model_path(args):
    print(os.path.abspath(args.model))


def print_all_models(args):
    print('Available:', ', '.join(allowed_models))
    print('Default:', default_model)


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
    fparser.add_argument('--truth_haplotag', help='Two-letter tag defining haplotype of alignments for polyploidy labels.')
    fparser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel execution.')

    # Training program
    tparser = subparsers.add_parser('train',
        help='Train a model from features.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    tparser.set_defaults(func=medaka.inference.train)
    tparser.add_argument('features', nargs='+', help='Paths to training data.')
    tparser.add_argument('--train_name', type=str, default='keras_train', help='Name for training run.')
    tparser.add_argument('--model', action=ResolveModel, help='Model definition and initial weights .hdf, or .yml with kwargs to build model.')
    tparser.add_argument('--label_scheme', default='HaploidLabelScheme', help='Labelling scheme.',
                         choices=sorted(medaka.labels.label_schemes.keys()))
    tparser.add_argument('--max_label_len', type=int, default=1, help='Maximum label length.')
    tparser.add_argument('--epochs', type=int, default=5000, help='Maximum number of trainig epochs.')
    tparser.add_argument('--batch_size', type=int, default=200, help='Training batch size.')
    tparser.add_argument('--max_samples', type=int, default=np.inf, help='Only train on max_samples.')
    tparser.add_argument('--mini_epochs', type=int, default=1, help='Reduce fraction of data per epoch by this factor')
    tparser.add_argument('--balanced_weights', action='store_true', help='Balance label weights.')
    tparser.add_argument('--multi_label', action='store_true', help='Multi-classification training.')
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
    cparser.set_defaults(func=medaka.inference.predict)
    cparser.add_argument('output', help='Output file.')
    cparser.add_argument('--threads', type=int, default=1, help='Number of threads used by inference.')
    cparser.add_argument('--check_output', action='store_true', default=False,
            help='Verify integrity of output file after inference.')
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
    sparser.set_defaults(func=medaka.stitch.stitch)
    sparser.add_argument('inputs', nargs='+', help='Consensus .hdf files.')
    sparser.add_argument('output', help='Output .fasta.', default='consensus.fasta')
    sparser.add_argument('--regions', default=None, nargs='+', help='Limit stitching to these reference names')

    pparser = subparsers.add_parser('variant',
        help='Decode probabilities to VCF.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pparser.set_defaults(func=medaka.variant.variants_from_hdf)
    pparser.add_argument('ref_fasta', help='Reference sequence .fasta file.')
    pparser.add_argument('inputs', nargs='+', help='Consensus .hdf files.')
    pparser.add_argument('output', help='Output .vcf.', default='medaka.vcf')
    pparser.add_argument('--regions', default=None, nargs='+', help='Limit variant calling to these reference names')
    pparser.add_argument('--threshold', default=0.04, type=float,
                         help="""Threshold for considering secondary calls. Only used in SNPDecoder.
                         A value of 1 will result in haploid decoding.""")
    pparser.add_argument('--ref_vcf', default=None, help='Reference vcf to compare to, only used in SNPDecoder.')
    pparser.add_argument('--decoder', default='SNPDecoder', help='Variant decoder.',
                         choices=sorted(medaka.variant.variant_decoders.keys()))
    pparser.add_argument('--multi_label', action='store_true', help='Multi-classification decoding.')


    # Tools
    toolparser = subparsers.add_parser('tools',
        help='tools sub-command.',
        parents=[_log_level()],
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
    h2dparser = toolsubparsers.add_parser('haploid2diploid',
        help='Merge two haploid VCFs into a diploid VCF.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    h2dparser.set_defaults(func=medaka.vcf.haploid2diploid)
    h2dparser.add_argument('vcf1', help='Input .vcf file.')
    h2dparser.add_argument('vcf2', help='Input .vcf file.')
    h2dparser.add_argument('ref_fasta', help='Reference .fasta file.')
    h2dparser.add_argument('vcfout', help='Output .vcf.')
    h2dparser.add_argument('--adjacent', action='store_true',
                         help=('Merge adjacent variants (i.e. variants with non-overlapping genomic ranges) instead' +
                               ' of just overlapping ones. If set to True, all runs of adjacent variants will be' +
                               ' merged including those which appear in just one of the input VCFs.'))

    # split a diploid VCF into two haploid VCFs.
    d2hparser = toolsubparsers.add_parser('diploid2haploid',
        help='Split a diploid VCF into two haploid VCFs.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    d2hparser.set_defaults(func=medaka.vcf.diploid2haploid)
    d2hparser.add_argument('vcf', help='Input .vcf file.')
    d2hparser.add_argument('--notrim', action='store_true',
                         help='Do not trim variant to minimal ref and alt and update pos.')

    # classify variants in a vcf writing new vcfs containing subs, indels etc.
    clparser = toolsubparsers.add_parser('classify_variants',
        help='Classify variants and write vcf for each type and subtype.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    clparser.set_defaults(func=medaka.vcf.classify_variants)
    clparser.add_argument('vcf', help='Input .vcf file.')
    clparser.add_argument('--replace_info', action='store_true',
                         help='Replace info tag (useful for visual inspection of types).')

    # convert a vcf to tsv
    tsvparser = toolsubparsers.add_parser('vcf2tsv',
        help='convert vcf to tsv, unpacking info and sample columns.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    tsvparser.set_defaults(func=medaka.vcf.vcf2tsv)
    tsvparser.add_argument('vcf', help='Input .vcf file.')


    # find homozygous and heterozygous regions in a VCF
    hzregparser = toolsubparsers.add_parser('homozygous_regions',
        help='Find homozygous regions from a diploid vcf.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    hzregparser.set_defaults(func=medaka.vcf.get_homozygous_regions)
    hzregparser.add_argument('vcf', help='Input .vcf file.')
    hzregparser.add_argument('region', help='Genomic region within which to find homozygous sub-regions.')
    hzregparser.add_argument('--min_len', type=int, default=1000,
                             help='Minimum region length.')
    hzregparser.add_argument('--suffix', help='Output suffix.', default='regions.txt')

    # resolve model and print full model file path to screen
    rmparser = toolsubparsers.add_parser('resolve_model',
        help='Resolve model and print full file path.',
        parents=[_log_level(), _model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rmparser.set_defaults(func=print_model_path)

    # print all model tags followed by default
    lmparser = toolsubparsers.add_parser('list_models',
        help='List all models.',
        parents=[_log_level(), _model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    lmparser.set_defaults(func=print_all_models)

    args = parser.parse_args()

    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger(__package__)
    logger.setLevel(args.log_level)

    if args.command == 'tools' and not hasattr(args, 'func'):
        # display help if given `medaka tools (--help)`
        toolparser.print_help()
    else:
        #TODO: do common argument validation here: e.g. rle_ref being present if
        #      required by model
        args.func(args)


if __name__ == '__main__':
    main()
