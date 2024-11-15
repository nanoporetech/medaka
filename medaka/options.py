"""Stores miscellaneous data for medaka."""

# note this module is imported into setup.py, so do not use any
# non-stdlib packages
import os
import pathlib

import pkg_resources

# the models used by default for CLI entry points
default_models = {
    'consensus': 'r1041_e82_400bps_sup_v5.0.0',
    'variant': 'r1041_e82_400bps_sup_variant_v5.0.0'
}

# current models are those included in PyPI packages
current_models = [
    # r1041 e82 (kit14) consensus
    'r1041_e82_400bps_hac_v5.0.0',
    'r1041_e82_400bps_sup_v5.0.0',
    # r1041 variant calling
    'r1041_e82_400bps_hac_variant_v5.0.0',
    'r1041_e82_400bps_sup_variant_v5.0.0',
]

# mapping from basecaller model names to medaka models.
# name: (consensus, variant)
# models here are automatically added to full archived list
basecaller_models = {
    # R10.3
    'dna_r10.3_450bps_hac':
        ('r103_hac_g507', 'r103_hac_variant_g507'),
    'dna_r10.3_450bps_hac_prom':
    ('r103_hac_g507', 'r103_hac_variant_g507'),
    # R10.4.1 260bps
    'dna_r10.4.1_e8.2_260bps_hac':
    ('r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_hac_variant_g632'),
    'dna_r10.4.1_e8.2_260bps_hac@v4.0.0':
        ('r1041_e82_260bps_hac_v4.0.0', None),
    'dna_r10.4.1_e8.2_260bps_hac@v4.1.0':
        ('r1041_e82_260bps_hac_v4.1.0', 'r1041_e82_260bps_hac_variant_v4.1.0'),
    'dna_r10.4.1_e8.2_260bps_hac_prom':
        ('r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_hac_variant_g632'),
    'dna_r10.4.1_e8.2_260bps_sup@v4.0.0':
        ('r1041_e82_260bps_sup_v4.0.0', None),
    'dna_r10.4.1_e8.2_260bps_sup@v4.1.0':
        ('r1041_e82_260bps_sup_v4.1.0', 'r1041_e82_260bps_sup_variant_v4.1.0'),
    # R10.4.1 400bps
    'dna_r10.4.1_e8.2_400bps_hac':
        ('r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_hac_variant_g632'),
    'dna_r10.4.1_e8.2_400bps_hac@v3.5.1':
        ('r1041_e82_400bps_hac_g615', 'r1041_e82_400bps_hac_variant_g615'),
    'dna_r10.4.1_e8.2_400bps_hac@v3.5.2':
        ('r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_hac_variant_g632'),
    'dna_r10.4.1_e8.2_400bps_hac@v4.0.0':
        ('r1041_e82_400bps_hac_v4.0.0', None),
    'dna_r10.4.1_e8.2_400bps_hac@v4.1.0':
        ('r1041_e82_400bps_hac_v4.1.0', 'r1041_e82_400bps_hac_variant_v4.1.0'),
    'dna_r10.4.1_e8.2_400bps_hac@v4.2.0':
        ('r1041_e82_400bps_hac_v4.2.0', 'r1041_e82_400bps_hac_variant_v4.2.0'),
    'dna_r10.4.1_e8.2_5khz_400bps_hac@v4.2.0':  # alias for v4.2.0
        ('r1041_e82_400bps_hac_v4.2.0', 'r1041_e82_400bps_hac_variant_v4.2.0'),
    'dna_r10.4.1_e8.2_400bps_hac@v4.3.0':
        ('r1041_e82_400bps_hac_v4.3.0', 'r1041_e82_400bps_hac_variant_v4.3.0'),
    'dna_r10.4.1_e8.2_400bps_hac@v5.0.0':
        ('r1041_e82_400bps_hac_v5.0.0', 'r1041_e82_400bps_hac_variant_v5.0.0'),
    'dna_r10.4.1_e8.2_400bps_hac_prom':
        ('r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_hac_variant_g632'),
    'dna_r10.4.1_e8.2_400bps_sup@v3.5.2':
        ('r1041_e82_400bps_sup_g615', 'r1041_e82_400bps_sup_variant_g615'),
    'dna_r10.4.1_e8.2_400bps_sup@v3.5.1':
        ('r1041_e82_400bps_sup_g615', 'r1041_e82_400bps_sup_variant_g615'),
    'dna_r10.4.1_e8.2_400bps_sup@v4.0.0':
        ('r1041_e82_400bps_sup_v4.0.0', None),
    'dna_r10.4.1_e8.2_400bps_sup@v4.1.0':
        ('r1041_e82_400bps_sup_v4.1.0', 'r1041_e82_400bps_sup_variant_v4.1.0'),
    'dna_r10.4.1_e8.2_400bps_sup@v4.2.0':
        ('r1041_e82_400bps_sup_v4.2.0', 'r1041_e82_400bps_sup_variant_v4.2.0'),
    'dna_r10.4.1_e8.2_5khz_400bps_sup@v4.2.0':  # alias for v4.2.0
        ('r1041_e82_400bps_sup_v4.2.0', 'r1041_e82_400bps_sup_variant_v4.2.0'),
    'dna_r10.4.1_e8.2_400bps_sup@v4.3.0':
        ('r1041_e82_400bps_sup_v4.3.0', 'r1041_e82_400bps_sup_variant_v4.3.0'),
    'dna_r10.4.1_e8.2_400bps_sup@v5.0.0':
        ('r1041_e82_400bps_sup_v5.0.0', 'r1041_e82_400bps_sup_variant_v5.0.0'),
    # R9.4.1 This is a little dodgy
    # note: 'dna_r9.4.1_450bps_hac' is not present here as there is not a
    #       injective mapping.
    'dna_r9.4.1_e8_fast@v3.4':
        ('r941_min_fast_g507', 'r941_min_fast_variant_g507'),
    'dna_r9.4.1_e8_hac@v3.3':
        ('r941_min_hac_g507', 'r941_min_hac_variant_g507'),
    'dna_r9.4.1_e8_sup@v3.3':
        ('r941_min_sup_g507', 'r941_min_sup_variant_g507'),
    # APK joint polishing
    # Not adding this since simultaneously matching on multiple tags is not
    # implemented
    # 'dna_r10.4.1_e8.2_260bps_apk_sup@v5.0.0':
    #     ('r1041_e82_260bps_joint_apk_ulk_v5.0.0', None),
}

# archived models are not included in packages but can be downloaded on the fly
# this list SHOULD NOT be added to, any new models should go into the above
# basecaller_models structure
archived_models = [
    # r9 consensus
    'r941_sup_plant_g610',
    'r941_min_fast_g507', 'r941_prom_fast_g507',
    # r10 consensus
    'r103_fast_g507', 'r103_hac_g507', 'r103_sup_g507',
    # r104 e81 consensus
    'r104_e81_fast_g5015', 'r104_e81_sup_g5015', 'r104_e81_hac_g5015',
    'r104_e81_sup_g610',
    # r104 e81 variant calling
    'r104_e81_fast_variant_g5015', 'r104_e81_hac_variant_g5015',
    'r104_e81_sup_variant_g610',
    # r1041 e82 consensus
    'r1041_e82_400bps_hac_g615',  'r1041_e82_400bps_fast_g615',
    'r1041_e82_400bps_fast_g632', 'r1041_e82_260bps_fast_g632',
    'r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_sup_g615',
    'r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_sup_g632',
    'r1041_e82_400bps_hac_v4.0.0', 'r1041_e82_400bps_sup_v4.0.0',
    'r1041_e82_260bps_hac_v4.0.0', 'r1041_e82_260bps_sup_v4.0.0',
    'r1041_e82_260bps_hac_v4.1.0', 'r1041_e82_260bps_sup_v4.1.0',
    'r1041_e82_400bps_hac_v4.1.0', 'r1041_e82_400bps_sup_v4.1.0',
    'r1041_e82_400bps_hac_v4.2.0', 'r1041_e82_400bps_sup_v4.2.0',
    'r1041_e82_400bps_hac_v4.3.0', 'r1041_e82_400bps_sup_v4.3.0',
    # r1041 e82 variant calling
    'r1041_e82_400bps_hac_variant_g615',
    'r1041_e82_400bps_fast_variant_g615',
    'r1041_e82_400bps_fast_variant_g632',
    'r1041_e82_260bps_fast_variant_g632',
    'r1041_e82_400bps_hac_variant_g632',
    'r1041_e82_400bps_sup_variant_g615',
    'r1041_e82_260bps_hac_variant_g632',
    'r1041_e82_260bps_sup_variant_g632',
    'r1041_e82_260bps_hac_variant_v4.1.0',
    'r1041_e82_260bps_sup_variant_v4.1.0',
    'r1041_e82_400bps_hac_variant_v4.1.0',
    'r1041_e82_400bps_sup_variant_v4.1.0',
    'r1041_e82_400bps_hac_variant_v4.2.0',
    'r1041_e82_400bps_sup_variant_v4.2.0',
    'r1041_e82_400bps_hac_variant_v4.3.0',
    'r1041_e82_400bps_sup_variant_v4.3.0',
    # snp and variant - crf guppy507+
    'r941_sup_plant_variant_g610',
    'r941_min_fast_snp_g507', 'r941_min_fast_variant_g507',
    'r941_min_hac_snp_g507',
    'r941_min_sup_snp_g507', 'r941_min_sup_variant_g507',
    'r941_prom_fast_snp_g507', 'r941_prom_fast_variant_g507',
    'r941_prom_hac_snp_g507',
    'r941_prom_sup_snp_g507', 'r941_prom_sup_variant_g507',
    'r103_fast_snp_g507', 'r103_fast_variant_g507',
    'r103_hac_snp_g507', 'r103_hac_variant_g507',
    'r103_sup_snp_g507', 'r103_sup_variant_g507',
    # r9 consensus
    'r941_min_hac_g507', 'r941_min_sup_g507',
    'r941_prom_hac_g507', 'r941_prom_sup_g507',
    # r9 variant calling
    'r941_min_hac_variant_g507',
    'r941_prom_hac_variant_g507',
    # r941 e81 consensus
    'r941_e81_fast_g514', 'r941_e81_hac_g514', 'r941_e81_sup_g514',
    # r941 e81 variant calling
    'r941_e81_fast_variant_g514', 'r941_e81_hac_variant_g514',
    'r941_e81_sup_variant_g514',
    # APK joint polishing consensus
    'r1041_e82_260bps_joint_apk_ulk_v5.0.0',
    # bacterial methylation consensus
    'r1041_e82_400bps_bacterial_methylation',
]

bact_methyl_model = 'r1041_e82_400bps_bacterial_methylation'
bact_methyl_compatible_models = [
    'r1041_e82_400bps_hac_v4.2.0', 'r1041_e82_400bps_sup_v4.2.0',
    'r1041_e82_400bps_hac_v4.3.0', 'r1041_e82_400bps_sup_v4.3.0',
    'r1041_e82_400bps_hac_v5.0.0', 'r1041_e82_400bps_sup_v5.0.0',
]

# Previous models that have not been converted to v2.0 format. Attempting to
# load any of these sohuld lead to an error message instructing the user to
# use v1.x.
deprecated_models = [
    # r9 consensus
    'r941_min_fast_g303', 'r941_min_high_g303', 'r941_min_high_g330',
    'r941_prom_fast_g303', 'r941_prom_high_g303', 'r941_prom_high_g330',
    'r941_min_high_g344', 'r941_min_high_g351', 'r941_min_high_g360',
    'r941_prom_high_g344', 'r941_prom_high_g360', 'r941_prom_high_g4011',
    # r10 consensus
    'r10_min_high_g303', 'r10_min_high_g340',
    'r103_min_high_g345', 'r103_min_high_g360', 'r103_prom_high_g360',
    # snp and variant - flipflop
    'r941_prom_snp_g303', 'r941_prom_variant_g303',
    'r941_prom_snp_g322', 'r941_prom_variant_g322',
    'r941_prom_snp_g360', 'r941_prom_variant_g360',
    'r103_prom_snp_g3210', 'r103_prom_variant_g3210',
    # rle consensus
    'r941_min_high_g340_rle',
]

# add basecaller models, then deduplicate and sort
for models in basecaller_models.values():
    archived_models.extend((m for m in models if m is not None))
known_models = sorted(set(current_models + archived_models))
allowed_models = sorted(set(known_models) - set(deprecated_models))


# where we look for model files and store them
model_subdir = 'data'
model_stores = (
    pkg_resources.resource_filename(__package__, model_subdir),
    os.path.join(
        str(pathlib.Path.home()), '.{}'.format(__package__), model_subdir)
)
model_url_template = \
    'https://github.com/nanoporetech/{pkg}/raw/master/{pkg}/{subdir}/{fname}'

# suspect this isn't used anymore...
alignment_params = {
    'rle': "-M 5 -S 4 -O 2 -E 3",
    'non-rle': "-M 2 -S 4 -O 4,24 -E 2,1"}


class DeprecationError(ValueError):
    """Raised when trying to resolve a deprecated model."""

    def __init__(self, model):
        """Initialise and add model name to error message."""
        super().__init__(
            f"Model '{model}' is deprecated in medaka v2.0. "
            f"Please use medaka version 1.x.")
