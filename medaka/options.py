"""Stores miscellaneous data for medaka."""

# note this module is imported into setup.py, so do not use any
# non-stdlib packages
import os
import pathlib

import pkg_resources


model_subdir = 'data'
model_stores = (
    pkg_resources.resource_filename(__package__, model_subdir),
    os.path.join(
        str(pathlib.Path.home()), '.{}'.format(__package__), model_subdir)
)
model_url_template = \
    'https://github.com/nanoporetech/{pkg}/raw/master/{pkg}/{subdir}/{fname}'

current_models = [
    # r1041 variant calling
    'r1041_e82_400bps_hac_variant_g615', 'r1041_e82_400bps_sup_variant_g615',
    'r1041_e82_400bps_fast_variant_g615',
    # r9 consensus
    'r941_min_hac_g507', 'r941_min_sup_g507',
    'r941_prom_hac_g507', 'r941_prom_sup_g507',
    # r1041 e82 (kit14) consensus
    'r1041_e82_400bps_fast_g615', 'r1041_e82_400bps_hac_g615',
    'r1041_e82_400bps_sup_g615',
    # r9 variant calling
    'r941_min_hac_variant_g507',
    'r941_prom_hac_variant_g507',

]
archived_models = [
    # r9 consensus
    'r941_sup_plant_g610',
    'r941_min_fast_g507', 'r941_prom_fast_g507',
    'r941_min_fast_g303', 'r941_min_high_g303', 'r941_min_high_g330',
    'r941_prom_fast_g303', 'r941_prom_high_g303', 'r941_prom_high_g330',
    'r941_min_high_g344', 'r941_min_high_g351', 'r941_min_high_g360',
    'r941_prom_high_g344', 'r941_prom_high_g360', 'r941_prom_high_g4011',
    # r10 consensus
    'r10_min_high_g303', 'r10_min_high_g340',
    'r103_min_high_g345', 'r103_min_high_g360', 'r103_prom_high_g360',
    'r103_fast_g507', 'r103_hac_g507', 'r103_sup_g507',
    # r104 e81 consensus
    'r104_e81_fast_g5015', 'r104_e81_sup_g5015', 'r104_e81_hac_g5015',
    'r104_e81_sup_g610',
    # r104 e81 variant calling
    'r104_e81_fast_variant_g5015', 'r104_e81_hac_variant_g5015',
    'r104_e81_sup_variant_g610',
    # snp and variant - flipflop
    'r941_prom_snp_g303', 'r941_prom_variant_g303',
    'r941_prom_snp_g322', 'r941_prom_variant_g322',
    'r941_prom_snp_g360', 'r941_prom_variant_g360',
    'r103_prom_snp_g3210', 'r103_prom_variant_g3210',
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
    # rle consensus
    'r941_min_high_g340_rle',
    # r941 e81 consensus
    'r941_e81_fast_g514', 'r941_e81_hac_g514', 'r941_e81_sup_g514',
    # r941 e81 variant calling
    'r941_e81_fast_variant_g514', 'r941_e81_hac_variant_g514',
    'r941_e81_sup_variant_g514',
]
allowed_models = sorted(current_models + archived_models)
default_models = {
    'consensus': 'r941_min_hac_g507',
    'variant': 'r941_min_hac_variant_g507'}

alignment_params = {
    'rle': "-M 5 -S 4 -O 2 -E 3",
    'non-rle': "-M 2 -S 4 -O 4,24 -E 2,1"}
