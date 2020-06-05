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
    # r9 consensus
    'r941_min_high_g344', 'r941_min_high_g351', 'r941_min_high_g360',
    'r941_prom_high_g344', 'r941_prom_high_g360',
    # rle consensus
    'r941_min_high_g340_rle',
    # r10 consensus
    'r103_min_high_g345', 'r103_min_high_g360', 'r103_prom_high_g360',
    # snp and variant
    'r941_prom_snp_g360', 'r941_prom_variant_g360',
    'r103_prom_snp_g3210', 'r103_prom_variant_g3210']
archived_models = [
    # r9 consensus
    'r941_min_fast_g303', 'r941_min_high_g303', 'r941_min_high_g330',
    'r941_prom_fast_g303', 'r941_prom_high_g303', 'r941_prom_high_g330',
    # r10 consensus
    'r10_min_high_g303', 'r10_min_high_g340',
    # snp and variant
    'r941_prom_snp_g303', 'r941_prom_variant_g303',
    'r941_prom_snp_g322', 'r941_prom_variant_g322']
allowed_models = sorted(current_models + archived_models)
default_models = {
    'consensus': 'r941_min_high_g360',
    'snp': 'r941_prom_snp_g360',
    'variant': 'r941_prom_variant_g360'}

alignment_params = {
    'rle': "-M 5 -S 4 -O 2 -E 3",
    'non-rle': "-M 2 -S 4 -O 4,24 -E 2,1"}
