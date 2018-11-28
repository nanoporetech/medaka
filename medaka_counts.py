from medaka.features import FeatureEncoder, pileup_counts
from medaka.common import Region
import sys

import numpy as np

from timeit import default_timer as now
import logging
logging.basicConfig(level=logging.DEBUG)

np.set_printoptions(precision=4, linewidth=100)

reads_bam = sys.argv[1]
reg = sys.argv[2]
do_print = sys.argv[3] == 'print'

region = Region.from_string(reg)

kwargs={
    'log_min': None,
    'max_hp_len': 1,
    'is_compressed': False,
    'consensus_as_ref': False,
    'ref_mode': None,
    'with_depth': False,
}

def _print(samples):
   if do_print:
       for p, f in zip(samples.positions, samples.features):
           print('{}\t{}\t0\t{}\t{}'.format(p[0], p[1], '\t'.join('{:.3f}'.format(x) if x>0.0 else '-' for x in f), sum(f)))

dtype_options = ([''], ['r94', 'r10'])
norm_options = (None, 'total', 'fwd_rev')

#dtype_options=([''],)
#norm_options=(None,)

for dtypes in dtype_options:
    kwargs['dtypes'] = dtypes
    for norm in norm_options:
        kwargs['normalise'] = norm

        print("###########################################################")
        print(kwargs)
        encoder = FeatureEncoder(**kwargs)
    
        # py-style
        t0=now()
        samples = encoder.bam_to_sample(reads_bam, region, force_py=True)
        t1=now()
        _print(samples)
        print(samples.features.shape)
        print("---------------------")
    
        # C-style
        t2=now()
        samples = encoder.bam_to_sample(reads_bam, region)
        t3=now()
        _print(samples)
        print(samples.features.shape)
    
        print("pysam time:", t1 - t0)
        print("hts time:", t3 - t2)
