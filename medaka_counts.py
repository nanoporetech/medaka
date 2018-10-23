from medaka.compress import FeatureEncoder, pileup_counts
from medaka.common import parse_regions
import sys

import numpy as np
from collections import Counter

from timeit import default_timer as now
np.set_printoptions(precision=4, linewidth=100)

reads_bam = sys.argv[1]
reg = sys.argv[2]
do_print = sys.argv[3] == 'print'

region = parse_regions([reg])[0]
print(region)

# py-style
t0=now()
kwargs={'consensus_as_ref': False, 'is_compressed': False, 'log_min': None, 'max_hp_len': 1, 'normalise': 'total', 'ref_mode': None, 'with_depth': False}
kwargs['normalise'] = None   # change this just for simple comparison
encoder = FeatureEncoder(**kwargs)
samples = encoder.bam_to_sample(reads_bam, region, ref_rle_fq=None, read_fraction=None)
t1=now()
if do_print:
    for p, f in zip(samples.positions, samples.features):
        print('{}\t{}\t0\t{}\t{}'.format(p[0], p[1], '\t'.join(str(x) for x in f), sum(f)))
print(samples.features.shape)


# C-style
t2 = now()
counts, major_pos, minor_pos = pileup_counts(reg, reads_bam)
t3 = now()
if do_print:
    print("-----------------------------------------------------------")
    for maj, minor, f in zip(major_pos, minor_pos, counts):
        print('{}\t{}\t{}\t{}'.format(maj, minor, '\t'.join(str(x) for x in f), int(sum(f))))
print(counts.shape)

print("pysam time:", t1 - t0)
print("hts time:", t3 - t2)
