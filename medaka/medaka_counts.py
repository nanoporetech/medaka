from medaka.features import FeatureEncoder, pileup_counts
from medaka.common import Region
import sys

import numpy as np

import argparse
from timeit import default_timer as now
import logging

def main():
    # Entry point for testing/checking
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    np.set_printoptions(precision=4, linewidth=100)
    
    
    parser = argparse.ArgumentParser('medaka', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('bam', help='alignment file.')
    parser.add_argument('region', help='alignment region to sample.')
    parser.add_argument('--print', action='store_true', help='print counts.')
    parser.add_argument('--dtypes', nargs='+', help='perform a multi-datatype tests.')
    parser.add_argument('--norm', nargs='+', help='additional normalisation tests. (total, fwd_rev)')

    args = parser.parse_args()
    
    region = Region.from_string(args.region)
    
    kwargs={
        'log_min': None,
        'max_hp_len': 1,
        'is_compressed': False,
        'consensus_as_ref': False,
        'ref_mode': None,
        'with_depth': False,
    }
    
    def _print(samples):
       if args.print:
           for p, f in zip(samples.positions, samples.features):
               print('{}\t{}\t0\t{}\t{}'.format(p[0], p[1], '\t'.join('{:.3f}'.format(x) if x>0.0 else '-' for x in f), sum(f)))
    
    dtype_options = [('',)]
    if args.dtypes is not None:
        dtype_options.append(args.dtypes)
    norm_options = [None, ]
    if args.norm is not None:
        norm_options.extend(args.norm)
    
    for dtypes in dtype_options:
        kwargs['dtypes'] = dtypes
        for norm in norm_options:
            kwargs['normalise'] = norm
    
            print("###########################################################")
            print(kwargs)
            encoder = FeatureEncoder(**kwargs)
        
            # py-style
            t0=now()
            samples = encoder.bam_to_sample(args.bam, region, force_py=True)[0]
            t1=now()
            if not samples.is_empty:
                print(samples.features.shape)
                _print(samples)
            else:
                print("Samples is empty")
            print("---------------------")
        
            # C-style
            t2=now()
            samples = encoder.bam_to_sample(args.bam, region)[0]
            t3=now()
            if not samples.is_empty:
                print(samples.features.shape)
                _print(samples)
            else:
                print("Samples is empty")
        
            print("pysam time:", t1 - t0)
            print("hts time:", t3 - t2)


if __name__ == '__main__':
    main()

