# coding: utf-8
from copy import deepcopy
import numpy as np
import intervaltree
import medaka.vcf
bases = np.array(['A', 'C', 'G', 'T'])
seq = ''.join(np.random.choice(bases, size=2000, replace=True))

with open('test_ref.fasta','w') as fh:
    fh.write('>contig1\n')
    fh.write(''.join(seq) + '\n')


def make_variant(seq, tree, max_indel=5, chrom='contig1'):
    """Return a random variant (variable-length indel or single base sub) which
    does not overlap with any intervals in a supplied tree and update the tree.
    """
    attempts = 0
    while True:
        attempts += 1
        vpos = np.random.randint(max_indel + 1, len(seq)-max_indel - 1)
        indel_len = np.random.randint(-max_indel, max_indel + 1)
        if indel_len < 0:  # del
            ref = seq[vpos: vpos + abs(indel_len) + 1]
            alt = seq[vpos]
        elif indel_len > 1:  # ins
            ref = seq[vpos]
            alt = ref + ''.join(np.random.choice(bases, indel_len))
        else:  # single nucleotide sub
            ref = seq[vpos]
            alt = np.random.choice([b for b in bases if str(b) != str(ref)])
        v = medaka.vcf.Variant(chrom, vpos, ref=ref, alt=alt, qual=np.random.randint(1,10), genotype_data={'GT':'1/1'})
        vtrimmed = v.trim()
        end = vtrimmed.pos + len(vtrimmed.ref)
        interval = intervaltree.Interval(v.pos, end + 1)
        if not tree.overlaps(interval):
            tree.add(interval)
            # print('Variant {} succeeded after {}'.format(len(tree), attempts))
            break
    return vtrimmed


# make sure variants on one haplotype don't overlap with variants on the same
# haplotype
tree1 = intervaltree.IntervalTree()
max_indel=10

variants_common = [make_variant(seq, tree1) for p in range(100)]
tree2 = deepcopy(tree1)
variants1 = [make_variant(seq, tree1) for p in range(100)]
variants2 = [make_variant(seq, tree2) for p in range(100)]
with medaka.vcf.VCFWriter('test_hap1.vcf', 'w') as vcf1:
    vcf1.write_variants(variants1 + variants_common)

with medaka.vcf.VCFWriter('test_hap2.vcf', 'w') as vcf2:
    vcf2.write_variants(variants2 + variants_common)
