from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import re
import subprocess
import sys
import tempfile
from timeit import default_timer as now

import mappy
import numpy as np
import parasail
import pysam

import medaka.common

re_split_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)")
def first_cigar(cigar):
    """Extract details of the first operation in a cigar string.

    :param cigar: cigar string.

    :returns: op. length, op. type

    """
    m = re.search(re_split_cigar, cigar)
    return m.group('len'), m.group('op')


def parasail_to_sam(result, seq):
    """Extract reference start and sam compatible cigar string
    from a parasail result object.

    :param result: parasail alignment result.
    :param seq: query sequence.

    :returns: reference start coordinate, cigar string.

    """
    cigstr = result.cigar.decode.decode()

    first = first_cigar(cigstr)
    prefix = ''.join(first)
    rstart = 0
    cliplen = result.cigar.beg_query
    clip = '' if cliplen == 0 else '{}S'.format(cliplen)
    if first[1] == 'I':
        pre = '{}S'.format(int(first[0]) + cliplen)
    elif first[1] == 'D':
        pre = clip
        rstart = int(first[0])
    else:
        pre = '{}{}'.format(clip, prefix)

    mid = cigstr[len(prefix):]
    end_clip = len(seq) - result.end_query - 1
    suf = '{}S'.format(end_clip) if end_clip > 0 else ''
    new_cigstr = ''.join((pre, mid, suf))
    return rstart, new_cigstr


comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    #'-': '-'
}
comp_trans = str.maketrans(''.join(comp.keys()), ''.join(comp.values()))


def reverse_complement(seq):
    """Reverse complement sequence.

    :param: input sequence string.

    :returns: reverse-complemented string.

    """
    return seq.translate(comp_trans)[::-1]


spoa='spoa'
Subread = namedtuple('Subread', 'name seq')
Alignment = namedtuple('Alignment', 'rname qname flag rstart seq cigar')

class Read(object):

    def __init__(self, name, subreads, initialize=False):
        """Functionality to extract information from a read with subreads.

        :param name: read name.
        :param subreads: list of subreads.
        :param initialize: initialize subread alignments.

        """
        self.name = name
        self.subreads = subreads
        if len(self.subreads) == 0:
            raise ValueError("Cannot create a read with fewer than 0 subreads.")
        self.consensus = self.subreads[0].seq
        # SW alignments of subreads to consensus
        self._alignments = None
        self._alignments_valid = False
        # orientation of subreads wrt consensus
        self._orient = None
        # has initialize() been run (i.e. are the above filled in)
        self._initialized = False
        # has a consensus been run
        self.consensus_run = False


    def initialize(self):
        # we find initial alignments with parasail as mappy may not find some
        if not self._initialized:
            self._alignments = self.orient_subreads()
            self._alignments_valid = True
            self._initialized = True
        return None


    @classmethod
    def from_fastx(cls, fastx, name=None):
        """Create a Read from a fasta/q file.

        :param fastx: input file path.
        :param name: name of Read. If not given the filename will be used.

        """
        try:
            read = next(cls.multi_from_fastx(fastx, take_all=True))
        except Exception as e:
            raise IOError("Could not create Read from file {}.".format(fastx))
        return read


    @classmethod
    def multi_from_fastx(cls, fastx, take_all=False, read_id=None, depth_filter=1, length_filter=0):
        """Create multiple `Read` s from a fasta/q file, assuming subreads
        are grouped by read and named with <read_id>_<subread_id>.

        :param fastx: input file path.
        :param take_all: skip check on subread_ids, take all subreads in one
            `Read`.
        :param read_id: name of `Read`. Only used for `take_all == True`. If
            not given the basename of the input file is used.
        :param depth_filter: require reads to have at least this many subreads.
        :param length_filter: require reads to have a median subread length above
            this value.

        """
        depth_filter = max(1, depth_filter)
        if take_all and read_id is None:
            read_id = os.path.splitext(os.path.basename(fastx))[0]
        else:
            read_id = None
        subreads = []
        with pysam.FastxFile(fastx) as fh:
            for entry in fh:
                if not take_all:
                    cur_read_id = entry.name.split("_")[0]
                    if cur_read_id != read_id:
                        if len(subreads) >= depth_filter:
                            med_length = np.median([len(x.seq) for x in subreads])
                            if med_length > length_filter:
                                yield cls(read_id, subreads)
                        read_id = cur_read_id
                        subreads = []
                if len(entry.sequence) > 0:
                    subreads.append(Subread(entry.name, entry.sequence))

            if len(subreads) >= depth_filter:
                med_length = np.median([len(x.seq) for x in subreads])
                if med_length > length_filter:
                    yield cls(read_id, subreads)


    @property
    def seqs(self):
        """A list of the subread sequences."""
        return [x.seq for x in self.subreads]


    @property
    def nseqs(self):
        """The number of subreads contained in the read."""
        return len(self.subreads)


    def poa_consensus(self, additional_seq=None, method='racon'):
        """Create a consensus sequence for the read."""
        self.initialize()
        spoa_seq = None
        with tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False) as fh:
            if additional_seq is not None:
                fh.write(">{}\n{}\n".format('additional', additional_seq))
            for orient, subread in zip(self._orient, self.subreads):
                if method == 'spoa':
                    seq = subread.seq if orient else reverse_complement(subread.seq)
                else:
                    seq = subread.seq
                fh.write(">{}\n{}\n".format(subread.name, seq))
            fh.flush()

            if method == 'spoa':
                consensus_seq = self._run_spoa(fh.name)
            elif method == 'racon':
                consensus_seq = self._run_racon(fh.name)
            else:
                raise ValueError('Unrecognised method: {}.'.format(method))
        self.consensus = consensus_seq
        self._alignments_valid = False
        self.consensus_run = True
        return consensus_seq


    def _run_spoa(self, fasta):
        opts = []
        #opts=['-l', '1','-r', '0', '-m', '8', '-n', '-6', '-g', '-8', '-e', '-8']
        out = subprocess.check_output([spoa, fasta] + opts)
        spoa_seq = out.decode().split()[2]
        return spoa_seq


    def _run_racon(self, fasta):
        tname = 'consensus_{}'.format(self.name)
        file_ext = 'sam' # paf route is not enabled since mappy no faster
        source = self._alignments
        if not self._alignments_valid:
            self._alignments = self.align_to_template(self.consensus, tname)
            self._alignments_valid = True
            source = self._alignments

        header = {
            'HD': {'VN': 1.0},
            'SQ': [{
                'LN': len(self.consensus),
                'SN': tname,
            }]
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            ref_fasta = os.path.join(tmpdir, 'racon_ref.fasta')
            with open(ref_fasta, 'w') as fh:
                fh.write(">{}\n{}\n".format(tname, self.consensus))

            overlaps = os.path.join(tmpdir, 'racon_in.{}'.format(file_ext))
            if file_ext == 'sam':
                write_bam(overlaps, [self._alignments], header, bam=False)
            else: #paf
                with open(overlaps, 'w') as fh:
                    for src in source:
                        fh.write('{}\n'.format(src))

            opts = ['-m', '8', '-x', '-6', '-g', '-8']
            try:
                out = subprocess.check_output(
                    ['racon', fasta, overlaps, ref_fasta] + opts,
                    stderr=subprocess.PIPE
                )
            except subprocess.CalledProcessError as e:
                print("RACON FAILED", file_ext)
                print(e.stdout)
                print(e.stderr)
                print(e.cmd)
            racon_seq = out.decode().splitlines()[1]
        return racon_seq


    def orient_subreads(self):
        """Find orientation of subreads with respect to consensus sequence.

        :returns: `Alignment` s of subreads to consensus.

        """
        # TODO: use a profile here
        # TODO: refactor with align_to_template
        self._orient = []
        alignments = []
        for sr in self.subreads:
            rc_seq = reverse_complement(sr.seq)
            result_fwd = parasail.sw_trace_striped_16(sr.seq, self.consensus, 8, 4, parasail.pam100)
            result_rev = parasail.sw_trace_striped_16(rc_seq, self.consensus, 8, 4, parasail.pam100)
            is_fwd = result_fwd.score > result_rev.score
            self._orient.append(is_fwd)
            result = result_fwd if is_fwd else result_rev
            seq = sr.seq if is_fwd else rc_seq
            if result.cigar.beg_ref >= result.end_ref or result.cigar.beg_query >= result.end_query:
                # unsure why this can happen
                continue
            rstart, cigar = parasail_to_sam(result, seq)
            flag = 0 if is_fwd else 16
            aln = Alignment('consensus_{}'.format(self.name), sr.name, flag, rstart, seq, cigar)
            alignments.append(aln)
        return alignments


    def align_to_template(self, template, template_name):
        """Align subreads to a template sequence using Smith-Waterman.

        :param template: sequence to which to align subreads.
        :param template_name: name of template sequence.

        :returns: `Alignment` tuples.

        """
        alignments = []
        for orient, sr in zip(self._orient, self.subreads):
            seq = sr.seq if orient else reverse_complement(sr.seq)
            result = parasail.sw_trace_striped_16(seq, template, 8, 4, parasail.pam100)
            if result.cigar.beg_ref >= result.end_ref or result.cigar.beg_query >= result.end_query:
                # unsure why this can happen
                continue
            rstart, cigar = parasail_to_sam(result, seq)
            flag = 0 if orient else 16
            aln = Alignment(template_name, sr.name, flag, rstart, seq, cigar)
            alignments.append(aln)
        return alignments


    def mappy_to_template(self, template, template_name, align=True):
        """Align subreads to a template sequence using minimap.

        :param template: sequence to which to align subreads.
        :param template_name: name of template sequence.
        :param align: retrieve cigar string (else produce paf)

        :returns: `Alignment` tuples.

        """
        # align False requires forked minimap2, and isn't much faster for
        # a small number of sequences due to index construction time.
        align = True
        alignments = []
        aligner = mappy.Aligner(seq=template, preset='map-ont')
        for sr in self.subreads:
            try:
                hit = next(aligner.map(sr.seq))
            except StopIteration:
                continue
            else:
                flag = 0 if hit.strand == 1 else 16
                seq = sr.seq if hit.strand == 1 else reverse_complement(sr.seq)
                if align:
                    clip = ['' if x == 0 else '{}S'.format(x) for x in (hit.q_st, len(sr.seq) - hit.q_en)]
                    if hit.strand == -1:
                        clip = clip[::-1]
                    cigstr = ''.join((clip[0], hit.cigar_str, clip[1]))
                    aln = Alignment(template_name, sr.name, flag, hit.r_st, seq, cigstr)
                else:
                    # return paf string
                    aln = '\t'.join(str(x) for x in (
                        sr.name, len(sr.seq), hit.q_st, hit.q_en, '+' if hit.strand == +1 else '-', template_name, hit.ctg_len,
                        hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq,
                        'tp:A:P','ts:A:.','cg:Z:'+hit.cigar_str
                    ))
                alignments.append(aln)
                hit = None
        return alignments


def write_bam(fname, alignments, header, bam=True):
    mode = 'wb' if bam else 'w'
    with pysam.AlignmentFile(fname, mode, header=header) as fh:
        for ref_id, subreads in enumerate(alignments):
            for aln in sorted(subreads, key=lambda x: x.rstart):
                a = pysam.AlignedSegment()
                a.reference_id = ref_id
                a.query_name = aln.qname
                a.query_sequence = aln.seq
                a.reference_start = aln.rstart
                a.cigarstring = aln.cigar
                a.flag = aln.flag
                a.mapping_quality = 60
                fh.write(a)
    if mode == 'wb':
        pysam.index(fname)


def _read_worker(read, align=True):
    read.initialize()
    if read.nseqs > 2: #skip if theres only one subread
        for it in range(2):
            read.poa_consensus(method='racon')
    aligns = None
    if align:
        aligns = read.mappy_to_template(template=read.consensus, template_name=read.name)
    return read.name, read.consensus, aligns


def poa_workflow(reads, threads):
    #TODO: this is quite memory inefficient, but we can only build the header
    #      by seeing everything.
    logger = medaka.common.get_named_logger('POAManager')
    pool = ThreadPoolExecutor(max_workers=threads)
    futures = []
    for read in reads:
        logger.debug("Adding {} to queue.".format(read.name))
        futures.append(pool.submit(_read_worker, read))

    header = {
        'HD': {'VN': 1.0},
        'SQ': [],
    }
    consensuses = []
    alignments = []

    for fut in as_completed(futures):
        try:
            res = fut.result()
        except Exception as e:
            logger.warning(e)
            pass
        else:
            rname, consensus, aligns = res
            logger.debug('Finished {}.'.format(rname))
            if consensus is not None:
                header['SQ'].append({
                    'LN': len(consensus),
                    'SN': rname,
                })
                consensuses.append([rname, consensus])
                alignments.append(aligns)
    return header, consensuses, alignments


def main(args):
    # arg parser does not supply these
    args.tag_name = None
    args.tag_value = None
    args.tag_keep_missing = False

    logger = medaka.common.get_named_logger('Smolecule')
    medaka.common.mkdir_p(args.output, info='Results will be overwritten.')

    def _multi_file_reader():
        for fname in args.fasta:
            try:
                yield Read.from_fastx(fname)
            except:
                pass

    if len(args.fasta) > 1:
        logger.info("Given {} input files, assuming one read per file.".format(len(args.fasta)))
        reads = _multi_file_reader()
    else:
        logger.info("Given one input file, subreads are assumed to be grouped by read.")
        reads = Read.multi_from_fastx(args.fasta[0], depth_filter=args.depth, length_filter=args.length)

    logger.info("Running pre-medaka POA consensus for all reads.")
    t0 = now()
    header, consensuses, alignments = poa_workflow(reads, args.threads)
    t1 = now()

    logger.info("Writing medaka input bam for {} reads.".format(len(alignments)))
    bam_file = os.path.join(args.output, 'subreads_to_spoa.bam')
    write_bam(bam_file, alignments, header)

    spoa_file = os.path.join(args.output, 'poa.fasta')
    with open(spoa_file, 'w') as fh:
        for rname, cons in consensuses:
            fh.write('>{}\n{}\n'.format(rname, cons))

    logger.info("Running medaka consensus.")
    t2 = now()
    args.bam = bam_file
    out_dir = args.output
    args.output = os.path.join(out_dir, 'consensus.hdf')
    medaka.inference.predict(args)
    t3 = now()

    logger.info("Running medaka stitch.")
    args.inputs = [args.output]
    args.output = os.path.join(out_dir, 'consensus.fasta')
    args.regions = None
    medaka.stitch.stitch(args)
    logger.info("Single-molecule consensus sequences written to {}.".format(args.output))
    logger.info("POA time: {:.0f}s, medaka time: {:.0f}s".format(t1 - t0, t3 - t2))


if __name__ == '__main__':
    main()
