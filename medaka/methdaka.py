"""Functionality for aggregating single read methylation calls."""

from collections import Counter, namedtuple
import queue
import sys
import threading
import time

import mappy
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list
from ont_fast5_api.fast5_interface import get_fast5_file
import pysam

import medaka.common
from medaka.executor import ProcessPoolExecutor, ThreadPoolExecutor


MODBASEPATH = 'BaseCalled_template/ModBaseProbs'
BASES = 'A 6mA C 5mC G T'.split()
MODTYPE = [(b, 'uint8') for b in BASES]
MOTIFS = {
    'dcm': {
        # seq: [fwd offset, rev offset]
        # NB: code will assume rev offset is larger!
        'motifs': {
            'CCAGG': (1, 3),
            'CCTGG': (1, 3)},
        'tag': 'MC'},
    'dam': {
        'motifs': {
            'GATC': (1, 2)},
        'tag': 'MA'},
    'cpg': {
        'motifs': {
            'CG': (0, 1)},
        'tag': 'MC'}
    }
Read = namedtuple('Read', ['read_id', 'sequence', 'quality'])


def align_read(aligner, read, tags=[]):
    """Aligns a `Read` namedtuple to produce a sam record.

    :param aligner: a `mappy.Aligner` instance`
    :param read: a `Read` namedtuple.
    :param tags: list containing additional tags to add to sam record.

    :returns: a (string) sam record.
    """
    try:
        align = list(aligner.map(read.sequence, MD=True, cs=True))[0]
    except IndexError:
        return None
    if align.strand == +1:
        flag = '0'
        seq = read.sequence
        qstring = read.quality
    else:
        flag = '16'
        seq = medaka.common.reverse_complement(read.sequence)
        qstring = read.quality[::-1]
    rname = align.ctg
    pos = str(align.r_st + 1)
    mapq = str(align.mapq)
    clip = [
        '' if x == 0 else '{}S'.format(x)
        for x in (align.q_st, len(seq) - align.q_en)]
    if align.strand == -1:
        clip = clip[::-1]
    cigar = clip[0] + align.cigar_str + clip[1]
    NM = 'NM:i:' + str(align.NM)

    # NOTE: tags written without reversal, see below
    sam = '\t'.join((
        read.read_id, flag, rname, pos, mapq, cigar,
        '*', '0', '0', seq, qstring, NM, *tags))
    return sam


def unaligned_read(read, tags=[]):
    """Create an unaligned sam record for a read.

    :param read: a `Read` namedtuple.
    :param tags: list containing additional tags to add to sam record.

    :returns: a (string) sam record.
    """
    sam = '\t'.join((
        read.read_id, '4', '*', '0', '255', '*',
        '*', '0', '0', read.sequence, read.quality, *tags))
    return sam


def hdf_to_sam_worker(fname):
    """Extract and align basecall and methylation data from `.fast5`.

    :param reference: `.fasta` file containing reference sequence(s).
    :param fname: `.'fast5` file containing read data.
    """
    logger = medaka.common.get_named_logger('ModExtract')
    logger.info("Processing {}.".format(fname))
    results = list()
    with get_fast5_file(fname, mode="r") as f5:
        reads = list(f5.get_read_ids())
        logger.debug("Found {} reads for {}.".format(len(reads), fname))
        for read_id in reads:
            read = f5.get_read(read_id)
            tool = Basecall1DTools(read)

            latest = read.get_latest_analysis('Basecall_1D')
            mod_base = read.get_analysis_dataset(latest, MODBASEPATH)
            mod_base = mod_base.view(dtype=MODTYPE)
            mA = 'MA:B:C,{}'.format(','.join(
                mod_base['6mA'].reshape(-1).astype(str)))
            mC = 'MC:B:C,{}'.format(','.join(
                mod_base['5mC'].reshape(-1).astype(str)))
            header, sequence, qstring = tool.get_called_sequence(
                'template', fastq=False)
            read = Read(read_id, sequence, qstring)
            results.append((read, (mA, mC)))
    return results


class Aligner(object):
    """Mappy alignment interface."""

    def __init__(self, reference, **kwargs):
        """Initialize `mappy` interface.

        :param reference: .fasta reference file.
        :param kwargs: keyword arguments for `mappy.Aligner`.

        """
        self.reference = reference
        self.aligner = None
        self._loadlock = threading.Lock()
        self.kwargs = kwargs
        self.logger = medaka.common.get_named_logger('Aligner')

    def map(self, read, tags=[]):
        """Map/align a read returning a sam record.

        :param read: a `Read` namedtuple.
        :param tags: list containing additional tags to add to sam record.

        """
        if self.aligner is None:
            self._loadlock.acquire()
            if self.aligner is None:
                self.logger.info('Loading alignment index...')
                self.aligner = mappy.Aligner(self.reference, **self.kwargs)
                self.logger.info('Alignment index loaded.')
            self._loadlock.release()
        return align_read(self.aligner, read, tags)


class Extractor(object):
    """Extracts read data from .fast5 files."""

    def __init__(
            self, path, recursive=False, workers=2,
            extractor=hdf_to_sam_worker, max_size=16000):
        """Initialize read data extraction.

        :param path: .fast5 directory path.
        :param recursive: search directory recursively.
        :param workers: number of worker processes.
        :param extractor: callback function to apply to each file.
        :param max_size: approximate maximum size of internal queue.
            When internal queue is at this size. Further read extraction
            will block until results are consumed.

        The results of the extractor can be accessed from the .get() method
        or by iteration.

        """
        self.path = path
        self.recursive = recursive
        self.workers = workers
        self.extractor = extractor
        self.max_size = max_size
        # setting max size for queue has no effect on pausing workers
        # as the queue is filled from a callback in the main thread.
        # Instead we just don't submit new jobs if the queue is above
        # the requested size.
        self.queue = queue.Queue()
        self.files_processed = 0
        self.total_files = 0

        self.logger = medaka.common.get_named_logger("Extractor")
        self.logger.info("Starting worker processes.")
        self._thread = threading.Thread(target=self._run)
        self._thread.daemon = True
        self._thread.start()

    def _run(self):
        """Iterate over input files and stores results into internal queue."""
        fast5s = get_fast5_file_list(self.path, recursive=self.recursive)
        self.total_files = len(fast5s)
        self.logger.info("Found {} files to process.".format(self.total_files))
        with ProcessPoolExecutor(
                self.workers, max_workers=self.workers) as executor:
            for fname in fast5s:
                while True:
                    if self.queue.qsize() < self.max_size:
                        future = executor.submit(self.extractor, fname)
                        future.add_done_callback(self._store)
                        break
                    else:
                        time.sleep(1)
        self.queue.put(None)

    def _store(self, future):
        """Store read data from a future."""
        try:
            results = future.result()
            for r in results:
                self.queue.put(r)
        except Exception:
            pass
        # https://bugs.python.org/issue27144
        future._result = None
        self.files_processed += 1
        self.logger.info("Extracted {}/{} files.".format(
            self.files_processed, self.total_files))

    def get(self):
        """Get a read from the internal queue."""
        item = self.queue.get()
        if item is not None:
            self.queue.task_done()
        else:
            raise StopIteration('All items processed.')
        return item

    def __iter__(self):
        """Iterate over read data from files."""
        while True:
            try:
                yield self.get()
            except StopIteration:
                break


def hdf_to_sam(args):
    """Entry point for converting guppy methylcalled fast5s to sam."""
    logger = medaka.common.get_named_logger('ModExtract')
    logger.info(
        "NOTE: Mod. base scores are output w.r.t the sequencing direction, "
        "not the aligned read orientation.")
    extractor = Extractor(
        args.path,
        recursive=args.recursive, workers=args.io_workers)

    sys.stdout.write('\t'.join(('@HD', 'VN:1.5', 'SO:unsorted')))
    sys.stdout.write('\n')
    sys.stdout.write('\t'.join((
        '@CO', 'Guppy basecaller mod. base tags are stored w.r.t. the '
        'sequencing direction, they should be reversed for reads '
        'aligning to the reverse strand.\n')))
    if args.reference is None:
        # write unaligned sam
        for read, tags in extractor:
            sam = unaligned_read(read, tags)
            sys.stdout.write('{}\n'.format(sam))
    else:
        for name, seq, _ in mappy.fastx_read(
                args.reference, read_comment=False):
            sys.stdout.write('@SQ\tSN:{}\tLN:{}\n'.format(name, len(seq)))
        aligner = Aligner(
            args.reference, preset='map-ont', n_threads=args.workers)

        def _write(future):
            try:
                sam = future.result()
                if sam is not None:
                    sys.stdout.write('{}\n'.format(sam))
            except Exception:
                pass
            # https://bugs.python.org/issue27144
            future._result = None

        with ThreadPoolExecutor(
                max_items=args.workers, max_workers=args.workers) as executor:
            for read, tags in extractor:
                future = executor.submit(aligner.map, read, tags)
                future.add_done_callback(_write)


def call_methylation(args):
    """Entry point for calling methylation from bam file."""
    logger = medaka.common.get_named_logger('Mextract')
    logger.info("Processing {}.".format(args.bam))
    region = medaka.common.get_regions(args.bam, [args.region])[0]
    ref = pysam.FastaFile(args.reference)
    ref = ref.fetch(region.ref_name)
    motifs = MOTIFS[args.meth]['motifs']
    tag_name = MOTIFS[args.meth]['tag']

    def _find_loci(ref, motif, rel_pos, start, end):
        index = start
        while index < end:
            index = ref.find(motif, index)
            if index == -1:
                break
            yield index, index + rel_pos[0], index + rel_pos[1]
            index += 1

    with open(args.output, 'w') as fh:
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for motif, rel_pos in motifs.items():
                logger.info("Searching for motif '{}'.".format(motif))
                loci = _find_loci(
                    ref, motif, rel_pos, region.start, region.end)
                locus = next(loci)

                # we assume the reverse strand position is after
                # the forward strand
                is_rev = False
                next_pos = locus[1]
                meth = Counter()
                not_meth = Counter()
                for pileupcolumn in bam.pileup(
                        region.ref_name, region.start, region.end):
                    if pileupcolumn.pos < next_pos:
                        continue
                    for read in pileupcolumn.pileups:
                        aln = read.alignment
                        if read.is_del \
                                or read.is_refskip \
                                or aln.is_reverse != is_rev:
                            continue
                        qpos = read.query_position
                        tag = read.alignment.get_tag(tag_name)
                        # NOTE: tags are not reversed, see above
                        if not aln.is_reverse:
                            mscore = tag[qpos]
                        else:
                            mscore = tag[len(tag) - qpos - 1]

                        if mscore > args.filter[1]:
                            meth[aln.is_reverse] += 1
                        elif mscore < args.filter[0]:
                            not_meth[aln.is_reverse] += 1
                    if is_rev:
                        fh.write(
                            '\t'.join(str(x) for x in (
                                region.ref_name, locus[0], motif,
                                meth[False], meth[True],
                                not_meth[False], not_meth[True])))
                        fh.write('\n')
                        try:
                            locus = next(loci)
                        except StopIteration:
                            break
                        else:
                            is_rev = False
                            next_pos = locus[1]
                            meth = Counter()
                            not_meth = Counter()
                    else:
                        is_rev = True
                        next_pos = locus[2]
