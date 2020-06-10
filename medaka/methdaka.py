"""Functionality for aggregating single read methylation calls."""

from collections import Counter, namedtuple
import queue
import sys
import threading
import time

import mappy
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list
from ont_fast5_api.fast5_interface import get_fast5_file
import pysam

import libmedaka
import medaka.common
from medaka.executor import ProcessPoolExecutor, ThreadPoolExecutor

BASECALLANALYSIS = 'Basecall_1D'
MODBASEPATH = 'BaseCalled_template/ModBaseProbs'
FASTQPATH = 'BaseCalled_template/Fastq'
BASES = 'A 6mA C 5mC G T'.split()
MODTYPE = [(b, 'uint8') for b in BASES]
MOTIFS = {
    'dcm': {
        # seq (on fwd strand): [fwd offset, rev offset]
        'CCAGG': ((1, 3), 'MC'),
        'CCTGG': ((1, 3), 'MC')},
    'dam': {
        'GATC': ((1, 2), 'MA')},
    'cpg': {
        'CG': ((0, 1), 'MC')},
    'all': {
        # TODO: these will appear distinct in output, should fix
        'C': ((0, None), 'MC'),
        'G': ((None, 0), 'MC'),
        'A': ((0, None), 'MA'),
        'T': ((None, 0), 'MA')}
    }
Read = namedtuple('Read', ['read_id', 'sequence', 'quality'])


def format_uint8_list(array):
    """Serialize a numpy array to a comma separated list.

    :param array: a uint8 array.

    :returns: a comma separate string.

    """
    ffi = libmedaka.ffi
    lib = libmedaka.lib

    pointer = ffi.cast(
        "uint8_t *",
        array.flatten().ctypes.data)
    result = ffi.new("char[{}]".format(4 * len(array)))
    lib.format_uint8_array(pointer, len(array), result)

    string = ffi.string(result).decode()
    return string


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
            latest = read.get_latest_analysis(BASECALLANALYSIS)
            # get modified base data
            mod_base = read.get_analysis_dataset(latest, MODBASEPATH)
            mod_base = mod_base.view(dtype=MODTYPE)
            mA = 'MA:B:C,{}'.format(format_uint8_list(mod_base['6mA']))
            mC = 'MC:B:C,{}'.format(format_uint8_list(mod_base['5mC']))
            # get basecalling data
            fastq = read.get_analysis_dataset(latest, FASTQPATH)
            header, sequence, _, qstring = fastq.splitlines()
            # put everything together
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
        except Exception as e:
            sys.stderr.write(e)
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
        with pysam.FastaFile(args.reference) as fa:
            for name, length in zip(fa.references, fa.lengths):
                sys.stdout.write('@SQ\tSN:{}\tLN:{}\n'.format(name, length))
        sys.stdout.flush()
        aligner = Aligner(
            args.reference, preset='map-ont', n_threads=args.workers)

        with ThreadPoolExecutor(
                max_items=args.workers, max_workers=args.workers) as executor:
            for read, tags in extractor:
                future = executor.submit(aligner.map, read, tags)
                future.add_done_callback(_write_sam_future)


def _write_sam_future(future):
    try:
        sam = future.result()
        if sam is not None:
            sys.stdout.write('{}\n'.format(sam))
    except Exception:
        pass
    # https://bugs.python.org/issue27144
    future._result = None


class MotifTracker():
    """Track methylation counts of a motif while parsing a pileup."""

    def __init__(self, ref, region, seq, strand_offsets, tag):
        """Initialize the tracker.

        :param ref: reference sequence.
        :param region: region of reference to parse.
        :param seq: the motif search to examine.
        :param strand_offsets: positions of motif to examine on forward and
            reverse strands (`None` if a particular strand shouldn't be
            examined.)
        :param tag: bam tag to look up data.

        """
        self.motif = seq
        self.ref = ref
        self.region = region
        self.seq = seq
        self.tag = tag
        self.fwd_offset, self.rev_offset = strand_offsets
        self.has_fwd = self.fwd_offset is not None
        self.has_rev = self.rev_offset is not None
        self.loci = self._find_loci(
            self.ref, self.motif, self.region.start, self.region.end)
        self._queue = list()
        self.reset_counters()

    def reset_counters(self):
        """Reset methylation counters to zero."""
        self.meth = Counter()
        self.not_meth = Counter()

    def add(self, is_meth, is_rev=None):
        """Add a methylation/non methylation count to the tallies.

        :param is_meth: boolean, the count to which to add.
        :param is_rev: if given checked against internal state to ensure
            caller is adding to the correct count.

        """
        if is_rev is not None:
            assert is_rev == self.is_rev
        if is_meth:
            self.meth[self.is_rev] += 1
        else:
            self.not_meth[self.is_rev] += 1

    @property
    def summary(self):
        """Return a tuple of counts and meta data.

        The data returned is:

            * reference name
            * reference position
            * motif sequence
            * modified count fwd. strand
            * modified count rev. strand
            * non-modified count fwd. strand
            * non-modified count rev. strand

        """
        return (
            self.region.ref_name, self.index, self.motif,
            self.meth[False], self.meth[True],
            self.not_meth[False], self.not_meth[True])

    def __iter__(self):
        """Return self."""
        return self

    def __next__(self):
        """Iterate over the reference positions which need to be examined."""
        if len(self._queue) == 0:
            self.index = next(self.loci)
            if self.has_fwd:
                self._queue.append((False, self.index + self.fwd_offset))
            if self.has_rev:
                self._queue.append((True, self.index + self.rev_offset))
            if self.has_fwd and self.has_rev:
                if self.rev_offset < self.fwd_offset:
                    self._queue = self._queue[::-1]
        item = self._queue.pop(0)
        self.is_rev, self.pos = item
        self.taken_all = len(self._queue) == 0
        return item

    @staticmethod
    def _find_loci(ref, motif, start, end):
        index = start
        while index < end:
            index = ref.find(motif, index)
            if index == -1:
                break
            yield index
            index += 1


def call_methylation(args):
    """Entry point for calling methylation from bam file."""
    logger = medaka.common.get_named_logger('Mextract')
    logger.info("Processing {}.".format(args.bam))
    region = medaka.common.get_regions(args.bam, [args.region])[0]
    ref = pysam.FastaFile(args.reference)
    ref = ref.fetch(region.ref_name)
    motifs = MOTIFS[args.meth]

    additional_text = ''
    if args.meth == 'all':
        additional_text = (
            " Note motifs are stranded, such that both 'C' and 'G' (and 'A' "
            "and 'T') bases will be seen in the output for 5mC (6mA). "
            "Output with a 'G' ('T' for 6mA) correspond to the reverse "
            "strand. Users should post process the output to join/sum lines "
            "belonging to the same genomic loci.")
    elif args.meth == 'dcm':
        additional_text = (
            " The output is sorted by motif, users may wish to subsequently "
            "sorted by position.")
    logger.info(
        "Calling methylation at {} sites.{}".format(
            args.meth, additional_text))

    with open(args.output, 'w') as fh:
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            for seq, motif_info in motifs.items():
                tracker = MotifTracker(ref, region, seq, *motif_info)
                logger.info("Searching for motif '{}'.".format(seq))
                try:
                    next(tracker)
                except StopIteration:
                    continue

                for pileupcolumn in bam.pileup(
                        region.ref_name, region.start, region.end):
                    if (pileupcolumn.pos - region.start) % 1000000 == 0:
                        done = pileupcolumn.pos - region.start
                        pct_done = 100 * done // (region.end - region.start)
                        logger.info(
                            "Processed {} ref positions ({}%)".format(
                                done, pct_done))
                    if pileupcolumn.pos != tracker.pos:
                        continue
                    for read in pileupcolumn.pileups:
                        aln = read.alignment
                        if read.is_del \
                                or read.is_refskip \
                                or aln.is_reverse != tracker.is_rev:
                            continue
                        qpos = read.query_position
                        tag = read.alignment.get_tag(tracker.tag)
                        # NOTE: tags are not reversed, see above
                        if not aln.is_reverse:
                            mscore = tag[qpos]
                        else:
                            mscore = tag[len(tag) - qpos - 1]
                        if mscore > args.filter[1]:
                            tracker.add(True)
                        elif mscore < args.filter[0]:
                            tracker.add(False)

                    if tracker.taken_all:
                        fh.write('\t'.join(str(x) for x in tracker.summary))
                        fh.write('\n')
                        tracker.reset_counters()
                    try:
                        next(tracker)
                    except StopIteration:
                        break
