"""Module for generating consensus sequences from reads."""
from contextlib import ExitStack
import dataclasses
import math
import multiprocessing
import os
import tempfile
import time
from timeit import default_timer as now
from typing import Dict, List, Optional

import pysam

import medaka.common
import medaka.smolecule
from medaka.tandem.alignment import align_consensus_fx_to_ref
from medaka.tandem.io import SpanningReadsExtractor
from medaka.tandem.polisher import Polisher
from medaka.tandem.record_name import RecordName
from medaka.tandem.spanning_read_clusterer import SpanningReadClusterer


class InsufficientCoverage(Exception):
    """Exception for tracking cases of insufficient read coverage."""

    pass


@dataclasses.dataclass
class ConsensusResult:
    """Simple class for organizing consensus results and alignments."""

    rec: RecordName
    subreads: tuple
    consensus_seq: str = ""
    consensus_alignments: tuple = dataclasses.field(default_factory=tuple)
    ref_seq: str = ""
    exception: Optional[Exception] = None


class ConsensusGenerator:
    """A class to generate consensus sequences for specified regions."""

    PARASAIL_ALIGNER_NAME = "nw_trace_scan_32"

    def __init__(
        self,
        regions: List[RecordName],
        bam: str,
        ref: str,
        reads_clusterer: SpanningReadClusterer,
        min_depth: int,
        reads_filter: Dict[str, int],
        output_prefix: str,
        process_large_regions: bool = False,
        model=None,
    ):
        """Initialize the ConsensusGenerator with the specified parameters.

        Args:
            regions (List[RecordName]): list of genomic regions to process.
            bam (str): Path to the BAM file with sequencing reads.
            ref (str): Path to the reference genome file.
            reads_clusterer (SpanningReadClusterer): Clustering Method
            ('abpoa', 'prephased', 'hybrid').
            min_depth (int): Minimum depth threshold for filtering reads.
            reads_filter (Dict[str, int]): Filtering criteria for reads.
            process_large_regions (bool): Flag to process large regions.
            output_prefix (str): Prefix for output files.
            model (medaka.model.Model): Medaka model for polishing.

        """
        self.regions = regions
        self.bam_reader = SpanningReadsExtractor(bam, reads_filter)
        self.ref = ref
        self.reads_clusterer = reads_clusterer
        self.min_depth = min_depth
        self.process_large_regions = process_large_regions
        self.output_prefix = output_prefix
        self.output_files = []
        self.max_region_size = 10000
        self.logger = medaka.common.get_named_logger("ConsensusGenerator")
        self.min_mapq = reads_filter.get("min_mapq", 0)

        self.header = {"HD": {"VN": 1.0}, "SQ": []}
        self.all_consensus_alignments = []
        self.model = model

        op = output_prefix
        self.poa_file = os.path.join(op, "poa.fasta")
        self.trimmed_reads_file = os.path.join(op, "trimmed_reads.fasta")
        self.ref_chunks_file = os.path.join(op, "ref_chunks.fasta")
        self.skipped_bed_file = os.path.join(op, "skipped.bed")
        self.skipped_large_file = os.path.join(op, "skipped_large.bed")
        self.trimmed_to_poa_bam = os.path.join(op, "trimmed_reads_to_poa.bam")
        self.cons_to_ref_bam = os.path.join(op, "medaka_to_ref.bam")

    def _open_output_files(self, stack: ExitStack) -> None:
        """Open all output files under the supplied ExitStack."""
        for f in (
            self.poa_file,
            self.trimmed_reads_file,
            self.ref_chunks_file,
        ):
            medaka.common.remove_if_exist(f + ".fai")

        # region‑metrics
        self.metrics_fhs = {
            k: stack.enter_context(
                open(
                    os.path.join(self.output_prefix,
                                 f"{k}_region_metrics.txt"),
                    "w",
                )
            )
            for k in ("prephased", "abpoa", "unphased")
        }

        # main outputs
        self.trimmed_reads_fh = stack.enter_context(
            open(self.trimmed_reads_file, "w")
        )
        self.ref_chunks_fh = stack.enter_context(
            open(self.ref_chunks_file, "w")
        )
        self.skipped_bed_fh = stack.enter_context(
            open(self.skipped_bed_file, "w")
        )
        self.skipped_large_fh = stack.enter_context(
            open(self.skipped_large_file, "w")
        )
        self.poa_fh = stack.enter_context(open(self.poa_file, "w"))

    def write_metrics(self, rec: RecordName, m: dict):
        """Write the metrics for a given record.

        Args:
            rec (Record): The record to write metrics for.
            m (dict): The metrics to write.

        """
        metrics = rec.to_unpadded_region()._asdict()
        method = m.pop("phasing_method")
        metrics.update(m)
        if not self.metrics_fhs[method].tell():  # file empty
            self.metrics_fhs[method].write("\t".join(metrics.keys()) + "\n")
        self.metrics_fhs[method].write(
            "\t".join(str(v) for v in metrics.values()) + "\n"
        )

    def write_spanning_reads(
        self, spanning_reads: List[medaka.smolecule.Subread]
    ) -> None:
        """Write the spanning reads to the trimmed reads file.

        Args:
            spanning_reads (list): list of spanning reads.

        """
        for read in spanning_reads:
            self.trimmed_reads_fh.write(f">{read.name}\n{read.seq}\n")

    def write_consensus(self, consensus: ConsensusResult):
        """Write the consensus sequence and reference sequence."""
        self.poa_fh.write(
            ">{}\n{}\n".format(str(consensus.rec), consensus.consensus_seq)
        )
        self.ref_chunks_fh.write(
            f">{str(consensus.rec)}\n{consensus.ref_seq}\n"
        )
        self.header["SQ"].append(
            {"LN": len(consensus.consensus_seq), "SN": str(consensus.rec)}
        )
        self.all_consensus_alignments.append(consensus.consensus_alignments)

    def get_subreads(self, rec: RecordName) -> List[medaka.smolecule.Subread]:
        """Get the subreads for a given record.

        Args:
            rec (RecordName): The record to get subreads for.

        Returns:
            list: The list of subreads.

        """
        sub_reads = self.bam_reader.get_subreads(rec)
        if len(sub_reads) < self.min_depth:
            self.logger.info(
                f"{rec}: Retrieved too few reads "
                f"({len(sub_reads)} < {self.min_depth})"
            )
            self.skipped_bed_fh.write(
                f"{rec.ref_name}\t{rec.ref_start}\t{rec.ref_end}\t{rec}\n"
            )
            return []

        if not self.process_large_regions:
            max_trimmed_read_len = max(len(r.seq) for r in sub_reads)
            if max_trimmed_read_len > self.max_region_size:
                self.logger.info(
                    f"{rec}: Region of length({max_trimmed_read_len})"
                    f" > {self.max_region_size}"
                )
                self.skipped_large_fh.write(
                    f"{rec.ref_name}\t"
                    f"{rec.ref_start}\t"
                    f"{rec.ref_end}\t"
                    f"{rec}\n"
                )
                return []

        return sub_reads

    def generate_consensus(
        self,
        clustered_spanning_reads: Dict[RecordName,
                                       List[medaka.smolecule.Subread]]
    ):
        """Generate the consensus sequences for clustered subreads.

        Args:
            clustered_spanning_reads (Dict[RecordName,
                        List[medaka.smolecule.Subread]]):
                                        clustered subreads.

        """
        for record, sub_reads in clustered_spanning_reads.items():
            if (
                record.hap == 0
            ):  # dont generate consensus for ambigous reads
                self.write_spanning_reads(sub_reads)
                continue
            if record.hap == 2 and "_HOM" in record.query_name:
                continue
            if len(sub_reads) < self.min_depth:
                self.logger.info(
                    f"{record}: Retrieved too few reads "
                    f"({len(sub_reads)} < {self.min_depth})"
                )
                self.skipped_bed_fh.write(
                    f"{record.ref_name}\t"
                    f"{record.ref_start}\t"
                    f"{record.ref_end}\t{record}\n"
                )
                continue

            self.write_spanning_reads(sub_reads)
            consensus = self.consensus_pileup_from_reads(
                record, sub_reads
            )
            self.write_consensus(consensus)

    def polish(self):
        """Polish the consensus sequences."""
        self.polished_consensus = os.path.join(
            self.output_prefix, "consensus.fasta"
        )
        polisher = Polisher(
            input_fasta=self.poa_file,
            bam=self.trimmed_to_poa_bam,
            threads=1,
            output_dir=self.output_prefix,
            output_fasta=self.polished_consensus,
        )
        try:
            polisher.polish(
                batch_size=50,
                chunk_len=1000,
                chunk_ovlp=250,
                save_features=False,
                model=self.model,
                min_mapq=self.min_mapq,
            )
        except Exception as e:
            self.logger.error(
                f"Error while polishing {self.output_prefix}: {e}"
            )
            raise e

    def process(self) -> int:
        """Process regions and write polished consensus sequences."""
        with ExitStack() as stack:
            self._open_output_files(stack)

            for rec in self.regions:
                sub_reads = self.get_subreads(rec)
                if not sub_reads:
                    continue

                m, clustered = self.reads_clusterer.cluster_spanningreads(
                    rec, sub_reads
                )
                self.write_metrics(rec, m)
                self.generate_consensus(clustered)

        medaka.smolecule.write_bam(
            self.trimmed_to_poa_bam,
            self.all_consensus_alignments,
            self.header,
        )

        if self.all_consensus_alignments:
            # free space for no longer needed alignments
            self.all_consensus_alignments = []
            self.polish()
            align_consensus_fx_to_ref(
                self.polished_consensus,
                self.cons_to_ref_bam,
                self.ref,
            )

        return len(self.regions)

    def consensus_pileup_from_reads(
        self, rec: RecordName, subreads: List[medaka.smolecule.Subread]
    ) -> ConsensusResult:
        """Run consensus and generate read-consensus alignments."""
        if isinstance(rec, str):
            rec = RecordName.from_str(rec)
        non_empty_subreads = [s for s in subreads if s.seq != "N"]
        if len(non_empty_subreads) < self.min_depth:
            num_empty_subreads = len(subreads) - len(non_empty_subreads)
            self.logger.info(
                f"{rec}: {num_empty_subreads} out of {len(subreads)} reads "
                "support the deletion of the tandem repeat array."
            )
            # Refactor this part by having another stream for empty sequences
            # so that medaka won't crash on empty consensus
            res = ConsensusResult(rec, subreads=tuple(subreads))
            res.consensus_seq = "N"
            subread_names = [s.name for s in subreads]
            res.consensus_alignments = (
                ConsensusGenerator.create_empty_alignments(
                    subread_names, str(rec)
                )
            )
            return res
        # Sort reads for consistent results.
        # I am starting with shorter reads first
        # since we hypothesize that the basecaller is overcalling.
        non_empty_subreads.sort(
            key=lambda read: (len(read.seq), read.name), reverse=True
        )
        res = ConsensusResult(rec, subreads=tuple(non_empty_subreads))

        if (
            len(res.subreads) < self.min_depth
        ):  # definately insufficient coverage
            # though if partial == True could still be insufficient depth z
            # extra checks later
            res.exception = InsufficientCoverage(
                f"{rec}: Too few reads "
                f"({len(res.subreads)} < {self.min_depth})"
            )
            return res

        consensus_read = medaka.smolecule.Read(
            name=str(rec), subreads=res.subreads
        )
        # avoid realignment to discover strandedness
        consensus_read._orient = [
            RecordName.from_str(s.name).strand == "fwd" for s in res.subreads
        ]
        consensus_read._initialized = True
        res.consensus_seq = consensus_read.poa_consensus(method="abpoa")
        # estimate depth of coverage from number of bp and round to be generous
        depth = sum(len(s) for s in subreads)
        if round(depth) < self.min_depth:
            res.exception = InsufficientCoverage(
                f"{rec}: Estimated coverage too low ({depth:.2f} < "
                f"{self.min_depth})"
            )
            return res

        # use global aligner
        consensus_read.parasail_aligner_name = self.PARASAIL_ALIGNER_NAME

        # align sub-reads (i.e. trimmed pileup) to consensus
        res.consensus_alignments = consensus_read.align_to_template(
            template=res.consensus_seq,
            template_name=consensus_read.name,
        )
        return res

    @staticmethod
    def merge_results(out_dir, temp_outputs_prefix: List[str]):
        """Merge the results of the consensus generation process.

        Args:
            out_dir (str): Path to the output directory.
            temp_outputs_prefix (List[str]): list of temporary output folders.

        Returns:
            None

        """
        logger = medaka.common.get_named_logger("TR")
        trimmed_to_poa = os.path.join(out_dir, "trimmed_reads_to_poa.bam")
        # write bam of trimmed reads aligned to poa consensus
        logger.info(f"Writing trimmed reads to {trimmed_to_poa}.")
        temp_bams = [
            os.path.join(prefix, "trimmed_reads_to_poa.bam")
            for prefix in temp_outputs_prefix
        ]
        temp_bams = [bam for bam in temp_bams if os.path.exists(bam)]
        if temp_bams != []:
            pysam.merge("-f", trimmed_to_poa, *temp_bams)
            pysam.index(trimmed_to_poa)

        consensus_to_ref = os.path.join(out_dir, "medaka_to_ref.bam")
        # write bam of trimmed reads aligned to poa consensus
        logger.info(f"Writing consensus mapping to ref {consensus_to_ref}.")
        temp_bams = [
            os.path.join(prefix, "medaka_to_ref.bam")
            for prefix in temp_outputs_prefix
        ]
        temp_bams = [bam for bam in temp_bams if os.path.exists(bam)]
        if temp_bams != []:
            pysam.merge("-f", consensus_to_ref, *temp_bams)
            pysam.index(consensus_to_ref)

        outputFiles = [
            "trimmed_reads.fasta",
            "ref_chunks.fasta",
            "skipped.bed",
            "skipped_large.bed",
            "poa.fasta",
            "consensus.fasta",
        ]

        for output in outputFiles:
            out_file = os.path.join(out_dir, output)
            temp_files = [
                os.path.join(prefix, output) for prefix in temp_outputs_prefix
            ]
            temp_files = [f for f in temp_files if os.path.exists(f)]
            logger.info(f"Writing to {out_file}.")
            if ".fasta" in output:
                medaka.common.remove_if_exist(out_file + ".fai")
            medaka.common.concat_files(temp_files, out_file, has_header=False)

        for method in ("prephased", "abpoa", "unphased"):
            out_file_suffix = f"{method}_region_metrics.txt"
            out_file = os.path.join(out_dir, out_file_suffix)
            temp_files = [
                os.path.join(prefix, f"{out_file_suffix}")
                for prefix in temp_outputs_prefix
            ]
            temp_files = [f for f in temp_files if os.path.exists(f)]
            logger.info(f"Writing to {out_file}.")
            medaka.common.concat_files(temp_files, out_file, has_header=True)

    @staticmethod
    def create_empty_alignments(subread_names, template_name):
        """Create empty alignments when subreads or template contain 'N'.

        :param subread_names: List of subread names.
        :param template_name: Name of the template sequence.

        :returns: List of `Alignment` objects representing empty alignments.
        """
        alignments = []
        for sr in subread_names:
            rstart = 0
            cigar = "1M"
            flag = 0
            aln = medaka.smolecule.Alignment(
                template_name, sr, flag, rstart, "N", cigar
            )
            alignments.append(aln)
        return alignments


class ParallelConsensusGenerator:
    """A class to run ConsensusGenerator in parallel using multiprocessing."""

    PERCENT_OF_INPUT_TO_BE_PROCESSED_SIMULTANEOUSLY = 15
    BIG_REGION_SIZE = 2000

    def __init__(
        self,
        regions: List[RecordName],
        bam: str,
        ref: str,
        reads_clusterer: SpanningReadClusterer,
        min_depth: int,
        reads_filter: Dict[str, int],
        process_large_regions: bool,
        output_prefix: str,
        num_processes: int,
        model=None,
    ):
        """Initialize the ParallelConsensusGenerator class.

        Args:
            regions (List[RecordName]): list of genomic regions to process.
            bam (str): Path to the BAM file with sequencing reads.
            ref (str): Path to the reference genome file.
            reads_clusterer (SpanningReadClusterer): The method to use for
            consensus generation('abpoa', 'prephased', 'hybrid').
            min_depth (int): Minimum depth threshold for filtering reads.
            reads_filter (Dict[str, int]): Filtering criteria for reads.
            process_large_regions (bool): Flag to process large regions.
            output_prefix (str): Prefix for output files.
            num_processes (int): Number of processes to run in parallel.
            model (medaka.model.Model): Medaka model for polishing.

        """
        self.n_regions = len(regions)
        self.bam = bam
        self.ref = ref
        self.reads_clusterer = reads_clusterer
        self.min_depth = min_depth
        self.reads_filter = reads_filter
        self.process_large_regions = process_large_regions
        self.output_prefix = output_prefix
        self.num_processes = num_processes
        self.t0 = now()
        self.t1, self.tlast = self.t0, self.t0
        self.n_done = 0
        self.split_regions(regions)
        self.model = model

    def process(self):
        """Process the genomic regions in parallel using multiprocessing."""
        running_jobs = []
        finished_temp_results = []
        logger = medaka.common.get_named_logger("ParallelConsensusGenerator")
        error = False
        with multiprocessing.Pool(
            processes=self.num_processes, maxtasksperchild=1
        ) as pool:
            for sub_region in self.sub_regions:
                result = pool.apply_async(
                    self.run_consensus_generator,
                    (
                        sub_region,
                        self.bam,
                        self.ref,
                        self.reads_clusterer,
                        self.min_depth,
                        self.reads_filter,
                        self.process_large_regions,
                        self.output_prefix,
                        self.model,
                    ),
                )
                running_jobs.append(result)

            num_sub_regions = len(self.sub_regions)
            self.sub_regions = []

            while len(running_jobs) > 0 and not error:
                res = self._collect_results(running_jobs, logger)
                running_jobs, finished_results, error = res
                finished_temp_results.extend(finished_results)
                time.sleep(
                    1
                )  # Sleep for a short period before checking the result again

        if error or len(finished_temp_results) != num_sub_regions:
            logger.error(
                "Some processes failed. "
                "Please check the logs for more information."
            )
            return False

        logger.info(
            "All processes finished successfully. Ready to merge the results."
        )
        ConsensusGenerator.merge_results(
            self.output_prefix, finished_temp_results
        )
        self.clean_subregions_files()
        for temp_results in finished_temp_results:
            medaka.common.remove_directory(temp_results)

        total_time = now() - self.t0
        logger.info(f"Total time taken: {total_time:.1f}s")

        return True

    def _collect_results(self, jobs, logger):
        running_jobs = []
        error = False
        finished_temp_results = []
        for result in jobs:
            if result.ready():
                if result.successful():
                    temp_results, n_processed_regions = result.get()
                    finished_temp_results.append(temp_results)
                    self.update_progress(n_processed_regions, logger)
                else:
                    try:
                        error = True
                        result.get()
                    except Exception as e:
                        logger.error(f"Error encountered: {e}")
                        break
            else:
                running_jobs.append(result)
        return running_jobs, finished_temp_results, error

    def update_progress(self, n_processed, logger):
        """Update the progress of the consensus generation process."""
        self.t1 = now()
        self.n_done += n_processed
        if self.t1 - self.tlast > 10:
            self.tlast = self.t1
            logger.info(
                f"{self.n_done/self.n_regions:.1%} Done "
                f"({self.n_done}/{self.n_regions} regions) "
                f"in {self.t1-self.t0:.1f}s"
            )

    def clean_subregions_files(self):
        """Remove the temporary files associated with the sub-regions."""
        for sub_region in self.sub_regions:
            medaka.common.remove_if_exist(sub_region)

    def split_regions(self, regions):
        """Split regions into sub-regions for parallel processing.

        Args:
            regions (list): list of regions to be split.

        Returns:
            None

        """
        # Split regions into fixed number of sub-regions
        logger = medaka.common.get_named_logger("ParallelConsensusGenerator")

        num_regions_per_process = int(
            math.ceil(
                len(regions)
                / (
                    self.num_processes
                    * self.PERCENT_OF_INPUT_TO_BE_PROCESSED_SIMULTANEOUSLY
                )
            )
        )
        logger.info(
            f"Splitting {len(regions)} regions into "
            f"{num_regions_per_process} regions per process."
        )

        self.sub_regions = []
        big_regions = [
            r for r in regions
            if (r.ref_end - r.ref_start) > self.BIG_REGION_SIZE
        ]
        rest = [
            r for r in regions
            if (r.ref_end - r.ref_start) <= self.BIG_REGION_SIZE
        ]
        regions = big_regions + rest
        for i in range(0, len(regions), num_regions_per_process):
            region_file = tempfile.NamedTemporaryFile(
                dir=self.output_prefix, delete=False
            )
            region_file.write(
                "\n".join(
                    str(r) for r in regions[i: i + num_regions_per_process]
                ).encode()
            )
            region_file.close()
            self.sub_regions.append(region_file.name)

    @staticmethod
    def run_consensus_generator(
        sub_regions,
        bam,
        ref,
        reads_clusterer,
        min_depth,
        reads_filter,
        process_large_regions,
        output_prefix,
        model,
    ):
        """Run the consensus generator for a given set of sub-regions.

        Args:
            sub_regions (str): Path to the file containing sub-regions.
            bam (str): Path to the BAM file.
            ref (str): Path to the reference genome.
            reads_clusterer (SpanningReadClusterer): Clustering method:
            ('abpoa', 'prephased', 'hybrid').
            min_depth (int): Minimum depth threshold.
            reads_filter (Dict[str, int]): Filtering criteria for reads.
            process_large_regions (bool): Flag to process large regions.
            output_prefix (str): Prefix for the output files.
            model (str): Path to the model.

        Returns:
            tuple: A tuple containing the path to the temporary folder
                    and the number of processed regions.

        """
        # Create a temporary folder for each sub-region
        logger = medaka.common.get_named_logger("ParallelConsensusGenerator")
        processed = 0
        try:
            temp_folder = tempfile.mkdtemp(dir=output_prefix)
            with open(sub_regions, "r") as f:
                regions = [
                    RecordName.from_str(line.strip()) for line in f.readlines()
                ]
            medaka.common.remove_if_exist(sub_regions)
            consensus_generator = ConsensusGenerator(
                regions=regions,
                bam=bam,
                ref=ref,
                reads_clusterer=reads_clusterer,
                min_depth=min_depth,
                reads_filter=reads_filter,
                process_large_regions=process_large_regions,
                output_prefix=temp_folder,
                model=model,
            )
            processed = consensus_generator.process()
        except Exception as e:
            logger.error(f"Error encountered: {e}")
            raise e

        return temp_folder, processed
