"""Spanning read clustering module."""

from abc import ABC, abstractmethod
import collections
from copy import copy
from itertools import filterfalse
from typing import Dict, List, Tuple


import edlib
import numpy as np
import wurlitzer

from medaka import abpoa
import medaka.common
import medaka.smolecule
from medaka.tandem.record_name import RecordName


class SpanningReadClusterer(ABC):
    """Abstract base class for spanning read clustering."""

    def cluster_spanningreads(
        self, rec: RecordName, spanning_reads: List[medaka.smolecule.Subread]
    ) -> Tuple[
        Dict[str, int], Dict[RecordName, List[medaka.smolecule.Subread]]
    ]:
        """Clusters spanning reads based on the chosen method.

        Args:
            rec (RecordName): A record name object with an identifier.
            spanning_reads (List[medaka.smolecule.Subread]): spanning reads.

        Returns:
            Tuple[
            Dict[str, int],
            Dict[RecordName,
            List[medaka.smolecule.Subread]]]:
            1. Statistics Dictionary
            2. Clustered spanning reads dict with updated RecordName keys.

        """
        if rec.ploidy == 1:
            new_rec = copy(rec)
            new_rec.hap = 1
            d = self.summarize_reads(
                [s.name for s in spanning_reads],
                prefix="hap1_",
                bhp_counts=False,
            )
            d["phasing_method"] = "unphased"
            return d, {new_rec: spanning_reads}
        elif rec.ploidy == 2:
            return self._cluster_spanningreads(rec, spanning_reads)
        else:
            raise ValueError(f"Unsupported ploidy: {rec.ploidy}")

    @abstractmethod
    def _cluster_spanningreads(
        self,
        recordname: RecordName,
        spanningreads: List[medaka.smolecule.Subread],
    ) -> Tuple[
        Dict[str, int], Dict[RecordName, List[medaka.smolecule.Subread]]
    ]:
        """Clusters spanning reads based on the chosen method.

        Args:
            recordname (RecordName): A record name object with an identifier.
            spanningreads (List[medaka.smolecule.Subread]): spanning reads.

        Returns:
            Tuple[
            Dict[str, int],
            Dict[RecordName,
            List[medaka.smolecule.Subread]]]:
            1. Statistics Dictionary
            2. Clustered spanning reads dict with updated RecordName keys.

        """
        pass

    def summarize_reads(
        self, names, prefix="", bhp_counts=False
    ) -> Dict[str, int]:
        """Generate counts of reads by strand and optionally bam HP tags."""
        records = [RecordName.from_str(n) for n in names]
        counts = collections.Counter()

        if bhp_counts:
            # we want to have 0 counts so initialize
            for h in range(records[0].ploidy + 1):
                counts[f"{prefix}n_reads_BHP{h}"] = 0
            counts.update(f"{prefix}n_reads_BHP{r.hap}" for r in records)

        for strand in "fwd", "rev":
            counts[f"{prefix}n_reads_{strand}"] = 0
        counts.update(f"{prefix}n_reads_{r.strand}" for r in records)

        counts[f"{prefix}n_reads"] += len(names)

        return dict(counts)


class PrephasedClusterer(SpanningReadClusterer):
    """Clusters spanning reads by looking at HP and PS tags from the BAM."""

    def __init__(
        self, remove_outliers: bool = True, min_depth_for_outliers: int = 5
    ):
        """Initialize the PrephasedClusterer.

        Args:
            remove_outliers (bool): flag to remove outliers based on seq len.
            min_depth_for_outliers (int): Min depth threshold.
        """
        self.remove_outlier = remove_outliers
        self.min_depth_for_outliers = min_depth_for_outliers

    def _cluster_spanningreads(
        self, rec: RecordName, spanningreads: List[medaka.smolecule.Subread]
    ) -> Tuple[
        Dict[str, int], Dict[RecordName, List[medaka.smolecule.Subread]]
    ]:
        """Clusters spanning reads by looking at HP and PS tags from the BAM.

        Args:
            rec (RecordName): A record name object with an identifier.
            spanningreads (List[medaka.smolecule.Subread]): spanning reads.

        Returns:
            Tuple[
            Dict[str, int],
            Dict[RecordName,
            List[medaka.smolecule.Subread]]]:
            1. Statistics Dictionary
            2. Clustered spanning reads dict with updated RecordName keys.

        """
        res = PrephasedClusterer._filter_reads_by_dominant_phased_set(
            spanningreads
        )
        spanningreads, filtered = res
        subreads_by_hap = collections.defaultdict(list)
        phasedset_by_hap = collections.defaultdict(int)
        for s in spanningreads:
            rn = RecordName.from_str(s.name)
            subreads_by_hap[rn.hap].append(s)
            phasedset_by_hap[rn.hap] = rn.phased_set

        clustered_reads = {}
        d = {}
        filtered = subreads_by_hap[0]
        for h in range(1, rec.ploidy + 1):
            new_rec = copy(rec)
            new_rec.hap = h
            new_rec.phased_set = phasedset_by_hap[h]
            reads, outliers = self._remove_outliers_subreads(
                subreads_by_hap[h]
            )
            clustered_reads[new_rec] = reads
            filtered += outliers
            d.update(
                self.summarize_reads(
                    [s.name for s in clustered_reads[new_rec]],
                    prefix=f"hap{h}_",
                    bhp_counts=False,
                )
            )

        # handle filtered read as hap0
        # reads are filtered if they are untagged HP:0  or flagged as outliers.
        new_rec = copy(rec)
        new_rec.hap = 0
        clustered_reads[new_rec] = filtered
        d.update(
            self.summarize_reads(
                [s.name for s in clustered_reads[new_rec]],
                prefix="hap0_",
                bhp_counts=False,
            )
        )
        d["phasing_method"] = "prephased"
        return d, clustered_reads

    def _remove_outliers_subreads(self, subreads, multiplier=2):
        """Remove outliers from a list of subreads based on sequence length.

        Args:
            subreads (List[medaka.smolecule.Subread]): List of subreads.
            multiplier (int): Multiplier for the IQR.

        Returns:
            List[medaka.smolecule.Subread]: Filtered subreads without outliers.
            List[medaka.smolecule.Subread]: outliers.

        """
        # don't remove the outliers for few reads because it can be biased
        if (
            not self.remove_outlier
            or len(subreads) <= self.min_depth_for_outliers
        ):
            return subreads, []
        # Convert data to a NumPy array for easier calculations
        read_lengths = np.array([len(read.seq) for read in subreads])
        # Calculate Q1 (25th percentile) and Q3 (75th percentile)
        Q1 = np.percentile(read_lengths, 25)
        Q3 = np.percentile(read_lengths, 75)

        # Calculate the IQR
        IQR = Q3 - Q1
        # Define the lower and upper bounds
        lower_bound = Q1 - multiplier * IQR
        upper_bound = Q3 + multiplier * IQR

        def pred(read):
            return lower_bound <= len(read.seq) <= upper_bound

        reads_pass = list(filter(pred, subreads))
        read_filtered = list(filterfalse(pred, subreads))
        # Filter out the outliers
        return reads_pass, read_filtered

    @staticmethod
    def _filter_reads_by_dominant_phased_set(reads):
        """Filter a list of items based on the most common phased_set.

        Args:
            reads (List[medaka.smolecule.Subread]): subreads to be filtered.

        Returns:
            Tuple[
            List[medaka.smolecule.Subread],
            List[medaka.smolecule.Subread]
            ]:
            1. Reads belonging to the most common phased_set
            2. Rejected reads

        """
        parsed_reads = [RecordName.from_str(item.name) for item in reads]
        # Count occurrences of each phased_set within this hap group
        phased_set_counts = collections.Counter(
            read.phased_set for read in parsed_reads if read.hap != 0
        )

        if len(phased_set_counts) == 0:
            return [], []  # Return an empty list if phased_set_counts is empty

        # Find the most common phased_set
        most_common_phased_set = phased_set_counts.most_common(1)[0][0]

        def pred(read):
            record_name = RecordName.from_str(read.name)
            return record_name.phased_set == most_common_phased_set

        # Filter items to only include those with the most common phased_set
        reads_pass = list(filter(pred, reads))
        reads_filtered = list(filterfalse(pred, reads))

        return reads_pass, reads_filtered


class ABPOAClusterer(SpanningReadClusterer):
    """Cluster spanning reads using the ABPOA method."""

    def __init__(
        self, put_bam_hp_in_name: bool = True, orderings: str = "dsc"
    ):
        """Initialize the ABPOAClusterer.

        Args:
            put_bam_hp_in_name (bool): flag to include the hap identifier
            in the record name.
            orderings (str): The ordering of clusters; 'asc' for ascending,
            'desc' for descending.

        """
        self.put_bam_hp_in_name = put_bam_hp_in_name
        self.orderings = orderings

    @staticmethod
    def rle_seq(seq):
        """
        Run-length encodes a sequence.

        Args:
            seq (str): The input sequence to be encoded.

        Returns:
            str: The run-length encoded sequence.

        """
        return "".join(medaka.common.rle(seq)["value"])

    def _process_clusters(self, recordname, subreads, d):
        """Process clusters of subreads to expected output.

        Args:
            recordname (RecordName): The record name of the region.
            subreads (dict): A dictionary of subreads.
            d (dict): A dictionary containing read information.

        Returns:
            dict: A dictionary of clustered subreads.

        """
        clustered_subreads = {}
        subreads = {
            s.name: s for s in subreads
        }
        for h in range(recordname.ploidy + 1):
            d.update(
                self.summarize_reads(d[f"hap{h}_reads"], prefix=f"hap{h}_")
            )
            new_rec = copy(recordname)
            new_rec.hap = h
            new_rec.query_name += "_HOM" if d["is_homozygous"] else "_HET"
            clustered_subreads[new_rec] = []
            for name in d[f"hap{h}_reads"]:
                s = subreads[name]
                # optionally move SNP-based phasing HP tags into query name
                # and replace with cluster id with
                rn = RecordName.from_str(name)
                if self.put_bam_hp_in_name:
                    rn.query_name += f"_BHP{rn.hap}"
                rn.hap = h
                clustered_subreads[new_rec].append(
                    medaka.smolecule.Subread(str(rn), s.seq)
                )
            del d[f"hap{h}_reads"]
        return clustered_subreads

    def _run_abpoa_clustering(self, subreads, recordname):
        """Run the ABPOA clustering algorithm on a list of subreads.

        Args:
            subreads (list): A list of subreads to be clustered.
            recordname: The name of the record.

        Returns:
            Statistics Dictionary

        Raises:
            None

        """
        logger = medaka.common.get_named_logger("ABPOA")
        # run-length compress sequences
        subreads_rle = []
        for s in subreads:
            subreads_rle.append(
                medaka.smolecule.Subread(s.name, ABPOAClusterer.rle_seq(s.seq))
            )
        # sort reads by length, see https://github.com/yangao07/abPOA/issues/48
        subreads_rle_asc = sorted(subreads_rle, key=lambda s: len(s.seq))
        subreads_rle_dsc = subreads_rle_asc[::-1]
        aligner = abpoa.msa_aligner(aln_mode="g")
        if self.orderings in {"asc", "both"}:
            result_rle_asc, err_asc = ABPOAClusterer.abpoa_consensus(
                aligner, subreads_rle_asc, max_n_cons=recordname.ploidy
            )
            if self.orderings == "asc":
                result_rle_dsc, err_dsc = result_rle_asc, err_asc
                subreads_rle_dsc = subreads_rle_asc
        if self.orderings in {"dsc", "both"}:
            result_rle_dsc, err_dsc = ABPOAClusterer.abpoa_consensus(
                aligner, subreads_rle_dsc, max_n_cons=recordname.ploidy
            )
            if self.orderings == "dsc":
                result_rle_asc, err_asc = result_rle_dsc, err_dsc
                subreads_rle_asc = subreads_rle_dsc
        abpoa_stderr = err_asc + err_dsc

        if abpoa_stderr:
            logger.debug(
                f"Abpoa stderr for {recordname}:\n{abpoa_stderr}"
            )  # orig non-rle compressed reads

        return ABPOAClusterer.check_cluster_read_partitioning(
            result_rle_asc, result_rle_dsc, subreads_rle_asc, subreads_rle_dsc
        )

    def _cluster_spanningreads(
        self, recordname: RecordName, subreads: List[medaka.smolecule.Subread]
    ) -> Tuple[
        Dict[str, int], Dict[RecordName, List[medaka.smolecule.Subread]]
    ]:
        """Clusters spanning reads using the ABPOA method.

        Args:
            recordname (RecordName): A record name object with an identifier.
            spanningreads (List[medaka.smolecule.Subread]): spanning reads.

        Returns:
            Tuple[
            Dict[str, int],
            Dict[RecordName,
            List[medaka.smolecule.Subread]]]:
            1. Statistics Dictionary
            2. Clustered spanning reads dict with updated RecordName keys.

        """
        d = self._run_abpoa_clustering(subreads, recordname)

        clustered_subreads = self._process_clusters(recordname, subreads, d)

        d["phasing_method"] = "abpoa"
        return d, clustered_subreads

    @staticmethod
    def abpoa_consensus(aligner, subreads, max_n_cons=2):
        """Run abpoa consensus."""
        records = [RecordName.from_str(s.name) for s in subreads]
        seqs = [
            s.seq
            if r.strand == "fwd"
            else medaka.common.reverse_complement(s.seq)
            for s, r in zip(subreads, records)
        ]
        # capture abpoa c-code error messages written to sterr so we can report
        # which regions / records had errors.
        with wurlitzer.pipes(bufsize=0) as (out, err):
            result = aligner.msa(
                seqs,
                out_cons=True,
                out_msa=False,
                max_n_cons=max_n_cons,
                min_freq=0.3,
            )
        return result, err.read()

    @staticmethod
    def check_cluster_read_partitioning(res_a, res_b, subreads_a, subreads_b):
        """Assess dependence of read ordering on phasing and generate metrics.

        :param res_a: `abpoa.msa_result` using ordering of subreads subreads_a.
        :param res_b: `abpoa.msa_result` using ordering of subreads subreads_b.
        :param subreads_a: iterable of `medaka.smolecule.Subread` instances.
        :param subreads_b: same as `subreads_a` but with alternative ordering.

        :returns: dict
        """
        assert res_a.n_seq == res_b.n_seq
        # subreads should be in a different order but check we have same reads
        assert {s.name for s in subreads_a} == {s.name for s in subreads_b}
        assert all(r.n_cons <= 2 for r in (res_a, res_b))
        # check consensus seq edit distances to see
        # if cluster ordering has flipped
        cluster_edits = np.zeros(shape=(2, 2), dtype=int)
        cluster_ovlps = [[set(), set()], [set(), set()]]

        def edit_dist(seq1, seq2):
            return edlib.align(seq1, seq2, mode="NW")["editDistance"]

        for a, (a_inds, a_cons) in enumerate(
            zip(res_a.clu_read_ids, res_a.cons_seq)
        ):
            for b, (b_inds, b_cons) in enumerate(
                zip(res_b.clu_read_ids, res_b.cons_seq)
            ):
                a_reads = [subreads_a[i].name for i in a_inds]
                b_reads = [subreads_b[i].name for i in b_inds]
                cluster_ovlps[a][b].update(set(a_reads).intersection(b_reads))
                cluster_edits[a][b] = edit_dist(a_cons, b_cons)

        diag_edits = cluster_edits.trace()
        off_diag_edits = cluster_edits.sum() - diag_edits

        consensus_is_flipped = off_diag_edits < diag_edits
        if not consensus_is_flipped:
            cluster1_reads = cluster_ovlps[0][0]
            cluster2_reads = cluster_ovlps[1][1]
            ambig_reads = cluster_ovlps[0][1] | cluster_ovlps[1][0]
        else:
            cluster1_reads = cluster_ovlps[0][1]
            cluster2_reads = cluster_ovlps[1][0]
            ambig_reads = cluster_ovlps[0][0] | cluster_ovlps[1][1]
            diag_edits, off_diag_edits = off_diag_edits, diag_edits

        # check if clu_read_ids of each cluster are superset of subreads
        # abpoa reports that in come cases reads are not assigned a
        # cluster - these should be inluded in counts of ambigous reads
        all_names = {s.name for s in subreads_a}
        unassigned = set()
        for res, reads in ((res_a, subreads_a), (res_b, subreads_b)):
            assigned = {reads[i].name for c in res.clu_read_ids for i in c}
            # track unassigned reads with either ordering.
            unassigned.update(all_names - assigned)

        # in an ideal case of no POA order dependence, diagonal edits
        # (i.e. between a single haplotype with alternative orderings)
        # will be zero, so edits
        # ratio will be zero. Large values imply stronger order depedence of
        # the differences in consensus sequence between haplotypes.
        edits_ratio = diag_edits / off_diag_edits if diag_edits else 0

        is_homozygous = all(r.n_cons == 1 for r in (res_a, res_b))

        # if we have two clusters, but all reads in second cluster are ambigous
        # and hence removed  such that second cluster has no reads,
        #  make this hom
        # TODO - or should we skip such regions given such issues?
        empty_second_cluster = False
        if (
            max(res_a.n_cons, res_b.n_cons) > 1
            and min(len(cluster1_reads), len(cluster2_reads)) == 0
        ):
            is_homozygous = True
            empty_second_cluster = True
            cluster1_reads = cluster1_reads | cluster2_reads | ambig_reads
            cluster2_reads = set()
            ambig_reads = set()

        n_same_tags, n_switched = None, None
        if not is_homozygous:
            # if we have SNP-based bam HP phasing, use to phase clusters.
            cluster_bhp_ovlps = np.zeros(shape=(2, 2), dtype=int)
            subreads_by_bhp = {1: set(), 2: set()}
            for name in cluster1_reads | cluster2_reads:
                rn = RecordName.from_str(name)
                if rn.hap in subreads_by_bhp:  # don't want to keep hap0
                    subreads_by_bhp[rn.hap].add(name)
            for cid, cluster_names in enumerate(
                [cluster1_reads, cluster2_reads]
            ):
                for bhp_id, bhp_names in subreads_by_bhp.items():
                    cluster_bhp_ovlps[cid, bhp_id - 1] = len(
                        cluster_names & bhp_names
                    )

            n_same_tags = cluster_bhp_ovlps.trace()
            n_switched = cluster_bhp_ovlps.sum() - n_same_tags

            if n_switched > n_same_tags:  # switch cluster1/2 to match HP tags
                cluster1_reads, cluster2_reads = cluster2_reads, cluster1_reads
                n_same_tags, n_switched = n_switched, n_same_tags

        return {
            "n_reads": len(subreads_a),
            "hap1_reads": cluster1_reads,
            "hap2_reads": cluster2_reads,
            "hap0_reads": ambig_reads | unassigned,
            "is_homozygous": is_homozygous,
            "empty_second_cluster": empty_second_cluster,
            "n_ambig_reads": len(ambig_reads),
            "n_unasign_reads": len(unassigned),
            "edits_ratio": round(edits_ratio, 3),
            "diag_edits": diag_edits,
            "nreads_cluster_phasing_matches_bhp": n_same_tags,
            "nreads_cluster_phasing_switched_wrt_bhp": n_switched,
        }


class HybridClusterer(SpanningReadClusterer):
    """Clusters spanning reads using a hybrid method with minimum depth.

    Attributes:
        min_depth (int): Minimum depth threshold for clustering.

    """

    def __init__(self, min_depth: int, remove_outliers: bool = True):
        """Initialize the HybridClusterer with a minimum depth.

        Args:
            min_depth (int): Minimum depth threshold for clustering.
            remove_outliers (bool): flag to remove outliers based on seq len.

        """
        self.min_depth = min_depth
        self.prephased_clusterer = SpanningReadClusterFactory.create_clusterer(
            "prephased", remove_outliers=remove_outliers
        )
        self.abpoa_clusterer = SpanningReadClusterFactory.create_clusterer(
            "abpoa"
        )

    def _cluster_spanningreads(
        self,
        recordname: RecordName,
        spanningreads: List[medaka.smolecule.Subread],
    ) -> Tuple[
        Dict[str, int], Dict[RecordName, List[medaka.smolecule.Subread]]
    ]:
        """Clusters spanning reads using the hybrid method.

        Args:
            recordname (RecordName): A record name object with an identifier.
            spanningreads (List[medaka.smolecule.Subread]): spanning reads.

        Returns:
            Tuple[
            Dict[str, int],
            Dict[RecordName,
            List[medaka.smolecule.Subread]]]:
            1. Statistics Dictionary
            2. Clustered spanning reads dict with updated RecordName keys.

        """
        d, clusters = self.prephased_clusterer.cluster_spanningreads(
            recordname, spanningreads
        )
        for record, cluster in clusters.items():
            if record.hap != 0 and len(cluster) < self.min_depth:
                return self.abpoa_clusterer.cluster_spanningreads(
                    recordname, spanningreads
                )
        return d, clusters


class SpanningReadClusterFactory:
    """Factory to create instances of SpanningReadClusterer."""

    clustering_techniques = ["prephased", "hybrid", "abpoa", "unphased"]

    @staticmethod
    def create_clusterer(method: str, **kwargs) -> SpanningReadClusterer:
        """Create and returns a SpanningReadClusterer.

        Args:
            method (str): The clustering method to use.
            ('prephased', 'abpoa', 'hybrid').
            kwargs: Additional parameters required for clustering methods,
            such as 'min_depth' for the hybrid method.

        Returns:
            SpanningReadClusterer: An instance of SpanningReadClusterer.

        Raises:
            ValueError: If an unknown clustering method is specified
            or required parameters are missing.

        """
        if method in ["prephased", "unphased"]:
            remove_outliers = kwargs.get("remove_outliers", True)
            return PrephasedClusterer(remove_outliers=remove_outliers)
        elif method == "abpoa":
            return ABPOAClusterer(
                put_bam_hp_in_name=kwargs.get("put_bam_hp_in_name", True),
                orderings=kwargs.get("orderings", "dsc"),
            )
        elif method == "hybrid":
            remove_outliers = kwargs.get("remove_outliers", True)
            min_depth = kwargs.get("min_depth")
            if min_depth is None:
                raise ValueError(
                    "Hybrid clustering requires 'min_depth' " "as an argument."
                )
            return HybridClusterer(
                min_depth=min_depth,
                remove_outliers=remove_outliers
                )
        else:
            raise ValueError(f"Unknown clustering method: {method}")
