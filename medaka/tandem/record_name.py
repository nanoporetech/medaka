"""Recordname structure for encoding/decoding metadata to/from strings."""
import re

import medaka.common


class RecordName(object):
    """Simple class for encoding/decoding sequence metadata to/from strings."""

    def __init__(
        self,
        *,
        query_name,
        ref_name,
        ref_start,
        ref_end,
        hap=0,
        phased_set=0,
        ploidy=1,
        strand="fwd",
        ref_start_padded=None,
        ref_end_padded=None,
    ):
        """Initialize repeat read analysis.

        :param query_name: str
        :param ref_name: str
        :param ref_start: int
        :param ref_end: int
        :param hap: int
        :param phased_set: int
        :param ploidy: int
        :param strand: str, {fwd, rev}
        :param ref_start_pad: padded reference start, int
        :param ref_end_pad: padded reference end, int.
        """
        self.query_name = query_name
        self.ref_name = ref_name
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.hap = hap
        self.phased_set = phased_set
        self.ploidy = ploidy
        self.strand = strand
        self.ref_start_padded = (
            ref_start if ref_start_padded is None else ref_start_padded
        )
        self.ref_end_padded = (
            ref_end if ref_end_padded is None else ref_end_padded
        )

    def __str__(self):
        """Encode as a string."""
        return (
            f"{self.query_name}_{self.ref_name}_"
            f"{self.ref_start}_{self.ref_end}_"
            f"pad_{self.ref_start_padded}_{self.ref_end_padded}_"
            f"{self.strand}_"
            f"hap{self.hap}_"
            f"phased-set{self.phased_set}_"
            f"ploidy{self.ploidy}"
        )

    @classmethod
    def from_str(cls, name):
        """Decode from string."""
        pattern = (
            r"(?P<query_name>.+)_(?P<ref_name>.+)_"
            r"(?P<ref_start>\d+)_(?P<ref_end>\d+)_"
            r"pad_(?P<ref_start_padded>\d+)_(?P<ref_end_padded>\d+)_"
            r"(?P<strand>fwd|rev)_hap(?P<hap>\d+)_"
            r"phased-set(?P<phased_set>\d+)_ploidy(?P<ploidy>\d+)"
        )
        m = re.match(pattern, name)
        if m is None:
            raise ValueError(f"Could not parse {name} using {pattern}")
        d = m.groupdict()
        for field in (
            "ref_start",
            "ref_end",
            "hap",
            "ref_start_padded",
            "ref_end_padded",
            "phased_set",
            "ploidy",
        ):
            d[field] = int(d[field])
        return cls(**d)

    def sorter(self):
        """Get sorting keys."""
        return self.ref_name, self.ref_start

    def to_padded_region(self):
        """Convert to padded `medaka.common.Region`."""
        return medaka.common.Region(
            self.ref_name, self.ref_start_padded, self.ref_end_padded
        )

    def to_unpadded_region(self):
        """Convert to unpadded `medaka.common.Region`."""
        return medaka.common.Region(
            self.ref_name, self.ref_start, self.ref_end
        )
