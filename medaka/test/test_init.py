import unittest

from medaka import check_htslib_tool_version, report_binaries

class TestHTSLibVersion(unittest.TestCase):

    def test_000_version_return_type(self):
        check_htslib_tool_version('samtools')

    def test_001_report_binaries(self):
        report_binaries()

