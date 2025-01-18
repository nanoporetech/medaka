"""Consensus and variant calling of Nanopore sequencing data."""
from distutils.version import LooseVersion
import functools
import os
import subprocess
import sys

__version__ = "2.0.1"

try:
    import pyabpoa as abpoa
except ImportError:
    sys.stderr.write(
        "Cannot import pyabpoa, some features may not be available.\n")
    abpoa = None
try:
    import wurlitzer
except ImportError:
    sys.stderr.write(
        "Cannot import wurlitzer, some features may not be available.\n")
    wurlitzer = None


def check_minimap2_version():
    """Get minimap2 version."""
    try:
        proc = subprocess.run(
            ["minimap2", "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        # 2.11-r797
        version = LooseVersion(proc.stdout.decode().split('-')[0])
    except Exception:
        return None

    return version


def check_htslib_tool_version(tool, pos=2):
    """Get version of an htslib program.

    :param tool: program name.
    :param pos: the position index of the item containing the version
        information. e.g. `pos=2` extracts `1.3.1` from
        `tabix (htslib) 1.3.1`.

    :returns: the `LooseVersion` number.
    """
    try:
        proc = subprocess.run([tool, "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        # tabix (htslib) 1.3.1\n...
        first_line = proc.stdout.decode('utf-8',errors='replace').split("\n", 1)[0]
        version = first_line.split()[pos]
        version = LooseVersion(version)
    except Exception:
        return None

    return version


check_tabix_version = functools.partial(
    check_htslib_tool_version, 'tabix')
check_bgzip_version = functools.partial(
    check_htslib_tool_version, 'bgzip')
check_samtools_version = functools.partial(
    check_htslib_tool_version, 'samtools', 1)
check_bcftools_version = functools.partial(
    check_htslib_tool_version, 'bcftools', 1)

required_version = {
    'minimap2': LooseVersion('2.11'),
    'samtools': LooseVersion('1.11'),
    'tabix': LooseVersion('1.11'),
    'bgzip': LooseVersion('1.11'),
    'bcftools': LooseVersion('1.11'),
}

get_version = {
    'minimap2': check_minimap2_version,
    'samtools': check_samtools_version,
    'tabix': check_tabix_version,
    'bgzip': check_bgzip_version,
    'bcftools': check_bcftools_version,
}


def report_binaries():
    """Print a report of versions of required programs.

    The program will exit with status 1 if a bad version is found.

    """
    versions = {prog: get_version[prog]() for prog in required_version.keys()}

    width = 9
    cols = "Program Version Required Pass".split()
    print('  '.join([x.ljust(width) for x in cols]))
    all_good = True
    for prog in sorted(required_version.keys()):
        if versions[prog] is None:
            good = False
            versions[prog] = 'Not found'
        else:
            good = versions[prog] >= required_version[prog]
        print('  '.join((
            str(x).ljust(width)
            for x in (prog, versions[prog], required_version[prog], good))))
        all_good &= good
    if not all_good:
        os._exit(1)
