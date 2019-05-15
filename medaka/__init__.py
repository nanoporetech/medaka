from distutils.version import LooseVersion
import functools
import os
import subprocess

__version__ = '0.7.0'


def check_minimap2_version():
    """get minimap2 version."""
    try:
        proc = subprocess.run(["minimap2", "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        # 2.11-r797
        version = LooseVersion(proc.stdout.decode().split('-')[0])
    except:
        return None

    return version


def check_htslib_tool_version(tool, pos=2):
    try:
        proc = subprocess.run([tool, "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        # tabix (htslib) 1.3.1\n...
        first_line = proc.stdout.decode().split("\n", 1)[0]
        version = first_line.split()[pos]
        version = LooseVersion(version)
    except:
        return None

    return version

check_tabix_version = functools.partial(check_htslib_tool_version, 'tabix')
check_bgzip_version = functools.partial(check_htslib_tool_version, 'bgzip')
check_samtools_version = functools.partial(check_htslib_tool_version, 'samtools', 1)


required_version = {
    'minimap2':LooseVersion('2.11'),
    'samtools':LooseVersion('1.9'),
    'tabix':LooseVersion('1.9'),
    'bgzip':LooseVersion('1.9'),
}

get_version = {
    'minimap2':check_minimap2_version,
    'samtools':check_samtools_version,
    'tabix':check_tabix_version,
    'bgzip':check_bgzip_version,
}



def report_binaries():
    """Print a report of versions of required programs, exits with status 1
    if a bad version is found.

    """
    fail = False
    versions = {prog:get_version[prog]() for prog in required_version.keys()}

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
        print('  '.join((str(x).ljust(width) for x in (prog, versions[prog], required_version[prog], good))))
        all_good &= good
    if not all_good:
        os._exit(1)

