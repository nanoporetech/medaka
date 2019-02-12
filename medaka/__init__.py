from distutils.version import LooseVersion
import os
import subprocess

__version__ = '0.5.2'


def check_minimap2_version():
    """Checks minimap2 version is greater than or equal to that required."""
    try:
        proc = subprocess.run(["minimap2", "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        # 2.11-r797
        version = LooseVersion(proc.stdout.decode().split('-')[0])
    except:
        return None

    return version


def check_samtools_version():
    """Checks samtools version is greater than or equal to that required."""
    try:
        proc = subprocess.run(["samtools", "--version"], stdout=subprocess.PIPE)
        if proc.returncode != 0:
            return None
        # samtools 1.3.1\n...
        first_line = proc.stdout.decode().split("\n", 1)[0]
        version = first_line.split()[1]
        version = LooseVersion(version)
    except:
        return None

    return version


required_version = {
    'samtools':LooseVersion('1.3.1'),
    'minimap2':LooseVersion('2.11'),
}

get_version = {
    'samtools':check_samtools_version,
    'minimap2':check_minimap2_version,
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

