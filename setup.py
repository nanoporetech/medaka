import os
import re
import sys
from glob import glob
import shutil
from setuptools import setup, find_packages, Extension
from setuptools import Distribution, Command
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
import subprocess

#TODO: fill in these
__pkg_name__ = 'medaka'
__author__ = 'syoung'
__description__ = 'Neural network sequence error correction.'

# Use readme as long description and say its github-flavour markdown
from os import path
this_directory = path.abspath(path.dirname(__file__))
kwargs = {'encoding':'utf-8'} if sys.version_info.major == 3 else {}
with open(path.join(this_directory, 'README.md'), **kwargs) as f:
    __long_description__ = f.read()
__long_description_content_type__ = 'text/markdown'

__path__ = os.path.dirname(__file__)
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_name__))

verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))


dir_path = os.path.dirname(__file__)
install_requires = []
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.strip().startswith('#')
    )
    for req in reqs:
        if req.startswith('git+https'):
            req.split('/')[-1].split('@')[0]
        install_requires.append(req)


data_files = []
if os.environ.get('MEDAKA_BINARIES') is not None:
    exes = ['samtools', 'minimap2']
    data_files.append(
        ('exes', [
            'bincache/{}'.format(x, x) for x in exes
        ])
    )


class HTSBuild(build_ext):
    # uses the Makefile to build libhts.a, this will get done before the cffi extension
    def run(self):

        def compile_hts():
            subprocess.check_call(['make', 'libhts.a'])

        self.execute(compile_hts, [], 'Compiling htslib using Makefile')
        build_ext.run(self)


setup(
    name='medaka',
    version=__version__,
    url='https://github.com/nanoporetech/{}'.format(__pkg_name__),
    author=__author__,
    author_email='{}@nanoporetech.com'.format(__author__),
    description=__description__,
    long_description=__long_description__,
    long_description_content_type=__long_description_content_type__,
    python_requires='>=3.4.*,<3.7',
    packages=find_packages(),
    package_data={
        __pkg_name__:[os.path.join('data','*.hdf5')],
    },
    cffi_modules=["build.py:ffibuilder"],
    install_requires=install_requires,
    #place binaries as package data, below we'll copy them to standard path in dist
    data_files=data_files,
    entry_points = {
        'console_scripts': [
            '{0} = {0}.{0}:main'.format(__pkg_name__),
            'medaka_counts = {0}.medaka_counts:main'.format(__pkg_name__),
            'medaka_data_path = {0}.{1}:{2}'.format(__pkg_name__, 'common', 'print_data_path'),
        ]
    },
    scripts=['scripts/medaka_consensus', 'scripts/mini_align'],
    zip_safe=False,
    cmdclass={
        'build_ext': HTSBuild
    },
)


# Nasty hack to get binaries into bin path
class GetPaths(install):
    def run(self):
        self.distribution.install_scripts = self.install_scripts
        self.distribution.install_libbase = self.install_libbase

def get_setuptools_script_dir():
    # Run the above class just to get paths
    dist = Distribution({'cmdclass': {'install': GetPaths}})
    dist.dry_run = True
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()

    print(dist.install_libbase)
    src_dir = glob(os.path.join(dist.install_libbase, 'medaka-*', 'exes'))[0]
    for exe in (os.path.join(src_dir, x) for x in os.listdir(src_dir)):
        print("Copying", os.path.basename(exe), '->', dist.install_scripts)
        shutil.copy(exe, dist.install_scripts)
    return dist.install_libbase, dist.install_scripts

if os.environ.get('MEDAKA_BINARIES') is not None:
    print("\nCopying utility binaries to your path.")
    get_setuptools_script_dir()
