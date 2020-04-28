import os
import re
import sys
from glob import glob
import importlib.util
import shutil
from setuptools import setup, find_packages, Extension
from setuptools import Distribution, Command
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
import subprocess

__pkg_name__ = 'medaka'
__author__ = 'ont-research'
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

# find version
verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))

_options_path = os.path.join(__pkg_name__, 'options.py')
spec = importlib.util.spec_from_file_location("medaka.options", _options_path)
_options = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_options)

def check_model_lfs():
    # determine if data files look like LFS stubs, fail if they are
    model = _options.default_models['consensus']
    default_model = os.path.join(__pkg_path__, 'data', '{}_model.hdf5'.format(model))
    stub_signature = "ASCII text"
    if os.path.exists(default_model):
        stdout = subprocess.check_output(['file', default_model])
    else:
        print("Failed to determine filetype of {}.".format(default_model))
        sys.exit(1)
    github_link = "github.com/nanoporetech/{}".format(__pkg_name__)
    if stub_signature in stdout.decode():
        print(
            "\nModel files appear to be git lfs stubs. Please refer to README.md\n"
            "for instructions regarding git lfs. If you obtained {} through\n"
            "a method other than cloning the repository from github, please raise\n"
            "an issue on {}.".format(__pkg_name__, github_link))
        sys.exit(1)

# create requirements from requirements.txt
dir_path = os.path.dirname(__file__)
install_requires = []
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.strip().startswith('#')
    )
    for req in reqs:
        if req.startswith('git+https'):
            req = req.split('/')[-1].split('@')[0]
        install_requires.append(req)

# locate and any third party binaries
data_files = []
if os.environ.get("MEDAKA_BINARIES") is not None:
    with open(os.path.join(dir_path, 'Makefile')) as fh:
        for line in fh.readlines():
            tokens = line.split('=')
            if tokens[0] == 'BINARIES':
                exes = tokens[1].split()
                break
    #place binaries as package data, below we'll copy them to standard path in dist
    data_files.append(
        ('exes', [
            'bincache/{}'.format(x, x) for x in exes
        ])
    )

# to avoid pypi pacakges getting too big we only bundle some models
if os.environ.get('MEDAKA_DIST') is not None:
    bundled_models = _options.current_models
else:
    bundled_models = _options.allowed_models
print("Bundling models: {}".format(bundled_models))

class HTSBuild(build_ext):
    # uses the Makefile to build libhts.a, this will get done before the cffi extension
    def run(self):

        def compile_hts():
            subprocess.check_call(['make', 'libhts.a'])

        self.execute(compile_hts, [], 'Compiling htslib using Makefile')
        build_ext.run(self)

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


if __name__ == '__main__':
    check_model_lfs()

    pymajor, pyminor = sys.version_info[0:2]
    if (pymajor < 3) or (pyminor not in {5, 6}):
        raise RuntimeError(
            '`medaka` is unsupported on your version of python, '
            'please use python 3.5 or python 3.6.')

    setup(
        name='medaka',
        version=__version__,
        url='https://github.com/nanoporetech/{}'.format(__pkg_name__),
        author=__author__,
        description=__description__,
        long_description=__long_description__,
        long_description_content_type=__long_description_content_type__,
        python_requires='>=3.5.*,<3.7',
        packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
        package_data={
            __pkg_name__:[os.path.join('data', '{}_model.hdf5'.format(f))
                          for f in bundled_models],
        },
        cffi_modules=["build.py:ffibuilder"],
        install_requires=install_requires,
        #place binaries as package data, below we'll copy them to standard path in dist
        data_files=data_files,
        entry_points = {
            'console_scripts': [
                '{0} = {0}.{0}:main'.format(__pkg_name__),
                '{0}_counts = {0}.medaka_counts:main'.format(__pkg_name__),
                '{0}_data_path = {0}.{1}:{2}'.format(__pkg_name__, 'common', 'print_data_path'),
                '{0}_version_report = {0}:report_binaries'.format(__pkg_name__, )
            ]
        },
        scripts=['scripts/medaka_consensus', 'scripts/medaka_variant', 'scripts/mini_align'],
        zip_safe=False,
        cmdclass={
            'build_ext': HTSBuild
        },
    )

    if os.environ.get('MEDAKA_BINARIES') is not None:
        print("\nCopying utility binaries to your path.")
        get_setuptools_script_dir()
