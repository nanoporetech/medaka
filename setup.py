import os
import re
import sys
from glob import glob
import importlib.util
import platform
import shutil
from setuptools import setup, find_packages, Extension
from setuptools import Distribution, Command
from setuptools.command.build_ext import build_ext
import subprocess

__pkg_name__ = 'medaka'
__dist_name__ = 'medaka'
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

# load the options module so we can decide which models to use in package
_options_path = os.path.join(__pkg_name__, 'options.py')
spec = importlib.util.spec_from_file_location("medaka.options", _options_path)
_options = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_options)

def check_model_lfs():
    # determine if data files look like LFS stubs, fail if they are
    model = _options.default_models['consensus']
    default_model = os.path.join(__pkg_path__, 'data', '{}_model.tar.gz'.format(model))
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
skip_deps = {}
rename_deps = {}
if platform.machine() in {"aarch64", "arm64"}:
    # we used to set parasail and pyspoa here, but both can be fairly easily
    # be built on both macOS and linux nowadays with basic build tool depedencies
    skip_deps = {}
    if platform.system() == "Darwin":
        rename_deps['tensorflow'] = 'tensorflow-macos'
else:
    if os.environ.get('MEDAKA_CPU') is not None:
        rename_deps['tensorflow'] = 'tensorflow-cpu'
dir_path = os.path.dirname(__file__)
install_requires = []
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.strip().startswith('#')
    )
    for req in reqs:
        # possibly rename, then decide if we want it. Doesn't handle changing version
        base = re.split("[<>~=]", req)[0]
        if base in rename_deps:
            req = req.replace(base, rename_deps[base])
        if base not in skip_deps:
            install_requires.append(req)


# to avoid pypi packages getting too big we only bundle some models
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


if __name__ == '__main__':
    check_model_lfs()

    pymajor, pyminor = sys.version_info[0:2]
    if (pymajor < 3) or (pyminor not in {8, 9, 10}):
        raise RuntimeError(
            '`medaka` is unsupported on your version of python, '
            'please use python 3.8-3.10 (inclusive)')

    setup(
        name=__dist_name__,
        version=__version__,
        url='https://github.com/nanoporetech/{}'.format(__pkg_name__),
        author=__author__,
        description=__description__,
        long_description=__long_description__,
        long_description_content_type=__long_description_content_type__,
        python_requires='>=3.8,<3.11',
        packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
        package_data={
            __pkg_name__:[
                os.path.join('data', '{}_model.tar.gz'.format(f))
                for f in bundled_models],
        },
        cffi_modules=["build.py:ffibuilder"],
        install_requires=install_requires,
        entry_points = {
            'console_scripts': [
                '{0} = {0}.{0}:main'.format(__pkg_name__),
                '{0}_counts = {0}.medaka_counts:main'.format(__pkg_name__),
                '{0}_data_path = {0}.{1}:{2}'.format(__pkg_name__, 'common', 'print_data_path'),
                '{0}_version_report = {0}:report_binaries'.format(__pkg_name__, )
            ]
        },
        scripts=[
            'scripts/medaka_consensus',
            'scripts/medaka_haploid_variant',
            'scripts/mini_align', 'scripts/hdf2tf.py'],
        zip_safe=False,
        cmdclass={
            'build_ext': HTSBuild
        },
    )
