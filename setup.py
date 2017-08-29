import os
import re
from setuptools import setup, find_packages

__path__ = os.path.dirname(__file__)
__pkg_name__ = 'medaka'
__pkg_path__ = os.path.join(__path__, __pkg_name__)

verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))

dir_path = os.path.dirname(__file__)
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    install_requires = [
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.strip().startswith('#')
    ]

setup(
    name='medaka',
    version=__version__,
    url='https://github.com/nanoporetech/medaka',
    author='syoung',
    author_email='stephen.young@nanoporetech.com',
    description='Neural network sequence error correction.',
    packages=find_packages(),
    package_data={
        __pkg_name__:[os.path.join('data','*.h5')],
    },
    install_requires=install_requires,
    entry_points = {
        'console_scripts': [
            '{0}_prepare = {0}.prepare:main'.format(__pkg_name__),
            '{0}_train = {0}.train:main'.format(__pkg_name__),
            '{0}_correct = {0}.correct:main'.format(__pkg_name__),
        ]
    },
    zip_safe=False,
)
