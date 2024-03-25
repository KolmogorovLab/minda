import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from minda.__version__ import __version__


setup(name='minda',
      version=__version__,
      description='A tool for somatic structural variant calling using long reads',
      url='https://github.com/KolmogorovLab/minda',
      author='Asher Bryant',
      author_email = 'asher.bryant@nih.gov',
      license='BSD-3-Clause',
      packages=['minda'],
      entry_points={'console_scripts': ['minda = minda.main:main']},
      )