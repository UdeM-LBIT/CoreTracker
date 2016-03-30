#!/usr/bin/env python
#
# setup for CoreTracker library packages
#
# use the following to install:
#   python setup.py build
#   python setup.py install
#
from __future__ import division, absolute_import, print_function
import os, sys
from numpy.distutils.core  import Extension
import setuptools
from coretracker import __project__, __version__
sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "python")))

def configuration(top_path='') :
    fexact = Extension(name='statlib.fexact', sources=[os.path.join(top_path, 'coretracker/FisherExact/statlib', 'FEXACT.F90')])
    asa159 = Extension(name='statlib.asa159', sources=[os.path.join(top_path, 'coretracker/FisherExact/statlib', 'asa159.f90')])
    asa205 = Extension(name='statlib.asa205', sources=[os.path.join(top_path, 'coretracker/FisherExact/statlib', 'asa205.f90')])
    return [fexact, asa205, asa159]

def setup_package():
    if os.path.exists('README.md'):
        README = open('README.md').read()
    else:
        README = ""  # a placeholder, readme is generated on release
    print("\nVersion : %s\n"%__version__)

    from numpy.distutils.core import setup
    fortran_extnsion = configuration()

    setup(
        name=__project__,
        version=__version__,
        maintainer='UdeM-LBIT',
        description="CoreTracker, A codon reassignment tracker",
        url='https://github.com/UdeM-LBIT/CoreTracker',
        author='Emmanuel Noutahi',
        author_email='fmr.noutahi@umontreal.ca',
        scripts = ['bin/coretracker', 'bin/coretranslate', 'bin/corefusion'],

        entry_points={'console_scripts': []},

        keywords="bioinformatics codon reassignment tracker",

        long_description=(README + '\n'),
        license='GPL',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Natural Language :: English',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Operating System :: POSIX',
            'Programming Language :: Python',
        ],

        setup_requires=['numpy'],
        install_requires = [
            'ete3',
            'numpy >= 1.8.1',
            'pandas',
            'scikit-learn',
            'scipy',
            'biopython',
            'matplotlib',
            'WeasyPrint',
            'PyYAML',
        ],
        ext_modules = fortran_extnsion
    )

if __name__ == '__main__':
    setup_package()
