#!/usr/bin/env python
#
# setup for CoreTracker library packages

from __future__ import absolute_import
from __future__ import print_function
import os, sys
try:
    from setuptools import find_packages
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import find_packages

from numpy.distutils.core  import Extension as Ext
from numpy.distutils.core import setup

from coretracker import __project__, __version__
sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "python")))

def configuration(top_path='') :
    fexact_sources = [
        'coretracker/FisherExact/statlib/fexact.pyf',
        'coretracker/FisherExact/statlib/FEXACT.F90',
        'coretracker/FisherExact/statlib/prterr.f'
    ]

    fexact = Ext(name='coretracker.FisherExact.statlib.fexact', sources=[os.path.join(top_path,x) for x in fexact_sources])
    asa159 = Ext(name='coretracker.FisherExact.statlib.asa159', sources=[os.path.join(top_path, 'coretracker/FisherExact/statlib', 'asa159.f90')])
    asa205 = Ext(name='coretracker.FisherExact.statlib.asa205', sources=[os.path.join(top_path, 'coretracker/FisherExact/statlib', 'asa205.f90')])
    return [fexact, asa205, asa159]

def setup_package():
    if os.path.exists('README.md'):
        README = open('README.md').read()
    else:
        README = ""  # a placeholder, readme is generated on release
    print("\nVersion : %s\n"%__version__)

    fortran_extnsion = configuration()

    setup(
        name=__project__,
        version=__version__,
        maintainer='UdeM-LBIT',
        description="CoreTracker, A codon reassignment tracker",
        url='https://github.com/UdeM-LBIT/CoreTracker',
        download_url = 'https://github.com/UdeM-LBIT/CoreTracker/tarball/1.0.1',
        author='Emmanuel Noutahi',
        author_email='fmr.noutahi@umontreal.ca',
        scripts = ['bin/coretracker', 'bin/coretranslate', 'bin/corefusion'],
        packages=find_packages(),
        package_data = {
            'coretracker.classifier': ['models/new/*', 'models/old/*'],
            'coretracker.coreutils': ['templates/*'],
        },
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
            'Programming Language :: Python :: 2.7',
        ],

        setup_requires=['numpy'],
        install_requires = [
            'ete3',
            'numpy>= 1.8.1',
            'pandas',
            'scikit-learn',
            'scipy>=0.16.1',
            'biopython>=1.65',
            'matplotlib',
            'WeasyPrint',
            'PyYAML',
        ],
        ext_modules = fortran_extnsion
    )

if __name__ == '__main__':
    setup_package()
