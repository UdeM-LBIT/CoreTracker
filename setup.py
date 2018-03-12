#!/usr/bin/env python
#
# setup for CoreTracker library packages

from __future__ import absolute_import
from __future__ import print_function
import glob
import os
import sys
from distutils import spawn
try:
    from setuptools import find_packages
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import find_packages

from numpy.distutils.core import Extension as Ext
from numpy.distutils.core import setup

from coretracker import __project__, __version__
sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), "python")))


def configuration(top_path=''):
    fexact_sources = [
        'coretracker/FisherExact/statlib/fexact.pyf',
        'coretracker/FisherExact/statlib/FEXACT.F90',
        'coretracker/FisherExact/statlib/prterr.f'
    ]

    fexact = Ext(name='coretracker.FisherExact.statlib.fexact', sources=[
                 os.path.join(top_path, x) for x in fexact_sources])
    asa159 = Ext(name='coretracker.FisherExact.statlib.asa159', sources=[
                 os.path.join(top_path, 'coretracker/FisherExact/statlib', 'asa159.f90')])
    asa205 = Ext(name='coretracker.FisherExact.statlib.asa205', sources=[
                 os.path.join(top_path, 'coretracker/FisherExact/statlib', 'asa205.f90')])
    return [fexact, asa205, asa159]


def setup_package():
    if os.path.exists('README.rst'):
        README = open('README.rst').read()
    else:
        README = ""  # a placeholder, readme is generated on release
    print("\nVersion : %s\n" % __version__)

    fortran_extnsion = configuration()

    setup(
        name=__project__,
        version=__version__,
        maintainer='Emmanuel Noutahi',
        description="CoreTracker, A codon reassignment tracker",
        url='https://github.com/UdeM-LBIT/CoreTracker',
        download_url='https://github.com/UdeM-LBIT/CoreTracker/archive/%s.tar.gz' % __version__,
        author='Emmanuel Noutahi',
        author_email='fmr.noutahi@umontreal.ca',
        scripts=glob.glob('bin/*'),
        packages=find_packages(),
        package_data={
            'coretracker.classifier': ['models/*/*'],
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
        install_requires=[
            'ete3',
            'numpy>=1.8.2',
            'scipy>=0.16.1',
            'pandas',
            'scikit-learn==0.17',
            'matplotlib>=1.5.1',
            'biopython>=1.65',
            'WeasyPrint',
            'PyYAML',
            'psutil'
        ],
        ext_modules=fortran_extnsion
    )


def binaries_checker():
    """Check if some of the required binaries (alignment) are installed"""
    alignment = ['muscle', 'mafft']
    hmm = ['hmmbuild', 'hmmalign']
    total = [alignment, hmm]
    nfound = []
    for bintype in total:
        for exe in bintype:
            found = spawn.find_executable(exe)
            if not found:
                nfound.append(exe)

    if nfound:
        print("Some binaries where not found : \n-%s" % "\n-".join(nfound))
        answer = 'n'
        try:
            answer = raw_input(
                "Do you still want to continue the installation (y/n) ? ")
        except:
            answer = input(
                "Do you still want to continue the installation (y/n) ? ")

        return answer.lower().startswith('y')

    return True


if __name__ == '__main__':

    cn_continue = binaries_checker()
    if not cn_continue:
        sys.exit(0)

    setup_package()
