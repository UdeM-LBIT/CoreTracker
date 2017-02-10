[![Build Status](https://travis-ci.org/UdeM-LBIT/CoreTracker.svg?branch=master)](https://travis-ci.org/UdeM-LBIT/CoreTracker) [![PyPI version](https://badge.fury.io/py/CoreTracker.svg)](https://badge.fury.io/py/CoreTracker)
# CoreTracker
CoreTracker detects evidences of codon reassignment from the protein repertoire of a set
of genomes by successively applying different algorithms. It’s therefore a filtering approach that
explore all possible reassignments in every genomes from the input set, and retain only the most promising one.

# Installation

First install the system dependencies which include `gfortran`, `PyQt4` `muscle`, `mafft` and `hmmer`. `PyQt4` also require `Sip` and `qt`. It's easier to install those two using distribution specific packages. You can now download the github project and install using `python setup.py install` or pip (the package is also available on pypi). I recommend setting a virtual environment through `virtualenv`.

Alternatively, you can also install it with `conda`, which is the easiest way. 


# Basic Help
After installation, run `coretracker -h` for help.

An example of execution is :
``./coretracker.py -t speciestree.nw -p protein.ali -n nucsequences.core --gapfilter 0.4 --iccontent 0.3  --idfilter 0.5  --norefine --wdir outdir --params param.yml ``

Additionnal parameters could be set using the ``--params`` option. See the provided template (param.yml).

# Detailled information

Detailles information about the package and tutorials are available here ==> [http://udem-lbit.github.io/CoreTracker/](http://udem-lbit.github.io/CoreTracker/)

# TODO

    - consensus class to not ignore column with same aa frequency
    - full integration of rna editing part in the software
    - code quality improvement and cleaning classifier part 


