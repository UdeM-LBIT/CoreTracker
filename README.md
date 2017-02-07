[![Build Status](https://travis-ci.org/UdeM-LBIT/CoreTracker.svg?branch=master)](https://travis-ci.org/UdeM-LBIT/CoreTracker)
# CoreTracker
CoreTracker detects evidences of codon reassignment from the protein repertoire of a set
of genomes by successively applying different algorithms. It’s therefore a filtering approach that
explore all possible reassignments in every genomes from the input set, and retain only the most promising one.

# Installation

Download the github project and install using pip. I recommend setting a virtual environment through `virtualenv` or `conda`.

To use `CoreTracker`, you should have `PyQt4` installed, which also require `Sip` and `qt`. It's easier to install those two using distribution specific packages.

# Help
After installation, run `coretracker -h` for help.

An example of execution is :
``./coretracker.py -t speciestree.nw -p protein.ali -n nucsequences.core --gapfilter 0.4 --iccontent 0.3  --idfilter 0.5  --norefine --wdir outdir --params param.yml ``

Additionnal parameters could be set using the ``--params`` option. See the provided template (param.yml).
