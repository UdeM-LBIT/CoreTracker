#!/usr/bin/env python
#
# setup for CoreTracker library packages
#
# use the following to install:
#   python setup.py build
#   python setup.py install
#

import setuptools
import os, sys

if os.path.exists('README.md'):
	README = open('README.md').read()
else:
	README = ""  # a placeholder, readme is generated on release

sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "python")))
from coretracker import __project__, __version__

setuptools.setup(
	name=__project__,
	version=__version__,

	description="CoreTracker, A codon reassignment tracker",
	url='https://github.com/UdeM-LBIT/CoreTracker',
	author='Emmanuel Noutahi',
	author_email='fmr.noutahi@umontreal.ca',
	scripts = ['scripts/coretracker', 'scripts/translate'],

	packages=setuptools.find_packages(exclude=['tests']),

	entry_points={'console_scripts': []},

	keywords="bioinformatics codon reassignment tracker",

	long_description=(README + '\n'),
	license='MIT',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Environment :: Console',
		'Natural Language :: English',
		'Intended Audience :: Developers',
		'Intended Audience :: Education',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License (GPL)',
		'Operating System :: POSIX',
		'Programming Language :: Python',
	],
	install_requires=[
		'ete3',
		'numpy >= 1.8.1',
		'pandas',
		'scikit-learn',
		'scipy',
		'biopython',
		'matplotlib',
		'WeasyPrint',
		'PyYAML',
	]
)
