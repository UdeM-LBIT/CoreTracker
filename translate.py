#!/usr/bin/env python
from Bio import SeqIO
from utils import CoreFile, purge_directory
from Bio import AlignIO
import sys
import warnings
import os
warnings.filterwarnings("ignore") 
import argparse

# argument parser
parser = argparse.ArgumentParser(description='Translator, translate nucleotide alignment. You need to install hmmalign, muscle or mafft and biopython')
parser.add_argument('--wdir', '--outdir', dest="outdir", default="trans_output", help="Working directory")
parser.add_argument('--nuc', '--input', '-n', dest='dnafile', required=True, help="Dnafile input")
parser.add_argument('--gcode', type=int, default=1, dest='gcode', help="Genetic code to use for translation. Default value is 1")
parser.add_argument('--prog', default="mafft", dest='prog', choices=['mafft', 'muscle'] ,help="Genetic code to use for translation. Default value is 1")
parser.add_argument('--align', dest='align', action='store_true', help="Whether we should align or not")
parser.add_argument('--noclean', dest='noclean', action='store_false', help="Whether we should clean unnecessary file or not")

args = parser.parse_args()
purge_directory(args.outdir)

dnafile = args.dnafile
gcode = args.gcode
align = args.align

prog = "mafft --auto"
if args.prog =='muscle':
	prog = 'muscle'

# translate here
translated_prot = CoreFile.translate(dnafile, gcode)
alignment = {}

if align:
	for gene, seqs in translated_prot.items():
		al = CoreFile._align(seqs, prog, None, 1.0, args.outdir)
		refine =  CoreFile._refine(al, 9999, args.outdir, clean=args.noclean)
        alignment[gene] = al

	CoreFile.write_corefile(alignment, os.path.join(args.outdir, "prot_aligned.core"))
else :
	CoreFile.write_corefile(translated_prot, os.path.join(args.outdir, "prot.core"))