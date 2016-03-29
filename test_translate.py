#!/usr/bin/env python
from Bio import SeqIO
from coretracker.coreutils import SequenceLoader, CoreFile
import coretracker.coreutils.utils as utils
from Bio import AlignIO
import sys
import warnings
import os
import glob
import logging
warnings.filterwarnings("ignore")
import argparse


logging.basicConfig(level=logging.DEBUG)

# argument parser
parser = argparse.ArgumentParser(
    description='Translator, translate nucleotide alignment. You need to install hmmalign, muscle or mafft and biopython')
parser.add_argument('--wdir', '--outdir', dest="outdir",
                    default="trans_output", help="Working directory")
parser.add_argument('--nuc', '--input', '-n', dest='dnafile',
                    required=True, help="Dnafile input")
parser.add_argument('--gcode', type=int, default=1, dest='gcode',
                    help="Genetic code to use for translation. Default value is 1")
parser.add_argument('--prog', default="mafft", dest='prog', choices=[
                    'mafft', 'muscle'], help="Genetic code to use for translation. Default value is 1")
parser.add_argument('--align', dest='align', action='store_true',
                    help="Whether we should align or not")
parser.add_argument('--refine', dest='refine', action='store_true',
                    help="Whether we should refine the alignment with hmm")
parser.add_argument('--noclean', dest='noclean', action='store_false',
                    help="Whether we should clean unnecessary file or not")
parser.add_argument('--hmmdir', dest='hmmdir',
                    help="Link a directory with hmm files for alignment. Each hmmfile should be named in the following format : genename.hmm")


args = parser.parse_args()
utils.purge_directory(args.outdir)

dnafile = args.dnafile
gcode = args.gcode
align = args.align

prog = "mafft --auto"
if args.prog == 'muscle':
    prog = 'muscle'


hmmfiles = {}
if args.hmmdir:
    try:
        hfiles = glob.glob(os.path.join(args.hmmdir, '*'))
        for f in hfiles:
            genename = os.path.basename(f).split('.hmm')[0]
            hmmfiles[genename] = f
    except Exception as e:
        print e
        pass

# filter the sequences to remove genes where all the sequences do not have a len
# that is a multiple of 3
tmp_core_inst = CoreFile(dnafile, alphabet='nuc')
# do not treat as missing sequence just remove the entire gene
core_inst = {}
for gene, seqs in tmp_core_inst.items():
    good_data = []
    disp_data = []
    for seq in seqs:
        if len(seq) % 3 != 0:
            disp_data.append(seq)
        else:
            good_data.append(seq)
    if not disp_data:
        core_inst[gene] = good_data
        print("==> OK (%s)"%gene)
    else:
        print("==> Fail (%s), %d specs have frame-shifting : \n(%s)"%(gene, len(disp_data), ", ".join([x.name for x in disp_data])))

if not core_inst:
    sys.exit("After removing frameshifting, there isn't ny remaining genes")

# translate here
translated_prot = SequenceLoader.translate(core_inst, gcode)
alignment = {}

if align:
    for (gene, seqs) in translated_prot.items():
        al = SequenceLoader._align(seqs, prog, None, 1.0, args.outdir)
        if args.refine:
            al = SequenceLoader._refine(al, 9999, args.outdir, loop=10,
                                  clean=args.noclean, hmmfile=hmmfiles.get(gene, None))
        alignment[gene] = al

    CoreFile.write_corefile(alignment, os.path.join(
        args.outdir, "prot_aligned.core"))
else:
    CoreFile.write_corefile(
        translated_prot, os.path.join(args.outdir, "prot.core"))
