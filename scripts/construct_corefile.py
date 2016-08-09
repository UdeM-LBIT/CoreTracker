#!/usr/bin/env python
from Bio import SeqIO
from coretracker.coreutils import CoreFile
import os, sys
import glob
from collections import defaultdict as ddict
import argparse


def check_already_found(name_map, outkey, inkey, core_inst, new_seq):
	if name_map[outkey][inkey] > 1 :
		for old_seq in core_inst[outkey]:
			if old_seq.id == inkey and old_seq.seq == new_seq.seq:
				sys.stderr.write("Found and removed duplicated case : %s | %s "%(outkey,inkey))
				return False
			elif old_seq.id == inkey :
				raise ValueError('%s present multiple time in %s, with different sequences'%(outkey, inkey))
	return True


# argument parser
if __name__ == '__main__':

	informat = ('fasta', 'gb', 'nexus', 'stockholm', 'clustal')
	parser = argparse.ArgumentParser(
	    description='Convert sequences to corefile format. You need to install biopython')

	parser.add_argument('--informat', dest='inform', choices=informat, default='fasta', help="Dnafile input")
	parser.add_argument('-o', '--outfile',  dest="outfile", default="out.core", help="Outfile file to save files into")

	input_type = parser.add_mutually_exclusive_group(required=True)
	convert_type = parser.add_mutually_exclusive_group(required=True)

	input_type.add_argument('--inputdir', dest="indir", help="Input directory directory")
	input_type.add_argument('-i', '--infile', dest="infiles",  nargs='+', help="List of input files")

	convert_type.add_argument('--gene_to_core', dest='genetype', action='store_true')
	convert_type.add_argument('--genome_to_core', dest='genometype', action='store_true')

	args = parser.parse_args()

	files = []
	if args.infiles:
		files =  [x.strip() for x in args.infiles]

	elif args.indir:
		files = glob.glob(os.path.join(args.indir, "*"))


	if args.inform == 'genbank':
		raise NotImplementedError("Genbank not supported right now")

	name_map = ddict(dict)
	core_inst = ddict(list)
	for f in files:
		fname = os.path.basename(f).split('.')[0]
		fseq = None
		try:
			fseq = SeqIO.parse(f, args.inform)
			if args.genetype:
				core_inst[fname] = fseq
			elif args.genometype:
				for seq in fseq:
					gene = seq.name
					seq.name = seq.id = fname
					name_map[gene][fname] = name_map[gene].get(fname, 0) + 1
					if check_already_found(name_map, gene, fname, core_inst, seq):
						core_inst[gene].append(seq)
		except:
			sys.stderr.write("Wrong format, skipping this file : %s"%(f))

	if core_inst:
		CoreFile.write_corefile(core_inst, args.outfile)
