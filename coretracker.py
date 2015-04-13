#! /usr/bin/env python
from __future__ import division

import argparse
import collections
import itertools
import json, os, glob
import numpy as np
import pandas as pd
import subprocess
import sys, random, time

from collections import Counter

from Bio import AlignIO
from Bio import Alphabet
from Bio import SeqIO
from Bio import SubsMat
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from TreeLib import TreeClass
from cStringIO import StringIO
from ete2 import PhyloTree
from ete2 import TreeStyle
from settings import *


__author__ = "Emmanuel Noutahi"
__version__ = "0.1"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"



if not os.path.exists(TMP):
	os.makedirs(TMP)

# hardcorded settings variables
alpha = Alphabet.Gapped(IUPAC.protein)
DIST_THRES = 1e-10
aa_letters = "ACDEFGHIKLMNPQRSTVWY"
aa_letters_1to3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp',
    'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met',
    'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp',
    'Y': 'Tyr',
}
nuc_letters = "ACTG"

MAFFT_AUTO_COMMAND = ['linsi', 'ginsi', 'einsi', 'fftnsi', 'fftns', 'nwnsi', 'nwns']
MAFFT_DETAILLED_COMMAND = ['auto', 'retree', 'maxiterate', 'nofft', 'memsave', 'parttree']
# This is the multiple alignment used to estimate the accepted replacement matrix 

class Output(object):
	"""A simple Class to output results
	Either in a file or on stdout
	"""
	def __init__(self, file=None):
		self.file = file

		if(self.file):
			out = open(file,'w')
			self.out = out
		else:
			self.out = sys.stdout

	def write(self, line):
		self.out.write(line)
	
	def close(self):
		if self.out is not sys.stdout:
			self.out.close()

	@staticmethod
	def error(message, type="Error"):
		sys.stderr.write("%s: %s\n" % (type,message))


def get_argsname(command_dict, command_list, prefix=""):
	""" Get the mafft command from the argparse dict : command_dict
	and a command_list
	"""
	return ["%s%s"%(prefix, command) for command in command_list if command_dict[command]]


def timeit(func):
	"""Utility function to time command line execution"""
	def timed(*args, **kw):
		tstart = time.time()
		result = func(*args, **kw)
		tend = time.time()
		ttime = tend-tstart

		print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
		return ttime, result

	return timed


def executeCMD(cmd, prog):
	"""Execute a command line in the shell"""
	print "\n", cmd
	p=subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	out, err = p.communicate()
	Output.error(err)
	print "%s : \n----------------- %s"%(prog, out)
	return err


def filter_align_position(alignment, index_array):
	"""Remove columns specified by index_array from the alignment"""
	edited_alignment = alignment[:,slice(index_array[0], index_array[0]+1)]

	for i in index_array[1:]:
		edited_alignment += alignment[:, slice(i, i+1)]

	return edited_alignment


def clean_alignment(alignment, characs=['-'], threshold=0.5):
	"""Remove position of alignment which contain character from characs"""
	align_array = np.array([list(rec) for rec in alignment], np.character)
	indel_array = np.where((np.mean(np.in1d(align_array, \
		characs).reshape(align_array.shape), axis=0) > threshold)==False)[0].tolist()

	return filter_align_position(alignment, indel_array)


def filter_alignment(alignment, threshold=0.8, remove_identity=False):
	"""Filter an alignment using threshold as the minimum aa identity per columns"""
	aligninfo = AlignInfo.SummaryInfo(alignment)
	# Smart move : Make a consensus sequence with threshold
	# remove all gap position 
	consensus = aligninfo.gap_consensus(threshold=threshold, ambiguous='X')
	cons_array = [i for i, c in enumerate(consensus) if c != 'X']

	if(remove_identity):
		matched_consensus =  aligninfo.gap_consensus(threshold=1, ambiguous='X')
		cons_array = [i for i, c in enumerate(consensus) if c != 'X' and matched_consensus[i]=='X']

	return filter_align_position(alignment, cons_array)


def get_identity(seq1, seq2, ignore_gap=True, percent=True):
	"""Return percent of columns with same aa in a multiple alignment"""
	assert len(seq1)==len(seq2), "Sequence not aligned"
	identity = 0.0
	for i, character in enumerate(seq1):
		if not ignore_gap or character not in ['-', '.']:
			if(character == seq2[i]):
				identity += 1
	return identity*(99.0*percent + 1.0)/len(seq1)


def convert_tree_to_mafft(tree, seq_order, output):
	"""Convert a tree in a newick format to the mafft matrix format"""
	seqnames = [-1,-1]
	branchlens = [-1,-1]
	for node in tree.traverse("postorder"):
		if node.dist <= DIST_THRES:
			node.dist = 0.0
		
		if not node.is_leaf():
			left_branch = node.get_child_at(0)
			right_branch = node.get_child_at(1)
			left_branch_leaf = left_branch.get_leaf_name()
			right_branch_leaf = right_branch.get_leaf_name()

			seqnames = [min(map(lambda x : seq_order.index(x), left_branch_leaf)) + 1,  min(map(lambda x : seq_order.index(x), right_branch_leaf)) + 1,]
			#seqnames = [seq_order.index(left_branch_leaf.name)+1, seq_order.index(right_branch_leaf.name)+1]
			#seqnames = [left_branch_leaf.name, right_branch_leaf.name]
			branchlens = [left_branch.dist*args.scale, right_branch.dist*args.scale]

			if seqnames[1] < seqnames[0]:
				seqnames.reverse()
				branchlens.reverse()
			
			if filter(lambda x: x>10, branchlens):
				raise ValueError("Your branch length cannot be greater than 10.0. Mafft won't run.")

			output.write("%5d %5d %10.5f %10.5f\n"%(seqnames[0], seqnames[1], branchlens[0], branchlens[1]))


def get_aa_filtered_alignment(alignment, aa, threshold):
	"""Get the filtered alignment for an amino acid"""
	aa = aa.upper()
	assert aa in aa_letters, "Amino acid %s not found"%aa
	aligninfo = AlignInfo.SummaryInfo(alignment)
	consensus = aligninfo.gap_consensus(threshold=threshold, ambiguous='X')
	cons_array = [i for i, c in enumerate(consensus) if c.upper() == aa]
	return filter_align_position(alignment, cons_array)


def realign(msa, tree, enable_mafft):
	"""Realign a msa according to a tree"""
	tmpseq = TMP+"tmp0msa.fasta"
	inseq = TMP+"tmp0nodeseq.fasta"
	seq_ord = []
	for seqrec in msa :
		seqrec._set_seq(seqrec.seq.ungap())
		seq_ord.append(seqrec.description)
	AlignIO.write(msa, open(tmpseq, 'w'), 'fasta')
	out = Output(file=TMP+"tmp0treeout.mafft")
	convert_tree_to_mafft(TreeClass(tree.write(features=['name', 'support', 'dist'], format_root_node=True)),seq_ord, out)
	out.close()
	execute_mafft(enable_mafft+" --treein %s %s > %s"%(out.file, tmpseq, inseq))
	msa = AlignIO.read(inseq, 'fasta', alphabet=alpha)
	filelist = glob.glob(TMP+"tmp0*")
	for f in filelist :
		os.remove(f)
	return msa


@timeit
def execute_mafft(cmdline):
	"""Execute the mafft command line"""
	executeCMD(cmdline, 'mafft')


if __name__ == '__main__':
		
	# argument parser
	parser = argparse.ArgumentParser(description='Convert newick tree format to mafft format')
	parser.add_argument('-s', '--scale', type=float, default=1.0, dest='scale', help="Scale to compute the branch format")
	parser.add_argument('-o', '--output', type=argparse.FileType('w+'), dest="output", help="Output file" )
	parser.add_argument('-r', '--resample', type=int, help="For debug and memory purpose. Choose only x sequence")
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	parser.add_argument('--excludegap', type=float, dest='excludegap', help="Remove position with gap from the alignment, using excludegap as threshold. The absolute values are taken")
	parser.add_argument('--idfilter', type=float, default=0.8, dest='idfilter', help="Conserve only position with at least idfilter residue identity")
	parser.add_argument('--verbose', '-v', action='store_true', dest="verbose", help="Verbosity level")


	mafft_group = parser.add_mutually_exclusive_group()
	mafft_group.add_argument('--linsi', dest='linsi', action='store_true',help="L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)")
	mafft_group.add_argument('--ginsi', dest='ginsi', action='store_true', help="G-INS-i (suitable for sequences of similar lengths; recommended for <200 sequences; iterative refinement method incorporating global pairwise alignment information)")
	mafft_group.add_argument('--einsi', dest='einsi', action='store_true', help="E-INS-i (suitable for sequences containing large unalignable regions; recommended for <200 sequences)")
	mafft_group.add_argument('--fftnsi', dest='fftnsi', action='store_true', help="FFT-NS-i (iterative refinement method; two cycles only)")
	mafft_group.add_argument('--fftns', dest='fftns', action='store_true', help="FFT-NS-2 (fast; progressive method)")
	mafft_group.add_argument('--nwnsi', dest='nwnsi', action='store_true', help="NW-NS-i (iterative refinement method without FFT approximation; two cycles only)")
	mafft_group.add_argument('--nwns', dest='nwns', action='store_true', help="NW-NS-2 (fast; progressive method without the FFT approximation)")

	parser.add_argument('-t', '--intree', dest="tree", help='Input specietree in newick format', required=True)
	parser.add_argument('-a', '--alignment', dest='seq', help="The sequence input in fasta format", required=True)

	subparsers = parser.add_subparsers(title="Mafft options", description="Use mafft command to perform alignment")
	subparsers.required = False 
	subparsers.dest = 'command'

	parser_align = subparsers.add_parser('align', help='Perform alignment with mafft: The following additional argument will be mandatory. (See mafft doc) "--auto", "--retree", "--maxiterate", "--nofft", "--memsave", "--parttree"' )
	parser_align.add_argument('--auto', dest='auto', action='store_true', help="If unsure which option to use, try this option")
	parser_align.add_argument('--retree', dest='retree', type=int, help="Guide tree is built number times in the progressive stage")
	parser_align.add_argument('--maxiterate', dest='maxiterate', type=int, help="number cycles of iterative refinement are performed")
	parser_align.add_argument('--nofft', dest='nofft', action='store_true', help="Do not use FFT approximation in group-to-group alignment.")
	parser_align.add_argument('--memsave', dest='memsave', action='store_true', help="Use the Myers-Miller (1988) algorithm.")
	parser_align.add_argument('--parttree', dest='parttree', action='store_true', help="Use a fast tree-building method (PartTree, Katoh and Toh 2007) with the 6mer distance.")

	args= parser.parse_args()

	# Check mafft command input
	mafft_short = get_argsname(args.__dict__, MAFFT_AUTO_COMMAND)
	mafft_detail = get_argsname(args.__dict__, MAFFT_DETAILLED_COMMAND, prefix="--")
	enable_mafft = False

	if(mafft_detail and mafft_short):
		parser.error("You cannot use the shortcuts and the detailled parameters for mafft at the same time")

	if(mafft_detail or mafft_short):
		enable_mafft = "mafft " + " ".join(mafft_detail) if mafft_detail else "".join(mafft_short)
		args.output = args.output or TMP+"tree.mafft"

	# Output stream setting 
	output = Output(args.output)

	# Manage fasta sequence input
	seq_order = [] #sequence id in the multi-alignment fasta file
	fasta_sequences = SeqIO.parse(args.seq, "fasta")
	for fasta in fasta_sequences:
		seq_order.append(fasta.description)

	#check duplicated sequence name in the alignment file
	assert len(set(seq_order))==len(seq_order), \
		"Duplicate found for sequences :%s" %([x for x, y in collections.Counter(seq_order).items() if y > 1])

	# Use a portion of the sequence for test if asked by the user
	if args.resample:
		seq_order = random.sample(seq_order, args.resample)

	original_len = len(seq_order)

	# load tree
	specietree = TreeClass(args.tree)
	leave_names = specietree.get_leaf_name()

	#Check for sequence not in the tree and correct
	seq_order[:] = [seq for seq in seq_order if seq in leave_names]
	if len(seq_order) != original_len :
		Output.error("Sequence not in tree found, attempt to correct... ", "Warning")
		record_dict = SeqIO.index(args.seq, 'fasta')
		args.seq = args.seq+"_"
		with open(args.seq, 'w') as FASTAOUT:
			for seq in seq_order:
				FASTAOUT.write(record_dict[seq].format("fasta"))

	assert (set(seq_order) <= set(leave_names)), \
		"Some sequences id are not in your tree, Aborting... %s"% (set(seq_order) - set(leave_names))

	# prune tree to sequence list
	specietree.prune(seq_order)

	#Convert tree to mafft format
	convert_tree_to_mafft(specietree, seq_order, output)
	output.close()

	# execute mafft  
	if(enable_mafft):
		execute_mafft(enable_mafft + " --treein %s %s > %s"%(output.file, args.seq, MAFFT_OUTPUT))

	# Reload mafft alignment and filter alignment (remove gap positions and positions not conserved)
	alignment = AlignIO.read(MAFFT_OUTPUT, 'fasta', alphabet=alpha)

	if(args.excludegap):
		alignment = clean_alignment(alignment, threshold=(abs(args.excludegap)<=1 or 0.01)*abs(args.excludegap))
		AlignIO.write(alignment, open(MAFFT_OUTPUT+"_ungapped", 'w'), 'fasta')

	filtered_alignment = alignment

	# Already enabled by default in the arguments list
	if(args.idfilter):
		filtered_alignment = filter_alignment(alignment, threshold=(abs(args.idfilter)<=1 or 0.01)*abs(args.idfilter))
		AlignIO.write(filtered_alignment, open(MAFFT_OUTPUT+"_filtered", 'w'), 'fasta')
		identity_removed_alignment = filter_alignment(filtered_alignment, remove_identity=True, threshold=(abs(args.idfilter)<=1 or 0.01)*abs(args.idfilter))
		AlignIO.write(identity_removed_alignment, open(MAFFT_OUTPUT+"_matchremoved", 'w'), 'fasta')


	## Get the substitution matrice at each node
	sptree = PhyloTree(specietree.write(), alignment=filtered_alignment.format("fasta"), alg_format="fasta")
	
	# Compute expected frequency from the entire alignment and update at each node
	armseq = alignment
	
	# We could use another alignment to compute those frequency if it's set
	if(ARM_SEQUENCE):
		armseq = AlignIO.read(ARM_SEQUENCE, 'fasta', alphabet=alpha)

	summary_info = AlignInfo.SummaryInfo(armseq)
	acc_rep_mat = SubsMat.SeqMat(summary_info.replacement_dictionary())
	obs_freq_mat = SubsMat._build_obs_freq_mat(acc_rep_mat)
	exp_freq_table = SubsMat._exp_freq_table_from_obs_freq(obs_freq_mat)

	for node in sptree.traverse():
		if not node.is_leaf():
			spec_under_node = node.get_leaf_names()
			node_alignment = MultipleSeqAlignment([record[:] for record in alignment if record.id in spec_under_node], alphabet=alpha)
			
			if(REALIGN_AT_EACH_NODE and not node.is_root()):
				node_alignment = realign(node_alignment, node, enable_mafft)

			node_align_info =  AlignInfo.SummaryInfo(node_alignment)
			replace_info = node_align_info.replacement_dictionary()
			node_align_arm = SubsMat.SeqMat(replace_info)
			
			# add alignment content at root		
			if(node.is_root()):
				total_ic_content = node_align_info.information_content(chars_to_ignore=['X'])

				node.add_feature('ic', node_align_info.ic_vector)
		
			node_sub_matrix = SubsMat.make_log_odds_matrix(node_align_arm, exp_freq_table=exp_freq_table)

			node.add_feature('submatrix', node_sub_matrix)

			if(args.verbose):
				print "Sequences under node"
				print spec_under_node
				print "ARM matrix"
				print node_align_arm
				print "SUB matrix"
				print node_sub_matrix

	# Compute sequence identity in global and filtered alignment
	seq_names = [rec.id for rec in alignment]
	number_seq = len(seq_names)
	sim_json = collections.defaultdict(list)
	aa_json = collections.defaultdict(list)
	expected_freq = collections.defaultdict(float)
	#alignment length
	af_length = len(filtered_alignment[0])
	ag_length = len(alignment[0])
	count_max = 0
	id_matrix = np.zeros(shape=(number_seq, number_seq))
	global_paired_id = pd.DataFrame(id_matrix, index=seq_names, columns=seq_names)
	filter_paired_id = global_paired_id.copy()
	
	for i in xrange(number_seq):
		#Get using frequence of each aa per sequence
		count1 = Counter(alignment[i])
		count2 = Counter(filtered_alignment[i])
		for aa in aa_letters:
			expected_freq[aa_letters_1to3[aa]] = exp_freq_table[aa]
			global_val = count1[aa]/(ag_length*exp_freq_table[aa])
			filtered_val = count2[aa]/(af_length*exp_freq_table[aa])
			aa_json[aa_letters_1to3[aa]].append({'global': global_val, 'filtered': filtered_val, "gene":seq_names[i]})
			count_max = max(filtered_val, global_val, count_max)

		for j in xrange(i+1):
			if(j==i):
				global_paired_id.iloc[i,j] = filter_paired_id.iloc[i,j] = 1.0
			else :
				global_paired_id.iloc[i,j] = get_identity(alignment[i], alignment[j])
				filter_paired_id.iloc[i,j] = get_identity(filtered_alignment[i], filtered_alignment[j])
			sim_json[seq_names[i]].append({"global" : global_paired_id.iloc[i,j], "filtered" : filter_paired_id.iloc[i,j], "gene" : seq_names[j]})
			# do not add identity to itself twice
			if(i!=j):
				sim_json[seq_names[j]].append({"global" : global_paired_id.iloc[i,j], "filtered" : filter_paired_id.iloc[i,j], "gene" : seq_names[i]})
	
	sptree.write(outfile=TMP+"phylotree.nw")
	# dumping in json to reload with the web interface using d3.js
	with open(TMP+"similarity.json", "w") as outfile:
		json.dump(sim_json, outfile, indent=4)
	with open(TMP+"aafrequency.json", "w") as outfile:
		json.dump({"AA" : aa_json, "EXP": expected_freq, "MAX" : count_max}, outfile, indent=4)
	