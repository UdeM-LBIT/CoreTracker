
import argparse
import collections
import glob
import itertools
import json
import numpy as np
import os
import random
import subprocess
import sys
import time

from collections import Counter

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import codonalign
from Bio.codonalign.codonalphabet import default_codon_alphabet, get_codon_alphabet
from Bio.Data import CodonTable
from Bio.codonalign.codonseq import  _get_codon_list, CodonSeq
from Bio import SubsMat
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, generic_protein
from TreeLib import TreeClass
from cStringIO import StringIO
from ete2 import PhyloTree
from settings import *
from PPieChartFace import PPieChartFace, LineFace


__author__ = "Emmanuel Noutahi"
__version__ = "0.1"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"

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

dct = CodonTable.unambiguous_dna_by_id[abs(GENETIC_CODE)]
if GENETIC_CODE < 0:
	for codon in ['CTG', 'CTA', 'CTC', 'CTT']:
		dct.forward_table[codon] = 'L'
	dct.back_table = CodonTable.make_back_table(dct.forward_table, dct.stop_codons[0])

back_table = collections.defaultdict(list)
for aa, codon in zip(dct.forward_table.values(), dct.forward_table.keys()):
	back_table[aa].append(codon)

aa_letters_3to1 = dict((x[1], x[0]) for x in
							aa_letters_1to3.items())
nuc_letters = "ACTG"

class Output(object):

	"""A simple Class to output results
	Either in a file or on stdout
	"""

	def __init__(self, file=None):
		self.file = file

		if(self.file):
			out = open(file, 'w')
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
		sys.stderr.write("%s: %s\n" % (type, message))


class CodonReaData(object):
	"""A representation of a reassignment in a species"""
	
	def __init__(self, aas, consensus, codon_align_dict):
		
		self.aa1, self.aa2 = aas
		self.codon_alignment = codon_align_dict # 
		self.consensus = consensus
		self.reacodons = collections.defaultdict(collections.Counter)
		# find codon usage in the sequence
		self.usedcodons = collections.defaultdict(collections.Counter)
		self.mixtecodon = collections.defaultdict(collections.Counter)
		#print len(codon_align_dict.values()[0])
		self.specs_amino_count = collections.defaultdict(collections.Counter)
		self._spec_codon_usage()

	def _spec_codon_usage(self):
		for i, aa in enumerate(self.consensus):
			for spec in self.codon_alignment.keys():
				spec_codon = self.codon_alignment[spec].seq.get_codon(i)
				spec_aa = dct.forward_table[spec_codon] if '-' not in spec_codon else None
				if(spec_aa):
					# determine aa use in each specie
					if spec_aa == self.aa2 :
						self.specs_amino_count[spec][self.aa2]+=1
					elif spec_aa == self.aa1:
						self.specs_amino_count[spec][self.aa1]+=1

					# potential reassignment counter
					if aa==self.aa1 and spec_aa==self.aa2:
						self.reacodons[spec].update([spec_codon])
					elif aa==self.aa2 and spec_aa==self.aa2:
						self.usedcodons[spec].update([spec_codon])
					elif(spec_aa == self.aa2):
						self.mixtecodon[spec].update([spec_codon])

	
	def __getitem__(self, key):
		return self.get_reacodons(key)

	def get_reacodons(self, specie):
		return self.reacodons[specie]

	def get_mixtecodons(self, specie):
	   return self.mixtecodon[specie]

	def get_usedcodons(self, specie):
		return self.usedcodons[specie]

	def get_usedcodon_freq(self, specie):
		return self.get_usedcodons(specie)*1.0 / self.specs_amino_count[specie][self.aa2]

	def get_reacodon_freq(self,specie):
		return self.get_reacodons(specie)*1.0 / self.specs_amino_count[specie][self.aa2]        

	def get_mixtecodon_freq(self,specie):
		return self.get_mixtecodons(specie)*1.0 / self.specs_amino_count[specie][self.aa2]        

	def get_string(self, specie):
		""" Return a string representation for the specie"""
		return "%s\t%s"%(dict(self.get_reacodons(specie)), dict(self.get_usedcodons(specie)))

	def is_valid(self):
		pass


class NaiveFitch(object):
	"""A NaiveFitch algorithm for finding the most parcimonious solution"""
	def __init__(self, tree, reassigned, ori_aa="", dest_aa="", codon_rea=(None, None)):
		self.id = {}
		self.tree = tree
		self.corr = {'0':ori_aa, '1':dest_aa}
		self.codon_rea_global, self.codon_rea_filtered = codon_rea
		colors=['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f']
		codon_list = back_table[aa_letters_3to1[ori_aa]]
		self.colors = dict(zip(codon_list,colors))
		for leaf in tree:
			if leaf.name in reassigned:
				leaf.add_features(reassigned={1})
				leaf.add_features(rea=self.corr['1'])
				leaf.add_features(state=dest_aa)
			else:
				leaf.add_features(reassigned={0})
				leaf.add_features(rea=self.corr['0'])
				leaf.add_features(state=ori_aa)
		
		self.ori_aa = ori_aa
		self.dest_aa = dest_aa
		self.newick = tree.write(features=['name', 'dist', 'support', 'state'])
		self._bottomup(self.tree)
		self._topdown(self.tree)
	 
	def update_codon_data(codon_rea):
		"""Update codon reassignment data (global and filtered)
		"""
		self.codon_rea_global, self.codon_rea_filtered = codon_rea

	def write_tree(self, outfile):
		"""Export newick to file, if we are going to use R later to find
		ancestral character
		"""
		with open(outfile, 'w') as OUT:
			OUT.write(self.newick)
	
	def is_valid(self):
		tmptree = self.tree.copy()
		for l in tmptree.traverse():
			if not l.is_leaf():
				l.del_feature('reassigned')
				l.del_feature('rea')
			elif 'lost' in l.features and l.lost:
				l.add_features(reassigned={0})
				l.add_features(rea=self.corr['0'])
				l.add_features(state=self.ori_aa)
			elif 'lost' is l.features and not l.lost:
				l.add_features(reassigned={1})
				l.add_features(rea=self.corr['1'])
				l.add_features(state=self.dest_aa)              

		self._bottomup(tmptree)

		for node in tmptree.traverse():
			if node.rea == self.dest_aa:
				return True
		return False

	def _bottomup(self, tree):
		"""Fitch algorithm part 1 : bottom up"""
		for node in tree.traverse('postorder'):
			if not node.is_leaf() :
				intersect = []
				union = set()
				for c in node.get_children():
					union = union | c.reassigned
					intersect.append(c.reassigned)
				intersect = set.intersection(*intersect)
				
				if(intersect):
					node.add_features(reassigned=intersect)
					
					#if len(intersect) == 1:
						# if only one one element in the intersection
						# the children should have that elemenent
					#    for c in node.get_children():
					#        c.reassigned = intersect
				else:
					node.add_features(reassigned=union)

				node.add_features(rea="/".join([self.corr[str(r)] for r in node.reassigned]))

	def _topdown(self, tree):
		"""Fitch algorithm part 2 : top down"""
		pass

	def is_reassigned(cls, node):
		""" return True if a node has undergoned reassignment """
		if "reassigned" in node.features and node.reassigned == {1}:
			return True
		return False

	def get_species_list(self, limit_to_reassigned=True):
		"""Get the species list for which we are almost certain
		that there is a reassignment
		"""
		slist = set()
		if limit_to_reassigned :
			for node in self.tree.traverse():
				if 'reassigned' in node.features and (1 in node.reassigned):
					slist.update(node.get_leaf_names())
		else :
			slist = set(self.tree.get_tree_root().get_leaf_names())
		return slist

	
	def has_codon_data(self):
		return (self.codon_rea_global, self.codon_rea_filtered) != (None, None)


	def render_tree(self, output="", suffix="", pie_size=50, format=IMAGE_FORMAT):

		GRAPHICAL_ACCESS = True

		try:
			from ete2 import TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
		except ImportError, e:
			GRAPHICAL_ACCESS = False

		if not output:
			output = TMP+self.ori_aa+"_to_"+self.dest_aa+suffix+"."+format

		if(GRAPHICAL_ACCESS):
			ts = TreeStyle()
			ts.show_leaf_name = False
			ts.show_branch_length = False
			ts.show_branch_support = False

			rea_style = NodeStyle()
			rea_style["shape"] = "square"
			rea_style["size"] = 7
			rea_style["fgcolor"] = "crimson"
			rea_style["hz_line_type"] = 0

			other_style = NodeStyle()
			other_style["shape"] = "circle"
			other_style["size"] = 7
			other_style["fgcolor"] ="seagreen"
			other_style["node_bgcolor"] = "crimson"
			other_style["hz_line_type"] = 0

			prob_lost = NodeStyle()
			prob_lost["hz_line_type"] = 1
			prob_lost["hz_line_color"] = "#cccccc"
			
			def layout(node):
				N = AttrFace("rea", fsize=10)
				faces.add_face_to_node(N, node, 0, position="branch-right")
				has_count, has_fcount = 0, 0
				if('count' in node.features):
					faces.add_face_to_node(AttrFace("count", fsize=6, fgcolor="firebrick"), node, column=1, position="branch-bottom")
					has_count = node.count

				if('filter_count' in node.features):
					faces.add_face_to_node(AttrFace("filter_count", fsize=6, fgcolor="indigo"), node, column=1, position="branch-top")
					has_fcount = node.filter_count

				if 'lost' in node.features and node.lost:
					faces.add_face_to_node(AttrFace("name", fgcolor="#cccccc",text_suffix="_G"), node, 0, position="aligned")

				elif 'lost' in node.features and not node.lost:
					if(node.is_leaf() and node.reassigned == {0}):
						faces.add_face_to_node(AttrFace("name", fgcolor="seagreen", text_suffix="_V"), node, 0, position="aligned")
					else : 
						faces.add_face_to_node(AttrFace("name", fgcolor="red", fstyle ="italic",text_suffix="_R"), node, 0, position="aligned")

				elif has_count and not has_fcount:
					faces.add_face_to_node(AttrFace("name", fgcolor="gold", text_suffix="_O"), node, 0, position="aligned")

				else:
					faces.add_face_to_node(AttrFace("name"), node, 0, position="aligned")

				if(self.has_codon_data() and node.is_leaf()):

					node_colors = aa_letters_3to1[self.ori_aa]
					spec_codonrea_g = self.codon_rea_global.get_reacodons(node.name)
					spec_codonused_g = self.codon_rea_global.get_usedcodons(node.name)
					
					faces.add_face_to_node(PPieChartFace(spec_codonrea_g.values(), pie_size, pie_size, labels=[], \
						colors=[self.colors[k] for k in spec_codonrea_g.keys()]), node, column=1, position="aligned")
					
					faces.add_face_to_node(PPieChartFace(spec_codonused_g.values(), pie_size, pie_size, labels=[],\
					 colors=[self.colors[k] for k in spec_codonused_g.keys()]), node, column=2, position="aligned")
					
					if(SHOW_MIXTE_CODONS):
						spec_mixtecodon_g = self.codon_rea_global.get_mixtecodons(node.name)
						faces.add_face_to_node(PPieChartFace(spec_mixtecodon_g.values(), pie_size, pie_size, labels=[], \
							colors=[self.colors[k] for k in spec_mixtecodon_g.keys()]), node, column=3, position="aligned")
					
					if(SHOW_FILTERED_CODON_DATA): 

						spec_codonrea_f = self.codon_rea_filtered.get_reacodons(node.name)
						spec_codonused_f = self.codon_rea_filtered.get_usedcodons(node.name)
					
						# add separator 
						faces.add_face_to_node(LineFace(pie_size, pie_size, None), node, column=4, position="aligned")
						
						# add data
						faces.add_face_to_node(PPieChartFace(spec_codonrea_f.values(), pie_size, pie_size, labels=[],\
							colors=[self.colors[k] for k in spec_codonrea_f.keys()]), node, column=5, position="aligned")

						faces.add_face_to_node(PPieChartFace(spec_codonused_f.values(),  pie_size, pie_size, labels=[],\
							colors=[self.colors[k] for k in spec_codonused_f.keys()]), node, column=6, position="aligned")
		
						if(SHOW_MIXTE_CODONS):
							spec_mixtecodon_f = self.codon_rea_filtered.get_mixtecodons(node.name)
							faces.add_face_to_node(PPieChartFace(spec_mixtecodon_f.values(), pie_size, pie_size, labels=[], \
								colors=[self.colors[k] for k in spec_mixtecodon_f.keys()]), node, column=7, position="aligned")
			

			ts.layout_fn = layout
			ts.title.add_face(TextFace(self.ori_aa+" --> "+self.dest_aa, fsize=14), column=0)
			if(self.has_codon_data()):
				for cod, col in self.colors.items():
					ts.legend.add_face(CircleFace((pie_size/3), col), column=0) 
					ts.legend.add_face(TextFace("  "+cod+" ", fsize=8), column=1)
				ts.legend_position = 4

			# Apply node style
			for n in self.tree.traverse():
				if n.reassigned == {1}:
					n.set_style(rea_style)
				elif len(n.reassigned) > 1:
					n.set_style(other_style)
				if 'lost' in n.features and n.lost:
					n.set_style(prob_lost)

			# Save tree as figure
			self.tree.render(output, dpi=400, tree_style=ts)


def get_argsname(command_dict, command_list, prefix=""):
	""" Get the mafft command from the argparse dict : command_dict
	and a command_list
	"""
	return ["%s%s" % (prefix, command) for command in command_list if command_dict[command]]


def clean_spec_list(dna_dict, prot_dict):
	""" Restrict analyses to a list of species, 
	multi-alignment and tree will be pruned to that list"""
	
	if not isinstance(dna_dict, dict):
		dna_dict =  dict((d.id, d) for d in dna_dict)
	if not isinstance(prot_dict, dict):
		prot_dict = dict((p.id, p) for p in prot_dict)
	common_genome = set(dna_dict.keys()).intersection(prot_dict.keys())
	# remove every dna sequence that are not cds
	#print [(x, len(dna_dict[x]), len(prot_dict[x].seq.ungap('-'))+1, len(dna_dict[x])==(len(prot_dict[x].seq.ungap('-'))+1)*3) for x in common_genome ]
	#common_genome = [x for x in common_genome if len(dna_dict[x])==(len(prot_dict[x].seq.ungap('-'))+1)*3 ]
	return common_genome, dict((k, v) for k, v in dna_dict.items() if k in common_genome), prot_dict


def copy_codon_record(record, codon_table=dct, gap_char='-'):
	alphabet = get_codon_alphabet(codon_table, gap_char=gap_char)

	return SeqRecord(CodonSeq(record.seq._data, alphabet=alphabet), id=record.id)


def copy_codon_alignment(codon_alignment, alphabet=default_codon_alphabet):
	return codonalign.CodonAlignment((copy_codon_record(rec) for rec in codon_alignment._records), alphabet=alphabet)


def codon_align(dnaseq, prot_dict, keep_global_pos, keep_filtered_pos, get_dict=True):
	""" Perform a codon alignment based on a protein multiple alignment
	and remove gaped positions
	"""

	# dnaseq should be open with a Bio.SeqIO.parse
	# because duplicated key can cause problems

	# at this step, we are confident all genome in the prot alignement are in the tree

	common_genome, dna_seq, prot_dict = clean_spec_list(dnaseq, prot_dict)
	print common_genome
	#print common_genome
	prot_seq =  [s for s in prot_dict.values() if s.id in common_genome]
	for s in prot_seq:
		s.seq.alphabet = generic_protein

	prot_align = MultipleSeqAlignment(prot_seq, alphabet=alpha)
	# build codon alignment
	codon_alignment = codonalign.build(prot_align, dna_seq, codon_table=dct)
	
	keep_global_pos = sorted(keep_global_pos)
	
	codon_align =  copy_codon_alignment(codon_alignment)
	fcodon_align = copy_codon_alignment(codon_alignment)

	# alignment[:, slice(index_array[0], index_array[0] + 1)]
	for i in xrange(len(codon_align)) :
		codseq = CodonSeq('')
		fcodseq = CodonSeq('')

		for j in keep_global_pos:
			codseq += codon_align[i].seq.get_codon(j)
			if(j in keep_filtered_pos):
				fcodseq += fcodon_align[i].seq.get_codon(j)

		codon_align[i].seq = codseq
		fcodon_align[i].seq = fcodseq

	# fcodon_align does not necessary contains all species
	if get_dict:
		return SeqIO.to_dict(codon_align), SeqIO.to_dict(fcodon_align)
	return codon_align, fcodon_align


def timeit(func):
	"""Utility function to time command line execution"""
	def timed(*args, **kw):
		tstart = time.time()
		result = func(*args, **kw)
		tend = time.time()
		ttime = tend - tstart

		print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
		return ttime, result

	return timed


def executeCMD(cmd, prog):
	"""Execute a command line in the shell"""
	print "\n", cmd
	p = subprocess.Popen(
		cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	out, err = p.communicate()
	Output.error(err)
	print "%s : \n----------------- %s" % (prog, out)
	return err


def filter_align_position(alignment, index_array):
	"""Keep only columns specified by index_array from the alignment"""
	edited_alignment =  alignment[:,:]
	index_array = sorted(index_array)
	# alignment[:, slice(index_array[0], index_array[0] + 1)]
	for seqrec in edited_alignment :
		seq = Seq('', alpha)
		for i in index_array:
			seq += seqrec[i]
		seqrec.seq = seq

	return edited_alignment


def clean_alignment(alignment, characs=['-'], threshold=0.5):
	"""Remove position of alignment which contain character from characs"""
	align_array = np.array([list(rec) for rec in alignment], np.character)
	indel_array = np.where((np.mean(np.in1d(align_array,
											characs).reshape(align_array.shape), axis=0) > threshold) == False)[0].tolist()

	return filter_align_position(alignment, indel_array), indel_array


def filter_alignment(alignment, threshold=0.8, remove_identity=False):
	"""Filter an alignment using threshold as the minimum aa identity per columns"""
	aligninfo = AlignInfo.SummaryInfo(alignment)
	# Smart move : Make a consensus sequence with threshold
	# remove all gap position
	consensus = aligninfo.gap_consensus(threshold=threshold, ambiguous='X')
	cons_array = [i for i, c in enumerate(consensus) if c != 'X']

	if(remove_identity):
		matched_consensus = aligninfo.gap_consensus(threshold=1, ambiguous='X')
		cons_array = [i for i, c in enumerate(
			consensus) if c != 'X' and matched_consensus[i] == 'X']

	return filter_align_position(alignment, cons_array), cons_array


def get_identity(seq1, seq2, ignore_gap=True, percent=True):
	"""Return percent of columns with same aa in a multiple alignment"""
	assert len(seq1) == len(seq2), "Sequence not aligned"
	identity = 0.0
	for i, character in enumerate(seq1):
		if not ignore_gap or character not in ['-', '.']:
			if(character == seq2[i]):
				identity += 1
	return identity * (99.0 * percent + 1.0) / len(seq1)


def convert_tree_to_mafft(tree, seq_order, output, scale):
	"""Convert a tree in a newick format to the mafft matrix format"""
	seqnames = [-1, -1]
	branchlens = [-1, -1]
	for node in tree.traverse("postorder"):
		if node.dist <= DIST_THRES:
			node.dist = 0.0

		if not node.is_leaf():
			left_branch = node.get_child_at(0)
			right_branch = node.get_child_at(1)
			left_branch_leaf = left_branch.get_leaf_name()
			right_branch_leaf = right_branch.get_leaf_name()

			seqnames = [min(map(lambda x: seq_order.index(x), left_branch_leaf)) +
						1,  min(map(lambda x: seq_order.index(x), right_branch_leaf)) + 1, ]
			#seqnames = [seq_order.index(left_branch_leaf.name)+1, seq_order.index(right_branch_leaf.name)+1]
			#seqnames = [left_branch_leaf.name, right_branch_leaf.name]
			branchlens = [
				left_branch.dist * scale, right_branch.dist * scale]

			if seqnames[1] < seqnames[0]:
				seqnames.reverse()
				branchlens.reverse()

			if filter(lambda x: x > 10, branchlens):
				raise ValueError(
					"Your branch length cannot be greater than 10.0. Mafft won't run.")

			output.write("%5d %5d %10.5f %10.5f\n" %
						 (seqnames[0], seqnames[1], branchlens[0], branchlens[1]))


def get_consensus(alignment, threshold):
	"""return consensus, using a defined threshold"""
	aligninfo = AlignInfo.SummaryInfo(alignment)
	consensus = aligninfo.gap_consensus(threshold=threshold, ambiguous='X')
	return consensus


def get_aa_filtered_alignment(consensus, aa):
	"""Get the filtered alignment for an amino acid"""
	aa = aa.upper()
	assert aa in aa_letters, "Amino acid %s not found" % aa
	cons_array = [i for i, c in enumerate(consensus) if c.upper() == aa]
	return cons_array


def realign(msa, tree, enable_mafft):
	"""Realign a msa according to a tree"""
	tmpseq = TMP + "tmp0msa.fasta"
	inseq = TMP + "tmp0nodeseq.fasta"
	seq_ord = []
	for seqrec in msa:
		seqrec._set_seq(seqrec.seq.ungap())
		seq_ord.append(seqrec.description)
	AlignIO.write(msa, open(tmpseq, 'w'), 'fasta')
	out = Output(file=TMP + "tmp0treeout.mafft")
	convert_tree_to_mafft(TreeClass(tree.write(
		features=['name', 'support', 'dist'], format_root_node=True)), seq_ord, out)
	out.close()
	execute_mafft(enable_mafft + " --treein %s %s > %s" %
				  (out.file, tmpseq, inseq))
	msa = AlignIO.read(inseq, 'fasta', alphabet=alpha)
	filelist = glob.glob(TMP + "tmp0*")
	for f in filelist:
		os.remove(f)
	return msa

def is_aligned(seqfile, format):
	try:
		a = AlignIO.read(open(seqfile), format)
	except ValueError:
		return False
	return True

@timeit
def execute_mafft(cmdline):
	"""Execute the mafft command line"""
	executeCMD(cmdline, 'mafft')