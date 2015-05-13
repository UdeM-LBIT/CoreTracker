#! /usr/bin/env python
from __future__ import division

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

from multiprocessing import Pool

from collections import Counter

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SubsMat
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from TreeLib import TreeClass
from cStringIO import StringIO
from ete2 import PhyloTree
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

aa_letters_3to1 = dict((x[1], x[0]) for x in
                            aa_letters_1to3.items())
nuc_letters = "ACTG"

MAFFT_AUTO_COMMAND = [
    'linsi', 'ginsi', 'einsi', 'fftnsi', 'fftns', 'nwnsi', 'nwns']
MAFFT_DETAILLED_COMMAND = [
    'auto', 'retree', 'maxiterate', 'nofft', 'memsave', 'parttree', 'leavegappyregion']
# This is the multiple alignment used to estimate the accepted replacement
# matrix


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


class NaiveFitch(object):
    """A NaiveFitch algorithm for finding the most parcimonious solution"""
    def __init__(self, tree, reassigned, ori_aa="", dest_aa=""):
        self.id = {}
        self.tree = tree
        self.corr = {'0':ori_aa, '1':dest_aa}
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
     

    def write_tree(self, outfile):
        """Export newick to file, if we are going to use R later to find
        ancestral character
        """
        with open(outfile, 'w') as OUT:
            OUT.wrte(self.newick)
    
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

    def render_tree(self, output="", suffix=""):

        GRAPHICAL_ACCESS = True

        try:
            from ete2 import TreeStyle, NodeStyle, faces, AttrFace, TextFace
        except ImportError, e:
            GRAPHICAL_ACCESS = False

        if not output:
            output = TMP+self.ori_aa+"_to_"+self.dest_aa+suffix+".pdf"

        if(GRAPHICAL_ACCESS):
            ts = TreeStyle()
            ts.show_leaf_name = False

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
                if('count' in node.features):
                    faces.add_face_to_node(AttrFace("count", fsize=6), node, column=1, position="branch-bottom")
                if('filter_count' in node.features):
                    faces.add_face_to_node(AttrFace("filter_count", fsize=5, fgcolor="indigo"), node, column=1, position="branch-top")

                if 'lost' in node.features and node.lost:
                    faces.add_face_to_node(AttrFace("name", fgcolor="#cccccc",text_suffix="_G"), node, 0, position="aligned")
                elif 'lost' in node.features and not node.lost:
                    if(node.is_leaf() and node.reassigned == {0}):
                        faces.add_face_to_node(AttrFace("name", fgcolor="seagreen", text_suffix="_V"), node, 0, position="aligned")
                    else : 
                        faces.add_face_to_node(AttrFace("name", fgcolor="red", fstyle ="italic",text_suffix="_R"), node, 0, position="aligned")

                else:
                    faces.add_face_to_node(AttrFace("name"), node, 0, position="aligned")

            ts.layout_fn = layout
            ts.legend.add_face(TextFace(self.ori_aa+" --> "+self.dest_aa, fsize=14), column=1)
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
    """Remove columns specified by index_array from the alignment"""
    edited_alignment =  alignment[:,:]
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

    return filter_align_position(alignment, indel_array)


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

    return filter_align_position(alignment, cons_array)


def get_identity(seq1, seq2, ignore_gap=True, percent=True):
    """Return percent of columns with same aa in a multiple alignment"""
    assert len(seq1) == len(seq2), "Sequence not aligned"
    identity = 0.0
    for i, character in enumerate(seq1):
        if not ignore_gap or character not in ['-', '.']:
            if(character == seq2[i]):
                identity += 1
    return identity * (99.0 * percent + 1.0) / len(seq1)


def convert_tree_to_mafft(tree, seq_order, output):
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
                left_branch.dist * args.scale, right_branch.dist * args.scale]

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


if __name__ == '__main__':

    # argument parser
    parser = argparse.ArgumentParser(
        description='Convert newick tree format to mafft format')
    parser.add_argument('-s', '--scale', type=float, default=1.0,
                        dest='scale', help="Scale to compute the branch format")
    parser.add_argument(
        '-o', '--output', type=argparse.FileType('w+'), dest="output", help="Output file")
    parser.add_argument('-r', '--resample', type=int,
                        help="For debug and memory purpose. Choose only x sequence")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('--excludegap', type=float,  default=0.6, dest='excludegap',
                        help="Remove position with gap from the alignment, using excludegap as threshold. The absolute values are taken")
    parser.add_argument('--idfilter', type=float, default=0.8, dest='idfilter',
                        help="Conserve only position with at least idfilter residue identity")
    parser.add_argument(
        '--verbose', '-v', action='store_true', dest="verbose", help="Verbosity level")

    parser.add_argument(
        '--debug', action='store_true', dest="debug", help="Print debug infos")

    parser.add_argument(
        '--sfx', dest="sfx", default="", help="PDF rendering suffix to differentiate runs.")

    mafft_group = parser.add_mutually_exclusive_group()
    mafft_group.add_argument('--linsi', dest='linsi', action='store_true',
                             help="L-INS-i (probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information)")
    mafft_group.add_argument('--ginsi', dest='ginsi', action='store_true',
                             help="G-INS-i (suitable for sequences of similar lengths; recommended for <200 sequences; iterative refinement method incorporating global pairwise alignment information)")
    mafft_group.add_argument('--einsi', dest='einsi', action='store_true',
                             help="E-INS-i (suitable for sequences containing large unalignable regions; recommended for <200 sequences)")
    mafft_group.add_argument('--fftnsi', dest='fftnsi', action='store_true',
                             help="FFT-NS-i (iterative refinement method; two cycles only)")
    mafft_group.add_argument(
        '--fftns', dest='fftns', action='store_true', help="FFT-NS-2 (fast; progressive method)")
    mafft_group.add_argument('--nwnsi', dest='nwnsi', action='store_true',
                             help="NW-NS-i (iterative refinement method without FFT approximation; two cycles only)")
    mafft_group.add_argument('--nwns', dest='nwns', action='store_true',
                             help="NW-NS-2 (fast; progressive method without the FFT approximation)")

    parser.add_argument('-t', '--intree', dest="tree",
                        help='Input specietree in newick format', required=True)
    parser.add_argument('-a', '--alignment', dest='seq',
                        help="The sequence input in fasta format", required=True)

    subparsers = parser.add_subparsers(
        title="Mafft options", description="Use mafft command to perform alignment")
    subparsers.required = False
    subparsers.dest = 'command'

    parser_align = subparsers.add_parser(
        'align', help='Perform alignment with mafft: The following additional argument will be mandatory. (See mafft doc) "--auto", "--retree", "--maxiterate", "--nofft", "--memsave", "--parttree", "--leavegappyregion"')
    parser_align.add_argument(
        '--auto', dest='auto', action='store_true', help="If unsure which option to use, try this option")
    parser_align.add_argument('--retree', dest='retree', type=int,
                              help="Guide tree is built number times in the progressive stage")
    parser_align.add_argument('--maxiterate', dest='maxiterate', type=int,
                              help="number cycles of iterative refinement are performed")
    parser_align.add_argument('--nofft', dest='nofft', action='store_true',
                              help="Do not use FFT approximation in group-to-group alignment.")
    parser_align.add_argument(
        '--memsave', dest='memsave', action='store_true', help="Use the Myers-Miller (1988) algorithm.")
    parser_align.add_argument('--parttree', dest='parttree', action='store_true',
                              help="Use a fast tree-building method (PartTree, Katoh and Toh 2007) with the 6mer distance.")

    parser_align.add_argument('--leavegappyregion', dest='leavegappyregion', action='store_true',
                              help="Leave gappy region. The default gap scoring scheme has been changed in version 7.110 (2013 Oct). It tends to insert more gaps into gap-rich regions than previous versions. To disable this change, add the --leavegappyregion option")

    args = parser.parse_args()

    # Check mafft command input
    mafft_short = get_argsname(args.__dict__, MAFFT_AUTO_COMMAND)
    mafft_detail = get_argsname(
        args.__dict__, MAFFT_DETAILLED_COMMAND, prefix="--")
    enable_mafft = False

    if(mafft_detail and mafft_short):
        parser.error(
            "You cannot use the shortcuts and the detailled parameters for mafft at the same time")

    if(mafft_detail or mafft_short):
        enable_mafft = "mafft " + \
            " ".join(mafft_detail) if mafft_detail else "".join(mafft_short)
        args.output = args.output or TMP + "tree.mafft"

    # Output stream setting
    mtoutput = Output(args.output)

    # Manage fasta sequence input
    seq_order = []  # sequence id in the multi-alignment fasta file
    fasta_sequences = SeqIO.parse(args.seq, "fasta")
    for fasta in fasta_sequences:
        seq_order.append(fasta.description)

    # check duplicated sequence name in the alignment file
    assert len(set(seq_order)) == len(seq_order), \
        "Duplicate found for sequences :%s" % (
            [x for x, y in collections.Counter(seq_order).items() if y > 1])

    # Use a portion of the sequence for test if asked by the user
    if args.resample:
        seq_order = random.sample(seq_order, args.resample)

    original_len = len(seq_order)

    # load tree
    specietree = TreeClass(args.tree)
    leave_names = specietree.get_leaf_name()

    # Check for sequence not in the tree and correct
    seq_order[:] = [seq for seq in seq_order if seq in leave_names]
    if len(seq_order) != original_len:
        Output.error(
            "Sequence not in tree found, attempt to correct... ", "Warning")
        record_dict = SeqIO.index(args.seq, 'fasta')
        args.seq = args.seq + "_"
        with open(args.seq, 'w') as FASTAOUT:
            for seq in seq_order:
                FASTAOUT.write(record_dict[seq].format("fasta"))

    assert (set(seq_order) <= set(leave_names)), \
        "Some sequences id are not in your tree, Aborting... %s" % (
            set(seq_order) - set(leave_names))

    # prune tree to sequence list
    specietree.prune(seq_order)
    debug_infos = []

    # Convert tree to mafft format
    convert_tree_to_mafft(specietree, seq_order, mtoutput)
    mtoutput.close()

    # execute mafft
    is_already_aligned = is_aligned(args.seq, "fasta")
    if(enable_mafft and not (SKIPMAFFT and is_already_aligned)):
        execute_mafft(enable_mafft + " --treein %s %s > %s" %
                      (mtoutput.file, args.seq, MAFFT_OUTPUT))

    # check if is already aligned
    elif is_already_aligned:
        MAFFT_OUTPUT = args.seq

    # Reload mafft alignment and filter alignment (remove gap positions and
    # positions not conserved)
    alignment = AlignIO.read(MAFFT_OUTPUT, 'fasta', alphabet=alpha)
    debug_infos.append("Initial alignment length : %d"%alignment.get_alignment_length())

    if(args.excludegap):
        alignment = clean_alignment(
            alignment, threshold=(abs(args.excludegap) <= 1 or 0.01) * abs(args.excludegap))
        AlignIO.write(
            alignment, open(MAFFT_OUTPUT + "_ungapped", 'w'), 'fasta')

        debug_infos.append("Alignment length after removing gaps : %d"%alignment.get_alignment_length())


    # Get the substitution matrice at each node
    sptree = PhyloTree(specietree.write(), alignment=alignment.format(
        "fasta"), alg_format="fasta")

    # Compute expected frequency from the entire alignment and update at each
    # node
    armseq = alignment

    # We could use another alignment to compute those frequency if it's set
    if(ARM_SEQUENCE):
        armseq = AlignIO.read(ARM_SEQUENCE, 'fasta', alphabet=alpha)

    summary_info = AlignInfo.SummaryInfo(armseq)
    acc_rep_mat = SubsMat.SeqMat(summary_info.replacement_dictionary())
    # format of the acc_rep_mat:
    # {('A','C'): 10, ('C','H'): 12, ...}
    obs_freq_mat = SubsMat._build_obs_freq_mat(acc_rep_mat)
    # format of the obs_freq_matrix
    # same as the acc_rep_mat, but, we have frequency instead of count
    exp_freq_table = SubsMat._exp_freq_table_from_obs_freq(obs_freq_mat)
    # Not a dict of tuples anymore and it's obtained from the obs_freq_mat
    # {'A': 0.23, ...}
    # it's just the sum of replacement frequency for each aa.
    # if the two aa are differents, add half of the frequency
    # else add the frequency value
    # this unfortunaly assume that A-->C and C-->A have the same probability

    if (not SKIPSUBMATRIX):
        for node in sptree.traverse():
            if not node.is_leaf():
                spec_under_node = node.get_leaf_names()
                node_alignment = MultipleSeqAlignment(
                    [record[:] for record in alignment if record.id in spec_under_node], alphabet=alpha)

                if(REALIGN_AT_EACH_NODE and not node.is_root()):
                    node_alignment = realign(node_alignment, node, enable_mafft)

                node_align_info = AlignInfo.SummaryInfo(node_alignment)
                replace_info = node_align_info.replacement_dictionary()
                node_align_arm = SubsMat.SeqMat(replace_info)

                node_sub_matrix = SubsMat.make_log_odds_matrix(
                    node_align_arm, exp_freq_table=exp_freq_table)

                node.add_feature('submatrix', node_sub_matrix)

                if(args.verbose):
                    print "Sequences under node"
                    print spec_under_node
                    print "ARM matrix"
                    print node_align_arm
                    print "SUB matrix"
                    print node_sub_matrix

    
    # Filter using the ic content
    if IC_INFO_THRESHOLD:
        align_info = AlignInfo.SummaryInfo(alignment)
        ic_content = align_info.information_content(e_freq_table=(exp_freq_table if USE_EXPECTED_FREQ_FOR_IC else None))
        max_val = max(align_info.ic_vector.values())*IC_INFO_THRESHOLD
        ic_pos = (np.asarray(align_info.ic_vector.values())>max_val).nonzero()
        filtered_alignment = filter_align_position(alignment, ic_pos[0])
        AlignIO.write(filtered_alignment, open(MAFFT_OUTPUT + "_filtered_IC", 'w'), 'fasta')
        
        debug_infos.append("Alignment length after removing low IC columns : %d"%filtered_alignment.get_alignment_length())


    # Filter using the match percent per columns
    # Already enabled by default in the arguments list to filter 
    if(args.idfilter):
        filtered_alignment = filter_alignment(
            filtered_alignment, threshold=(abs(args.idfilter) <= 1 or 0.01) * abs(args.idfilter))
        AlignIO.write(
            filtered_alignment, open(MAFFT_OUTPUT + "_filtered", 'w'), 'fasta')

        debug_infos.append("Alignment length after removing columns less conserved than threshold (%f) : %d"%(args.idfilter, filtered_alignment.get_alignment_length()))
        
        # no keeping a variable for this
        filtered_alignment = filter_alignment(filtered_alignment, remove_identity=True, threshold=(
            abs(args.idfilter) <= 1 or 0.01) * abs(args.idfilter))
        AlignIO.write(filtered_alignment, open(
            MAFFT_OUTPUT + "_matchremoved", 'w'), 'fasta')

        debug_infos.append("Alignment length after removing columns less conserved than threshold (%f) and 100%% identical : %d"%(args.idfilter, filtered_alignment.get_alignment_length()))


    #################################################################################################

    # Compute sequence identity in global and filtered alignment
    seq_names = [rec.id for rec in alignment]
    number_seq = len(seq_names)
    sim_json = collections.defaultdict(list)
    aa_json = collections.defaultdict(list)
    expected_freq = collections.defaultdict(float)
    # alignment length
    af_length = len(filtered_alignment[0])
    ag_length = len(alignment[0])
    count_max = -np.inf
    count_min = np.inf
    matCalc = DistanceCalculator('identity')
    global_paired_distance = matCalc.get_distance(alignment)
    filtered_paired_distance = matCalc.get_distance(filtered_alignment)
    suspect_species = collections.defaultdict(Counter)

    if(EXCLUDE_AA):
        aa_letters = "".join([aa for aa in aa_letters if aa not in EXCLUDE_AA])

    genome_aa_freq = collections.defaultdict(dict)
    for i in xrange(number_seq):
        # Get aa count for each sequence
        count1 = Counter(alignment[i])
        count2 = Counter(filtered_alignment[i])
        for aa in aa_letters:
            expected_freq[aa_letters_1to3[aa]] = exp_freq_table[aa]
            global_val = count1[aa] / (ag_length * exp_freq_table[aa])
            filtered_val = count2[aa] / (af_length * exp_freq_table[aa])
            aa_json[aa_letters_1to3[aa]].append(
                {'global': global_val, 'filtered': filtered_val, "species": seq_names[i]})
            count_max = max(filtered_val, global_val, count_max)
            count_min = min(filtered_val, global_val, count_min)
            genome_aa_freq[seq_names[i]][aa_letters_1to3[aa]] = global_val

        for j in xrange(i + 1):

            sim_json[seq_names[i]].append({"global": (1-global_paired_distance[seq_names[i], seq_names[j]])
                                          , "filtered":(1-filtered_paired_distance[seq_names[i], seq_names[j]]), "species": seq_names[j]})
            if i != j:
                sim_json[seq_names[j]].append({"global": (1-global_paired_distance[seq_names[i], seq_names[j]])
                                          , "filtered": (1-filtered_paired_distance[seq_names[i], seq_names[j]]), "species": seq_names[i]})
                
            # do not add identity to itself twice
         
    sptree.write(outfile=TMP + "phylotree.nw")

    # for performance issues, better make another loop to
    # get the data to plot the conservation of aa in each column
    # of the alignment

    consensus = get_consensus(filtered_alignment, AA_MAJORITY_THRESH)
    
    # A little pretraitment to speed access to each record later
    global_consensus = get_consensus(alignment, AA_MAJORITY_THRESH)
    record2seq = {}
    for record in alignment:
        record2seq[record.id] = record

    if(args.debug):
        print "Filtered alignment consensus : ", consensus

    aa2alignment = {}
    aa2identy_dict = {}
    for aa in aa_letters:
        cons_array = get_aa_filtered_alignment(consensus, aa)
        #print aa, "\n", cons_array, 
        if(cons_array):
            aa_filtered_alignment = filter_align_position(filtered_alignment, cons_array)
            aa2alignment[aa_letters_1to3[aa]] = aa_filtered_alignment
            aa2identy_dict[aa] =  matCalc.get_distance(aa_filtered_alignment)


    def get_aa_maj(mapinput):
        aa_shift_json = collections.defaultdict(list)
        i, j = mapinput
        global global_paired_distance
        global aa_letters
        global seq_names
        global aa2identy_dict
        gpaired = 1-global_paired_distance[seq_names[i], seq_names[j]]

        for aa in aa_letters:
            if(aa in aa2identy_dict.keys()):
                fpaired = 1 - aa2identy_dict[aa][seq_names[i], seq_names[j]]
                aa_shift_json[aa_letters_1to3[aa]].append({'global':gpaired , 'filtered': fpaired, "species": "%s_%s" % (seq_names[i], seq_names[j])})
        return aa_shift_json


    p = Pool(PROCESS_ENABLED)
    result = p.map(get_aa_maj, [(max(ind), min(ind)) for ind in itertools.combinations(xrange(number_seq), r=2)])

    aa_shift_json = Counter()
    for r in result:
        aa_shift_json += Counter(r)

 
    # dumping in json to reload with the web interface using d3.js
    with open(TMP + "similarity.json", "w") as outfile1:
        json.dump(sim_json, outfile1, indent=4)
    with open(TMP + "aafrequency.json", "w") as outfile2:
        json.dump(
            {"AA": aa_json, "EXP": expected_freq, "MAX": count_max, "MIN" : count_min}, outfile2, indent=4)
    with open(TMP + "aause.json", "w") as outfile3:
        json.dump(aa_shift_json, outfile3, indent=4)


    # transposing counter to dict and finding potential species 
    # that had a codon reassignment
    most_common = collections.defaultdict(list)
    for key, value in aa_shift_json.iteritems():
        for v in value:
            specs = v['species'].split('_')
            if(v['global']> v['filtered']):
                suspect_species[key][specs[0]] = (suspect_species[key][specs[0]] or 0) + 1
                suspect_species[key][specs[1]] = (suspect_species[key][specs[1]] or 0) + 1

        common_list =  suspect_species[key].most_common()
        i = 0
        while i < len(common_list) and common_list[i][1] > COUNTER_THRESHOLD*len(seq_names):
            global_freq_use = genome_aa_freq[common_list[i][0]][key]
            # we are going to use a radius around our expected values
            #if(abs(1-global_freq_use) <= 1):
            most_common[key].append(common_list[i]+(global_freq_use,))
            i += 1

    # Let's say that at this step we have the most suspected species for each aa.
    # for each aa let's find the targeted aa
    # FIX ME: this part of the code is now almost useless

    aa2aa_rea = collections.defaultdict(dict)
    for key, values in most_common.iteritems():
        susspeclist = [val[0] for val in values]
        aa_alignment = aa2alignment[key]

        for s in aa_alignment:
            if s.id in susspeclist:
                suspected_aa = []
                for i in range(len(s)):
                    if(s[i]!='-' and aa_letters_1to3[s[i]] != key):
                        suspected_aa.append(aa_letters_1to3[s[i]])
                        try :
                            aa2aa_rea[key][aa_letters_1to3[s[i]]].add(s.id)
                        except KeyError:
                            aa2aa_rea[key] = collections.defaultdict(set)
                            aa2aa_rea[key][aa_letters_1to3[s[i]]].add(s.id)

                pos = susspeclist.index(s.id)
                most_common[key][pos] += (suspected_aa,)


    fitch_tree = []
    counts = []
    for key1, dict2 in aa2aa_rea.iteritems():
        
        key1_alignment = aa2alignment[key1]
        for key2, val in dict2.iteritems():
            t = sptree.copy("newick")
            n = NaiveFitch(t, val, key2, key1)
            slist = n.get_species_list(LIMIT_TO_SUSPECTED_SPECIES)
            sslist = n.get_species_list()
            for s in slist:
                rec = record2seq[s]
                leaf = (n.tree&s)
                ori_count = 0
                try :
                                 
                    ori_count = len([y for x in key1_alignment for y in x if x.id==s and y!='-' and aa_letters_1to3[y]==key2])
                except Exception:
                    #wtver happen, do nothing
                    pass

                leaf.add_features(count=0)
                leaf.add_features(filter_count=ori_count)
                for position in range(len(rec)):
                    if global_consensus[position] == aa_letters_3to1[key1] \
                        and rec[position] == aa_letters_3to1[key2]:
                        leaf.count+=1
                
                if s in sslist:
                    leaf.add_features(lost=False)
                    if(leaf.count < COUNT_THRE):
                        leaf.lost = True

                counts.append((s, str(leaf.count), str(leaf.filter_count)))
            if(args.debug):
                print "\n\n", key1, " to ", key2 , ": "
                print "(Suspected species, global_count, filtered_count) : \n", "\n".join("\t".join(c) for c in counts)

            if(n.is_valid()):
                n.render_tree(suffix=args.sfx)
                fitch_tree.append(n)

    if(args.debug):
        for line in debug_infos:
            print line
        print "After validating the ancestral state and checking in the global alignment, %d cases were found interesting"%len(fitch_tree)

    # Use a dayhoff matrix to determine common substitution and filter
    # and filter result based on that
    # with that, we can remove false positive