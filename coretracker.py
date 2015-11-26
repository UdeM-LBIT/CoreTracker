#!/usr/bin/env python

import argparse
import os, sys
import traceback

import warnings
warnings.filterwarnings("ignore") 
sys.path.append(os.path.dirname(__file__))
from Bio import AlignIO
from Bio import SeqIO

from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import Seq
from ete3 import Tree
from utils import *
from shutil import copyfile
import random

__author__ = "Emmanuel Noutahi"
__version__ = "0.1"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"


MAFFT_AUTO_COMMAND = ['linsi', 'ginsi', 'einsi', 'fftnsi', 'fftns', 'nwnsi', 'nwns', 'auto']

def convert_tree_to_mafft(tree, seq_order, output, scale):
    """Convert a tree in a newick format to the mafft matrix format"""
    seqnames = [-1, -1]
    branchlens = [-1, -1]
    for node in tree.traverse("postorder"):
        if node.dist <= DIST_THRES:
            node.dist = 0.0

        if not node.is_leaf():
            left_branch = node.children[0]
            right_branch = node.children[0]
            left_branch_leaf = left_branch.get_leaf_names()
            right_branch_leaf = right_branch.get_leaf_names()

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


def get_argsname(command_dict, command_list, prefix=""):
    """ Get the mafft command from the argparse dict : command_dict
    and a command_list
    """
    return ["%s%s" % (prefix, command) for command in command_list if command_dict[command]]


def testing_a_lot(args, settings):
    t = random.randint(30, 90)
    time.sleep(t)
    return [settings.EXCLUDE_AA, settings.AA_MAJORITY_THRESH, settings.FREQUENCY_THRESHOLD, 
            settings.GENETIC_CODE, settings.COUNT_THRESHOLD, settings.LIMIT_TO_SUSPECTED_SPECIES]


def is_aligned(seqfile, format):
    try:
        a = AlignIO.read(open(seqfile), format)
    except ValueError:
        return False
    return True


def coretracker(args, settings):
    """Run coretracker on the argument list"""
      # Check mafft command input
    run_alignment = False
    mafft_detail = get_argsname(args.__dict__, MAFFT_AUTO_COMMAND, prefix="--")
    mafft_cmd = "mafft "
    if(mafft_detail):
        run_alignment = True
        mafft_cmd += mafft_detail[0]
    
    if args.outdir:
        settings.OUTDIR = args.outdir

    input_alignment = args.seq


    # Manage fasta sequence input
    protseq = SeqIO.parse(input_alignment, "fasta")
    seq_names = set([seq.id for seq in protseq])
    # load tree
    specietree = Tree(args.tree)
    leaf_names = specietree.get_leaf_names()
    c_genome = seq_names.intersection(leaf_names)
    
    # get dna to perform codon alignment
    dnaseq = SeqIO.parse(args.dnaseq, 'fasta', generic_nucleotide)
    dnadict = dict((seq.id, seq) for seq in dnaseq)
    c_genome.intersection_update(dnadict.keys())

    if len(c_genome) <  0.5 * (len(leaf_names) + len(dnadict.keys()) + len(seq_names))/3.0 :
        raise ValueError("Too many ID do not match in protein sequence, nucleic sequence and species tree.") 

    # execute mafft    
    is_already_aligned = is_aligned(input_alignment, 'fasta')

    if not (settings.SKIP_ALIGNMENT and is_already_aligned):
        if not run_alignment:
            raise ValueError("Could not align sequences\nAlignment parameters not provided")
        else :
            # this is the case where we need to convert the tree to mafft tree before
            # Output stream setting
            mtoutput = Output(os.path.join(settings.OUTDIR,"tree_mafft.nw"))
            # Convert tree to mafft format
            convert_tree_to_mafft(specietree, seq_names, mtoutput, args.scale)
            mtoutput.close()
            unaligned_seqs = input_alignment+"_unaligned"
            copyfile(input_alignment, unaligned_seqs)
            execute_mafft(mafft_cmd + " --treein %s %s > %s" %
                          (mtoutput.file, unaligned_seqs, input_alignment))



    # create sequence set
    alignment = AlignIO.read(input_alignment, 'fasta', alphabet=alpha)
    setseq = SequenceSet(dnadict, alignment, specietree, settings.GENETIC_CODE, args.stopcodon)
    setseq.prot_filtering(args.idfilter, args.gapfilter, args.iccontent)
    reafinder =  ReaGenomeFinder(setseq, settings)
    reafinder.get_genomes()
    reafinder.save_all()


if __name__ == '__main__':

    # argument parser
    parser = argparse.ArgumentParser(
        description='CoreTracker, A codon reassignment tracker newick tree format to mafft format')

    parser.add_argument(
        '--wdir', '--outdir', dest="outdir", help="Working directory")

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('--gapfilter', '--gap', type=float,  default=0.6, dest='gapfilter',
                        help="Remove position with gap from the alignment, using gapfilter as threshold. The absolute values are taken")
    
    parser.add_argument('--idfilter', '--id', type=float, default=0.8, dest='idfilter',
                        help="Conserve only position with at least idfilter residue identity")
    
    parser.add_argument('--iccontent', '--ic', type=float, default=0.5, dest='iccontent',
                        help="Shannon entropy threshold (default : 0.5 ). This will be used to discard column where IC < max(IC_INFO_VECTOR)*(IC_INFO_THRESHOLD/ 100.0)))")
    

    parser.add_argument(
        '--verbose', '-v', choices=[0,1,2], type=int, default=0, dest="verbose", help="Verbosity level")


    parser.add_argument(
        '--sfx', dest="sfx", default="", help="PDF rendering suffix to differentiate runs.")

    parser.add_argument('-t', '--intree', dest="tree",
                        help='Input specietree in newick format', required=True)
    
    parser.add_argument('-s', '--scale', type=float, default=1.0,
                        dest='scale', help="Scale to compute the branch format")

    parser.add_argument('--protseq', '--prot', '-p', dest='seq',
                        help="Protein sequence input in fasta format", required=True)

    parser.add_argument('--dnaseq', '--dna', '-n', dest='dnaseq',
                        help="Nucleotides sequences input in fasta format", required=True)
   
    parser.add_argument('--stopcodon', dest='stopcodon', action='store_true',
                        help="Whether or not stop present in protein alignment and dna sequences.")

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
   
    mafft_group.add_argument(
        '--auto', dest='auto', action='store_true', help="If unsure which option to use, try this option")


    args = parser.parse_args()
    setting =  Settings()
    setting.set()
    coretracker(args, setting)