#! /usr/bin/env python
from __future__ import division

import argparse
import collections
import itertools
import json
import numpy as np
import os
import random
import subprocess

from multiprocessing import Pool

from collections import Counter

from Bio import AlignIO
from Bio import SeqIO
from Bio import SubsMat
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_nucleotide
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.codonalign.codonseq import _get_codon_list
from TreeLib import TreeClass
from ete2 import PhyloTree
from settings import *
from utils import *
from shutil import copyfile


__author__ = "Emmanuel Noutahi"
__version__ = "0.1"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"


if not os.path.exists(TMP):
    os.makedirs(TMP)

MAFFT_AUTO_COMMAND = [
    'linsi', 'ginsi', 'einsi', 'fftnsi', 'fftns', 'nwnsi', 'nwns']

MAFFT_DETAILLED_COMMAND = [
    'auto', 'retree', 'maxiterate', 'nofft', 'memsave', 'parttree', 'leavegappyregion']
# This is the multiple alignment used to estimate the accepted replacement
# matrix


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
    record_dict = {}
    for fasta in fasta_sequences:
        seq_order.append(fasta.description)
        record_dict[fasta.description] = fasta

    # load tree
    specietree = TreeClass(args.tree)
    leave_names = set(specietree.get_leaf_name())
    
    # Use a portion of the sequence for test if asked by the user
    if args.resample:
        seq_order = random.sample(seq_order, args.resample)
    
    # check duplicated sequence name in the alignment file
    original_len = len(seq_order)

    # get dna to perform codon alignment
    dnaseq = SeqIO.parse(DNA_SEQ_PATH, 'fasta',generic_nucleotide)

    if len(set(seq_order)) != original_len or (set(seq_order) - leave_names):
        Output.error("Sequence not matching found, attempt to correct... ", "Warning")
        seq_order[:] = list(set(seq_order) & leave_names)

        args.seq = args.seq + "_"
        with open(args.seq, 'w') as FASTAOUT:
            for seq in seq_order:
                FASTAOUT.write(record_dict[seq].format("fasta"))


    # prune tree to sequence list
    specietree.prune(seq_order)
    debug_infos = []
    dnaseq = [dna for dna in dnaseq if dna.id in seq_order]

    # Convert tree to mafft format
    convert_tree_to_mafft(specietree, seq_order, mtoutput, args.scale)
    mtoutput.close()

    # execute mafft
    is_already_aligned = is_aligned(args.seq, "fasta")
    if(enable_mafft and not (SKIPMAFFT and is_already_aligned)):
        execute_mafft(enable_mafft + " --treein %s %s > %s" %
                      (mtoutput.file, args.seq, MAFFT_OUTPUT))

    # check if is already aligned
    elif is_already_aligned:
        copyfile(args.seq, MAFFT_OUTPUT)

    # Reload mafft alignment and filter alignment (remove gap positions and
    # positions not conserved)
    alignment = AlignIO.read(MAFFT_OUTPUT, 'fasta', alphabet=alpha)
    record2seq = SeqIO.to_dict(alignment)
    debug_infos.append("Initial alignment length : %d"%alignment.get_alignment_length())
    
    # add this to keep trace of the gap filtered position 
    gap_filtered_position  = []
    tt_filter_position =  np.asarray(xrange(alignment.get_alignment_length()))

    if(args.excludegap): 
        alignment, gap_filtered_position = clean_alignment(
            alignment, threshold=(abs(args.excludegap) <= 1 or 0.01) * abs(args.excludegap))
        AlignIO.write(
            alignment, open(MAFFT_OUTPUT + "_ungapped", 'w'), 'fasta')

        debug_infos.append("Alignment length after removing gaps : %d"%alignment.get_alignment_length())

    # update list of position from the original alignment
    tt_filter_position = tt_filter_position[gap_filtered_position]

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
        # update list of position based on change in ic_position
        tt_filter_position = tt_filter_position[ic_pos]
        AlignIO.write(filtered_alignment, open(MAFFT_OUTPUT + "_filtered_IC", 'w'), 'fasta')
        
        debug_infos.append("Alignment length after removing low IC columns : %d"%filtered_alignment.get_alignment_length())


    # Filter using the match percent per columns
    # Already enabled by default in the arguments list to filter 
    if(args.idfilter):
        filtered_alignment, position = filter_alignment(
            filtered_alignment, threshold=(abs(args.idfilter) <= 1 or 0.01) * abs(args.idfilter))
        AlignIO.write(
            filtered_alignment, open(MAFFT_OUTPUT + "_filtered", 'w'), 'fasta')

        # update : remove position not conserved
        tt_filter_position = tt_filter_position[position]
        debug_infos.append("Alignment length after removing columns less conserved than threshold (%f) : %d"%(args.idfilter, filtered_alignment.get_alignment_length()))
        
        # no keeping a variable for this
        filtered_alignment, position = filter_alignment(filtered_alignment, remove_identity=True, threshold=(
            abs(args.idfilter) <= 1 or 0.01) * abs(args.idfilter))
        AlignIO.write(filtered_alignment, open(
            MAFFT_OUTPUT + "_matchremoved", 'w'), 'fasta')

        # update remove of identical column
        tt_filter_position = tt_filter_position[position]
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
    record2ungapedseq = {}
    for record in alignment:
        record2ungapedseq[record.id] = record

    debug_infos.append("Filtered alignment consensus : \n%s\n"%consensus)

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
    if(JSON_DUMP):
        with open(TMP + "similarity.json", "w") as outfile1:
            json.dump(sim_json, outfile1, indent=4)
        with open(TMP + "aafrequency.json", "w") as outfile2:
            json.dump({"AA": aa_json, "EXP": expected_freq, "MAX": count_max, "MIN" : count_min}, outfile2, indent=4)
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
        while i < len(common_list) and common_list[i][1] > FREQUENCY_THRESHOLD*len(seq_names):
            global_freq_use = genome_aa_freq[common_list[i][0]][key]
            most_common[key].append(common_list[i]+(global_freq_use,))
            i += 1

    # At this step we have the most suspected species for each aa.
    # for each aa let's find the targeted aa

    aa2aa_rea = collections.defaultdict(dict)
    for key, values in most_common.iteritems():
        susspeclist = [val[0] for val in values]
        aa_alignment = aa2alignment[key]

        for s in aa_alignment:
            if s.id in susspeclist:
                suspected_aa = []
                for i in range(len(s)):
                    one_letter_key = aa_letters_3to1[key]
                    if(s[i]!='-' and s[i] != one_letter_key):
                        suspected_aa.append(s[i])
                        try :
                            aa2aa_rea[one_letter_key][s[i]].add(s.id)
                        except KeyError:
                            aa2aa_rea[one_letter_key] = collections.defaultdict(set)
                            aa2aa_rea[one_letter_key][s[i]].add(s.id)

                #pos = susspeclist.index(s.id)
                #most_common[key][pos] += (suspected_aa,)

    # Parsimony fitch tree list
    fitch_tree = []

    codon_alignment, fcodon_alignment = codon_align(dnaseq, record2seq, gap_filtered_position, tt_filter_position)
    #codon_lst = []
    #for codon_aln in codon_alignment.values():
        #codon_lst.append(_get_codon_list(codon_aln.seq))

    for key1, dict2 in aa2aa_rea.iteritems():
    
        key1_alignment = aa2alignment[aa_letters_1to3[key1]]
        for key2, val in dict2.iteritems():
            gcodon_rea = CodonReaData((key1, key2), global_consensus, codon_alignment)
            fcodon_rea = CodonReaData((key1, key2), consensus, fcodon_alignment)
            counts = []
            t = sptree.copy("newick")
            n = NaiveFitch(t, val, aa_letters_1to3[key2], aa_letters_1to3[key1], (gcodon_rea, fcodon_rea))
            slist = n.get_species_list(LIMIT_TO_SUSPECTED_SPECIES)
            
            for s in slist:
                rec = record2ungapedseq[s]
                leaf = (n.tree&s)
                ori_count = 0
                try :
                    ori_count = len([y for x in key1_alignment for y in x if x.id==s and y!='-' and y==key2])
                except Exception:
                    #wtver happen, do nothing
                    pass

                leaf.add_features(count=0)
                leaf.add_features(filter_count=ori_count)
                leaf.add_features(lost=False)
                for position in range(len(rec)):
                    if global_consensus[position] == key1 \
                        and rec[position] == key2:
                        leaf.count+=1
                
                if('lost' in leaf.features and leaf.count < COUNT_THRESHOLD):
                    leaf.lost = True

                counts.append((s, str(leaf.count), str(leaf.filter_count), gcodon_rea.get_string(s, SHOW_MIXTE_CODONS), fcodon_rea.get_string(s, SHOW_MIXTE_CODONS)))
            debug_infos.append("\n\n Substitutions : " + key2 + " to "+ key1 + ": ")
            if(SHOW_MIXTE_CODONS):
                debug_infos.append("species\tglob_AA_count\tfilt_AA_count\tglob_reacodon_count\tglob_usedcodon_count\tglob_mixtecodon_count\tfilt_reacodon_count\tfilt_usedcodon_count\tfilt_mixtecodon_count\n" + "\n".join("\t".join(c) for c in counts))
            else:
                debug_infos.append("species\tglob_AA_count\tfilt_AA_count\tglob_reacodon_count\tglob_othercodon_count\tfilt_reacodon_count\tfilt_othercodon_count\n" + "\n".join("\t".join(c) for c in counts))

            if(n.is_valid()):
                n.render_tree(suffix=args.sfx)
                fitch_tree.append(n)

    if(args.debug):
        for line in debug_infos:
            print line
        print "After validating the ancestral state and checking in the global alignment, %d cases were found interesting"%len(fitch_tree)