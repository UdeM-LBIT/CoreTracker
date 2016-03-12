#!/usr/bin/env python

import argparse
import glob
import logging
import os
import random
import sys
import traceback
import warnings
from collections import defaultdict as ddict

from Bio import SeqIO
from ete3 import Tree

from utils.classifier import Classifier, getDataFromFeatures, read_from_json
from settings import Settings, parameters
from coreutils import *
import coreutils.utils as utils

warnings.filterwarnings("ignore")
sys.path.append(os.path.dirname(__file__))


__author__ = "Emmanuel Noutahi"
__version__ = "1.2"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"


def testing_a_lot(args, settings):
    t = random.randint(30, 90)
    time.sleep(t)
    return [settings.EXCLUDE_AA, settings.AA_MAJORITY_THRESH, settings.FREQUENCY_THRESHOLD,
            settings.GENETIC_CODE, settings.COUNT_THRESHOLD, settings.LIMIT_TO_SUSPECTED_SPECIES]

def coretracker(args, settings):
    """Run coretracker on the argument list"""
    # Check mafft command input
    progcmd = lambda x: x + ' --auto' if x == 'mafft' else x

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    if args.outdir:
        settings.OUTDIR = args.outdir

    msaprg = progcmd(args.align)

    input_alignment = args.seq
    use_tree = None
    specietree = Tree(args.tree)

    if args.usetree:
        use_tree = args.tree

    settings.scale = args.scale
    settings.showall = args.showall
    hmmfiles = {}
    if args.hmmdir:
        try:
            hfiles = glob.glob(os.path.join(args.hmmdir, '*'))
            for f in hfiles:
                genename = os.path.basename(f).split('.hmm')[0]
                hmmfiles[genename] = f
        except:
            pass

    seqloader = SequenceLoader(input_alignment, args.dnaseq, settings, args.gapfilter, has_stop=args.hasstop,
                    use_tree=args.usetree, refine_alignment=args.refine, msaprog=msaprg,  hmmdict=hmmfiles)

    # create sequence set
    setseq = SequenceSet(seqloader, specietree, settings.GENETIC_CODE)
    setseq.prot_filtering(args.idfilter, args.gapfilter,
                          args.iccontent, args.rmconst)

    reafinder = ReaGenomeFinder(setseq, settings)
    reafinder.get_genomes()
    reafinder.possible_aa_reassignation()
    codon_align, fcodon_align = reafinder.seqset.get_codon_alignment()
    cod_align = SeqIO.to_dict(fcodon_align)
    reafinder.set_rea_mapper()
    clf = Classifier.load_from_file(parameters.MODELPATH)
    etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac",
                    "N. rea", "N. used", "Cod. count", "Sub. count",
                    "G. len", "codon_lik", "N. mixte" ,"id"] #, 'total_aa']
    selected_feats = [2,3,4,5,6,7,8,9,11]

    if clf is None or not clf.trained:
        logging.error("Classifier not found or not trained!")

    for (fitch, data) in reafinder.run_analysis(codon_align, fcodon_align):
        s_complete_data = makehash()
        s_complete_data['aa'][fitch.ori_aa1][fitch.dest_aa1] = data
        s_complete_data['genome'] = reafinder.reassignment_mapper['genome']
        X_data, X_labels, _ = classifier.read_from_json(s_complete_data, None, use_global=False)
        # extract usefull features
        X_data, _ = classifier.getDataFromFeatures(X_data, etiquette, feats=selected_feats)
        pred_prob = clf.predict_proba(X_data)
        pred =  clf.predict(X_data)
        utils.get_report(fitch, data, reafinder, cod_align, (X_data, X_labels, pred_prob, pred))

    # Print list of interesting cases
    logging.debug("After validation, %d cases were found interesting" % len(reafinder.interesting_case))
    for case in reafinder.interesting_case:
        logging.debug(case)

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
        '--verbose', '-v', choices=[0, 1, 2], type=int, default=0, dest="verbose", help="Verbosity level")

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

    parser.add_argument('--stopcodon', dest='hasstop', action='store_true',
                        help="Whether or not stop are present in protein alignment and dna sequences.")

    parser.add_argument('--debug', dest='debug', action='store_true',
                        help="Enable debug printing")

    parser.add_argument('--rmconst', dest='rmconst', action='store_true',
                        help="Remove constant site from filtered alignment. ")

    parser.add_argument('--showall', dest='showall', action='store_true',
                        help="Save all result. Do not attempt to filter based on significance")

    parser.add_argument('--norefine', dest='refine', action='store_false',
                        help="Do not refine the alignment. By default the alignment will be refined. This option should never be used if you have made your own multiple alignment and concatenate it. Else you will have absurd alignment (TO FIX)")

    parser.add_argument('--align', dest='align', choices=[
                        'muscle', 'mafft'], default="muscle", help="Choose a program to align your sequences")

    parser.add_argument('--use_tree', dest='usetree', action="store_true",
                        help="This is helpfull only if the mafft alignment is selected. Perform multiple alignment, using species tree as guide tree.")
    parser.add_argument('--submat', dest='submat', choices=AVAILABLE_MAT, default="blosum62",
                        help="Choose a substitution matrix to compute codon alignment to amino acid likelihood, Default value is blosum62")

    parser.add_argument('--hmmdir', dest='hmmdir',
                        help="Link a directory with hmm files for alignment. Each hmmfile should be named in the following format : genename.hmm")

    args = parser.parse_args()
    setting = Settings()
    setting.set()
    setting.update_params(SUBMAT=args.submat)
    coretracker(args, setting)
