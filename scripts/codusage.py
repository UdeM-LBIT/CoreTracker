#!/usr/bin/env python

# CoreTracker Copyright (C) 2016  Emmanuel Noutahi
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import glob
import logging
import os
import random
import sys
from yaml import load
import traceback
from collections import defaultdict as ddict
from collections import Counter

from math import log10
from Bio import SeqIO
from ete3 import Tree

from coretracker.coreutils import *
from coretracker.classifier.models import ModelType
from coretracker.classifier import *
from coretracker.classifier import MODELPATH
from coretracker.settings import *
from coretracker import  __version__, __author__, date

ENABLE_PAR = True
CPU_COUNT = 0
try:
    from multiprocessing import cpu_count
    CPU_COUNT = cpu_count()
    from joblib import Parallel, delayed
except ImportError:
    try:
        from sklearn.externals.joblib import Parallel, delayed
    except:
        ENABLE_PAR = False

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac",
                    "N. rea", "N. used", "Cod. count", "Sub. count",
                    "G. len", "codon_lik", "N. mixte" ,"id"] #, 'total_aa']

def testing_a_lot(args, settings):
    t = random.randint(30, 90)
    time.sleep(t)
    return [settings.EXCLUDE_AA, settings.AA_MAJORITY_THRESH, settings.FREQUENCY_THRESHOLD,
            settings.GENETIC_CODE, settings.COUNT_THRESHOLD, settings.LIMIT_TO_SUSPECTED_SPECIES]

def set_coretracker(args, settings):
    """Set all data for coretracker from the argument list"""
    # Check mafft command input
    progcmd = lambda x: x + ' --auto' if x == 'mafft' else x

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
            # let original error handling
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

    clf = Classifier.load_from_file(MODELPATH%settings.MODEL_TYPE)
    model = ModelType(settings.MODEL_TYPE, etiquette)


    if clf is None or not clf.trained:
        raise ValueError("Classifier not found or not trained!")

    seqloader = SequenceLoader(input_alignment, args.dnaseq, settings, args.gapfilter, has_stop=args.hasstop,
                    use_tree=args.usetree, refine_alignment=args.refine, msaprog=msaprg,  hmmdict=hmmfiles)

    # create sequence set
    setseq = SequenceSet(seqloader, specietree, settings.GENETIC_CODE)
    setseq.prot_filtering(args.idfilter, args.gapfilter,
                          args.iccontent, args.rmconst)

    reafinder = ReaGenomeFinder(setseq, settings)
    reafinder.get_genomes()
    reafinder.possible_aa_reassignation()
    return reafinder, clf, model

def compile_result(x, clf, cod_align, model):
    """compile result from analysis"""
    reafinder, fitch, data = x
    s_complete_data = utils.makehash()
    s_complete_data['aa'][fitch.ori_aa1][fitch.dest_aa1] = data
    s_complete_data['genome'] = reafinder.reassignment_mapper['genome']
    X_data, X_labels, _ = read_from_json(s_complete_data, None, use_global=False)
    # extract usefull features
    X_data, X_dataprint, selected_et =  model.format_data(X_data)
    pred_prob = clf.predict_proba(X_data)
    pred =  clf.predict(X_data)
    sppval, outdir, rkp = utils.get_report(fitch, data, reafinder, cod_align, (X_data, X_labels, pred_prob, pred))
    utils.print_data_to_txt(os.path.join(outdir, fitch.ori_aa+"_to_"+fitch.dest_aa+"_data.txt"), selected_et, X_dataprint, X_labels, pred, pred_prob, sppval, fitch.dest_aa)
    return rkp

def codon_counter(codonseq, l):
    codlist = [str(codonseq[i*3:(i+1)*3].seq) for i in range(l)]
    return Counter(codlist)


def filter_to_pos_with_aa(alignment, codonalign, X, codtable):
    # identify position 
    cons = SequenceSet.get_consensus(alignment, 0.5, ambiguous='X')
    ind_array = [i for (i, aa) in enumerate(cons) if aa==X]
    codon_al = SequenceSet.filter_align_position(codonalign, ind_array, codontable=codtable)
    return SeqIO.to_dict(codon_al), codon_al.get_aln_length()


if __name__ == '__main__':

    # argument parser
    parser = argparse.ArgumentParser(
        description='CoreTracker, A codon reassignment tracker')

    parser.add_argument(
        '--wdir', '--outdir', dest="outdir", default="output", help="Working directory")

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('--gapfilter', '--gap', type=float,  default=0.6, dest='gapfilter',
                        help="Remove position with gap from the alignment, using gapfilter as threshold. The absolute values are taken")

    parser.add_argument('--idfilter', '--id', type=float, default=0.8, dest='idfilter',
                        help="Conserve only position with at least idfilter residue identity")

    parser.add_argument('--iccontent', '--ic', type=float, default=0.5, dest='iccontent',
                        help="Shannon entropy threshold (default : 0.5 ). This will be used to discard column where IC < max(IC_INFO_VECTOR)*(IC_INFO_THRESHOLD/ 100.0)))")

    parser.add_argument(
        '--verbose', '-v', choices=[0, 1, 2], type=int, default=0, dest="verbose", help="Verbosity level")

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

    parser.add_argument('--expos','--export_position', dest='expos', action="store_true",
                    help="Export a json file with the position of each reassignment in the corresponding genome.")

    parser.add_argument('--aalist', dest='aalist',
                    help="List of aa to consider in order to find their codon usage in both conserved and not conserved position")


    parser.add_argument('--submat', dest='submat', choices=utils.AVAILABLE_MAT, default="blosum62",
                        help="Choose a substitution matrix to compute codon alignment to amino acid likelihood, Default value is blosum62")

    parser.add_argument('--hmmdir', dest='hmmdir',
                        help="Link a directory with hmm files for alignment. Each hmmfile should be named in the following format : genename.hmm")

    parser.add_argument('--params', dest='params',
                        help="Use A parameter file to load parameters. If a parameter is not set, the default will be used")

    parser.add_argument('--parallel', dest='parallel', nargs='?', const=CPU_COUNT, type=int, default=0,
                    help="Use Parallelization during execution for each reassignment. This does not guarantee an increase in speed. CPU count will be used if no argument is provided")

    print("CoreTracker v:%s Copyright (C) %s %s"%(__version__, date, __author__))

    args = parser.parse_args()
    setting = Settings()

    paramfile = args.params
    ad_params = {}
    if paramfile :
        with open(paramfile) as f:
            ad_params = load(f, Loader=Loader)

    # setting.set(OUTDIR=args.outdir)  # this does nothing
    # added for position computation
    setting.COMPUTE_POS = args.expos
    setting.fill(ad_params)
    setting.update_params(SUBMAT=args.submat)
    parallel = args.parallel
    reafinder, clf, model = set_coretracker(args, setting)
    codon_align, fcodon_align = reafinder.seqset.get_codon_alignment()
    fcod_align = SeqIO.to_dict(fcodon_align)
    cod_align = SeqIO.to_dict(codon_align)
    cod_aa, codaalen = filter_to_pos_with_aa(reafinder.seqset.prot_align, codon_align, args.aalist, reafinder.seqset.codontable)
    #reafinder.set_rea_mapper()

    if args.expos:
        rea_pos_keeper = ddict(dict)
        allen = fcodon_align.get_aln_length()
        for spec, seq in fcod_align.items():
            rea_pos_keeper[spec] = codon_counter(seq, allen)
        exp_outfile = os.path.join(reafinder.settings.OUTDIR, "fpositions.json" )
        reafinder.export_position(rea_pos_keeper, exp_outfile)

        rea_pos_keeper = ddict(dict)
        allen = codon_align.get_aln_length()
        for spec, seq in cod_align.items():
            rea_pos_keeper[spec] = codon_counter(seq, allen)
        exp_outfile = os.path.join(reafinder.settings.OUTDIR, "gpositions.json" )
        reafinder.export_position(rea_pos_keeper, exp_outfile)

        rea_pos_keeper = ddict(dict)
        for spec, seq in cod_aa.items():
            rea_pos_keeper[spec] = codon_counter(seq, codaalen)
        exp_outfile = os.path.join(reafinder.settings.OUTDIR, "aapositions.json" )
        reafinder.export_position(rea_pos_keeper, exp_outfile)

