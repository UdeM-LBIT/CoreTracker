from __future__ import division
import argparse
import glob
import itertools
import json
import numpy as np
import os
import random
import subprocess
import sys
import shutil
import time

from collections import Counter, defaultdict

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Alphabet
from distutils import spawn
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio import SeqIO
from Bio import codonalign
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.codonalign.codonalphabet import default_codon_alphabet, get_codon_alphabet
from Bio.Data import CodonTable
from Bio.codonalign.codonseq import _get_codon_list, CodonSeq
from Bio import SubsMat
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet import IUPAC, generic_protein, generic_nucleotide
from cStringIO import StringIO
from ete3 import Tree
from PPieChartFace import PPieChartFace, LineFace
import Bio.SubsMat.MatrixInfo as MatrixInfo

import settings as default_settings
from functools import partial
import scipy.stats as ss
from FisherExact import fisher_exact
from scipy.cluster.vq import kmeans2
import logging
__author__ = "Emmanuel Noutahi"
__version__ = "0.1"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"


aa_letters_1to3 = {

    'A': 'Ala', 'C': 'Cys', 'D': 'Asp',
    'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met',
    'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp',
    'Y': 'Tyr',
}

aa_letters_3to1 = dict((x[1], x[0]) for x in aa_letters_1to3.items())

AVAILABLE_MAT = MatrixInfo.available_matrices

# hardcorded settings variables
alpha = Alphabet.Gapped(IUPAC.protein)
DIST_THRES = 1e-10

# if settings.GENETIC_CODE < 0:
# for codon in ['CTG', 'CTA', 'CTC', 'CTT']:
#dct.forward_table[codon] = 'L'
#dct.back_table = CodonTable.make_back_table(dct.forward_table, dct.stop_codons[0])

nuc_letters = "ACTG"


class CodonSeq(CodonSeq):

    def __init__(self, data='', alphabet=default_codon_alphabet,
                 gap_char="-", rf_table=None, enable_undef=False):

        Seq.__init__(self, data.upper(), alphabet=alphabet)
        self.gap_char = gap_char
        self.enable_undef =  enable_undef

        # check the length of the alignment to be a triple
        if rf_table is None:
            seq_ungapped = self._data.replace(gap_char, "")
            assert len(self) % 3 == 0, "Sequence length is not a triple number"
            self.rf_table = list(filter(lambda x: x % 3 == 0,
                                        range(len(seq_ungapped))))
            # check alphabet
            # Not use Alphabet._verify_alphabet function because it
            # only works for single alphabet
            for i in self.rf_table:
                if self._data[i:i + 3] not in alphabet.letters and not self.is_ambiguous_codon(self._data[i:i + 3]):
                    raise ValueError("Sequence contain undefined letters from"
                                     " alphabet "
                                     "({0})! ".format(self._data[i:i + 3]))
        else:
            # if gap_char in self._data:
            #    assert  len(self) % 3 == 0, \
            #            "Gapped sequence length is not a triple number"
            assert isinstance(rf_table, (tuple, list)), \
                    "rf_table should be a tuple or list object"
            assert all(isinstance(i, int) for i in rf_table), \
                    "elements in rf_table should be int that specify " \
                  + "the codon positions of the sequence"
            seq_ungapped = self._data.replace(gap_char, "")
            for i in rf_table:
                if seq_ungapped[i:i + 3] not in alphabet.letters:
                    raise ValueError("Sequence contain undefined letters "
                                     "from alphabet "
                                     "({0})!".format(seq_ungapped[i:i + 3]))
            self.rf_table = rf_table

    def is_ambiguous_codon(self, codon):
        return ('N' in codon.upper()) and self.enable_undef


class Settings():

    def __init__(self):
        pass

    def set(self, **kwargs):
        EXCLUDE_AA = kwargs.get('EXCLUDE_AA', default_settings.EXCLUDE_AA)
        self.EXCLUDE_AA_FROM = kwargs.get(
            'EXCLUDE_AA_FROM', default_settings.EXCLUDE_AA_FROM)
        self.AA_LETTERS = "".join(
            [aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa not in EXCLUDE_AA])
        self.OUTDIR = kwargs.get('OUTDIR', default_settings.OUTDIR)
        self.AA_MAJORITY_THRESH = kwargs.get(
            'AA_MAJORITY_THRESH', default_settings.AA_MAJORITY_THRESH)
        self.LIMIT_TO_SUSPECTED_SPECIES = kwargs.get(
            'LIMIT_TO_SUSPECTED_SPECIES', default_settings.LIMIT_TO_SUSPECTED_SPECIES)
        self.FREQUENCY_THRESHOLD = kwargs.get(
            'FREQUENCY_THRESHOLD', default_settings.FREQUENCY_THRESHOLD)
        self.COUNT_THRESHOLD = kwargs.get(
            'COUNT_THRESHOLD', default_settings.COUNT_THRESHOLD)
        self.GENETIC_CODE = kwargs.get(
            'GENETIC_CODE', default_settings.GENETIC_CODE)
        self.SHOW_MIXTE_CODONS = kwargs.get(
            'SHOW_MIXTE_CODONS', default_settings.SHOW_MIXTE_CODONS)
        self.SHOW_FILTERED_CODON_DATA = kwargs.get(
            'SHOW_FILTERED_CODON_DATA', default_settings.SHOW_FILTERED_CODON_DATA)
        self.IMAGE_FORMAT = kwargs.get(
            'IMAGE_FORMAT', default_settings.IMAGE_FORMAT)
        self.ADD_LABEL_TO_LEAF = kwargs.get(
            'ADD_LABEL_TO_LEAF', default_settings.ADD_LABEL_TO_LEAF)
        self.USE_CONSENSUS_FOR_LIKELIHOOD = kwargs.get(
            'USE_CONSENSUS_FOR_LIKELIHOOD', default_settings.USE_CONSENSUS_FOR_LIKELIHOOD)

        self.USE_GLOBAL = kwargs.get(
            'USE_GLOBAL', default_settings.USE_GLOBAL)
        self.mode = kwargs.get('MODE', default_settings.MODE)
        self.method = kwargs.get('MATRIX', default_settings.MATRIX)
        self.hmmbuild = kwargs.get('hmmbuild', 'hmmbuild')
        self.eslalimask = kwargs.get('eslalimask', 'esl-alimask')
        self.eslalimanip = kwargs.get('eslalimanip', 'esl-alimanip')
        self.hmmalign = kwargs.get('hmmalign', 'hmmalign')
        self.conf = kwargs.get('CONF', default_settings.CONF)
        # matrix used to compute score (close to likelihood score)
        self.SUBMAT = getattr(MatrixInfo, 'blosum62')

    def fill(self, params):
        self.__dict__.update(params.__dict__)

    def update_params(self, **kwargs):
        for k, v in kwargs.items():
            if k == 'SUBMAT':
                self.__dict__[k] = getattr(MatrixInfo, v)
            else:
                self.__dict__[k] = v


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

    def __init__(self, aas, alignment, consensus, codon_align_dict, dct, positions, genelimit, settings):

        # change aa2 is a list now
        self.aa1, aareas = aas
        self.aa2_list = aareas.keys()
        self.codon_alignment = codon_align_dict
        self.dct = dct
        self.init_back_table()

        self.alignment = alignment
        if not isinstance(self.codon_alignment, dict):
            self.codon_alignment = SeqIO.to_dict(self.codon_alignment)
        self.consensus = consensus

        self.reacodons = defaultdict(partial(defaultdict, Counter))
        self.usedcodons = defaultdict(partial(defaultdict, Counter))
        self.mixtecodon = defaultdict(partial(defaultdict, Counter))
        # find codon usage in the sequence
        self.settings = settings
        self.subsmat = settings.SUBMAT
        self.cons_for_lik = settings.USE_CONSENSUS_FOR_LIKELIHOOD

        # print len(codon_align_dict.values()[0])
        self.specs_amino_count = defaultdict(
            partial(defaultdict, Counter))  # this one is ok
        self.positions = positions
        self.genelimit = genelimit
        self.t_genelimit = len(genelimit)
        self.codons_distribution = defaultdict(
            partial(defaultdict, partial(defaultdict, set)))
        self.total_codons_distribution = defaultdict(
            partial(defaultdict, partial(defaultdict, set)))
        self.codon_map_in_genome = defaultdict(partial(defaultdict, Counter))

        self._spec_codon_usage()

    def init_back_table(self):
        self.back_table = defaultdict(list)
        for aa, codon in zip(self.dct.forward_table.values(), self.dct.forward_table.keys()):
            self.back_table[aa].append(codon)

    def _spec_codon_usage(self):
        for i, aa in enumerate(self.consensus):
            for spec in self.codon_alignment.keys():
                spec_codon = self.codon_alignment[spec].seq.get_codon(i)
                spec_aa = self.dct.forward_table[
                    spec_codon] if '-' not in spec_codon else None
                cur_pos = self.positions[i]
                if(spec_aa):
                    self.specs_amino_count[spec][spec_aa].update([spec_codon])
                    aa_set_list = aa
                    if not self.cons_for_lik:
                        aa_set_list = self.alignment[:, i]
                    self.codon_map_in_genome[spec][
                        spec_codon].update(aa_set_list)

                    # only check amino acid that are suspected
                    if spec_aa in self.aa2_list:
                        lim_start = 0
                        while lim_start < self.t_genelimit and cur_pos > self.genelimit[lim_start][2]:
                            lim_start += 1

                        if lim_start < self.t_genelimit:
                            self.total_codons_distribution[spec][spec_aa][
                                    spec_codon].add(self.genelimit[lim_start][0])
                        # potential reassignment counter
                        # species use aa2 while aa1 is prevalent
                        if aa == self.aa1:
                            self.reacodons[spec][spec_aa].update([spec_codon])
                            if lim_start < self.t_genelimit:
                                self.codons_distribution[spec][spec_aa][
                                    spec_codon].add(self.genelimit[lim_start][0])
                    
                        # species use aa2 with aa2 being the prevalent aa
                        elif aa == spec_aa:
                            self.usedcodons[spec][spec_aa].update([spec_codon])

                        # other position where aa2 is used in species 
                        # this mean, species use aa2, while we don't care about the major aa
                        else : 
                            self.mixtecodon[spec][spec_aa].update([spec_codon])


    def get_score(self, spec, aa_ori, aa_rea):
        codon_score = {}
        for codon in self.back_table[aa_ori]:
            amino_counters = self.codon_map_in_genome[spec][codon]
            total = 0.0
            numerator = 0
            for k, v in amino_counters.items():
                if k != '-':
                    total += v
                    try:
                        cost = self.subsmat[(aa_rea, k)]
                    except:
                        cost = self.subsmat[(k, aa_rea)]
                    numerator += cost * v

            if total > 0:
                codon_score[codon] = numerator / total
            else:
                codon_score[codon] = np.inf
        return codon_score

    def __getitem__(self, key):
        return self.get_reacodons(key)

    def get_reacodons(self, specie, aa):
        return self.reacodons[specie][aa]

    def get_mixtecodons(self, specie, aa):
        return self.mixtecodon[specie][aa]

    def get_usedcodons(self, specie, aa):
        return self.usedcodons[specie][aa]

    def get_aa_usage(self, specie, aa):
        return self.specs_amino_count[specie][aa]

    def get_all_aas_usage(self, specie):
        return self.specs_amino_count[specie]

    def get_rea_aa_codon_distribution(self, specie, aa):
        return dict((k, len(v)) for k, v in self.codons_distribution[specie][aa].items())

    def get_total_rea_aa_codon_distribution(self, specie, aa):
        return dict((k, len(v)) for k, v in self.total_codons_distribution[specie][aa].items())

    
class NaiveFitch(object):
    """A NaiveFitch algorithm for finding the most parcimonious solution"""

    def __init__(self, tree, reassigned, ori_aa, dest_aa, settings, dct, codon_rea=(None, None)):
        self.id = {}
        self.tree = tree
        self.settings = settings
        self.corr = {'0': ori_aa, '1': dest_aa}
        self.codon_rea_global, self.codon_rea_filtered = codon_rea
        colors = ['#a6cee3', '#1f78b4', '#b2df8a',
                  '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f']
        self.dct = dct
        self.init_back_table()
        codon_list = self.back_table[aa_letters_3to1[ori_aa]]
        self.colors = dict(zip(codon_list, colors))
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

    def init_back_table(self):
        self.back_table = defaultdict(list)
        for aa, codon in zip(self.dct.forward_table.values(), self.dct.forward_table.keys()):
            self.back_table[aa].append(codon)

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
            elif l.lost and (l.state == self.dest_aa or l.count < self.settings.COUNT_THRESHOLD):
                l.add_features(reassigned={0})
                l.add_features(rea=self.corr['0'])
                l.add_features(state=self.ori_aa)
            elif not l.lost and l.state == self.ori_aa:
                l.add_features(reassigned={1})
                l.add_features(rea=self.corr['1'])
                l.add_features(state=self.dest_aa)
        for node in tmptree:
            if node.rea == self.dest_aa:
                return True
        return False

    def _bottomup(self, tree):
        """Fitch algorithm part 1 : bottom up"""
        for node in tree.traverse('postorder'):
            if not node.is_leaf():
                intersect = []
                union = set()
                for c in node.get_children():
                    union = union | c.reassigned
                    intersect.append(c.reassigned)
                intersect = set.intersection(*intersect)
                if(intersect):
                    node.add_features(reassigned=intersect)

                else:
                    node.add_features(reassigned=union)

                node.add_features(
                    rea="/".join([self.corr[str(r)] for r in node.reassigned]))

    def _topdown(self, tree):
        """Fitch algorithm part 2 : top down"""
        pass

    def is_reassigned(cls, node, strict=True):
        """ return True if a node has undergoned reassignment """
        if "reassigned" in node.features and (node.reassigned == {1} or (not strict and 1 in node.reassigned)):
            return True
        return False

    def get_species_list(self, limit_to_reassigned=False):
        """Get the species list for which we are almost certain
        that there is a reassignment
        """
        slist = set()
        if limit_to_reassigned:
            for node in self.tree.traverse():
                if 'reassigned' in node.features and node.reassigned == {1}:
                    slist.update(node.get_leaf_names())
        else:
            slist = set(self.tree.get_tree_root().get_leaf_names())
        return slist

    def get_distance_to_rea_node(self, node):
        cur_node = node
        if isinstance(node, str):
            cur_node = self.tree & node
        dist = 0
        while cur_node != None and not self.is_reassigned(cur_node, strict=False):
            dist += 1
            cur_node = cur_node.up
        return dist

    def has_codon_data(self):
        return (self.codon_rea_global, self.codon_rea_filtered) != (None, None)

    def render_tree(self, output="", pie_size=45):

        GRAPHICAL_ACCESS = True
        try:
            from ete3 import TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
        except ImportError, e:
            GRAPHICAL_ACCESS = False

        if not output:
            output = os.path.join(self.settings.OUTDIR, self.ori_aa +
                                  "_to_" + self.dest_aa + "." + self.settings.IMAGE_FORMAT)

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
            other_style["fgcolor"] = "seagreen"
            other_style["node_bgcolor"] = "crimson"
            other_style["hz_line_type"] = 0

            prob_lost = NodeStyle()
            prob_lost["hz_line_type"] = 1
            prob_lost["hz_line_color"] = "#cccccc"

            def layout(node):

                get_suffix = lambda x: x if self.settings.ADD_LABEL_TO_LEAF else ""
                faces.add_face_to_node(
                    AttrFace("rea", fsize=10), node, column=0, position="branch-right")

                has_count, has_fcount = 0, 0
                if node.is_leaf() and ('count' in node.features):
                    faces.add_face_to_node(AttrFace(
                        "count", fsize=6, fgcolor="firebrick"), node, column=0, position="branch-bottom")
                    has_count = node.count

                if node.is_leaf() and ('filter_count' in node.features):
                    faces.add_face_to_node(AttrFace(
                        "filter_count", fsize=6, fgcolor="indigo"), node, column=0, position="branch-top")
                    has_fcount = node.filter_count

                if node.is_leaf():
                    if not (has_fcount or has_count):
                        faces.add_face_to_node(
                            AttrFace("name"), node, 0, position="aligned")

                    else:
                        # if lost in node.features, then node is already a leaf
                        # color change only concern leaf with count
                        if(self.is_reassigned(node)):
                            # was suspected and passed fisher test
                            if 'lost' in node.features and not node.lost:
                                faces.add_face_to_node(AttrFace("name", fgcolor="#ff1111", text_suffix=get_suffix(
                                    "_R")), node, column=0, position="aligned")
                            # was suspected and failed fisher test
                            else:
                                faces.add_face_to_node(AttrFace("name", fgcolor="#1a9850", text_suffix=get_suffix(
                                    "_G")), node, column=0, position="aligned")

                        else:
                            # was not suspected and passed fisher test
                            if 'lost' in node.features and not node.lost:
                                faces.add_face_to_node(AttrFace("name", fgcolor="#fc8d59", fstyle="italic", text_suffix=get_suffix(
                                    "_O")), node, column=0, position="aligned")
                            # was not suspected and failed fisher test
                            else:
                                faces.add_face_to_node(AttrFace("name", fgcolor="#cccccc", text_suffix=get_suffix(
                                    "_Gr")), node, column=0, position="aligned")

                if(self.has_codon_data() and node.is_leaf()):
                    node_colors = aa_letters_3to1[self.ori_aa]
                    spec_codonrea_g = self.codon_rea_global.get_reacodons(
                        node.name, aa_letters_3to1[self.ori_aa])
                    spec_codonused_g = self.codon_rea_global.get_usedcodons(
                        node.name, aa_letters_3to1[self.ori_aa])

                    faces.add_face_to_node(PPieChartFace(spec_codonrea_g.values(), pie_size, pie_size, labels=[],
                                                         colors=[self.colors[k] for k in spec_codonrea_g.keys()]), node, column=1, position="aligned")

                    faces.add_face_to_node(PPieChartFace(spec_codonused_g.values(), pie_size, pie_size, labels=[],
                                                         colors=[self.colors[k] for k in spec_codonused_g.keys()]), node, column=2, position="aligned")
                    next_column = 3
                    if(self.settings.SHOW_MIXTE_CODONS):
                        spec_mixtecodon_g = self.codon_rea_global.get_mixtecodons(
                            node.name, aa_letters_3to1[self.ori_aa])
                        faces.add_face_to_node(PPieChartFace(spec_mixtecodon_g.values(), pie_size, pie_size, labels=[],
                                                             colors=[self.colors[k] for k in spec_mixtecodon_g.keys()]), node, column=next_column, position="aligned")
                        next_column = 4

                    if(self.settings.SHOW_FILTERED_CODON_DATA):

                        spec_codonrea_f = self.codon_rea_filtered.get_reacodons(
                            node.name, aa_letters_3to1[self.ori_aa])
                        spec_codonused_f = self.codon_rea_filtered.get_usedcodons(
                            node.name, aa_letters_3to1[self.ori_aa])

                        # add separator
                        faces.add_face_to_node(LineFace(
                            pie_size, pie_size, None), node, column=next_column, position="aligned")

                        # add data
                        faces.add_face_to_node(PPieChartFace(spec_codonrea_f.values(), pie_size, pie_size, labels=[],
                                                             colors=[self.colors[k] for k in spec_codonrea_f.keys()]), node, column=next_column + 1, position="aligned")

                        faces.add_face_to_node(PPieChartFace(spec_codonused_f.values(),  pie_size, pie_size, labels=[],
                                                             colors=[self.colors[k] for k in spec_codonused_f.keys()]), node, column=next_column + 2, position="aligned")

                        if(self.settings.SHOW_MIXTE_CODONS):
                            spec_mixtecodon_f = self.codon_rea_filtered.get_mixtecodons(
                                node.name, aa_letters_3to1[self.ori_aa])
                            faces.add_face_to_node(PPieChartFace(spec_mixtecodon_f.values(), pie_size, pie_size, labels=[],
                                                                 colors=[self.colors[k] for k in spec_mixtecodon_f.keys()]), node, column=next_column + 3, position="aligned")

            ts.layout_fn = layout
            ts.title.add_face(
                TextFace(self.ori_aa + " --> " + self.dest_aa, fsize=14), column=0)
            if(self.has_codon_data()):
                for cod, col in self.colors.items():
                    ts.legend.add_face(CircleFace(
                        (pie_size / 3), col), column=0)
                    ts.legend.add_face(
                        TextFace("  " + cod + " ", fsize=8), column=1)
                ts.legend_position = 4

            # Apply node style
            for n in self.tree.traverse():
                if not n.is_leaf():
                    if n.reassigned == {1}:
                        n.set_style(rea_style)
                    elif len(n.reassigned) > 1:
                        n.set_style(other_style)
                if 'lost' in n.features and n.lost:
                    n.set_style(prob_lost)

            # Save tree as figure

            self.tree.render(output, dpi=400, tree_style=ts)


class SequenceSet(object):

    def __init__(self, coreinstance, phylotree, table_num):
        # using sequence set suppose that the alignment and the protein
        # concatenation was already done
        self.codontable = CodonTable.unambiguous_dna_by_id[abs(table_num)]

        self.prot_dict, self.dna_dict, self.gene_limits = coreinstance.concat()
        self.prot_align = MultipleSeqAlignment(self.prot_dict.values())
        self.compute_current_mat()
        self.core = coreinstance
        self.phylotree = phylotree
        self.common_genome = []
        self.restrict_to_common()
        self.codon_align()


    def compute_current_mat(self):
        suminfo = AlignInfo.SummaryInfo(self.prot_align)
        replace_info = suminfo.replacement_dictionary()
        self.protarm = SubsMat.SeqMat(replace_info)
        self.obs_freq_mat = SubsMat. _build_obs_freq_mat(self.protarm)
        self.exp_freq_mat = SubsMat._build_exp_freq_mat(
            SubsMat. _exp_freq_table_from_obs_freq(self.obs_freq_mat))
        self.subs_mat = SubsMat._build_subs_mat(
            self.obs_freq_mat, self.exp_freq_mat)
        self.lo_mat = SubsMat._build_log_odds_mat(self.subs_mat)

    def get_genes_per_species(self):
        """ Return the number of genes in each species"""
        gene_in_spec = {}
        for spec in self.common_genome:
            gene_in_spec[spec] = np.sum(
                [1 for gene, speclist in self.core.common_spec_per_gene.items() if spec in speclist])
        return gene_in_spec

    def restrict_to_common(self):
        """ Keep only the species present in the dna, prot seq and in the tree """

        dna_set = set(self.dna_dict.keys())

        species_set = set(self.phylotree.get_leaf_names())
        common_genome = dna_set.intersection(species_set)
        speclen = len(species_set)
        # we know that prot and dna dict contains the same species (keys)
        # remove every dna sequence that are not cds
        # If you put stop codon in protein, it should be in nucleotide also
        common_genome = set([x for x in common_genome if (
            len(self.dna_dict[x]) == len(self.prot_dict[x].seq.ungap('-')) * 3)])
        # print [(x, len(dna_dict[x]), len(prot_dict[x].seq.ungap('-'))+1, len(dna_dict[x])==(len(prot_dict[x].seq.ungap('-'))+1)*3) for x in common_genome ]
        # common_genome = [x for x in common_genome if len(dna_dict[x])==(len(prot_dict[x].seq.ungap('-'))+1)*3 ]
        # return common_genome, dict((k, v) for k, v in dna_dict.items() if k
        # in common_genome)
        if not common_genome:
            raise ValueError(
                'ID intersection for dna, prot and species tree is empty')
        if len(dna_set) != speclen or len(species_set) != speclen:
            logging.debug(
                'Non-uniform sequence in sequences and tree. The tree will be pruned.')

        # prune tree to sequence list
        self.phylotree.prune(common_genome)
        self.common_genome = common_genome
        logging.debug("List of common genome")
        logging.debug(common_genome)
        # return common_genome, dict((k, v) for k, v in dna_dict.items() if k
        # in common_genome), prot_dict
        self.dna_dict, self.prot_dict = dict((k, v) for k, v in self.dna_dict.items(
        ) if k in common_genome), dict((k, v) for k, v in self.prot_dict.items() if k in common_genome)

        for k, v in self.prot_dict.items():
            v.seq.alphabet = generic_protein
            v.id = k

        self.prot_align = MultipleSeqAlignment(
            self.prot_dict.values(), alphabet=alpha)

    @classmethod
    def copy_codon_record(clc, record, codontable, gap_char='-'):
        """Return a Codon seq sequence from a dna sequence"""
        alphabet = get_codon_alphabet(codontable, gap_char=gap_char)
        return SeqRecord(CodonSeq(record.seq._data, alphabet=alphabet), id=record.id)

    def copy_codon_alignment(self, codon_alignment, alphabet=default_codon_alphabet):
        """Return a codon alignment object from a list of codon seq sequences"""
        return codonalign.CodonAlignment((self.copy_codon_record(rec, self.codontable) for rec in codon_alignment._records), alphabet=alphabet)

    def codon_align(self, alphabet=default_codon_alphabet, gap_char='-'):
        """ Perform a codon alignment based on a protein multiple alignment
        and remove gaped positions
        """
        alphabet = get_codon_alphabet(self.codontable, gap_char=gap_char)
        if self.common_genome is None:
            self.restrict_to_common()
        # build codon alignment and return it
        codon_aln = []
        undef_codon = set([])
        for g in self.dna_dict.keys():
            codon_rec, undef_c = self._get_codon_record(self.dna_dict[g], self.prot_dict[g], self.codontable, alphabet)
            undef_codon.update(undef_c)
            codon_aln.append(codon_rec)

        self.codon_alignment = codonalign.CodonAlignment(
            codon_aln, alphabet=alphabet)
        if undef_codon:
            undef_codon = sorted(list(undef_codon))
            R_aa_position = np.asarray(xrange(current_alignment.get_alignment_length()))
            R_aa_position = np.setxor1d(tt_filter_position, np.asarray(undef_codon))
            self.prot_align = self.filter_align_position(self.prot_align, R_aa_position, alphabet=generic_protein)
            self.prot_dict = SeqIO.to_dict(self.prot_align)
            self.codon_alignment = next(self.filter_codon_alignment((R_aa_position,)))
            # remove all the position with undef codon from the dna_dict
            for k in dna_dict.keys():
                cur_seq = dna_dict[k].seq
                dna_dict[k] = SeqRecord(Seq("".join([cur_seq[x*3:(x+1)*3] for x in undef_codon]), generic_nucleotide), id=k, name=k)
            
        #self.codon_alignment = codonalign.build(self.prot_align, self.dna_dict, codon_table=self.codontable)

    def _get_codon_record(self, dnarec, protrec, codontable, alphabet, gap_char='-'):
        nuc_seq = dnarec.seq.ungap(gap_char)
        codon_seq = ""
        aa_num = 0
        max_error = 1000
        x_undecoded = []
        for aa in protrec.seq:
            if aa == gap_char:
                codon_seq += gap_char * 3
            else:
                next_codon = nuc_seq._data[aa_num * 3: (aa_num + 1) * 3]

                if 'N' in next_codon.upper():
                    x_undecoded.append(aa_num)
                    if aa.upper() != 'X':
                        logging.warn("%s(%s %d) decoded by undefined codon %s(%s)"
                                 % (protrec.id, aa, aa_num, dnarec.id, next_codon))
                        max_error -= 1
                elif not str(Seq(next_codon.upper()).translate(table=codontable)) == aa.upper():
                    logging.warn("%s(%s %d) does not correspond to %s(%s)"
                                 % (protrec.id, aa, aa_num, dnarec.id, next_codon))
                    max_error -= 1
                if max_error <= 0:
                    raise ValueError(
                        "You're obviously using the wrong genetic code or you have too much undefined codon coding for Amino acid")
                codon_seq += next_codon
                aa_num += 1
        return SeqRecord(CodonSeq(codon_seq, alphabet, enable_undef=True), id=dnarec.id), x_undecoded

    def filter_codon_alignment(self, ind_array=None, get_dict=False, alphabet=default_codon_alphabet):
        """Return the codon aligment from a list of in array"""
        if not ind_array:
            ind_array = (self.filt_position, )

        if not isinstance(ind_array, tuple):
            raise ValueError('Expect tuple for ind_array, got %s' %
                             (type(ind_array)))
        else:
            filt_codon_align = self.copy_codon_alignment(self.codon_alignment)

            for indexes in ind_array:
                indexes = sorted(indexes)
                for i in xrange(len(filt_codon_align)):
                    codseq = CodonSeq('')
                    for pos in indexes:
                        codseq += filt_codon_align[i].seq.get_codon(pos)
                    filt_codon_align[i].seq = codseq
                if get_dict:
                    filt_codon_align = SeqIO.to_dict(filt_codon_align)
                yield filt_codon_align

    def get_codon_alignment(self):
        r = self.filter_codon_alignment()
        self.fcodon_alignment = next(r)
        return self.codon_alignment, self.fcodon_alignment

    def prot_filtering(self, id_thresh=None, gap_thresh=None, ic_thresh=None, rmcnst=True):

        current_alignment = self.prot_align
        tt_filter_position = np.asarray(
            xrange(current_alignment.get_alignment_length()))
        self.position = np.asarray(
            xrange(current_alignment.get_alignment_length()))

        if(gap_thresh):
            self._gap_alignment, self._gap_filtered_position = self.clean_alignment(
                current_alignment, threshold=(abs(gap_thresh) <= 1 or 0.01) * abs(gap_thresh))
            current_alignment = self._gap_alignment
            logging.debug("Alignment length after removing gaps : %d" %
                          current_alignment.get_alignment_length())
            tt_filter_position = tt_filter_position[
                self._gap_filtered_position]

        if(ic_thresh):
            align_info = AlignInfo.SummaryInfo(current_alignment)
            ic_content = align_info.information_content()
            max_val = max(align_info.ic_vector.values()) * \
                ((abs(ic_thresh) <= 1 or 0.01) * abs(ic_thresh))
            ic_pos = (np.asarray(align_info.ic_vector.values())
                      >= max_val).nonzero()
            logging.debug(
                "Filtering with ic_content, vector to discard is %s" % str(ic_pos[0]))
            # ic_pos here is a tuple, so we should access the first element
            self._ic_alignment = self.filter_align_position(
                current_alignment, ic_pos[0])
            self._ic_filtered_positions = ic_pos
            current_alignment = self._ic_alignment
            tt_filter_position = tt_filter_position[
                self._ic_filtered_positions]
            logging.debug("Alignment length filtering with information content: %d" %
                          current_alignment.get_alignment_length())

        if(id_thresh):
            self._id_alignment, self._id_filtered_position = self.filter_alignment(
                current_alignment, remove_identity=rmcnst, threshold=(abs(id_thresh) <= 1 or 0.01) * abs(id_thresh))
            current_alignment = self._id_alignment
            tt_filter_position = tt_filter_position[self._id_filtered_position]
            logging.debug("Alignment length after filtering by sequence identity : %d" %
                          current_alignment.get_alignment_length())

        self.filt_prot_align = current_alignment
        self.filt_position = tt_filter_position

    def write_data(self, id_filtered=None, gap_filtered=None, ic_filtered=None, tree=None):
        """Save file depending on the file use as input"""
        if id_filtered:
            AlignIO.write(self._id_alignment, open(id_filtered, 'w'), 'fasta')
        if gap_filtered:
            AlignIO.write(self._gap_alignment, open(
                gap_filtered, 'w'), 'fasta')
        if ic_filtered:
            AlignIO.write(self._ic_alignment, open(ic_filtered, 'w'), 'fasta')
        if tree:
            self.phylotree.write(outfile=tree)

    @classmethod
    def filter_align_position(clc, alignment, index_array, alphabet=generic_protein):
        """Keep only columns specified by index_array from the alignment"""
        edited_alignment = alignment[:, :]
        index_array = sorted(index_array)
        # alignment[:, slice(index_array[0], index_array[0] + 1)]
        for seqrec in edited_alignment:
            seq = Seq('', alphabet)
            for i in index_array:
                seq += seqrec[i]
            seqrec.letter_annotations = _RestrictedDict(length=len(seq))
            seqrec.seq = seq.upper()

        return edited_alignment

    @classmethod
    def clean_alignment(clc, alignment=None, characs=['-'], threshold=0.5):
        """Remove position of alignment which contain character from characs"""
        align_array = np.array([list(rec) for rec in alignment], np.character)
        indel_array = np.where((np.mean(np.in1d(align_array,
                                                characs).reshape(align_array.shape), axis=0) >= threshold) == False)[0].tolist()

        return clc.filter_align_position(alignment, indel_array), indel_array

    @classmethod
    def filter_alignment(clc, alignment, threshold=0.8, remove_identity=False, ambiguous='X', alphabet=generic_protein):
        """Filter an alignment using threshold as the minimum aa identity per columns"""
        aligninfo = AlignInfo.SummaryInfo(alignment)
        # Smart move : Make a consensus sequence with threshold
        # remove all gap position
        consensus = aligninfo.gap_consensus(
            threshold=threshold, ambiguous=ambiguous, consensus_alpha=alphabet)
        cons_array = [i for i, c in enumerate(consensus) if c != ambiguous]

        if(remove_identity):
            matched_consensus = aligninfo.gap_consensus(
                threshold=1, ambiguous=ambiguous, consensus_alpha=alphabet)
            cons_array = [i for i, c in enumerate(
                consensus) if c != ambiguous and matched_consensus[i] == ambiguous]

        return clc.filter_align_position(alignment, cons_array), cons_array

    @classmethod
    def get_identity(clc, seq1, seq2, ignore_gap=True, percent=True):
        """Return percent of columns with same aa in a multiple alignment"""
        assert len(seq1) == len(seq2), "Sequence not aligned"
        identity = 0.0
        for i, character in enumerate(seq1):
            if not ignore_gap or character not in ['-', '.']:
                if(character == seq2[i]):
                    identity += 1
        return identity * (99.0 * percent + 1.0) / len(seq1)

    @classmethod
    def get_consensus(clc, alignment, threshold, ambiguous='X', alphabet=generic_protein):
        """return consensus, using a defined threshold"""
        aligninfo = AlignInfo.SummaryInfo(alignment)
        consensus = aligninfo.gap_consensus(
            threshold=threshold, ambiguous=ambiguous, consensus_alpha=alphabet)
        return consensus

    @classmethod
    def get_aa_filtered_alignment(clc, consensus, aa):
        """Get the filtered alignment for an amino acid"""
        aa = aa.upper()
        cons_array = [i for i, c in enumerate(consensus) if c.upper() == aa]
        return cons_array

    def get_total_genomes(self):
        """Return the total number of genomes. This should be genome present in dna, prot and tree """
        return len(self.common_genome)

    def get_genome_size(self, aligntype='global'):
        use_alignment = self.prot_align
        if aligntype == 'filtered':
            use_alignment =  self.filt_prot_align
        gsize = dict((k, len(v.seq.ungap())) for k, v in SeqIO.to_dict(use_alignment))
        return gsize

class ReaGenomeFinder:

    def __init__(self, seqset, settings):
        self.seqset = seqset
        self.mode = getattr(settings, 'mode', "count")
        self.method = getattr(settings, 'method', "identity")
        self.confd = getattr(settings, 'conf', 0.05)
        self.suspected_species = defaultdict(dict)
        self.aa2aa_rea = defaultdict(dict)
        self.settings = settings

        def makehash():
            return defaultdict(makehash)
        self.reassignment_mapper = makehash()

    def compute_sequence_identity(self, matCalc=None):
        if not matCalc:
            matCalc = DistanceCalculator(self.method)
        self.global_paired_distance = matCalc.get_distance(
            self.seqset.prot_align)
        self.filtered_paired_distance = matCalc.get_distance(
            self.seqset.filt_prot_align)

    def get_genomes(self, use_similarity=1):
        matCalc = DistanceCalculator(self.method)
        self.compute_sequence_identity(matCalc)

        self.filtered_consensus = self.seqset.get_consensus(
            self.seqset.filt_prot_align, self.settings.AA_MAJORITY_THRESH)
        self.global_consensus = self.seqset.get_consensus(
            self.seqset.prot_align, self.settings.AA_MAJORITY_THRESH)

        logging.debug("Filtered alignment consensus : \n%s\n" %
                      self.filtered_consensus)
        number_seq = self.seqset.get_total_genomes()
        self.seq_names = list(self.seqset.common_genome)
        self.sim_json = defaultdict(list)

        self.seqset.aa_filt_prot_align = {}
        self.aa_paired_distance = {}
        for aa in self.settings.AA_LETTERS:
            cons_array = self.seqset.get_aa_filtered_alignment(
                self.filtered_consensus, aa)
            if(cons_array):
                self.seqset.aa_filt_prot_align[aa_letters_1to3[
                    aa]] = SequenceSet.filter_align_position(self.seqset.filt_prot_align, cons_array)
                self.aa_paired_distance[aa] = matCalc.get_distance(
                    self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]])
                self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]] = SeqIO.to_dict(
                    self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]])

        self.aa_sim_json = defaultdict(list)
        aa2suspect = defaultdict(Counter)
        aa2suspect_dist = defaultdict(partial(defaultdict, list))
        for (j, i) in itertools.combinations(xrange(number_seq), r=2):
            gpaired = abs(
                use_similarity - self.global_paired_distance[self.seq_names[i], self.seq_names[j]])
            fpaired = abs(
                use_similarity - self.filtered_paired_distance[self.seq_names[i], self.seq_names[j]])
            self.sim_json[self.seq_names[i]].append(
                {"global": gpaired, "filtered": fpaired, "species": self.seq_names[j]})
            paired = fpaired
            if self.settings.USE_GLOBAL:
                paired = gpaired

            for aa in self.aa_paired_distance.keys():
                aapaired = abs(
                    use_similarity - self.aa_paired_distance[aa][self.seq_names[i], self.seq_names[j]])
                self.aa_sim_json[aa_letters_1to3[aa]].append(
                    {'global': gpaired, 'filtered': aapaired, "species": "%s||%s" % (self.seq_names[i], self.seq_names[j])})
                if (paired > aapaired and use_similarity) or (paired < aapaired and not use_similarity):
                    aa2suspect[aa][self.seq_names[i]] = aa2suspect[
                        aa].get(self.seq_names[i], 0) + 1
                    aa2suspect[aa][self.seq_names[j]] = aa2suspect[
                        aa].get(self.seq_names[j], 0) + 1

                aa2suspect_dist[aa][self.seq_names[
                    i]].append([paired, aapaired])
                aa2suspect_dist[aa][self.seq_names[
                    j]].append([paired, aapaired])

        if self.mode == 'count':
            self.get_suspect_by_count(aa2suspect, number_seq)
        elif self.mode in ['wilcoxon', 'mannwhitney', 'ttest']:
            self.get_suspect_by_stat(
                aa2suspect_dist, number_seq, use_similarity, test=self.mode, confd=self.confd)
        else:
            self.get_suspect_by_clustering(aa2suspect_dist, number_seq)
        logging.debug('The list of suspected species per amino acid is:')
        logging.debug(self.suspected_species)

    def get_suspect_by_count(self,  aa2suspect, seq_num):
        for aa in aa2suspect.keys():
            tmp_list = aa2suspect[aa].most_common()
            i = 0
            while i < len(tmp_list) and tmp_list[i][1] > self.settings.FREQUENCY_THRESHOLD * seq_num:
                self.suspected_species[aa][tmp_list[i][0]] = 1 # tmp_list[i][1]
                i += 1

    def get_suspect_by_stat(self, aa2suspect_dist, seq_num, use_similarity=1, test='wilcoxon', confd=0.05):
        for aa in aa2suspect_dist.keys():
            tmp_list = aa2suspect_dist[aa]
            for seq, val in tmp_list.items():
                val = np.asarray(val)
                rank, pval = self.get_paired_test(
                    val[:, 0], val[:, 1], use_similarity, test)
                if pval <= confd:
                    self.suspected_species[aa][seq] = 1# pval

    def get_suspect_by_clustering(self, aa2suspect_dist, number_seq, use_similarity=1):
        logging.debug('Clustering chosen, labels for each species')
        for aa in aa2suspect_dist.keys():
            aafilt2dict = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            tmp_list = aa2suspect_dist[aa]
            cluster_list = {}
            aa_set = set([])
            for seq, val in tmp_list.items():
                val = np.mean(val, axis=0)
                if (val[0] > val[1] and use_similarity) or (val[0] < val[1] and not use_similarity):
                    aa_set.add(seq)
                aa_freq = np.mean([x == aa for x in aafilt2dict[seq]])
                cluster_list[seq] = np.append(val, aa_freq)

            centroid, labels = self.get_cluster(cluster_list)
            logging.debug("======> Amino acid %s :" % aa_letters_1to3[aa])
            for index, elem in enumerate(cluster_list):
                logging.debug("%s\t%s\t%d" %
                              (elem, cluster_list[elem], labels[index]))

            lab1 = set(np.asarray(cluster_list.keys())[labels == 1])
            lab2 = set(np.asarray(cluster_list.keys())[labels == 0])
            size1 = aa_set.intersection(lab1)
            size2 = aa_set.intersection(lab2)
            lab = lab1 if len(size1) > len(size2) else lab2

            for specie in lab:
                self.suspected_species[aa][specie] = 1

    def get_cluster(self, cluster_list):
        features = np.asarray(cluster_list.values())
        return kmeans2(features, 2, iter=100)

    @classmethod
    def get_paired_test(clc, y1, y2, use_similarity, test="wilcoxon"):
        if test == 'wilcoxon':
            r, pval = ss.wilcoxon(y1, y2, correction=True)
        elif test == 'ttest':
            r, pval = ss.ttest_rel(y1, y2)
        else:
            r, pval = ss.mannwhitneyu(y1, y2, use_continuity=True)
        
        if use_similarity:
            return r, pval / 2.0
        else:
            return r, 1 - pval / 2.0

    # TODO :  CHANGE THIS FUNCTION
    def possible_aa_reassignation(self):
        """Find to which amino acid the reassignation are probable"""
        for aa, suspect in self.suspected_species.items():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            suspected_aa = []
            for spec in suspect.keys():
                for cur_aa in aa_alignment[spec]:
                    if cur_aa != '-' and cur_aa != aa and cur_aa not in self.settings.EXCLUDE_AA_FROM:
                        suspected_aa.append(cur_aa)
                        try:
                            self.aa2aa_rea[aa][cur_aa].add(spec)
                        except KeyError:
                            self.aa2aa_rea[aa] = defaultdict(set)
                            self.aa2aa_rea[aa][cur_aa].add(spec)

    def save_json(self):
        with open(os.path.join(self.settings.OUTDIR, "aause.json"), "w") as outfile1:
            json.dump(self.aa_sim_json, outfile1, indent=4)
        with open(os.path.join(self.settings.OUTDIR, "similarity.json"), "w") as outfile2:
            json.dump(self.sim_json, outfile2, indent=4)
        with open(os.path.join(self.settings.OUTDIR, "reassignment.json"), "w") as outfile3:
            json.dump(self.reassignment_mapper, outfile3, indent=4)

    def save_all(self):
        self.save_json()
        id_filtfile = os.path.join(self.settings.OUTDIR, "id_filt.fasta")
        ic_filtfile = os.path.join(self.settings.OUTDIR, "ic_filt.fasta")
        gap_filtfile = os.path.join(self.settings.OUTDIR, "gap_filt.fasta")
        newick = os.path.join(self.settings.OUTDIR, "tree.nwk")
        self.seqset.write_data(
            id_filtered=id_filtfile, gap_filtered=gap_filtfile, ic_filtered=ic_filtfile, tree=newick)

    def run_analysis(self):
        """ Run the filtering analysis of the current dataset in sequenceset"""
        self.get_genomes()
        self.interesting_case = []
        self.possible_aa_reassignation()
        codon_align, fcodon_align = self.seqset.get_codon_alignment()
        self.reassignment_mapper['genes'] = self.seqset.get_genes_per_species()
        self.reassignment_mapper['genome']['global'] = self.seqset.get_genome_size()
        self.reassignment_mapper['genome']['filtered'] = self.seqset.get_genome_size("filtered")
        for aa1, aarea in self.aa2aa_rea.items():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa1]]
            gcodon_rea = CodonReaData((aa1, aarea), self.seqset.prot_align, self.global_consensus, codon_align,
                                      self.seqset.codontable, self.seqset.position, self.seqset.gene_limits, self.settings)
            fcodon_rea = CodonReaData((aa1, aarea), self.seqset.filt_prot_align, self.filtered_consensus, fcodon_align,
                                      self.seqset.codontable, self.seqset.filt_position, self.seqset.gene_limits, self.settings)

            for aa2, species in aarea.items():
                logging.debug("%s to %s" % (aa2, aa1))
                counts = []
                t = self.seqset.phylotree.copy("newick")
                fitch = NaiveFitch(t, species, aa_letters_1to3[aa2], aa_letters_1to3[
                                   aa1], self.settings,  self.seqset.codontable, (gcodon_rea, fcodon_rea))
                slist = fitch.get_species_list(
                    self.settings.LIMIT_TO_SUSPECTED_SPECIES)

                alldata = {}
                for genome in slist:
                    rec = self.seqset.prot_dict[genome]
                    leaf = fitch.tree & genome
                    filt_count = 0
                    try:
                        filt_count = len(
                            [y for y in aa_alignment[genome] if y != '-' and y == aa2])
                    except Exception:
                        # wtver happen, do nothing
                        pass
                    leaf.add_features(count=0)
                    leaf.add_features(filter_count=filt_count)
                    leaf.add_features(lost=True)
                    for pos in xrange(len(self.global_consensus)):
                        if self.global_consensus[pos] == aa1 and rec[pos] == aa2:
                            leaf.count += 1
                    
                    reacodon = gcodon_rea.get_reacodons(genome, aa2)
                    usedcodon = gcodon_rea.get_usedcodons(genome, aa2)
                    #if not self.settings.USE_GLOBAL:
                    #    reacodon = fcodon_rea.get_reacodons(genome, aa2)
                    #    usedcodon = fcodon_rea.get_usedcodons(genome, aa2)
                    
                    fisher_passed, pval = independance_test(
                        reacodon, usedcodon, confd=0.05)
                    # print "%s\t%s\t%s ==> %s" %(aa2, aa1, genome,
                    # str(column_comparision) )
                    # CHANGED :  we do not need count anymore 
                    # if('lost' in leaf.features and fisher_passed and leaf.count > self.settings.COUNT_THRESHOLD
                    if('lost' in leaf.features and fisher_passed):
                        leaf.lost = False

                    #fcodon_rea.check_gain
                    # this is computed many time
                    gdata = {}
                    greacodon = gcodon_rea.get_reacodons(genome, aa2)
                    gusedcodon = gcodon_rea.get_usedcodons(genome, aa2)
                    gmixtecodon = gcodon_rea.get_mixtecodons(genome, aa2)
                    g_rea_dist = gcodon_rea.get_rea_aa_codon_distribution(
                        genome, aa2)
                    g_total_rea_dist = gcodon_rea.get_total_rea_aa_codon_distribution(genome, aa2)
                    gdata['global'] = {'rea_codon': greacodon, 'used_codon': gusedcodon, 'mixte_codon': gmixtecodon,
                                     'count': leaf.count, 'rea_distribution': g_rea_dist, 'total_rea_distribution': g_total_rea_dist}
                    freacodon = fcodon_rea.get_reacodons(genome, aa2)
                    fusedcodon = fcodon_rea.get_usedcodons(genome, aa2)
                    fmixtecodon = fcodon_rea.get_mixtecodons(genome, aa2)
                    f_rea_dist = fcodon_rea.get_rea_aa_codon_distribution(genome, aa2)
                    f_total_rea_dist = fcodon_rea.get_total_rea_aa_codon_distribution(genome, aa2)
                    gdata['filtered'] = {'rea_codon': freacodon, 'used_codon': fusedcodon, 'mixte_codon': fmixtecodon, 
                                         'count': leaf.filter_count, 'rea_distribution': f_rea_dist, 'total_rea_distribution': f_total_rea_dist}
                    gdata['suspected'] = self.suspected_species[aa1].get(genome, 0)
                    gdata['score'] = {'global' : gcodon_rea.get_score(genome, aa2, aa1), 'filtered' : fcodon_rea.get_score(genome, aa2, aa1)}

                    gdata['lost'] = {'pval': pval, 'pass' : 1 if leaf.lost else 0}
                    gdata['fitch'] = fitch.get_distance_to_rea_node(leaf)
                    gdata['codons'] = {'global' : gcodon_rea.get_aa_usage(genome, aa2), 'filtered': fcodon_rea.get_aa_usage(genome, aa2)}
                    alldata[genome] = gdata

                if(fitch.is_valid()):
                    fitch.render_tree()
                    self.reassignment_mapper['aa'][aa2][aa1] = alldata
                    self.interesting_case.append("%s to %s" % (aa2, aa1))


        logging.debug("After validation, %d cases were found interesting" % len(
            self.interesting_case))
        for case in self.interesting_case:
            logging.debug(case)


def executeCMD(cmd, prog):
    """Execute a command line in the shell"""
    logging.debug("The following will be executed : \n%s\n" % cmd)
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    if err:
        Output.error(err)
    if out:
        logging.debug("%s : \n----------------- %s" % (prog, out))
    return err


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


class CoreFile:
    """Parse corefile format:
        >>geneID
        >specie1_ID
        ATCGTGATCCTAAGCTGATAGAAT
        >specie2_ID
        ATCGTGCAGATAGCTGATAGATGA
        ...
        >specie10_ID
        ATCTTAGTCGTAGCATAGTACTCG

        >>gene2_ID
        >specie1_ID
        ATGTGTCGTGCAGTCGT
        >specie1_ID
        ATGTGTCGTGCAGTCGT
        ...
    """

    def __init__(self, infile, dnafile, settings, gap_thresh, has_stop=True, use_tree=None, refine_alignment=True, msaprog="muscle"):

        # try parsing the protein file as a corefile and if it fail
        # parse it as a fasta file
        self.infile = infile
        self.sequences = {}
        stop_removed = False
        self.dnasequences = {}
        gap_thresh = (abs(gap_thresh) <= 1 or 0.01) * abs(gap_thresh)
        try:
            self.sequences = self.parse_corefile(infile, generic_protein)
            # we got to the end without finding any sequences
            if not self.sequences:
                raise ValueError(
                    'Corefile detected but protein sequence is in the wrong format')
        except:
            seqset = SeqIO.parse(open(infile, 'r'), 'fasta', generic_protein)
            seqset = list(seqset)
            is_aligned = False
            try:
                a = MultipleSeqAlignment(seqset)
                is_aligned = True
            except:
                pass
            self.sequences, stop_removed = self._split_at_stop(
                seqset, is_aligned)
            # if we got another error here, do not attempt to correct it

        logging.debug('Protein sequence was read')

        # here we are parsing the dnafile :
        # dnafile should be a corefile, if it's not, then we raise an exception
        try:
            self.dnasequences = self.parse_corefile(
                dnafile, generic_nucleotide)
            # we got to the end without finding any sequences
            if not self.dnasequences:
                raise ValueError(
                    'Nucleotide file empty or not in the corefile format')
        except:
            raise ValueError('Nucleotide is not in the corefile format')

        logging.debug('DNA sequence was read')
        self.genes = set(self.dnasequences.keys()).intersection(
            self.sequences.keys())

        common_spec_per_gene = {}
        common_spec_per_gene_len = {}
        max_spec = -np.inf
        self.genome_list = set([])
        for gene in self.genes:
            g_dnaspec = [x.id for x in self.dnasequences[gene]]
            g_protspec = [x.id for x in self.sequences[gene]]
            g_common_spec = set(g_dnaspec).intersection(set(g_protspec))
            common_spec_per_gene[gene] = g_common_spec
            self.genome_list.update(g_common_spec)
            cur_size = len(g_common_spec)
            max_spec = max_spec if max_spec > cur_size else cur_size
            common_spec_per_gene_len[gene] = cur_size

        # dropping gene with a low fraction of species representing them
        self.genes = [x for x in self.genes if common_spec_per_gene_len[
            x] > max_spec * gap_thresh]
        self.common_spec_per_gene = common_spec_per_gene
        if not self.genes:
            raise ValueError(
                "Can't find genes present in protein and nucleotide file.")
        logging.debug('List of selected genes : %s ' % str(self.genes))

        self.true_spec_list = set().union(*self.common_spec_per_gene.values())
        if (stop_removed and has_stop):
            self.clean_stop()
            logging.debug('Stop codon removed for sequences')

        is_aligned, alignment = self._prot_is_aligned()

        if not is_aligned:
            alignment = self.align(
                settings, refine=True, tree=use_tree, msaprog=msaprog, scale=settings.scale)
            logging.debug('Sequence is not aligned, align and refine it.')

        elif is_aligned and refine_alignment:
            logging.debug(
                'Sequence is already aligned but refine was requested')
            for g in self.genes:
                alignment[g] = self.__class__._refine(alignment[
                                                     g], 9999, settings.OUTDIR, settings.hmmbuild, settings.hmmalign, settings.eslalimanip, settings.eslalimask)

        self.alignment = alignment
        logging.debug('Sequence alignment done')

    def save_alignment(self, outfile):
        """save alignment in a file """
        AlignIO.write(self.alignment, open(outfile, 'w'), 'fasta')

    def concat(self, missing='-', alpha=generic_protein):
        """Concatenate alignment use as input, a list of all the gene, 
        a list of all the speclist, and a dict of alignment
        """
        align_dict = {}
        genepos = []
        lastpos = 0
        dna_dict = {}
        for gene in self.genes:
            cur_al = self.alignment[gene]
            al_len = cur_al.get_alignment_length()
            genepos.append((gene, lastpos, al_len + lastpos))
            lastpos += al_len
            spec_dict = dict((x.id, x) for x in cur_al)
            dna_spec_dict = dict((x.id, x) for x in self.dnasequences[gene])
            for spec in self.true_spec_list:
                # default , put missing data
                adding = SeqRecord(
                    Seq(missing * al_len, alpha), id=spec, name=spec)
                dna_adding =  SeqRecord(
                    Seq("", generic_nucleotide), id=spec, name=spec)
                if(spec in self.common_spec_per_gene[gene]):
                    adding = spec_dict[spec]
                    dna_adding = dna_spec_dict[spec]
                try:
                    align_dict[spec] += adding
                except (KeyError, ValueError):
                    align_dict[spec] = adding
                try:
                    dna_dict[spec] += dna_adding
                except (KeyError, ValueError):
                    dna_dict[spec] = dna_adding

        return align_dict, dna_dict, genepos

    def align(self, settings, msaprog, refine=True, tree=None, scale=1.0, alpha=generic_protein):
        alignment = {}
        outdir = settings.OUTDIR
        for gene in self.genes:
            # in order to have the best possible alignment, we are going to
            # keep the
            al = self.__class__._align(
                self.sequences[gene], msaprog, tree, scale, outdir, alpha)
            if refine:
                al = self.__class__._refine(
                    al, 9999, outdir, settings.hmmbuild, settings.hmmalign, settings.eslalimanip, settings.eslalimask)
            alignment[gene] = al
        return alignment

    @classmethod
    def _align(clc, msa, msaprog, tree, scale, outdir, alpha=generic_protein, is_aligned=False):
        """Realign a msa according to a tree"""
        tmpseq = os.path.join(outdir, "tmp_seq.fasta")
        align_seq = os.path.join(outdir, "tmp_msa.fasta")
        seqs = msa
        if is_aligned:
            seqs = []
            for seqrec in msa:
                seqs.append(seqrec.seq.ungap())

        SeqIO.write(seqs, open(tmpseq, 'w'), 'fasta')

        if tree and 'mafft' in msaprog:
            seq_order = [seqrec.id for seqrec in seqs]
            out = Output(file=os.path.join(outdir, "tmp_tree.mafft"))
            convert_tree_to_mafft(Tree(tree), seq_order, out, scale)
            out.close()
            msaprog += " --treein %s" % out.file

        execute_alignment(msaprog, tmpseq, align_seq)
        msa = AlignIO.read(align_seq, 'fasta', alphabet=alpha)
        filelist = glob.glob(os.path.join(outdir + "tmp_*"))
        for f in filelist:
            os.remove(f)
        return msa

    @classmethod
    def _refine(clc, alignment, timeout, outdir, hmmbuild="hmmbuild", hmmalign="hmmalign",
               eslalimanip="esl-alimanip", eslalimask="esl-alimask", loop=10, minqual=9, clean=True):
        """Align and refine at the same time with hmmalign and muscle"""

        success, not_found = check_binaries(
            hmmbuild, hmmalign, eslalimask, eslalimanip)
        if not success:
            raise RuntimeError(
                "Could not refine alignment, Executable not found: %s!!!" % (not_found))

        def accessQuality(alignfile, minqual):
            ppcons = ''
            with open(alignfile, 'r') as text:
                for line in text:
                    if '#=GC' in line and 'PP_cons' in line:
                        ppcons += line.split()[-1].strip()
            highpos = [1 for x in ppcons if x ==
                       '*' or (isInt(x) and int(x) >= minqual)]
            return sum(highpos)

        cur_dir = []
        for i in xrange(loop):
            d = os.path.join(outdir, '%d' % i)
            purge_directory(d)
            cur_dir.append(d)

        bestiter = {}
        file_created = False
        inputFile = os.path.join(outdir, 'tmp_align.fasta')
        if isinstance(alignment, MultipleSeqAlignment):
            file_created = True
            AlignIO.write(alignment, open(inputFile, 'w'), 'fasta')
        else:
            inputFile = alignment

        logging.debug(
            '... trying to run hmmbuild and hmmalign ' + str(loop) + " times!")
        quality = []
        outlist = []
        for i in xrange(loop):
            hmmfile = os.path.join(cur_dir[i], "alignment.hmm")
            outputFile = os.path.join(cur_dir[i], "alignment.sto")
            buildline = hmmbuild + " --amino %s %s" % (hmmfile, inputFile)
            executeCMD(buildline, 'hmmbuild')
            # will continue if not exception is found
            alignline = hmmalign + \
                " -o %s %s %s" % (outputFile, hmmfile, inputFile)
            executeCMD(alignline, 'hmmalign')
            # finding quality
            quality.append(accessQuality(outputFile,  minqual))
            inputFile = remove_gap_position(outputFile, 'stockholm')
            outlist.append(inputFile)
            # next input for hmm is current output

        # find the iteration with the greatest number of aligned position
        bestiter = outlist[np.asarray(quality).argmax()]
        alignment = AlignIO.read(bestiter, format="fasta")
        # clean by removing anything we had
        if clean:
            for i in xrange(loop):
                d = os.path.join(outdir, '%d' % i)
                shutil.rmtree(d, ignore_errors=True)
            if file_created:
                os.remove(os.path.join(outdir, 'tmp_align.fasta'))
        return alignment

    def clean_stop(self, stop='*'):
        """Clean Stop from the sequence"""
        # this remove all the stop codon at the ends
        for gene in self.genes:
            seq_len_per_spec = {}
            for seqrec in self.sequences[gene]:
                common_sequence_per_genes[gene].add(seqrec)
                if seqrec.seq.endswith(stop):
                    seqrec.seq = seqrec.seq[0:-1]

                seq_len_per_spec[gene] = len(seqrec)

            for dseqrec in self.dnasequences[gene]:
                if dseqrec.id in self.common_spec_per_gene[gene]:
                    if len(dseqrec) == (seq_len_per_spec[gene] + 1) * 3:
                        dseqrec.seq = dseqrec.seq[0:-3]

    def _prot_is_aligned(self):
        alignment = {}
        try:
            for gene in self.genes:
                alignment[gene] = MultipleSeqAlignment(self.sequences[gene])
        except ValueError:
            return False, None

        return True, alignment

    def _split_at_stop(self, seqlist, is_aligned):
        """Split the sequence where we have a stop """
        stop_checked, has_stop, stopmapping = check_stop(seqlist, is_aligned)
        stop_removed = True
        if not has_stop:
            # in this case, there isn't any stop codon (*)
            return {'0': seqlist}, stop_removed

        # stop codon are strip from this representation
        elif stop_checked:
            stop_removed = False
            core_dict = {}
            gene_pos = stopmapping.values()[0]
            start = 0
            for i in xrange(len(gene_pos)):
                seq_rec_list = []
                for seqrec in seqlist:
                    s_rec = SeqRecord(seqrec.seq[start:stopmapping[seqrec.id][i]],
                                      id=seqrec.id, name=seqrec.name,
                                      description=seqrec.description,
                                      features=seqrec.features[:],
                                      annotations=seqrec.annotations.copy())
                    seq_rec_list.append(s_rec)
                start = stopmapping[seqrec.id][i] + 1
                core_dict[str(i)] = seq_rec_list

            return core_dict, stop_removed
        else:
            raise ValueError(
                "Stop codon no uniformely placed in your alignment. Expect them to be places at the same position in all sequences")

    @classmethod
    def write_corefile(clc, sequences, infile):
        with open(infile, 'w') as OUT:
            for gene, sequence in sequences.items():
                OUT.write('>>%s\n' % gene)
                for seq in sequence:
                    OUT.write('>%s\n' % seq.name)
                    OUT.write('%s\n' % seq.seq._data)

    @classmethod
    def parse_corefile(clc, infile, alphabet=generic_protein):
        core_dict = {}
        with open(infile, 'r') as handle:
            for gene, sequences in clc._internal_coreparser(handle, alphabet):
                core_dict[gene] = sequences
        return core_dict

    @classmethod
    def translate(clc, core_inst, gcode=1):
        codontable = CodonTable.unambiguous_dna_by_id[1]
        try :
            codontable = CodonTable.unambiguous_dna_by_id[abs(gcode)]
        except:
            logging.warn("Wrong genetic code, resetting it to 1")

        if not isinstance(core_inst, dict):
            core_inst =  clc.parse_corefile(core_inst, alphabet=generic_nucleotide)

        core_prot_inst = {}
        for gene, seqs in core_inst.items():
            translated = []
            for s in seqs:
                if len(s)%3 != 0:
                    raise ValueError("Frame-shifting detected in %s : [%s], current version does not supported it."%(s.id, gene))
                else:
                    trans_s = Seq("".join([codontable.forward_table.get(s[i:i+3], 'X') for i in xrange(0, len(s), 3)]), generic_protein)
                    translated.append(SeqRecord(trans_s, id=s.id, name=s.name))
            core_prot_inst[gene] = translated

        return core_prot_inst

    @classmethod
    def _internal_coreparser(clc, handle, alphabet):
        while True:
            line = handle.readline()
            if line == "":
                return  # file end abruptly
            if line.startswith('>>'):
                # start looking only if we have the first line starting with
                # '>>'
                break

        genecounter = 0
        while True:
            if not line.startswith('>>'):
                raise ValueError(
                    'Records in corefile format should start with an >> followed by the name of the gene')
            genename = line[2:].rstrip()
            if not genename:
                genename = str(genecounter)
                genecounter += 1

            sequences = []
            line = handle.readline()
            while True:
                if not line:
                    break
                if line.startswith('>>'):
                    break
                if line.startswith('>'):
                    title = line[1:].rstrip()
                    lines = []
                    line = handle.readline()
                    while True:
                        if not line:
                            break
                        if line.startswith('>'):
                            break
                        lines.append(line.rstrip())
                        line = handle.readline()

                    sequence = "".join(lines).replace(
                        " ", "").replace("\r", "")
                    try:
                        first_word = title.split(None, 1)[0]
                    except IndexError:
                        assert not title, repr(title)
                        first_word = ""

                    sequences.append(SeqRecord(Seq(sequence, alphabet),
                                               id=first_word, name=first_word, description=title))
                else:
                    line = handle.readline()

            if sequences:
                yield genename, sequences
            if not line:
                return

        assert False, "We should return before this point"


def check_stop(seq_list, is_aligned, stop='*'):
    result_map = {}
    for seq in seq_list:
        last_pos = len(seq.seq) - 1
        result_map[seq.id] = [pos for pos,
                              char in enumerate(seq.seq) if char == stop]
    result = result_map.values()
    length_list = sorted([len(x) for x in result])
    length = length_list[-1]
    if is_aligned:
        result = np.array([r + [None] * (length - len(r)) for r in result])
        first_stops = result[0]
        position_match = np.all(result == first_stops, axis=0)
        if not position_match:
            return False, False, result_map
        elif np.all(position_match):
            return True, True, result_map
        return False, True, result_map
    else:
        return np.all(np.asarray(length_list) == length), length != 0, result_map


def purge_directory(dirname):
    """Empty a repertory"""
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname)


def check_binaries(*args):
    not_found = []
    for exe in args:
        if not spawn.find_executable(exe):
            not_found.append(exe)

    return not_found == [], not_found


def timeit(func):

    def timed(*args, **kw):
        tstart = time.time()
        result = func(*args, **kw)
        tend = time.time()
        ttime = tend - tstart
        # print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
        return ttime, result

    return timed


def isInt(value):
    """ Verify if a string
    can be casted to an int"""
    try:
        int(value)
        return True
    except ValueError:
        return False


def remove_gap_position(alignfile, curformat, alimanip):
    align = AlignIO.read(alignfile, curformat)
    align, positions = SequenceSet.clean_alignment(align, threshold=1)
    fastafile = alignfile.split('.')[0] + ".fasta"
    AlignIO.write(align, open(fastafile, 'w'), 'fasta')
    return fastafile


def independance_test(rea, ori, confd=0.05):

    codon_list = set(rea.keys())
    codon_list.update(ori.keys())
    codon_list = [x for x in codon_list if not (rea.get(x,0) == 0 and ori.get(x,0) == 0)]
    nelmt = len(codon_list)
    # line 1 is rea, line 2 is observed
    if nelmt > 1:
        obs = np.zeros((nelmt, 2))
        for i in xrange(nelmt):
            obs[i, 0] = rea.get(codon_list[i], 0)
            obs[i, 1] = ori.get(codon_list[i], 0)

        pval = fisher_exact(obs, midP=True)
        return pval <= confd, pval
    # strangely, codon is used only in rea column
    elif len(rea.values()) != 0 and len(ori.values()) == 0:
        return True, 0
    # codon is used in both column
    elif len(rea.values()) != 0 and len(ori.values()) != 0:
        return rea.values()[0] >= ori.values()[0], 0
    else :
        return False, 1


def chi2_is_possible(obs):
    exp = ss.contingency.expected_freq(obs)
    limit = 0.8
    count = 0
    size = len(exp.ravel())
    for i in sorted(exp.ravel()):
        if i == 0:
            return False, False
        if i < 5:
            count += 1
        if count * 1.0 / size > 1 - limit:
            return False, True
    return True, count != 0


@timeit
def execute_alignment(cmdline, inp, out):
    """Execute the aligner command line"""
    prog = 'mafft'
    if 'muscle' in cmdline:
        prog = 'muscle'
        cmdline += " -in %s -out %s" % (inp, out)
    elif 'mafft' in cmdline:
        cmdline += " %s > %s" % (inp, out)

    executeCMD(cmdline, prog)


def check_align_upgrade(al1, al2, settings, method="wilcoxon", column_list=None):
    # al1 is the new alignement after changing the Gcode, and it's a dict
    # what we are testing is mean(al1_qual) >  mean(al2_qual)
    scoring_method = getattr(settings, 'method', "identity")
    if scoring_method != "identity" and isinstance(scoring_method, basestring):
        scoring_method = getattr(MatrixInfo, scoring_method)
    
    if isinstance(al1, dict):
        al1 =  al1.values()
    if isinstance(al2, dict):
        al2 =  al2.values()
    
    al_len = len(al1[0])
    nspec = len(al1)
    if not column_list:
        column_list = range(al_len)
    
    assert (nspec == len(al2) and al_len == len(al2[0])), "The alignement are different"
    al1_sim, al2_sim = compute_SP_per_col(al1, al2, column_list, nspec, scoring_method)
    rank, pval = ReaGenomeFinder.get_paired_test(al1_sim, al2_sim, 1, method)
    return pval

def compute_SP_per_col(al1, al2, columns, nspec,  scoring_matrix):
    def scoring_function(aa1, aa2, scoring_matrix):
        if scoring_matrix == 'identity':
            return (aa1 == aa2)*1
        else: 
            return scoring_matrix[aa1, aa2]

    al1_score = []
    al2_score = []
    for col in columns:
        al1_s = 0
        al2_s = 0
        for i in xrange(nspec-1):
            for j in xrange(i+1, nspec):
                al1_s += scoring_function(al1[i][col], al1[j][col], scoring_matrix)
                al2_s += scoring_function(al2[i][col], al2[j][col], scoring_matrix)
        al1_score.append(al1_s)
        al2_score.append(al2_s)
    return al1_score, al2_score


def check_gain(codon, cible_aa, speclist, table, codon_alignment, settings, alignment=None,  method="wilcoxon"):

    if not isinstance(codon_alignment, dict):
        codon_alignment = SeqIO.to_dict(codon_alignment)
    
    codontable = table.forward_table       
   
    def translate(codon_alignment, codontable, changes=None):
        if changes is None:
            changes = {}
        #return dict((spec, "".join(map(lambda x: codontable[x] if changes.get(spec, (None,))[0] != x else changes[spec][1],
        #                [prot[i:i + 3] for i in xrange(0, len(prot), 3)]))) for spec, prot in codon_alignment.items())
        translated_al = {}
        position = []
        for spec, prot in codon_alignment.items():
            translated = ""
            for i in xrange(0, len(prot), 3):
                codon = prot[i:i+3]
                if changes.get(spec, (None,))[0] != codon :
                    translated += codontable[codon]
                else :
                    # in this case, there were reassignment
                    translated += changes[spec][1]
                    position.append(int(i/3))
            translated_al[spec] = translated
        return translated_al, position

    if not alignment:
        alignment, _ = translate(codon_alignment, codontable)

    changes = {}
    for spec in speclist:
        changes[spec] = (codon, cible_aa)
    
    cor_alignment, position = translate(codon_alignment, codontable, changes)
    score_improve = check_align_upgrade(cor_alignment, alignment, settings, method, position)
    return score_improve, codontable