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
import time

from collections import Counter, defaultdict

from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import codonalign
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.codonalign.codonalphabet import default_codon_alphabet, get_codon_alphabet
from Bio.Data import CodonTable
from Bio.codonalign.codonseq import  _get_codon_list, CodonSeq
from Bio import SubsMat
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, generic_protein
from cStringIO import StringIO
from ete3 import PhyloTree, Tree
from PPieChartFace import PPieChartFace, LineFace
import settings as default_settings
from functools import partial
import scipy.stats as ss
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

# hardcorded settings variables
alpha = Alphabet.Gapped(IUPAC.protein)
DIST_THRES = 1e-10

#if settings.GENETIC_CODE < 0:
    #for codon in ['CTG', 'CTA', 'CTC', 'CTT']:
        #dct.forward_table[codon] = 'L'
    #dct.back_table = CodonTable.make_back_table(dct.forward_table, dct.stop_codons[0])

nuc_letters = "ACTG"

class Settings():
    def __init__(self, mode=0):
        # mode will set the type of construction we want
        self.mode = mode
   
    def set(self, **kwargs):
        self.PROCESS_ENABLED = kwargs.get('PROCESS_ENABLED', default_settings.PROCESS_ENABLED)
        EXCLUDE_AA = kwargs.get('EXCLUDE_AA', default_settings.EXCLUDE_AA)
        self.AA_LETTERS = "".join([aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa not in EXCLUDE_AA])
        self.OUTDIR = kwargs.get('OUTDIR', default_settings.OUTDIR)
        self.SKIP_ALIGNMENT = kwargs.get('SKIP_ALIGNMENT', default_settings.SKIP_ALIGNMENT)
        self.AA_MAJORITY_THRESH = kwargs.get('AA_MAJORITY_THRESH', default_settings.AA_MAJORITY_THRESH)
        self.LIMIT_TO_SUSPECTED_SPECIES = kwargs.get('LIMIT_TO_SUSPECTED_SPECIES', default_settings.LIMIT_TO_SUSPECTED_SPECIES)
        self.FREQUENCY_THRESHOLD = kwargs.get('FREQUENCY_THRESHOLD', default_settings.FREQUENCY_THRESHOLD)
        self.COUNT_THRESHOLD = kwargs.get('COUNT_THRESHOLD', default_settings.COUNT_THRESHOLD)
        self.GENETIC_CODE = kwargs.get('GENETIC_CODE', default_settings.GENETIC_CODE)
        self.SHOW_FILTERED_CODON_DATA = kwargs.get('SHOW_FILTERED_CODON_DATA', default_settings.SHOW_FILTERED_CODON_DATA)
        self.IMAGE_FORMAT = kwargs.get('IMAGE_FORMAT', default_settings.IMAGE_FORMAT)
        self.ADD_LABEL_TO_LEAF = kwargs.get('ADD_LABEL_TO_LEAF', default_settings.ADD_LABEL_TO_LEAF)
     
    def fill(self, setting):
        self.__dict__.update(setting.__dict__)


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
    
    def __init__(self, aas, consensus, codon_align_dict, dct):
        
        self.aa1, self.aa2 = aas
        self.codon_alignment = codon_align_dict # 
        self.consensus = consensus
        self.reacodons = defaultdict(Counter)
        # find codon usage in the sequence
        self.usedcodons = defaultdict(Counter)
        self.mixtecodon = defaultdict(Counter)
        #print len(codon_align_dict.values()[0])
        self.specs_amino_count = defaultdict(partial(defaultdict, Counter))
        self._spec_codon_usage()

    def _spec_codon_usage(self):
        for i, aa in enumerate(self.consensus):
            for spec in self.codon_alignment.keys():
                spec_codon = self.codon_alignment[spec].seq.get_codon(i)
                spec_aa = dct.forward_table[spec_codon] if '-' not in spec_codon else None
                if(spec_aa):
                    # determine aa use in each specie
                    self.specs_amino_count[spec][spec_aa].update([spec_codon])

                    # potential reassignment counter
                    # species use aa2 while aa1 is prevalent
                    if spec_aa==self.aa2 and aa ==self.aa1 :
                        self.reacodons[spec].update([spec_codon])
                    # species use aa2 while aa2 is prevalent
                    elif spec_aa==self.aa2 and aa==self.aa2 :
                        self.usedcodons[spec].update([spec_codon])
                    # other position where aa2 is used in species
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

    def get_aa_usage(self, specie):
        return self.specs_amino_count[specie]

    def get_string(self, specie, show_mixte=False):
        """ Return a string representation for the specie"""

        seq = "%s\t%s"%(dict(self.get_reacodons(specie)), dict(self.get_usedcodons(specie)))
        if show_mixte:
            seq += "\t%s" %(dict(self.get_mixtecodons(specie)))

        return seq


class NaiveFitch(object):
    """A NaiveFitch algorithm for finding the most parcimonious solution"""
    def __init__(self, tree, reassigned, ori_aa, dest_aa, settings, dct, codon_rea=(None, None)):
        self.id = {}
        self.tree = tree
        self.settings =  settings
        self.corr = {'0':ori_aa, '1':dest_aa}
        self.codon_rea_global, self.codon_rea_filtered = codon_rea
        colors=['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f']
        self.init_back_table()
        codon_list = self.back_table[aa_letters_3to1[ori_aa]]
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
    
    def init_back_table(self):
        self.back_table = defaultdict(list)
        for aa, codon in zip(dct.forward_table.values(), dct.forward_table.keys()):
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

    def get_species_list(self, limit_to_reassigned=False):
        """Get the species list for which we are almost certain
        that there is a reassignment
        """
        slist = set()
        if limit_to_reassigned :
            for node in self.tree.traverse():
                if 'reassigned' in node.features and node.reassigned=={1}:
                    slist.update(node.get_leaf_names())
        else :
            slist = set(self.tree.get_tree_root().get_leaf_names())
        return slist

    def get_distance_to_rea_node(self, node):
        if isinstance(node, str):
            node = self.tree&node
        dist = 0
        while not self.is_reassigned(node):
            dist += 1
            node = node.up
        return 0
    
    def has_codon_data(self):
        return (self.codon_rea_global, self.codon_rea_filtered) != (None, None)


    def render_tree(self, output="", suffix="", pie_size=50):

        GRAPHICAL_ACCESS = True
        try:
            from ete3 import TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
        except ImportError, e:
            GRAPHICAL_ACCESS = False

        if not output:
            output = TMP+self.ori_aa+"_to_"+self.dest_aa+suffix+"."+settings.IMAGE_FORMAT

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
            
            add_label = self.settings.ADD_LABEL_TO_LEAF
            def layout(node):

                get_suffix = lambda x : x if add_label else ""

                N = AttrFace("rea", fsize=10)
                faces.add_face_to_node(N, node, 0, position="branch-right")
                has_count, has_fcount = 0, 0
                if('count' in node.features):
                    faces.add_face_to_node(AttrFace("count", fsize=6, fgcolor="firebrick"), node, column=1, position="branch-bottom")
                    has_count = node.count

                if('filter_count' in node.features):
                    faces.add_face_to_node(AttrFace("filter_count", fsize=6, fgcolor="indigo"), node, column=1, position="branch-top")
                    has_fcount = node.filter_count

                if not (has_fcount or has_count):
                    faces.add_face_to_node(AttrFace("name"), node, 0, position="aligned")

                else:
                    # if lost in node.features, then node is already a leaf
                    if 'lost' in node.features and node.lost:
                        if(self.is_reassigned(node)):
                            faces.add_face_to_node(AttrFace("name", fgcolor="seagreen",text_suffix=get_suffix("_V")), node, 0, position="aligned")
                        else :
                            faces.add_face_to_node(AttrFace("name", fgcolor="#cccccc",text_suffix=get_suffix("_G")), node, 0, position="aligned")

                    elif 'lost' in node.features and not node.lost:
                        if(self.is_reassigned(node)):
                            faces.add_face_to_node(AttrFace("name", fgcolor="red", fstyle ="italic",text_suffix=get_suffix("_R")), node, 0, position="aligned")
                        else : 
                            faces.add_face_to_node(AttrFace("name", fgcolor="gold", text_suffix=get_suffix("_O")), node, 0, position="aligned")


                if(self.has_codon_data() and node.is_leaf()):
                    node_colors = aa_letters_3to1[self.ori_aa]
                    spec_codonrea_g = self.codon_rea_global.get_reacodons(node.name)
                    spec_codonused_g = self.codon_rea_global.get_usedcodons(node.name)
                    
                    faces.add_face_to_node(PPieChartFace(spec_codonrea_g.values(), pie_size, pie_size, labels=[], \
                        colors=[self.colors[k] for k in spec_codonrea_g.keys()]), node, column=1, position="aligned")
                    
                    faces.add_face_to_node(PPieChartFace(spec_codonused_g.values(), pie_size, pie_size, labels=[],\
                     colors=[self.colors[k] for k in spec_codonused_g.keys()]), node, column=2, position="aligned")
                    next_column = 3
                    if(settings.SHOW_MIXTE_CODONS):
                        spec_mixtecodon_g = self.codon_rea_global.get_mixtecodons(node.name)
                        faces.add_face_to_node(PPieChartFace(spec_mixtecodon_g.values(), pie_size, pie_size, labels=[], \
                            colors=[self.colors[k] for k in spec_mixtecodon_g.keys()]), node, column=next_column, position="aligned")
                        next_column = 4

                    if(settings.SHOW_FILTERED_CODON_DATA): 

                        spec_codonrea_f = self.codon_rea_filtered.get_reacodons(node.name)
                        spec_codonused_f = self.codon_rea_filtered.get_usedcodons(node.name)
                    
                        # add separator 
                        faces.add_face_to_node(LineFace(pie_size, pie_size, None), node, column=next_column, position="aligned")
                        
                        # add data
                        faces.add_face_to_node(PPieChartFace(spec_codonrea_f.values(), pie_size, pie_size, labels=[],\
                            colors=[self.colors[k] for k in spec_codonrea_f.keys()]), node, column=next_column+1, position="aligned")

                        faces.add_face_to_node(PPieChartFace(spec_codonused_f.values(),  pie_size, pie_size, labels=[],\
                            colors=[self.colors[k] for k in spec_codonused_f.keys()]), node, column=next_column+2, position="aligned")
        
                        if(SHOW_MIXTE_CODONS):
                            spec_mixtecodon_f = self.codon_rea_filtered.get_mixtecodons(node.name)
                            faces.add_face_to_node(PPieChartFace(spec_mixtecodon_f.values(), pie_size, pie_size, labels=[], \
                                colors=[self.colors[k] for k in spec_mixtecodon_f.keys()]), node, column=next_column+3, position="aligned")
            

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


class SequenceSet(object):

    def __init__(self, dnadict, prot_align, phylotree, table_num, stopcodon=False):
        # using sequence set suppose that the alignment and the protein concatenation was already done
        self.codontable = CodonTable.unambiguous_dna_by_id[abs(table_num)]
        self.prot_align = prot_align
        self.phylotree = phylotree
        self.common_genome = None
        self.stop_positions = []
        self.dna_dict = dnadict
        if stopcodon:
            self.split_and_concat()
        self.prot_dict = SeqIO.to_dict(prot_align)
        self.restrict_to_common()
        self.codon_align()

    def split_and_concat(self):
        # we expect all stop to be at the same position in all the aligned sequence
        # this is not checked, user should check that instead
        for i, c in enumerate(self.prot_align[0]):
            if c == '*':
                self.stop_positions.append(i)

        # we don't need to remove stop codons, since this will be done later in the filtering

    def restrict_to_common(self):
        """ Keep only the species present in the dna, prot seq and in the tree """

        dna_set = set(self.dna_dict.keys())
        prot_set = set(self.prot_dict.keys())
        species_set = set(self.phylotree.get_leaf_names())
        common_genome = dna_set.intersection(prot_set)
        # remove every dna sequence that are not cds
        common_genome = set([x for x in common_genome if (len(self.dna_dict[x]) in [len(self.prot_dict[x].seq.ungap('-'))*3, len(self.prot_dict[x].seq.ungap('-'))*3 +3])])
        # print [(x, len(dna_dict[x]), len(prot_dict[x].seq.ungap('-'))+1, len(dna_dict[x])==(len(prot_dict[x].seq.ungap('-'))+1)*3) for x in common_genome ]
        # common_genome = [x for x in common_genome if len(dna_dict[x])==(len(prot_dict[x].seq.ungap('-'))+1)*3 ]
        # return common_genome, dict((k, v) for k, v in dna_dict.items() if k in common_genome)
        speclen =  len(common_genome)
        common_genome.intersection_update(species_set)

        if not common_genome:
            raise ValueError('ID intersection for dna, prot and species tree is empty')
        if len(prot_set) != speclen or len(dna_set) != speclen or len(species_set) != common_genome:
            logging.debug('Non-uniform sequence in dna sequences, prot sequences and trees')
            logging.debug('Recheck the following id %s' % (set.union(species_set , prot_set , dna_set) - common_genome))

        # prune tree to sequence list
        self.phylotree.prune(common_genome)
        self.common_genome = common_genome
        logging.debug("List of common genome")
        logging.debug(common_genome)
        #return common_genome, dict((k, v) for k, v in dna_dict.items() if k in common_genome), prot_dict
        self.dna_dict, self.prot_dict = dict((k, v) for k, v in self.dna_dict.items() if k in common_genome), dict((k, v) for k, v in self.prot_dict.items() if k in common_genome)
        
        for k,v in self.prot_dict.items():
            v.seq.alphabet = generic_protein
            v.id = k

        self.prot_align = MultipleSeqAlignment(self.prot_dict.values(), alphabet=alpha)


    @staticmethod
    def copy_codon_record(record, codontable, gap_char='-'):
        """Return a Codon seq sequence from a dna sequence"""
        alphabet = get_codon_alphabet(codontable, gap_char=gap_char)
        return SeqRecord(CodonSeq(record.seq._data, alphabet=alphabet), id=record.id)

    
    def copy_codon_alignment(self, codon_alignment, alphabet=default_codon_alphabet):
        """Return a codon alignment object from a list of codon seq sequences"""
        return codonalign.CodonAlignment((copy_codon_record(rec, self.codontable) for rec in codon_alignment._records), alphabet=alphabet)


    def codon_align(self, alphabet=default_codon_alphabet, gap_char='-'):
        """ Perform a codon alignment based on a protein multiple alignment
        and remove gaped positions
        """
        alphabet = get_codon_alphabet(self.codontable, gap_char=gap_char)
        if self.common_genome is None:
            self.restrict_to_common()
        # build codon alignment and return it
        codon_aln = []
        for g in self.dna_dict.keys():
            codon_rec = self._get_codon_record(self.dna_dict[g], self.prot_dict[g], self.codontable, alphabet)
            codon_aln.append(codon_rec)

        self.codon_alignment = codonalign.CodonAlignment(codon_aln, alphabet=alphabet)
        #self.codon_alignment = codonalign.build(self.prot_align, self.dna_dict, codon_table=self.codontable)
    

    def _get_codon_record(self, dnarec, protrec, codontable, alphabet, gap_char='-'):
        nuc_seq =  dnarec.seq.ungap(gap_char)
        codon_seq = ""
        aa_num = 0
        max_error = 10
        for aa in protrec.seq:
            if aa == gap_char:
                codon_seq += gap_char*3
            else : 
                next_codon =  nuc_seq._data[aa_num*3 : (aa_num+1)*3]
                if not str(Seq(next_codon.upper()).translate(table=codontable)) == aa.upper():
                    logging.warn("%s(%s %d) does not correspond to %s(%s)"
                                  % (protrec.id, aa, aa_num, dnarec.id, next_codon))
                    max_error -= 1
                if max_error <= 0:
                    raise ValueError("You're obviously using the wrong genetic code")
                codon_seq += next_codon
                aa_num += 1
        return SeqRecord(CodonSeq(codon_seq, alphabet), id=dnarec.id)


    def filter_codon_alignment(self, ind_array=None, get_dict=False, alphabet=default_codon_alphabet):

        """Return the codon aligment from a list of in array"""
        if not ind_array:
            ind_array = (self.filt_position, )

        elif isinstance(ind_array, tuple):
            raise ValueError('Expect tuple for ind_array, got %s'%(type(ind_array)))
        else : 
            codon_align =  self.copy_codon_alignment(codon_alignment)

            for indexes in ind_array:
                indexes = sorted(indexes)
                filt_codon_align = codonalign.CodonAlignment([])
                for i in xrange(len(codon_align)) :
                    codseq = CodonSeq('') 
                    for pos in indexes:
                        codseq += codon_align[i].seq.get_codon(pos)
                    filt_codon_align[i].seq = codseq
                if get_dict:
                    filt_codon_align = SeqIO.to_dict(filt_codon_align)
                yield filt_codon_align


    def get_codon_alignment(self):
        self.fcodon_alignment = next(self.filter_codon_alignment())
        return self.codon_alignment, self.fcodon_alignment


    def prot_filtering(self, id_thresh=None, gap_thresh=None, ic_thresh=None, rmcnst=True):
        
        current_alignment = self.prot_align
        tt_filter_position = np.asarray(xrange(current_alignment.get_alignment_length()))

        if(gap_thresh): 
            self._gap_alignment, self._gap_filtered_position = self.clean_alignment(current_alignment, threshold=(abs(gap_thresh) <= 1 or 0.01) * abs(gap_thresh))
            current_alignment = self._gap_alignment
            logging.debug("Alignment length after removing gaps : %d"%current_alignment.get_alignment_length())
            tt_filter_position = tt_filter_position[self._gap_filtered_position]      
        
        if(ic_thresh):
            align_info = AlignInfo.SummaryInfo(current_alignment)
            ic_content = align_info.information_content()
            max_val = max(align_info.ic_vector.values()) * ((abs(ic_thresh) <= 1 or 0.01) * abs(ic_thresh))
            ic_pos = (np.asarray(align_info.ic_vector.values())>max_val).nonzero()
            logging.debug("Filtering with ic_content, vector to discard is %s"%str(ic_pos[0]))
            # ic_pos here is a tuple, so we should access the first element
            self._ic_alignment = self.filter_align_position(current_alignment, ic_pos[0])
            self._ic_filtered_positions = ic_pos
            current_alignment = self._ic_alignment
            tt_filter_position = tt_filter_position[self._ic_filtered_positions]    
            logging.debug("Alignment length filtering with information content: %d"%current_alignment.get_alignment_length())

        if(id_thresh):
            self._id_alignment, self._id_filtered_position = self.filter_alignment(current_alignment, remove_identity=rmcnst, threshold=(abs(id_thresh) <= 1 or 0.01) * abs(id_thresh))
            current_alignment = self._id_alignment
            tt_filter_position = tt_filter_position[self._id_filtered_position]
            logging.debug("Alignment length after filtering by sequence identity : %d"%current_alignment.get_alignment_length())
     
        self.filt_prot_align = current_alignment
        self.filt_position = tt_filter_position


    def write_data(self, id_filtered=None, gap_filtered=None, ic_filtered=None, tree=None):
        """Save file depending on the file use as input"""
        if id_filtered: 
            AlignIO.write(self._id_alignment, open(id_filtered, 'w'), 'fasta')
        if gap_filtered: 
            AlignIO.write(self._gap_alignment, open(gap_filtered, 'w'), 'fasta')
        if ic_filtered: 
            AlignIO.write(self._ic_alignment, open(ic_filtered, 'w'), 'fasta')
        if tree:
            self.phylotree.write(outfile=tree)


    @classmethod
    def filter_align_position(clc, alignment, index_array, alphabet=generic_protein):
        """Keep only columns specified by index_array from the alignment"""
        edited_alignment =  alignment[:,:]
        index_array = sorted(index_array)
        # alignment[:, slice(index_array[0], index_array[0] + 1)]
        for seqrec in edited_alignment :
            seq = Seq('', alphabet)
            for i in index_array:
                seq += seqrec[i]
            seqrec.seq = seq

        return edited_alignment
        

    @classmethod
    def clean_alignment(clc, alignment=None, characs=['-'], threshold=0.5):
        """Remove position of alignment which contain character from characs"""
        align_array = np.array([list(rec) for rec in alignment], np.character)
        indel_array = np.where((np.mean(np.in1d(align_array,
                                                characs).reshape(align_array.shape), axis=0) > threshold) == False)[0].tolist()

        return clc.filter_align_position(alignment, indel_array), indel_array


    @classmethod
    def filter_alignment(clc, alignment, threshold=0.8, remove_identity=False, ambiguous='X', alphabet=generic_protein):
        """Filter an alignment using threshold as the minimum aa identity per columns"""
        aligninfo = AlignInfo.SummaryInfo(alignment)
        # Smart move : Make a consensus sequence with threshold
        # remove all gap position
        consensus = aligninfo.gap_consensus(threshold=threshold, ambiguous=ambiguous, consensus_alpha=alphabet)
        cons_array = [i for i, c in enumerate(consensus) if c != ambiguous]

        if(remove_identity):
            matched_consensus = aligninfo.gap_consensus(threshold=1, ambiguous=ambiguous, consensus_alpha=alphabet)
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
    def get_consensus(clc, alignment, threshold):
        """return consensus, using a defined threshold"""
        aligninfo = AlignInfo.SummaryInfo(alignment)
        consensus = aligninfo.gap_consensus(threshold=threshold, ambiguous='X')
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


class ReaGenomeFinder:

    def __init__(self, seqset, settings, mode="count", method="identity"):
        self.seqset = seqset
        self.mode =  getattr(settings, 'mode', None)
        self.method =  getattr(settings, 'method', None)
        if self.mode is None:
            self.mode = mode
        if self.method is None:
            self.method = method
        self.suspected_species = defaultdict(dict)
        aa2aa_rea = defaultdict(dict)
        
        def makehash():
            return defaultdict(makehash)
        self.reassignment_mapper = makehash()
        

    def compute_sequence_identity(self, matCalc=None):
        if not matCalc:
            matCalc = DistanceCalculator(self.method)
        self.global_paired_distance = matCalc.get_distance(self.seqset.prot_align)
        self.filtered_paired_distance = matCalc.get_distance(self.seqset.filt_prot_align)


    def get_genomes(self, use_similarity=1):
        matCalc = DistanceCalculator(self.method)
        self.compute_sequence_identity(matCalc)

        self.filtered_consensus = self.seqset.get_consensus(self.seqset.filt_prot_align, self.settings.AA_MAJORITY_THRESH)  
        self.global_consensus = self.seqset.get_consensus(self.seqset.prot_align, self.settings.AA_MAJORITY_THRESH)  
       
        logging.debug("Filtered alignment consensus : \n%s\n"%filtered_consensus)    
        number_seq = self.seqset.get_total_genomes()
        seq_names =  self.seqset.common_genome
        self.sim_json = defaultdict(list)

        self.seqset.aa_filt_prot_align = {}
        self.aa_paired_distance = {}
        for aa in self.settings.AA_LETTERS:
            cons_array = get_aa_filtered_alignment(consensus, aa)
            if(cons_array):
                self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]] = SequenceSet.filter_align_position(self.filt_prot_align, cons_array)
                self.aa_paired_distance[aa] =  matCalc.get_distance(self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]])
                self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]] = SeqIO.to_dict(self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]])

        self.aa_sim_json = defaultdict(list)
        aa2suspect = defaultdict(Counter)
        aa2suspect_dist = defaultdict(partial(defaultdict, list))
        for (j, i) in itertools.combinations(xrange(number_seq), r=2):
            gpaired = abs(use_similarity - self.global_paired_distance[seq_names[i], seq_names[j]])
            fpaired = abs(use_similarity - self.filtered_paired_distance[seq_names[i], seq_names[j]])
            self.sim_json[seq_names[i]].append({"global": gpaired, "filtered": fpaired, "species": seq_names[j]})
            for aa in self.aa_paired_distance[aa].keys():
                aapaired = abs(use_similarity - self.aa_paired_distance[aa][seq_names[i], seq_names[j]])
                self.aa_sim_json[aa_letters_1to3[aa]].append({'global':gpaired , 'filtered': aapaired, "species": "%s||%s" % (seq_names[i], seq_names[j])})
                if (gpaired > aapaired and use_similarity) or (gpaired < aapaired and  not use_similarity):
                    aa2suspect[aa][seq_names[i]] = aa2suspect[aa].get(seq_names[i], 0) + 1
                    aa2suspect[aa][seq_names[j]] = aa2suspect[aa].get(seq_names[j], 0) + 1
                
                aa2suspect_dist[aa][seq_name[i]].append([gpaired, aapaired])
                aa2suspect_dist[aa][seq_name[j]].append([gpaired, aapaired])

        if self.mode == 'count':
            self.get_suspect_by_count(aa2suspect, number_seq)
        elif self.mode == 'wilcoxon':
            self.get_suspect_by_wilcoxon(aa2suspect_dist, number_seq, use_similarity)
        else :
            self.get_suspect_by_clustering(aa2suspect_dist, number_seq)


    def get_suspect_by_count(self,  aa2suspect, seq_num):
        for aa in aa2suspect.keys():
            tmp_list = aa2suspect[aa].most_common()
            i = 0
            while i < len(tmp_list) and tmp_list[i][1] > self.settings.FREQUENCY_THRESHOLD*seq_num:
                self.suspected_species[aa][tmp_list[i][0]] = tmp_list[i][1]
                i += 1
            
    def get_suspect_by_wilcoxon(self, aa2suspect_dist, seq_num, use_similarity=1):
        for aa in aa2suspect_dist.keys():
            tmp_list = aa2suspect_dist[aa]
            for seq, val in tmp_list.iteritems():
                val =  np.asarray(val)
                rank, pval = self.get_wilcoxon(val[:, 0], val[:, 1], use_similarity)
                if pval < self.settings.pval:
                    self.suspected_species[aa][seq] = pval


    def get_suspect_by_clustering(aa2suspect_dist, number_seq, use_similarity=1):
        for aa in aa2suspect_dist.keys():
            aafilt2dict =  self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            tmp_list = aa2suspect_dist[aa]
            cluster_list = {}
            aa_set = set([])
            for seq, val in tmp_list.iteritems():
                val =  np.mean(val, axis=0)
                if (val[0] > val[1] and use_similarity) or (val[0] < val[1] and not use_similarity):
                    aa_set.add(seq)
                aa_freq = np.mean([x == aa for x in aafilt2dict[seq]])
                cluster_list[seq] = np.append(val, aa_freq)

            centroid, labels = self.get_cluster(cluster_list)
            lab1 = set(np.asarray(cluster_list.keys())[labels==True])
            lab2 = set(np.asarray(cluster_list.keys())[labels==False])
            size1 = aa_set.intersection(lab1)
            size2 = aa_set.intersection(lab2)
            lab = lab1 if size1> size2 else lab2
            for specie in lab:
                self.suspected_species[aa][seq] = 1


    def get_cluster(self, cluster_list):
        features =  np.asarray(cluster_list.values())
        centroid, labels = kmeans2(features, 2, iter=100)


    def get_wilcoxon(self, y1, y2, use_similarity):
        rank, pval = ss.wilcoxon(y1, y2, correction=True)
        if use_similarity :
            return rank, pval/2.0
        else :
            return rank, 1 - pval/2.0


    def possible_aa_reassignation(self):
        """Find to which amino acid the reassignation are probable"""
        for aa, suspect in self.suspected_species.iteritems():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            suspected_aa = []
            for spec in suspect.keys():
                for cur_aa in aa_alignment[spec]:
                    if cur_aa !='-' and cur_aa != aa:
                        suspected_aa.append(cur_aa)
                        try : 
                            self.aa2aa_rea[aa][cur_aa].add(spec)
                        except KeyError:
                            self.aa2aa_rea[aa] = defaultdict(set)
                            self.aa2aa_rea[aa][cur_aa].add(spec)

    
    def save_json(self):
        with open(os.path.join(self.settings.OUTDIR,"aause.json"), "w") as outfile1:
            json.dump(self.aa_sim_json, outfile1, indent=4)
        with open(os.path.join(self.settings.OUTDIR,"similarity.json"), "w") as outfile2:
            json.dump(self.sim_json, outfile2, indent=4)
        with open(os.path.join(self.settings.OUTDIR,"reassignment.json"), "w") as outfile3:
            json.dump(self.reassignment_mapper, outfile3, indent=4)


    def save_all(self):
        self.save_json()
        id_filtfile = os.path.join(self.settings.OUTDIR,"id_filt.fasta")
        ic_filtfile = os.path.join(self.settings.OUTDIR,"ic_filt.fasta")
        gap_filtfile = os.path.join(self.settings.OUTDIR,"gap_filt.fasta")
        newick = os.path.join(self.settings.OUTDIR,"tree.nwk")
        self.seqset.write_data(id_filtered=id_filtfile, gap_filtered=gap_filtfile, ic_filtered=ic_filtfile, tree=newick)


    def run_analysis(self):
        """ Run the filtering analysis of the current dataset in sequenceset"""
        self.get_genomes()
        self.interesting_case = []
        self.possible_aa_reassignation()
        codon_align, fcodon_align = self.seqset.get_codon_alignment()
        for aa1, aarea  in self.aa2aa_rea.iteritems():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            for aa2, species in aarea.iteritems():
                gcodon_rea = CodonReaData((aa1, aa2), self.global_consensus, codon_align, self.seqset.codontable)
                fcodon_rea = CodonReaData((aa1, aa2), self.filtered_consensus, fcodon_align, self.seqset.codontable)
                counts = []
                t = self.seqset.phylotree.copy("newick")
                fitch = NaiveFitch(t, val, aa_letters_1to3[aa2], aa_letters_1to3[aa1], self.settings,  self.seqset.codontable, (gcodon_rea, fcodon_rea))
                slist = fitch.get_species_list(self.settings.LIMIT_TO_SUSPECTED_SPECIES)

                for genome in slist :
                    rec = self.seqset.prot_dict[genome]
                    leaf = fitch.tree&genome
                    filt_count = 0
                    try :
                        filt_count = len([y for y in aa_alignment[genome] if y !='-' and y==aa2])
                    except Exception:
                        #wtver happen, do nothing
                        pass
                    leaf.add_features(count=0)
                    leaf.add_features(filter_count=filt_count)
                    leaf.add_features(lost=False)
                    for pos in xrange(len(self.global_consensus)):
                        if self.global_consensus[pos] == aa1 and rec[pos] == aa2:
                            leaf.count += 1
                    if('lost' in leaf.features and leaf.count < settings.COUNT_THRESHOLD):
                        leaf.lost = True

                    self.reassignment_mapper['genome'][genome] = gcodon_rea.get_aa_usage(genome) 
                    gdata = {}
                    greacodon = gcodon_rea.get_reacodons(genome)
                    gusedcodon = gcodon_rea.get_usedcodons(genome)
                    gmixtecodon = gcodon_rea.get_mixtecodons(genome)
                    gdata['global'] = {'rea_codon' : greacodon, 'use_codon': gsusedcodon, 'mixte_codon' : gmixtecodon, 'count' : leaf.count}
                    freacodon = fcodon_rea.get_reacodons(genome)
                    fusedcodon = fcodon_rea.get_usedcodons(genome)
                    fmixtecodon = fcodon_rea.get_mixtecodons(genome)
                    gdata['filtered'] = {'rea_codon' : freacodon, 'use_codon': fsusedcodon, 'mixte_codon' : fmixtecodon, 'count' : leaf.filter_count}
                    gdata['suspected'] = self.suspected_species[aa1][genome]
                    gdata['fitch'] = fitch.get_distance_to_rea_node(leaf)
                    gdata['codons'] = self.reassignment_mapper['genome'][genome][aa2]
                    self.reassignment_mapper[aa2][aa1][genome] = gdata
                
            if(fitch.is_valid()):
                fitch.render_tree(self.settings.suffix)
                self.interesting_case.append(n)

        logging.debug("After validating the ancestral state and checking in the global alignment, %d cases were found interesting"%len(self.interesting_case))


def executeCMD(cmd, prog):
    """Execute a command line in the shell"""
    print "\n", cmd
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    Output.error(err)
    print "%s : \n----------------- %s" % (prog, out)
    return err

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
    convert_tree_to_mafft(Tree(tree.write(
        features=['name', 'support', 'dist'], format_root_node=True)), seq_ord, out)
    out.close()
    execute_mafft(enable_mafft + " --treein %s %s > %s" %
                  (out.file, tmpseq, inseq))
    msa = AlignIO.read(inseq, 'fasta', alphabet=alpha)
    filelist = glob.glob(TMP + "tmp0*")
    for f in filelist:
        os.remove(f)
    return msa


def align(seqlist):
    pass

def timeit(func):

    def timed(*args, **kw):
        tstart = time.time()
        result = func(*args, **kw)
        tend = time.time()
        ttime = tend - tstart
        #print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
        return ttime, result

    return timed

@timeit
def execute_mafft(cmdline):
    """Execute the mafft command line"""
    executeCMD(cmdline, 'mafft')