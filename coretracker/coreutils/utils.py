from __future__ import division

import argparse
import glob
import itertools
import json
import logging
import os
import re
import random
import shutil
import subprocess
import sys
import time
import pandas as pd
from collections import Counter, defaultdict
from cStringIO import StringIO
from distutils import spawn
from functools import partial

import Bio.SubsMat.MatrixInfo as MatrixInfo
import numpy as np
import scipy.stats as ss
from Bio import AlignIO, Alphabet, SeqIO, SubsMat, codonalign
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet import IUPAC, generic_nucleotide, generic_protein
from Bio.codonalign.codonalphabet import (default_codon_alphabet,
                                          get_codon_alphabet)
from Bio.codonalign.codonseq import CodonSeq, _get_codon_list
from Bio.Data import CodonTable
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from ete3 import Tree
from scipy.cluster.vq import kmeans2

from coretracker.FisherExact import fisher_exact
from Faces import LineFace, PPieChartFace, SequenceFace
from corefile import CoreFile
from output import Output
from pdfutils import *
import matplotlib.pyplot as plt

__author__ = "Emmanuel Noutahi"
__version__ = "1.01"
__email__ = "fmr.noutahi@umontreal.ca"
__license__ = "The MIT License (MIT)"

# conditional import
GRAPHICAL_ACCESS = True
try:
    from ete3 import TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
except ImportError, e:
    GRAPHICAL_ACCESS = False

# global declaration
aa_letters_1to3 = {

    'A': 'Ala', 'C': 'Cys', 'D': 'Asp',
    'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met',
    'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp',
    'Y': 'Tyr',
}

# Liste of available matrices
AVAILABLE_MAT = MatrixInfo.available_matrices
# define array to contains temp values
glob_purge = []
# hmm strange id
hmmidpattern = re.compile("^\d+\|\w+")
# define lowest pvalue
eps = np.finfo(np.float).eps
aa_letters_3to1 = dict((x[1], x[0]) for x in aa_letters_1to3.items())
alpha = Alphabet.Gapped(IUPAC.protein)
nuc_letters = "ACTG"

class SequenceLoader:
    """Load and format data for analysis"""

    def __init__(self, infile, dnafile, settings, gap_thresh, has_stop=True,
                 use_tree=None, refine_alignment=True, seqformat="core", hmmloop=10, msaprog="muscle", hmmdict={}):
        # try parsing the protein file as a corefile and if it failed
        # parse it as a fasta file
        self.infile = infile
        self.sequences = {}
        stop_removed = False
        self.dnasequences = {}
        self.hmmdict = hmmdict
        self.hmmloop = 10
        gap_thresh = (abs(gap_thresh) <= 1 or 0.01) * abs(gap_thresh)
        try:
            self.sequences = self.get_sequences(
                infile, seqformat, generic_protein)
            # we got to the end without finding any sequences
            if not self.sequences:
                raise ValueError('Unable to load sequences')
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
            self.dnasequences = self.get_sequences(
                dnafile, seqformat, generic_nucleotide)
            # we got to the end without finding any sequences
        except:
            raise ValueError('Nucleotide is not in the correct format')

        logging.debug('DNA sequence was read')
        self.genes = set(self.dnasequences.keys()).intersection(
            self.sequences.keys())

        common_spec_per_gene = {}
        common_spec_per_gene_len = {}
        max_spec = -np.inf
        self.genome_list = set([])
        # Only interested in gene in both nucleotide and aa sequences
        for gene in self.genes:
            g_dnaspec = [x.id for x in self.dnasequences[gene]]
            g_protspec = [x.id for x in self.sequences[gene]]
            g_common_spec = set(g_dnaspec).intersection(set(g_protspec))
            common_spec_per_gene[gene] = g_common_spec
            self.genome_list.update(g_common_spec)
            cur_size = len(g_common_spec)
            max_spec = max_spec if max_spec > cur_size else cur_size
            common_spec_per_gene_len[gene] = cur_size

        # dropping gene with present in not enough species
        # using gap_thresh to perform this
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

        # align sequences if it's not already done
        if not is_aligned:
            alignment = self.align(
                settings, refine=refine_alignment, tree=use_tree, msaprog=msaprog, scale=settings.scale)
            logging.debug('Sequence is not aligned, align and refine it.')

        # case where sequence is not aligned but refinement is required
        elif is_aligned and refine_alignment:
            logging.debug(
                'Sequence is already aligned but refine was requested')
            for g in self.genes:
                alignment[g] = self.__class__._refine(alignment[g], 9999, settings.OUTDIR,
                                                      settings.hmmbuild, settings.hmmalign, settings.eslalimanip,
                                                      settings.eslalimask, loop=self.hmmloop, hmmfile=self.hmmdict.get(g, None))
        self.alignment = alignment
        logging.debug('Sequence alignment done')

    def get_sequences(self, infile, fileformat, alphabet):
        """Get sequence from file"""
        if fileformat == "core":
            corefile = CoreFile(infile, alphabet)
            return corefile.get_sequences()
        else:
            raise NotImplementedError("Others format are not yet supported")

    def save_alignment(self, outfile):
        """save alignment in a file """
        AlignIO.write(self.alignment, open(outfile, 'w'), 'fasta')

    def concat(self, missing='-', alpha=generic_protein):
        """Concatenate alignment into one global alignment"""
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
                dna_adding = SeqRecord(
                    Seq("", generic_nucleotide), id=spec, name=spec)
                if(spec in self.common_spec_per_gene[gene]):
                    adding = spec_dict[spec]
                    dna_adding = dna_spec_dict[spec]
                    # guess in this case there is a stop at the end of the dna
                    if len(adding.seq.ungap(missing))*3 +3 == len(dna_adding.seq.ungap(missing)):
                        dna_adding = dna_adding[:-3]

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
        """Align sequences"""
        alignment = {}
        outdir = settings.OUTDIR
        for gene in self.genes:
            # in order to have the best possible alignment, we are going to
            # keep the
            al = self.__class__._align(
                self.sequences[gene], msaprog, tree, scale, outdir, alpha)
            if refine:
                al = self.__class__._refine(
                    al, 9999, outdir, settings.hmmbuild, settings.hmmalign, settings.eslalimanip,
                    settings.eslalimask, loop=self.hmmloop, hmmfile=self.hmmdict.get(gene, None))
            alignment[gene] = al
        return alignment

    @classmethod
    def _align(clc, msa, msaprog, tree, scale, outdir, alpha=generic_protein, is_aligned=False):
        """Align sequences using muscle of mafft"""
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
    def _refine(clc, alignment, timeout, outdir, hmmbuild="hmmbuild", hmmalign="hmmalign", eslalimanip="esl-alimanip",
                eslalimask="esl-alimask", hmmfile=None, loop=10, minqual=8,  strategie="greater", clean=True):
        """Align and refine at the same time with hmmalign and muscle"""

        success, not_found = check_binaries(
            hmmbuild, hmmalign, eslalimask, eslalimanip)
        if not success:
            raise RuntimeError(
                "Could not refine alignment, Executable not found: %s!!!" % (not_found))
        strategie = (strategie == "greater")

        def accessQuality(alignfile, minqual, greater=True):
            """Check quality of the hmm alignment"""
            ppcons = ''
            with open(alignfile, 'r') as text:
                for line in text:
                    if '#=GC' in line and 'PP_cons' in line:
                        ppcons += line.split()[-1].strip()
            highpos = []
            if greater:
                highpos = [1 for x in ppcons if x ==
                           '*' or (isInt(x) and int(x) >= minqual)]
            else:
                highpos = [(10 if x == '*' else int(x))
                           for x in ppcons if (x == '*' or isInt(x))]

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
            outputFile = os.path.join(cur_dir[i], "alignment.sto")
            if hmmfile:
                # in this case copy it to the current dir
                tmpfile = os.path.join(cur_dir[i], "alignment.hmm")
                shutil.copyfile(hmmfile, tmpfile)
                hmmfile = tmpfile
            else:
                hmmfile = os.path.join(cur_dir[i], "alignment.hmm")
                buildline = hmmbuild + " --amino %s %s" % (hmmfile, inputFile)
                executeCMD(buildline, 'hmmbuild')
            # will continue if not exception is found
            alignline = hmmalign + \
                " -o %s %s %s" % (outputFile, hmmfile, inputFile)
            executeCMD(alignline, 'hmmalign')
            # finding quality
            quality.append(accessQuality(outputFile,  minqual, strategie))
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
        """Check if protein sequence are aligned"""
        alignment = {}
        try:
            for gene in self.genes:
                alignment[gene] = MultipleSeqAlignment(self.sequences[gene])
        except ValueError:
            return False, None

        return True, alignment

    def _split_at_stop(self, seqlist, is_aligned):
        """Split the sequence at stop position in order to simulate genes"""
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
    def translate(clc, core_inst, gcode=1):
        """Translate nucleotide sequences into protein sequence given a genetic code"""
        codontable = CodonTable.unambiguous_dna_by_id[1]
        try:
            codontable = CodonTable.unambiguous_dna_by_id[abs(gcode)]
        except:
            logging.warn("Wrong genetic code, resetting it to 1")

        if isinstance(core_inst, basestring):
            core_inst = CoreFile(core_inst, alphabet=generic_nucleotide)

        core_prot_inst = {}
        stop_checker = lambda x : '*' if x.upper() in codontable.stop_codons else 'X'
        for gene, seqs in core_inst.items():
            translated = []
            for s in seqs:
                if len(s) % 3 != 0:
                    raise ValueError(
                        "Frame-shifting detected in %s : [%s], current version does not supported it." % (s.id, gene))
                else:
                    aa_list = [codontable.forward_table.get(
                        s.seq[i:i + 3].upper(), stop_checker(s.seq[i:i + 3])) for i in xrange(0, len(s), 3)]
                    trans_s = Seq("".join(aa_list), generic_protein)
                    translated.append(SeqRecord(trans_s, id=s.id, name=s.name))
            core_prot_inst[gene] = translated

        return core_prot_inst


class CodonSeq(CodonSeq):
    """Codon seq representation, inherit from Bio-python codon seq"""

    def __init__(self, data='', alphabet=default_codon_alphabet,
                 gap_char="-", rf_table=None, enable_undef=False):
        Seq.__init__(self, data.upper(), alphabet=alphabet)
        self.gap_char = gap_char
        self.enable_undef = enable_undef

        # check the length of the alignment to be a triple
        if rf_table is None:
            seq_ungapped = self._data.replace(gap_char, "")
            assert len(self) % 3 == 0, "Sequence length is not divisible by 3"
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
        """Check whether or not this codon is ambiguous"""
        return ('N' in codon.upper()) and self.enable_undef


class CodonReaData(object):
    """A representation of a reassignment in a species"""

    def __init__(self, aas, alignment, consensus, codon_align_dict, dct, positions, genelimit, settings):

        # change aa2 is a list now
        self.aa1, aareas = aas
        self.aa2_list = aareas.keys()
        self.codon_alignment = codon_align_dict
        self.dct = dct
        self.back_table = init_back_table(dct)

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

    def _spec_codon_usage(self):
        """Check codon usage in each species"""
        for i, aa in enumerate(self.consensus):
            for spec in self.codon_alignment.keys():
                spec_codon = self.codon_alignment[spec].seq.get_codon(i)
                spec_aa = self.dct.forward_table.get(spec_codon, None)
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
                        # this mean, species use aa2, while we don't care about
                        # the major aa
                        else:
                            self.mixtecodon[spec][spec_aa].update([spec_codon])

    def get_score(self, spec, aa_ori, aa_rea):
        """Get Telford score for each codon"""
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
        """Get the list of rea codons"""
        return self.reacodons[specie][aa]

    def get_mixtecodons(self, specie, aa):
        """Get the list of mixte codons"""
        return self.mixtecodon[specie][aa]

    def get_usedcodons(self, specie, aa):
        """Get the list of normally used codons"""
        return self.usedcodons[specie][aa]

    def get_aa_usage(self, specie, aa):
        """Get aa usage in a specific genome"""
        return self.specs_amino_count[specie][aa]

    def get_all_aas_usage(self, specie):
        """Get all aa usage in a specific genome"""
        return self.specs_amino_count[specie]

    def get_rea_aa_codon_distribution(self, specie, aa):
        """Get codon distribution in potentially reassigned positions"""
        return dict((k, list(v)) for k, v in self.codons_distribution[specie][aa].items())

    def get_total_rea_aa_codon_distribution(self, specie, aa):
        """Get total codon distribution"""
        return dict((k, list(v)) for k, v in self.total_codons_distribution[specie][aa].items())


class NaiveFitch(object):
    """A NaiveFitch algorithm for finding the most parcimonious solution"""

    def __init__(self, tree, reassigned, ori_aa, dest_aa, dct, codon_rea=(None, None)):
        self.id = {}
        self.tree = tree
        self.corr = {'0': ori_aa, '1': dest_aa}
        self.codon_rea_global, self.codon_rea_filtered = codon_rea
        colors = ['#a6cee3', '#1f78b4', '#b2df8a',
                  '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f']
        self.dct = dct
        self.back_table = init_back_table(dct)
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
        self.ori_aa1 = aa_letters_3to1[ori_aa]
        self.dest_aa1 = aa_letters_3to1[dest_aa]
        self.newick = tree.write(features=['name', 'dist', 'support', 'state'])
        self._bottomup(self.tree)
        self._topdown(self.tree)

    def update_codon_data(codon_rea):
        """Update codon reassignment data (global and filtered)
        """
        self.codon_rea_global, self.codon_rea_filtered = codon_rea

    def write_tree(self, outfile):
        """Export newick to  a file"""
        with open(outfile, 'w') as OUT:
            OUT.write(self.newick)

    def is_valid(self, thresh=1):
        """Naive function to test whether or not if the current reassignment is valid"""
        tmptree = self.tree.copy()

        for l in tmptree.traverse():
            if not l.is_leaf():
                l.del_feature('reassigned')
                l.del_feature('rea')
            elif (l.lost and l.state == self.dest_aa) or l.count < thresh:
                l.add_features(reassigned={0})
                l.add_features(rea=self.corr['0'])
                l.add_features(state=self.ori_aa)
            elif not l.lost and l.state == self.ori_aa and l.count >= thresh:
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
        # Since it's not need, just pass
        pass

    def is_reassigned(cls, node, strict=True):
        """ return True if a node has undergoned reassignment """
        if "reassigned" in node.features and (node.reassigned == {1} or (not strict and 1 in node.reassigned)):
            return True
        return False

    def get_species_list(self, limit_to_suspected=False):
        """Get the current species list for this tree"""
        slist = set()
        if limit_to_suspected:
            for node in self.tree.traverse():
                if 'reassigned' in node.features and node.reassigned == {1}:
                    slist.update(node.get_leaf_names())
        else:
            slist = set(self.tree.get_tree_root().get_leaf_names())
        return slist

    def get_distance_to_rea_node(self, node):
        """Get the distance of node to the closest reassigned lca"""
        cur_node = node
        if isinstance(node, str):
            cur_node = self.tree & node
        dist = 0
        while cur_node != None and not self.is_reassigned(cur_node, strict=False):
            dist += 1
            cur_node = cur_node.up
        return dist

    def has_codon_data(self):
        # base codon_data on filtered position only
        return self.codon_rea_filtered != None


class SequenceSet(object):
    """Representation of the current genome set and their sequences"""

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
        """ Get log matrix for current alignment"""
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
        c_genome = set([x for x in common_genome if (
                len(self.dna_dict[x].seq.ungap('-')) == len(self.prot_dict[x].seq.ungap('-')) * 3)])

        #print [(x, len(self.dna_dict[x]), len(self.prot_dict[x].seq.ungap('-'))*3, len(self.dna_dict[x])==(len(self.prot_dict[x].seq.ungap('-'))+1)*3) for x in common_genome ]
        if not c_genome:
            help1 = "1) Protein length do not match with dna and stop already checked! Look for mistranslation of Frame-shifting"
            help2 = "2) Wrong tree or sequences name not matching!"

            raise ValueError(
                    'ID intersection for dna, prot and species tree is empty. Possible cause:\n%s\n%s\n'%(help1, help2))
        common_genome = c_genome
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

    @classmethod
    def copy_codon_alignment(clc, codon_alignment, codontable, alphabet=default_codon_alphabet):
        """Return a codon alignment object from a list of codon seq sequences"""
        return codonalign.CodonAlignment((clc.copy_codon_record(rec, codontable) for rec in codon_alignment._records), alphabet=alphabet)

    def codon_align(self, alphabet=default_codon_alphabet, gap_char='-'):
        """ Perform a codon alignment based on a protein multiple alignment
        and remove gaped positions
        """
        alphabet = get_codon_alphabet(self.codontable, gap_char=gap_char)
        if self.common_genome is None:
            self.restrict_to_common()
        # build codon alignment and return it
        codon_aln = []
        all_undef_codon = {}
        for g in self.dna_dict.keys():
            codon_rec, undef_c = self._get_codon_record(
                self.dna_dict[g], self.prot_dict[g], self.codontable, alphabet)
            all_undef_codon[g] = undef_c
            codon_aln.append(codon_rec)

        self.codon_alignment = codonalign.CodonAlignment(
            codon_aln, alphabet=alphabet)

        if all_undef_codon:
            undef_codon = set().union(*all_undef_codon.values())
            logging.debug("Total position lost due to undef codon : %d" %
                          len(undef_codon))
            undef_codon = sorted(list(undef_codon))
            #R_aa_position = np.asarray(xrange(self.prot_align.get_alignment_length()))
            #R_aa_position = np.setxor1d(tt_filter_position, np.asarray(undef_codon))
            for gene, seqrec in self.prot_dict.items():
                seqrec.seq._data = seqrec.seq._data.replace('X', gap_char)
                #self.prot_dict[gene] = seqrec
            self.prot_align = MultipleSeqAlignment(
                self.prot_dict.values(), alphabet=alpha)

            # remove all the position with undef codon from the dna_dict
            for codseqrec in self.codon_alignment:
                k = codseqrec.id
                self.dna_dict[k] = SeqRecord(
                    codseqrec.seq.toSeq().ungap(gap_char), id=k, name=k)

        #self.codon_alignment = codonalign.build(self.prot_align, self.dna_dict, codon_table=self.codontable)

    def _get_codon_record(self, dnarec, protrec, codontable, alphabet, gap_char='-'):
        """Get a codon seq record from dna and prot seqrecord"""
        nuc_seq = dnarec.seq.ungap(gap_char)
        codon_seq = ""
        aa_num = 0
        max_error = 1000
        x_undecoded = []
        undef = {"R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
        for aa in protrec.seq:
            if aa == gap_char:
                codon_seq += gap_char * 3
            else:
                next_codon = nuc_seq._data[aa_num * 3: (aa_num + 1) * 3]

                if undef.intersection(set(next_codon.upper())):
                    x_undecoded.append(aa_num)
                    if aa.upper() != 'X':
                        logging.warn("%s(%s %d) decoded by undefined codon %s(%s)"
                                     % (protrec.id, aa, aa_num, dnarec.id, next_codon))
                        max_error -= 1

                    next_codon = gap_char * 3

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

    def filter_codon_alignment(self, codon_alignment=None, ind_array=None, get_dict=False, alphabet=default_codon_alphabet):
        """Return the codon aligment from a list of in array"""
        if not ind_array:
            ind_array = (self.filt_position, )

        if codon_alignment is None:
            codon_alignment = self.codon_alignment

        if not isinstance(ind_array, tuple):
            raise ValueError('Expect tuple for ind_array, got %s' %
                             (type(ind_array)))
        else:
            for indexes in ind_array:
                filt_codon_align = self.filter_align_position(codon_alignment, indexes,
                                                    alphabet=alphabet, codontable=self.codontable)
                if get_dict:
                    filt_codon_align = SeqIO.to_dict(filt_codon_align)
                yield filt_codon_align

    def get_codon_alignment(self):
        """Get codon alignment"""
        r = self.filter_codon_alignment()
        self.fcodon_alignment = next(r)
        return self.codon_alignment, self.fcodon_alignment

    def prot_filtering(self, id_thresh=None, gap_thresh=None, ic_thresh=None, rmcnst=True):
        """Filter protein alignment"""
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
    def filter_align_position(clc, alignment, index_array, alphabet=generic_protein, codontable=None):
        """Keep only columns specified by index_array from the alignment"""
        if not codontable:
            edited_alignment = alignment[:, :]
            index_array = sorted(index_array)
            # alignment[:, slice(index_array[0], index_array[0] + 1)]
            for seqrec in edited_alignment:
                seq = Seq('', alphabet)
                for i in index_array:
                    # do not know why this raise an error on pyhton3
                    seq += seqrec[int(i)]
                seqrec.letter_annotations = _RestrictedDict(length=len(seq))
                seqrec.seq = seq.upper()

        else:
            edited_alignment = clc.copy_codon_alignment(alignment, codontable)
            index_array = sorted(index_array)
            for i in xrange(len(edited_alignment)):
                codseq = CodonSeq('')
                for pos in index_array:
                    codseq += edited_alignment[i].seq.get_codon(pos)
                edited_alignment[i].seq = codseq
        return edited_alignment

    @classmethod
    def clean_alignment(clc, alignment=None, characs=['-'], threshold=0.5):
        """Remove position of alignment which contain character from characs"""
        align_array = np.array([list(rec) for rec in alignment], dtype=str)
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

    def get_genome_size(self, aligntype='global', gap_char='-'):
        use_alignment = self.prot_align
        if aligntype == 'filtered':
            use_alignment = self.filt_prot_align
        gsize = dict((srec.id, len(srec.seq.ungap(gap_char)))
                     for srec in use_alignment)
        return gsize


class ReaGenomeFinder:
    """Find genome reassignment"""

    def __init__(self, seqset, settings):
        self.seqset = seqset
        self.mode = getattr(settings, 'mode', "count")
        self.method = getattr(settings, 'method', "identity")
        self.confd = getattr(settings, 'conf', 0.05)
        self.suspected_species = defaultdict(dict)
        self.aa2aa_rea = defaultdict(dict)
        self.settings = settings
        self.interesting_case = []
        self.reassignment_mapper = makehash()

    def compute_sequence_identity(self, matCalc=None):
        """Compute a distance matrix from the alignment"""
        if not matCalc:
            matCalc = DistanceCalculator(self.method)
        self.global_paired_distance = matCalc.get_distance(
            self.seqset.prot_align)
        self.filtered_paired_distance = matCalc.get_distance(
            self.seqset.filt_prot_align)

    def get_genomes(self, use_similarity=1):
        """ Get suspected genomes """
        matCalc = DistanceCalculator(self.method)
        self.compute_sequence_identity(matCalc)
        #logging.debug("Distance matrix : ")
        # logging.debug(self.filtered_paired_distance)
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
        #
        for aa in self.settings.AA_LETTERS:
            cons_array = self.seqset.get_aa_filtered_alignment(
                self.filtered_consensus, aa)
            not_uniq = len(set(cons_array)) > 1
            if(cons_array and not_uniq and len(cons_array) > self.settings.COUNT_THRESHOLD):
                self.seqset.aa_filt_prot_align[aa_letters_1to3[
                    aa]] = SequenceSet.filter_align_position(self.seqset.filt_prot_align, cons_array)
                self.aa_paired_distance[aa] = matCalc.get_distance(
                    self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]])
                self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]] = SeqIO.to_dict(
                    self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]])
                # logging.debug("Distance matrix for amino acid %s: " % aa)
                # logging.debug(self.aa_paired_distance[aa])

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
                    aa2suspect[aa][self.seq_names[i]] = aa2suspect[aa].get(self.seq_names[i], 0) + 1
                    aa2suspect[aa][self.seq_names[j]] = aa2suspect[aa].get(self.seq_names[j], 0) + 1

                aa2suspect_dist[aa][self.seq_names[i]].append([paired, aapaired])
                aa2suspect_dist[aa][self.seq_names[j]].append([paired, aapaired])

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
        """Find suspected species by simple counting"""
        for aa in aa2suspect.keys():
            tmp_list = aa2suspect[aa].most_common()
            i = 0
            while i < len(tmp_list) and tmp_list[i][1] > self.settings.FREQUENCY_THRESHOLD * seq_num:
                self.suspected_species[aa][
                    tmp_list[i][0]] = 1  # tmp_list[i][1]
                i += 1

    def get_suspect_by_stat(self, aa2suspect_dist, seq_num, use_similarity=1, test='wilcoxon', confd=0.05):
        """Use a statistic test to find suspected species"""
        for aa in aa2suspect_dist.keys():
            tmp_list = aa2suspect_dist[aa]
            for seq, val in tmp_list.items():
                val = np.asarray(val)
                rank, pval = self.get_paired_test(
                    val[:, 0], val[:, 1], use_similarity, test)
                if pval <= confd:
                    self.suspected_species[aa][seq] = pval

    def get_suspect_by_clustering(self, aa2suspect_dist, number_seq, use_similarity=1):
        """Get list of suspected genome using clustering"""

        def get_cluster(cluster_list):
            """Perform kmean in order to find clusters"""
            features = np.asarray(cluster_list.values())
            return kmeans2(features, 2, iter=100)

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

            centroid, labels = get_cluster(cluster_list)
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

    @classmethod
    def get_paired_test(clc, y1, y2, use_similarity, test="wilcoxon"):
        """Return a paired-test pvalue for the hypothesis u(y1)=u(y2)"""
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

    def get_aa_count_in_alignment(self):
        """Get the aa count for each aa in the alignment"""
        aa_count = defaultdict(Counter)
        f_prot_align = SeqIO.to_dict(self.seqset.filt_prot_align)
        for spec in f_prot_align:
            aa_count[spec] = Counter(f_prot_align[spec])
        return aa_count

    def get_prob_thresh(self, c1, c2, l, prob):
        """Get virtual prob threshold. This is temp"""
        logodd = np.log((c1 * c2) / (prob * (l**2)))
        return np.exp(logodd) / (1 + np.exp(logodd))

    def set_rea_mapper(self):
        """Set common values of reassignment_mapper for all reassignment"""
        self.reassignment_mapper['genes'] = self.seqset.get_genes_per_species()
        self.reassignment_mapper['genome'][
                'global'] = self.seqset.get_genome_size()
        self.reassignment_mapper['genome'][
                'filtered'] = self.seqset.get_genome_size("filtered")

    def possible_aa_reassignation(self):
        """Find to which amino acid the reassignation are probable"""
        # TODO :  CHANGE THIS FUNCTION TO A MORE REALISTIC ONE

        aa_count_spec = self.get_aa_count_in_alignment()
        aa_count_cons = Counter(self.filtered_consensus)
        for aa, suspect in self.suspected_species.items():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            suspected_list = sorted(suspect.keys())
            for spec in self.seqset.common_genome:
                spec_aa_counter = Counter(aa_alignment[spec])
                for (cur_aa, aacounter) in spec_aa_counter.most_common():
                    aa_c1 = aa_count_spec[spec][aa] + aa_count_cons[aa]
                    cur_aa_c2 = aa_count_spec[spec][
                        cur_aa] + aa_count_cons[cur_aa]
                    filt_size = sum([aa_count_spec[spec][x]
                                     for x in aa_count_spec[spec].keys() if x != '-'])
                    cons_size = sum(
                        [aa_count_cons[x] for x in aa_count_spec.keys() if x not in ('-', 'X')])
                    tot = sum([spec_aa_counter[x]
                               for x in spec_aa_counter.keys() if x != '-'])
                    prob = spec_aa_counter[cur_aa] / tot if tot>0 else 0
                    if prob>0 and cur_aa != '-' and cur_aa != aa and cur_aa not in self.settings.EXCLUDE_AA_FROM:

                        sc = self.get_prob_thresh(
                            aa_c1 + aa_count_cons[aa], cur_aa_c2, filt_size + cons_size, prob)
                        # print '---------------------------------------------'
                        # print cur_aa_c2, aa_c1, filt_size, spec_aa_counter[cur_aa], tot
                        # print '%s : %s ==> %s = %f (%d)' %(spec, cur_aa, aa, sc, aacounter)
                        # our list of suspected is the genome that either pass sc test
                        # or were suspected with a total count for amino acid
                        was_suspected = (spec in suspected_list)
                        if sc < self.confd or (was_suspected and aa_count_spec[spec][cur_aa] > 1):
                            try:
                                self.aa2aa_rea[aa][cur_aa].add(spec)
                            except KeyError:
                                self.aa2aa_rea[aa] = defaultdict(set)
                                self.aa2aa_rea[aa][cur_aa].add(spec)

    def save_json(self):
        """Save result into a json file"""
        with open(os.path.join(self.settings.OUTDIR, "aause.json"), "w") as outfile1:
            json.dump(self.aa_sim_json, outfile1, indent=4)
        with open(os.path.join(self.settings.OUTDIR, "similarity.json"), "w") as outfile2:
            json.dump(self.sim_json, outfile2, indent=4)
        with open(os.path.join(self.settings.OUTDIR, "reassignment.json"), "w") as outfile3:
            json.dump(self.reassignment_mapper, outfile3, indent=4)

    def save_all(self):
        """Save everything"""
        self.save_json()
        id_filtfile = os.path.join(self.settings.OUTDIR, "id_filt.fasta")
        ic_filtfile = os.path.join(self.settings.OUTDIR, "ic_filt.fasta")
        gap_filtfile = os.path.join(self.settings.OUTDIR, "gap_filt.fasta")
        newick = os.path.join(self.settings.OUTDIR, "tree.nwk")
        self.seqset.write_data(
            id_filtered=id_filtfile, gap_filtered=gap_filtfile, ic_filtered=ic_filtfile, tree=newick)

    def run_analysis(self, codon_align, fcodon_align):
        """ Run the filtering analysis of the current dataset in sequenceset"""

        for aa1, aarea in self.aa2aa_rea.items():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa1]]
            gcodon_rea = CodonReaData((aa1, aarea), self.seqset.prot_align, self.global_consensus, codon_align,
                                      self.seqset.codontable, self.seqset.position, self.seqset.gene_limits, self.settings)
            fcodon_rea = CodonReaData((aa1, aarea), self.seqset.filt_prot_align, self.filtered_consensus, fcodon_align,
                                      self.seqset.codontable, self.seqset.filt_position, self.seqset.gene_limits, self.settings)

            for aa2, species in aarea.items():
                # logging.debug("%s to %s" % (aa2, aa1))
                counts = []
                t = self.seqset.phylotree.copy("newick")
                fitch = NaiveFitch(t, species, aa_letters_1to3[aa2], aa_letters_1to3[
                                   aa1], self.seqset.codontable, (gcodon_rea, fcodon_rea))
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

                    reacodon = fcodon_rea.get_reacodons(genome, aa2)
                    usedcodon = fcodon_rea.get_usedcodons(genome, aa2)
                    # if not self.settings.USE_GLOBAL:
                    #    reacodon = fcodon_rea.get_reacodons(genome, aa2)
                    #    usedcodon = fcodon_rea.get_usedcodons(genome, aa2)

                    fisher_passed, pval = independance_test(
                        reacodon, usedcodon, confd=self.confd, tot_size=len(slist))
                    # print "%s\t%s\t%s ==> %s" %(aa2, aa1, genome,
                    # str(column_comparision) )
                    # if('lost' in leaf.features and fisher_passed and
                    # leaf.count > self.settings.COUNT_THRESHOLD
                    if('lost' in leaf.features and fisher_passed):
                        leaf.lost = False

                    gdata = {}
                    greacodon = gcodon_rea.get_reacodons(genome, aa2)
                    gusedcodon = gcodon_rea.get_usedcodons(genome, aa2)
                    gmixtecodon = gcodon_rea.get_mixtecodons(genome, aa2)
                    g_rea_dist = gcodon_rea.get_rea_aa_codon_distribution(
                        genome, aa2)
                    g_total_rea_dist = gcodon_rea.get_total_rea_aa_codon_distribution(
                        genome, aa2)
                    gdata['global'] = {'rea_codon': greacodon, 'used_codon': gusedcodon, 'mixte_codon': gmixtecodon,
                                       'count': leaf.count, 'rea_distribution': g_rea_dist, 'total_rea_distribution': g_total_rea_dist}
                    freacodon = fcodon_rea.get_reacodons(genome, aa2)
                    fusedcodon = fcodon_rea.get_usedcodons(genome, aa2)
                    fmixtecodon = fcodon_rea.get_mixtecodons(genome, aa2)
                    f_rea_dist = fcodon_rea.get_rea_aa_codon_distribution(
                        genome, aa2)
                    f_total_rea_dist = fcodon_rea.get_total_rea_aa_codon_distribution(
                        genome, aa2)
                    gdata['filtered'] = {'rea_codon': freacodon, 'used_codon': fusedcodon, 'mixte_codon': fmixtecodon,
                                         'count': leaf.filter_count, 'rea_distribution': f_rea_dist, 'total_rea_distribution': f_total_rea_dist}
                    gdata['suspected'] = self.suspected_species[
                        aa1].get(genome, 0)
                    gdata['score'] = {'global': gcodon_rea.get_score(
                        genome, aa2, aa1), 'filtered': fcodon_rea.get_score(genome, aa2, aa1)}

                    gdata['lost'] = {'pval': pval,
                                     'lost': 1 if leaf.lost else 0}
                    gdata['fitch'] = fitch.get_distance_to_rea_node(leaf)
                    gdata['codons'] = {'global': gcodon_rea.get_aa_usage(
                        genome, aa2), 'filtered': fcodon_rea.get_aa_usage(genome, aa2)}
                    alldata[genome] = gdata

                if(fitch.is_valid(self.settings.COUNT_THRESHOLD) or self.settings.showall):
                    self.reassignment_mapper['aa'][aa2][aa1] = alldata
                    self.interesting_case.append("%s to %s" % (aa2, aa1))
                    yield (self, fitch, alldata)


def init_back_table(dct):
    """Get back table for the current genetic code"""
    back_table = defaultdict(list)
    for aa, codon in zip(dct.forward_table.values(), dct.forward_table.keys()):
        back_table[aa].append(codon)
    return back_table


def executeCMD(cmd, prog):
    """Execute a command line in the shell"""
    logging.debug("The following will be executed : \n%s\n" % cmd)
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    if err:
        logging.debug(err)
    if out:
        logging.debug("%s : \n----------------- %s" % (prog, out))
    return err


def convert_tree_to_mafft(tree, seq_order, output, scale, dist_thresh=1e-10):
    """Convert a tree in a newick format to the mafft matrix format"""
    seqnames = [-1, -1]
    branchlens = [-1, -1]
    for node in tree.traverse("postorder"):
        if node.dist <= dist_thresh:
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


def check_stop(seq_list, is_aligned, stop='*'):
    """ Check all position for stop in the protein sequences"""
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
    """Clean a directory"""
    shutil.rmtree(dirname, ignore_errors=True)
    os.makedirs(dirname)
    return dirname


def check_binaries(*args):
    """check if binaries exist"""
    not_found = []
    for exe in args:
        if not spawn.find_executable(exe):
            not_found.append(exe)

    return not_found == [], not_found


def timeit(func):
    """ Time a function """
    def timed(*args, **kw):
        tstart = time.time()
        result = func(*args, **kw)
        tend = time.time()
        ttime = tend - tstart
        # print '%r (%r, %r) %2.2f sec' % (func.__name__, args, kw, ttime)
        return ttime, result
    return timed


def isInt(value):
    """ Verify if a string can be casted into an int"""
    try:
        int(value)
        return True
    except ValueError:
        return False


def remove_gap_position(alignfile, curformat):
    """Remove all gap position from a file and return a new file"""
    align = AlignIO.read(alignfile, curformat)
    align, positions = SequenceSet.clean_alignment(align, threshold=1)
    for i in xrange(len(align)):
        seqname = align._records[i].name
        if hmmidpattern.match(seqname):
            seqname = seqname.split('|')[1]
            align._records[i].name = seqname
            align._records[i].id = seqname

    fastafile = alignfile.split('.')[0] + ".fasta"
    AlignIO.write(align, open(fastafile, 'w'), 'fasta')
    return fastafile


def independance_test(rea, ori, confd=0.05, tot_size=1):
    """Perform a Fisher's Exact test"""
    codon_list = set(rea.keys())
    codon_list.update(ori.keys())
    codon_list = [x for x in codon_list if not (
        rea.get(x, 0) == 0 and ori.get(x, 0) == 0)]
    nelmt = len(codon_list)
    # line 1 is rea, line 2 is observed
    # in this case, we can perform a fisher test
    if nelmt > 1:
        obs = np.zeros((nelmt, 2))
        for i in xrange(nelmt):
            obs[i, 0] = rea.get(codon_list[i], 0)
            obs[i, 1] = ori.get(codon_list[i], 0)
        try:
            pval = fisher_exact(obs, midP=True, attempt=3)
            # fallback to chi2 test if fisher is impossible
        except:
            logging.debug("Using ch2 instead of FISHEREXACT")
            c, pval, dof, t =  ss.chi2_contingency(obs)
        return pval <= confd, pval


    # strangely, codon is used only in rea column
    # complete reassignment ??
    # to avoid returning 0, we return the smallest positive int
    elif len(rea.values()) > 0 and len(ori.values()) == 0:
        return True, eps
    # codon is used in both column
    elif len(rea.values()) > 0 and len(ori.values()) > 0:
        fpval = abs(rea.values()[0] - ori.values()[0]) / \
            (rea.values()[0] + ori.values()[0])
        return rea.values()[0] >= ori.values()[0], fpval / tot_size

    # In this case, the codon is neither in the rea column nor in the
    # second column
    else:
        return False, 1


def chi2_is_possible(obs):
    """Check if we can perfom chi2 instead of fisher exact test"""
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
    """Construct alignment command line and execute it """
    prog = 'mafft'
    if 'muscle' in cmdline:
        prog = 'muscle'
        cmdline += " -in %s -out %s" % (inp, out)
    elif 'mafft' in cmdline:
        cmdline += " %s > %s" % (inp, out)

    executeCMD(cmdline, prog)


def compute_SP_per_col(al1, al2, columns, nspec,  scoring_matrix):
    """Compute a SP score per column"""
    def scoring_function(aa1, aa2, scoring_matrix):
        if scoring_matrix == 'identity':
            return (aa1 == aa2) * 1
        else:
            return scoring_matrix[aa1, aa2]
    al1_score = []
    al2_score = []
    for col in columns:
        al1_s = 0
        al2_s = 0
        for i in xrange(nspec - 1):
            for j in xrange(i + 1, nspec):
                al1_s += scoring_function(al1[i][col],
                                          al1[j][col], scoring_matrix)
                al2_s += scoring_function(al2[i][col],
                                          al2[j][col], scoring_matrix)
        al1_score.append(al1_s)
        al2_score.append(al2_s)
    # to obtain the SP score, just sum over each vector
    return al1_score, al2_score


def check_align_upgrade(al1, al2, scoring_method, method="wilcoxon", column_list=None):
    """Check if reassignment actually improve the alignment"""
    # al1 is the new alignement after changing the Gcode, and it's a dict
    # what we are testing is mean(al1_qual) >  mean(al2_qual)
    if scoring_method != "identity" and isinstance(scoring_method, basestring):
        scoring_method = getattr(MatrixInfo, scoring_method)
    if isinstance(al1, dict):
        al1 = al1.values()
    if isinstance(al2, dict):
        al2 = al2.values()
    al_len = len(al1[0])
    nspec = len(al1)
    if not column_list:
        column_list = range(al_len)
    assert (nspec == len(al2) and al_len == len(
        al2[0])), "The alignment are different"
    al1_sim, al2_sim = compute_SP_per_col(
        al1, al2, column_list, nspec, scoring_method)
    rank, pval = ReaGenomeFinder.get_paired_test(al1_sim, al2_sim, 1, method)
    return pval, al1_sim, al2_sim


def check_gain(codon, cible_aa, speclist, codontable, codon_alignment,
            scoring_method="identity", alignment=None,  method="wilcoxon"):
    """Check if there is an actuall gain in global sequence quality after applying reassignment"""
    if not isinstance(codon_alignment, dict):
        codon_alignment = SeqIO.to_dict(codon_alignment)

    def translate(codon_alignment, codontable, changes=None):
        if changes is None:
            changes = {}
        # return dict((spec, "".join(map(lambda x: codontable[x] \
        # if changes.get(spec, (None,))[0] != x else changes[spec][1],
        # [prot[i:i + 3] for i in xrange(0, len(prot), 3)]))) for spec, prot in codon_alignment.items())
        translated_al = {}
        position = set([])
        for spec, nucseq in codon_alignment.items():
            translated = ""
            nucseq_len = int(len(nucseq.seq)/3)
            for i in xrange(0, nucseq_len):
                cod = nucseq[i*3:(i+1)* 3].seq.tostring()
                if changes.get(spec, (None,))[0] != cod:
                    # use X for amino acid when stop codon is found
                    translated += codontable.get(cod, 'X')
                else:
                    # in this case, there were reassignment
                    translated += changes[spec][1]
                    position.add(i)
            translated_al[spec] = SeqRecord(Seq(translated, generic_protein), id=spec, name=spec)
        return translated_al, sorted(list(position))

    if not alignment:
        alignment, _ = translate(codon_alignment, codontable)

    changes = {}
    for spec in speclist:
        changes[spec] = (codon, cible_aa)

    cor_alignment, position = translate(codon_alignment, codontable, changes)
    score_improve, cor_al_sp, al_sp = check_align_upgrade(
        cor_alignment, alignment, scoring_method, method, position)

    align_info = AlignInfo.SummaryInfo(MultipleSeqAlignment(cor_alignment.values()))
    align_info.information_content()
    cor_ic_cont = align_info.ic_vector.values()
    align_info = AlignInfo.SummaryInfo(MultipleSeqAlignment(alignment.values()))
    align_info.information_content()
    ic_cont = align_info.ic_vector.values()

    return score_improve, (al_sp, cor_al_sp), (ic_cont, cor_ic_cont), (alignment, cor_alignment), position


def get_rea_genome_list(xlabel, y):
    """Return a dict containing each directory
    and if whether or it has a reassignation"""
    rea_genome = {}
    all_pos = xlabel[y==1, 0]
    for g in all_pos:
        rea_genome[g] = True
    return rea_genome


def get_codon_set_and_genome(y, xlabel, trueval=None):
    """Get for each codon the list of reassigned genome"""
    all_pos = xlabel[:, 0:2]
    if trueval:
        all_pos = xlabel[y==trueval, 0:2]
    result = defaultdict(list)
    for (g, cod) in all_pos:
        result[cod].append(g)
    return result


def get_report(fitchtree, gdata, reafinder, codon_align, prediction, output="", pie_size=45):
    """ Render tree to a pdf"""

    settings = reafinder.settings
    OUTDIR =  purge_directory(os.path.join(settings.OUTDIR, fitchtree.ori_aa +
                          "_to_" + fitchtree.dest_aa ))
    c_rea = get_rea_genome_list(prediction[1], prediction[3])
    if not output:
        output = os.path.join(OUTDIR, "Codon_data." + settings.IMAGE_FORMAT)
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

        default_style = NodeStyle()
        default_style["size"] = 0
        default_style["fgcolor"] = "black"
        show_n = settings.ADD_NUMBER_PIE

        def layout(node):

            get_suffix = lambda x: x if settings.ADD_LABEL_TO_LEAF else ""
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
                        AttrFace("name", text_suffix=get_suffix("_B")), node, 0, position="aligned")

                else:
                    # if lost in node.features, then node is already a leaf
                    # color change only concern leaf with count
                    if(c_rea.get(node.name, False)):
                        # Reassigned and passed fisher test
                        if 'lost' in node.features and not node.lost:
                            faces.add_face_to_node(AttrFace("name", fgcolor="#ff1111",
                                             text_suffix=get_suffix("_R")), node, column=0, position="aligned")
                        # Reassigned and failed fisher test
                        else:
                            faces.add_face_to_node(AttrFace("name", fgcolor="#fc8d59",
                                             text_suffix=get_suffix("_O")), node, column=0, position="aligned")

                    elif 'lost' in node.features and not node.lost:
                        faces.add_face_to_node(AttrFace("name", fgcolor="#1a9850",
                                        text_suffix=get_suffix("_G")), node, column=0, position="aligned")
                    else:
                        faces.add_face_to_node(
                            AttrFace("name", text_suffix=get_suffix("_B")), node, 0, position="aligned")



            if(fitchtree.has_codon_data() and node.is_leaf()):
                spec_codonrea_f = gdata[node.name]['filtered']['rea_codon']
                spec_codonused_f = gdata[node.name]['filtered']['used_codon']

                # add data
                faces.add_face_to_node(PPieChartFace(spec_codonrea_f.values(), pie_size, pie_size, show_label=show_n,
                                                     colors=[fitchtree.colors[k] for k in spec_codonrea_f.keys()]),
                                                     node, column=1, position="aligned")

                faces.add_face_to_node(PPieChartFace(spec_codonused_f.values(), pie_size, pie_size, show_label=show_n,
                                                     colors=[fitchtree.colors[k] for k in spec_codonused_f.keys()]),
                                                     node, column=2, position="aligned")

                next_column = 3

                if(settings.SHOW_MIXTE_CODONS):
                    spec_mixtecodon_f = gdata[node.name]['filtered']['mixte_codon']
                    faces.add_face_to_node(PPieChartFace(spec_mixtecodon_f.values(), pie_size, pie_size, show_label=show_n,
                                                         colors=[fitchtree.colors[k] for k in spec_mixtecodon_f.keys()]),
                                                         node, column=3, position="aligned")
                    next_column = 4

                if(settings.SHOW_GLOBAL_CODON_DATA):
                    spec_codonrea_g = gdata[node.name]['global']['rea_codon']
                    spec_codonused_g = gdata[node.name]['global']['used_codon']

                    # add separator
                    faces.add_face_to_node(LineFace(
                        pie_size, pie_size, None), node, column=next_column, position="aligned")

                    faces.add_face_to_node(PPieChartFace(spec_codonrea_g.values(), pie_size, pie_size, show_label=show_n,
                                    colors=[fitchtree.colors[k] for k in spec_codonrea_g.keys()]),
                                    node, column=next_column+1, position="aligned")

                    faces.add_face_to_node(PPieChartFace(spec_codonused_g.values(), pie_size, pie_size, show_label=show_n,
                                    colors=[fitchtree.colors[k] for k in spec_codonused_g.keys()]),
                                    node, column=next_column+2, position="aligned")

                    next_column += 3
                    if(settings.SHOW_MIXTE_CODONS):
                        spec_mixtecodon_g = gdata[node.name]['global']['mixte_codon']
                        faces.add_face_to_node(PPieChartFace(spec_mixtecodon_g.values(), pie_size, pie_size,  show_label=show_n,
                                               colors=[fitchtree.colors[k] for k in spec_mixtecodon_g.keys()]),
                                               node, column=next_column, position="aligned")
                        next_column += 1

                faces.add_face_to_node(LineFace(pie_size, pie_size, None), node, column=next_column, position="aligned")

                faces.add_face_to_node(TextFace("{:.2e}".format(gdata[node.name]['lost']['pval']), fsize=10, fgcolor="#000"),
                                           node, column=next_column+1, position="aligned")

        ts.layout_fn = layout

        # header declaration
        h1 = TextFace(fitchtree.dest_aa, fsize=10, fgcolor="#aa0000")
        h2 = TextFace(fitchtree.ori_aa, fsize=10, fgcolor="#aa0000")
        h3 = TextFace(fitchtree.ori_aa+" O.", fsize=10, fgcolor="#aa0000")
        h1g = TextFace(fitchtree.dest_aa, fsize=10, fgcolor="#aa0000")
        h2g = TextFace(fitchtree.ori_aa, fsize=10, fgcolor="#aa0000")
        h3g = TextFace(fitchtree.ori_aa+" O.", fsize=10, fgcolor="#aa0000")
        # center vertically and horizontally
        h1.vt_align, h2.vt_align, h3.vt_align = 1, 1, 1
        h1.hz_align, h2.hz_align, h3.hz_align = 1, 1, 2
        h1g.hz_align, h2g.hz_align, h3g.hz_align = 1, 1, 1
        h1g.hz_align, h2g.hz_align, h3g.hz_align = 1, 1, 2

        ts.aligned_header.add_face(h1, column=1)
        ts.aligned_header.add_face(h2, column=2)
        next_column = 3
        if(settings.SHOW_MIXTE_CODONS):
            ts.aligned_header.add_face(h3, column=3)
            next_column = 4

        if settings.SHOW_GLOBAL_CODON_DATA:
            ts.aligned_header.add_face(h1g, column=next_column+1)
            ts.aligned_header.add_face(h2g, column=next_column+2)
            next_column += 3
            if(settings.SHOW_MIXTE_CODONS):
                ts.aligned_header.add_face(h3g, column=next_column)
                next_column += 1

        f_pval_h = TextFace("FE. pval", fsize=10, fgcolor="#aa0000")
        f_pval_h.vt_align = 1
        f_pval_h.hz_align = 1
        ts.aligned_header.add_face(f_pval_h, column=next_column+1)

        ts.title.add_face(TextFace(fitchtree.ori_aa + " --> " + fitchtree.dest_aa, fsize=14), column=0)
        if(fitchtree.has_codon_data()):
            for cod, col in fitchtree.colors.items():
                ts.legend.add_face(CircleFace(
                    (pie_size / 3), col), column=0)
                ts.legend.add_face(
                    TextFace("  " + cod + " ", fsize=8), column=0)
            ts.legend_position = 4

        # Apply node style
        for n in fitchtree.tree.traverse():
            n.set_style(default_style)
            if n.reassigned == {1}:
                n.set_style(rea_style)
            elif len(n.reassigned) > 1:
                n.set_style(other_style)
            if 'lost' in n.features and n.lost:
                n.set_style(prob_lost)

        fitchtree.tree.render(output, dpi=800, tree_style=ts)
        glob_purge.append(output)

        # get report output
        rep_out = os.path.join(OUTDIR, "Report_"+fitchtree.ori_aa + "_to_" + fitchtree.dest_aa)
        # glob_purge.append(rep_out)
        pdf_format_data(fitchtree.ori_aa1, fitchtree.dest_aa1, gdata, prediction,
                                        'filtered', rep_out+".pdf")
        # now get sequence retranslation improvement
        table = fitchtree.dct.forward_table
        # add the gap to codon_table
        table['---'] = '-'
        table['...'] = '.'
        data_present, data_var = codon_adjust_improve(fitchtree, reafinder, codon_align, table, prediction, outdir=OUTDIR)
        return data_var, OUTDIR

def format_tree(tree, codon, cible, alignment, SP_score, ic_contents, pos=[],
                limits=(None, None, None), dtype="aa", codontable={}, codon_col={}):
    """Format the rendering of tree data for alignment"""
    t = tree.copy('newick')
    s = 0

    gpos, limiter, start_holder = limits

    def _get_genes(seq, jump=1):
        """ Return a dict of seq from the sequence"""
        gseq = {}
        ppos = 0
        if limiter is None :
            gseq['All'] = "".join([seq[i] for i in pos])
        else:
            for (l, name) in limiter:
                gseq[name] = "".join([seq[i*jump:jump*(i+1)] for i in pos[ppos:l]])
                ppos = l
        return gseq

    for node in t:
        node_seq = alignment[node.name].seq.tostring()
        if pos and dtype=="aa":
            node_seq = _get_genes(node_seq)
            s = len(node_seq)
        elif pos and dtype=="codon":
            node_seq = _get_genes(node_seq, 3)
        node.add_feature('sequence', node_seq)

    ts = TreeStyle()
    ts.branch_vertical_margin = 15
    ts.scale = 15
    ts.allow_face_overlap = False
    ts.show_scale = False
    ts.show_leaf_name = False

    if limiter is None:
        limiter = [(len(pos), 'All')]

    ns = NodeStyle()
    ns['shape'] = 'square'
    ns['fgcolor'] = 'black'
    ns['size'] = 0
    ns['node_bgcolor'] = 'black'

    def layout(node):
        node.img_style = ns
        if node.is_leaf():
            faces.add_face_to_node(AttrFace('name', fsize=14), node, 0, position="aligned")
            if hasattr (node, "sequence"):
                ind = 0
                for (k, name) in limiter:
                    seq = node.sequence[name]
                    seqface = faces.SequenceFace(seq, fsize=13)
                    if dtype != 'aa':
                        seqface =  SequenceFace(seq, cible, seqtype=dtype, fsize=13,
                                            codontable=codontable, spec_codon_col=codon_col)
                    faces.add_face_to_node(seqface, node, 2+ind, aligned=True)
                    ind += 2

    ts.layout_fn = layout

    if dtype == 'aa':
        gfunc = lambda x : '*' if x==1 else ' '
        if pos:

            ts.title.add_face(TextFace('(%s) - SP score : %.0f | IC = %.2f'%(codon, sum(SP_score), sum(ic_contents)),
                                                            fsize=14, fgcolor='red'), 0)
            ts.aligned_header.add_face(faces.RectFace(14, 14, 'white', 'white'), 1)

            ts.aligned_foot.add_face(faces.RectFace(14, 14, 'white', 'white'), 1)
            start = 0
            ind = 3
            for (l, name) in limiter:
                ic_content = np.asarray(ic_contents)[pos[start:l]]
                #sp_score = np.asarray(SP_score)[pos]
                footer_seq = SequenceFace("".join([gfunc(st) for st in start_holder[start:l]]), None, dtype, fsize=13)
                start = l
                ic_plot = faces.SequencePlotFace(ic_content, hlines=[(int(min(ic_content)-0.5))],
                                         hlines_col=['white'], ylim=(int(min(ic_content)-0.5),
                                         int(max(ic_content)+1)), fsize=10, col_width=14,
                                         header="IC", kind='bar')
                ts.aligned_header.add_face(ic_plot, ind)
                ts.aligned_foot.add_face(footer_seq, ind)
                ts.aligned_foot.add_face(TextFace(name), ind)
                ind += 2

    else:
        for (cod, col) in codon_col.items():
            ts.legend.add_face(faces.RectFace(50, 25, col, col), column=0)
            ts.legend.add_face(TextFace("  %s "%cod, fsize=8), column=1)
        ts.legend.add_face(faces.RectFace(50, 25, 'black', 'black'), column=0)
        ts.legend.add_face(TextFace("  %s\'s codons "%(cible), fsize=8), column=1)
        ts.legend_position = 3
        ind = 1
        lprev = 0
        for (l, name) in limiter:
            ts.aligned_foot.add_face(faces.RectFace(14*(l-lprev)*3, 14, '#fefefe', 'white'), ind+1)
            ts.aligned_foot.add_face(TextFace(' ', fsize=14), ind)
            ts.aligned_foot.add_face(TextFace(name, fsize=14), ind+1)
            lprev = l
            ind +=2

    t.dist = 0
    return t, ts


def identify_position_with_codon(fcodal, codon, spec_to_check):
    """Get all positions where a codon is used"""
    positions = []
    try:
        al_len = int(len(fcodal.values()[0])/3)
    except:
        al_len = fcodal.get_aln_length()
    for i in xrange(al_len):
        for spec in spec_to_check:
            if fcodal[spec].seq.get_codon(i) == codon:
                positions.append(i)
                break
    return positions


def violin_plot(vals, output, score, codon, cible, imformat="pdf"):
    """Return violin plot of SP score """

    figure = plt.figure()
    axes = figure.add_subplot(1,1,1)
    data = vals.values()
    pos = [y+1 for y in range(len(data))]
    fig = plt.violinplot(data, pos, showmedians=True, showextrema=True)
    for bod in fig['bodies']:
        bod.set_facecolor('lightblue')
        bod.set_edgecolor('darkblue')
        bod.set_linewidth(2)
    fig['cbars'].set_color('#555555')
    fig['cmaxes'].set_color('#555555')
    fig['cmedians'].set_color('#555555')
    fig['cmins'].set_color('#555555')
    axes.yaxis.grid(True)
    axes.set_xticks(pos)
    plt.setp(axes, xticks=pos, xticklabels=vals.keys())
    output = output+"."+imformat
    axes.set_title('p-values({} --> {}) : {:.2e}'.format(codon, cible, score))
    figure.savefig(output, format=imformat)
    logging.info('{} --> {}) : {:.2e}'.format(codon, cible, score))
    return output


def codon_adjust_improve(fitchtree, reafinder, codon_align, codontable, prediction, outdir=""):
    """Get a representation of the improvement after translation"""

    cible_aa = fitchtree.dest_aa1
    X_data, X_labels, pred_prob, pred = prediction
    true_codon_set = get_codon_set_and_genome(pred, X_labels, 1)
    settings = reafinder.settings
    genelimit =  reafinder.seqset.gene_limits
    filt_position = reafinder.seqset.filt_position
    sc_meth = settings.method
    method = settings.mode if settings.mode in ['wilcoxon', 'mannwhitney', 'ttest'] else 'wilcoxon'
    outputs = []
    violinout = None
    data_var = {}
    ori_al = None
    tree = fitchtree.tree

    def _limit_finder(codon_pos, gene_limit=genelimit, proche=settings.startdist):
        """ Return the gene for each position of codon pos"""
        gene_pos = []
        gene_break = []
        gll = len(gene_limit)
        entering, start = 0, 0
        brokn = False
        close_to_start = []
        for i in range(len(codon_pos)):
            while codon_pos[i] > gene_limit[start][2] and start < gll:
                start += 1
                brokn = True
            if brokn and i>0:
                gene_break.append((i, gene_limit[entering][0]))
                entering = start
            if codon_pos[i] - gene_limit[start][1] <= proche:
                close_to_start.append(1)
            else:
                close_to_start.append(0)
            brokn = False
            gene_pos.append(gene_limit[start][0])

        gene_break.append((len(codon_pos), gene_limit[start][0]))
        return gene_pos, gene_break, close_to_start

    for codon in true_codon_set.keys():
        speclist = true_codon_set[codon]
        if len(speclist) > 0:
            #pos = identify_position_with_codon(codon_align, codon, speclist)
            score_improve, alsp, alic, als, pos = check_gain(codon, cible_aa, speclist, codontable,
                                                    codon_align, scoring_method=sc_meth,
                                                    alignment=ori_al, method=method)
            codon_pos = [filt_position[i] for i in pos]
            limits =  _limit_finder(codon_pos)

            sp, cor_sp = alsp
            ic, cor_ic = alic
            ori_al, new_al = als
            ori_t, ori_ts = format_tree(tree, codon, cible_aa, ori_al, sp, ic, pos, limits)
            rea_t, rea_ts = format_tree(tree, codon, cible_aa, new_al, cor_sp, cor_ic, pos, limits)
            cod_t, cod_ts = format_tree(tree, codon, cible_aa, codon_align, None, None, pos, limits,
                                            codontable=codontable, dtype="codon", codon_col=fitchtree.colors)

            cod_ts.title.add_face(TextFace("Prediction validation for "+codon+" to "+fitchtree.dest_aa, fsize=14), column=0)

            violinout = violin_plot({'Original':sp, 'Corrected':cor_sp},
                                    os.path.join(outdir, "%s_violin"%codon),
                                    score_improve, codon, fitchtree.dest_aa, imformat=settings.IMAGE_FORMAT)

            cod_out = os.path.join(outdir, "%s_codons.%s"%(codon, settings.IMAGE_FORMAT))
            ori_out = os.path.join(outdir, "%s_ori.%s"%(codon, settings.IMAGE_FORMAT))
            rea_out = os.path.join(outdir, "%s_rea.%s"%(codon, settings.IMAGE_FORMAT))

            ori_t.render(ori_out, dpi=800, tree_style=ori_ts)
            rea_t.render(rea_out, dpi=800, tree_style=rea_ts)
            cod_t.render(cod_out, dpi=800, tree_style=cod_ts)

            data_var[codon] = score_improve

            glob_purge.append(cod_out)
            glob_purge.append(rea_out)
            glob_purge.append(ori_out)
            glob_purge.append(violinout)
            outputs.append(True)

    return len(outputs)>0, data_var


def pdf_format_data(ori_aa, dest_aa, gdata, prediction, dtype, output=""):
    """Render supplemental data into a pdf"""
    X_data, X_labels, pred_prob, pred = prediction
    codon_set = get_codon_set_and_genome(pred, X_labels)
    pd_data = defaultdict(dict)
    pred_dict = defaultdict(dict)
    val_list = X_labels.shape[0]
    for i in range(val_list):
        # This is in order to access data quicker later
        # Xlabels format : [genome, codon, aa2, aa1]
        gm = X_labels[i, 0]
        cd = X_labels[i, 1]
        pred_dict[gm][cd] = i

    for codon in codon_set.keys():
        genome_list =  codon_set[codon]
        for g in genome_list:
            pos = pred_dict[g][codon]
            dt_by_att = {
                            'Rea. count': gdata[g][dtype]['rea_codon'].get(codon,0),
                            'Use. count': gdata[g][dtype]['used_codon'].get(codon,0),
                            'Tot. count': gdata[g]['codons'][dtype].get(codon,0),
                            'Rea. Genes': len(gdata[g][dtype]['rea_distribution'].get(codon,[])),
                            'Tot. Genes': len(gdata[g][dtype]['total_rea_distribution'].get(codon,[])),
                            'Telford': gdata[g]['score'][dtype].get(codon,0),
                            'Reassigned': pred[pos] == 1, #position 1 is true prob
                            #'Delta Sim.' :
                            'Prob': pred_prob[pos, 1],  #position 1 is true prob
                         }

            pd_data[g][codon] = dt_by_att

    frames = []
    glist = []
    for (g, dt) in pd_data.items():
        glist.append(g)
        frames.append(pd.DataFrame.from_dict(dt, orient='index'))
    fdata = pd.concat(frames, keys=glist)
    html_data = fdata.to_html()
    data_var = {'title':'Codon reassignment prediction in all genome',
                'report_table':html_data,
                }
    return export_from_html(data_var, output)


def print_data_to_txt(outputfile, header, X, X_label, Y, codon_data, cible):
    """Print result in a text file"""
    out = Output(file=outputfile)
    out.write("### Random Forest prediction\n")
    out.write("\t".join(["genome", "codon",
                            "ori_aa", "rea_aa"] + header + ["prediction"]))

    total_elm = len(Y)
    for i in xrange(total_elm):
        out.write("\n"+"\t".join(list(X_label[i]) + [str(x) for x in X[i]] + [str(Y[i])]))

    out.write("\n\n### Alignment improvement:\n")
    for (codon,score) in codon_data.items():
        out.write("{}\t{}\t{:.2e}\n".format(codon, cible, score))
    out.close()


def makehash():
    """Utility method to make a multilevel dict"""
    return defaultdict(makehash)
