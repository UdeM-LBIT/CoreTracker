from __future__ import division

import glob
import itertools
import json
import logging
import os
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from distutils import spawn
from functools import partial

import Bio.SubsMat.MatrixInfo as MatrixInfo
import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt

from Bio import AlignIO, Alphabet, SeqIO, SubsMat, codonalign
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet import IUPAC, generic_nucleotide, generic_protein
from Bio.codonalign.codonalphabet import (default_codon_alphabet,
                                          get_codon_alphabet)
from Bio.codonalign.codonseq import CodonSeq
from Bio.Data import CodonTable
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from ete3 import Tree
from scipy.cluster.vq import kmeans2

from AncestralRecon import SingleNaiveRec, init_back_table
from corefile import CoreFile
from coretracker.FisherExact import fisher_exact
from Faces import LineFace, List90Face, PPieChartFace, SequenceFace
from letterconfig import *
from output import Output
from pdfutils import *

SEABORN = False
try:
    import seaborn as sns
    sns.set(context="paper", palette="deep")
    sns.set_style("whitegrid")
    SEABORN = True
except ImportError as e:
    # quiet fallback to matplotlib
    vplot = plt.violinplot

# conditional import
GRAPHICAL_ACCESS = True
try:
    from ete3 import TreeStyle, NodeStyle, faces, AttrFace, TextFace, CircleFace
except ImportError, e:
    GRAPHICAL_ACCESS = False

# hmm strange id
hmmidpattern = re.compile("^\d+\|\w+")
# define lowest pvalue
SYSMINVAL = np.finfo(np.float).tiny

alpha = Alphabet.Gapped(IUPAC.protein)


class ConsensusKeeper:

    def __init__(self, alignment, threshold, ambiguous='X', alphabet=generic_protein):
        self.alignment = alignment
        self.threshold = threshold
        self.ambiguous = ambiguous
        self.allen = alignment.get_alignment_length()
        self.number_seq = len(alignment)
        self.aligninfo = AlignInfo.SummaryInfo(alignment)
        self.consensus = self.aligninfo.gap_consensus(
            threshold=threshold, ambiguous=ambiguous, consensus_alpha=alphabet)
        self.pssm = self.aligninfo.pos_specific_score_matrix(self.consensus)

    def get_cons(self):
        return self.consensus

    def get_pssm(self):
        return self.pssm

    def get_alt_cons(self):
        rep_per_pos = defaultdict(list)
        for k in range(self.allen):
            added = False
            for aa, aacount in self.pssm[k].items():
                aa_freq_count = (aacount * 1.0) / self.number_seq
                if aa_freq_count >= self.threshold:
                    rep_per_pos[k].append((aa, aa_freq_count))
                    added = True
            if not added:
                rep_per_pos[k].append((self.ambiguous, 0))
        return rep_per_pos

    def get_maj_aa_pos(self, aa):
        rep_per_pos = self.get_alt_cons()
        filt_pos = []
        for pos in range(self.allen):
            aalist = rep_per_pos[pos]
            aalist = sorted(aalist, key=lambda x: x[1])
            ind = len(aalist) - 1
            sup_aa, score = aalist[ind]
            sc = score
            while ind >= 0 and sc == score:
                sup_aa, sc = aalist[ind]
                if sup_aa == aa:
                    filt_pos.append(pos)
                    break
                ind -= 1

        return filt_pos


class SequenceLoader:
    """Load and format data for analysis"""

    def __init__(self, infile, dnafile, settings, gap_thresh, use_tree=None, refine_alignment=True, msaprog="muscle", hmmdict={}, clean_stop=False):
        prog = msaprog if msaprog else "muscle"
        self.infile = infile
        self.sequences = {}
        self.dnasequences = {}
        self.hmmdict = hmmdict
        self.hmmloop = settings.HMMLOOP
        # try parsing the protein file as a corefile
        self.sequences = self.get_sequences(
            infile, settings.PROTFORMAT, generic_protein)
        # we got to the end without finding any sequences
        if not self.sequences:
            raise ValueError(
                'Unable to load protein sequences, Unknown format')
        logging.debug('Protein sequence was read')

        # here we are parsing the dnafile :
        # dnafile should be in corefile too
        self.dnasequences = self.get_sequences(
            dnafile, settings.DNAFORMAT, generic_nucleotide)
        if not self.dnasequences:
            raise ValueError(
                'Unable to load nucleotide sequences, Unknown format')
        logging.debug('DNA sequence was read')
        self.genes = set(self.dnasequences.keys()).intersection(
            self.sequences.keys())

        if not self.genes:
            raise ValueError(
                "Gene not found ! None of the genes in %s are found in %s" % (dnafile, infile))
        common_spec_per_gene = {}
        common_spec_per_gene_len = {}  # len of a gene in each genome
        max_spec = -np.inf
        # Only interested in gene in both nucleotide and aa sequences
        for gene in self.genes:
            g_dnaspec = [x.id for x in self.dnasequences[gene]]
            g_protspec = [x.id for x in self.sequences[gene]]
            g_common_spec = set(g_dnaspec).intersection(set(g_protspec))
            common_spec_per_gene[gene] = g_common_spec
            common_spec_per_gene_len[gene] = len(g_common_spec)
            max_spec = max(max_spec, common_spec_per_gene_len[gene])

        # dropping gene with present in not enough species
        # using gap_thresh to perform this
        if gap_thresh:
            gap_thresh = (abs(gap_thresh) <= 1 or 0.01) * abs(gap_thresh)
            self.genes = [x for x in self.genes if common_spec_per_gene_len[
                x] >= max_spec * (1 - gap_thresh)]
        self.common_spec_per_gene = dict(
            (g, common_spec_per_gene[g]) for g in self.genes)
        if not self.genes:
            raise ValueError(
                "Not enough conserved gene for the gap threshold provided.")
        logging.debug('List of selected genes :\n\t- %s ' %
                      "\n \t- ".join(self.genes))
        self.true_spec_list = set().union(*self.common_spec_per_gene.values())
        is_aligned, alignment = self._prot_is_aligned()

        # align sequences if it's not already done
        if not is_aligned or msaprog:
            logging.debug('Sequences alignment to be performed')
            alignment = self.align(
                settings, refine=refine_alignment, tree=use_tree, msaprog=prog, scale=settings.SCALE, is_aligned=is_aligned)

        # case where sequence is not aligned but refinement is required
        elif is_aligned and refine_alignment:
            logging.debug(
                'Refinement of alignment to be performed')
            for g in self.genes:
                alignment[g] = self.__class__._refine(alignment[g], 9999, settings.OUTDIR,
                                                      settings.hmmbuild, settings.hmmalign,
                                                      loop=self.hmmloop, hmmfile=self.hmmdict.get(g, None))
        self.alignment = alignment
        logging.debug('Sequence alignment done')

    def get_sequences(self, infile, fileformat, alphabet):
        """Get sequence from file"""
        if fileformat == "core":
            corefile = CoreFile(infile, alphabet)
            return corefile.get_sequences()
        else:
            return {'0': list(SeqIO.parse(infile, fileformat, alphabet))}

    def concat(self, missing='-', alpha=generic_protein, stopchar='*'):
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
                    # check if there is a stop at the end of the dna
                    ungap_prot = adding.seq.ungap(missing)
                    if ungap_prot[-1] != stopchar and len(ungap_prot) * 3 + 3 == len(dna_adding.seq.ungap(missing)):
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

    def align(self, settings, msaprog, refine=True, tree=None, scale=1.0, alpha=generic_protein, is_aligned=False):
        """Align sequences"""
        alignment = {}
        outdir = settings.OUTDIR
        for gene in self.genes:
            # in order to have the best possible alignment, we are going to
            # keep the
            if len(self.sequences[gene]) > 1:
                al = self.__class__._align(
                    self.sequences[gene], msaprog, tree, scale, outdir, alpha, is_aligned)
                if refine:
                    al = self.__class__._refine(
                        al, 9999, outdir, settings.hmmbuild, settings.hmmalign,
                        loop=self.hmmloop, hmmfile=self.hmmdict.get(gene, None))
                alignment[gene] = al
            else:
                logging.debug('%s dropped because only %s has it' %
                              (gene, self.sequences[gene][0].id))

        return alignment

    @classmethod
    def _align(clc, msa, msaprog, tree, scale, outdir, alpha=generic_protein, is_aligned=False, gap_char='-'):
        """Align sequences using muscle of mafft"""
        tmpseq = os.path.join(outdir, "tmp_seq.fasta")
        align_seq = os.path.join(outdir, "tmp_msa.fasta")
        seqs = msa
        if is_aligned:
            for seqrec in msa:
                seqrec.seq = seqrec.seq.ungap(gap_char)
        SeqIO.write(seqs, open(tmpseq, 'w'), 'fasta')

        if tree and 'mafft' in msaprog:
            seq_order = [seqrec.id for seqrec in seqs]
            out = Output(file=os.path.join(outdir, "tmp_tree.mafft"))
            tmpt = Tree(tree)
            tmpt.prune(seq_order, preserve_branch_length=False)
            # tmpt.resolve_polytomy(recursive=True)
            is_multifurcated = False
            for node in tmpt.traverse():
                if len(node.children) > 2:
                    is_multifurcated = True
                    break
            if is_multifurcated:
                raise ValueError(
                    "Input tree should be rooted and binary for mafft to work.")
            convert_tree_to_mafft(tmpt, seq_order, out, scale)
            out.close()
            msaprog += " --treein %s" % out.file

        execute_alignment(msaprog, tmpseq, align_seq)
        msa = AlignIO.read(align_seq, 'fasta', alphabet=alpha)
        filelist = glob.glob(os.path.join(outdir, "tmp_*"))
        for f in filelist:
            os.remove(f)
        return msa

    @classmethod
    def _refine(clc, alignment, timeout, outdir, hmmbuild="hmmbuild", hmmalign="hmmalign",
                hmmfile=None, loop=10, minqual=8, strategie="lesser", clean=True):
        """Align and refine at the same time with hmmalign and muscle"""

        success, not_found = check_binaries(hmmbuild, hmmalign)
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
        ungapedInputFile = os.path.join(outdir, 'tmp_seq.fasta')
        if isinstance(alignment, MultipleSeqAlignment):
            file_created = True
            AlignIO.write(alignment, open(inputFile, 'w'), 'fasta')
        else:
            inputFile = alignment

        ungapedseqrec = list(SeqIO.parse(inputFile, format="fasta"))
        for srec in ungapedseqrec:
            srec.seq._data = srec.seq._data.replace('-', "").replace('.', "")
        SeqIO.write(ungapedseqrec, ungapedInputFile, "fasta")

        logging.debug(
            '... trying to run hmmbuild and hmmalign ' + str(loop) + " times!")
        quality = []
        outlist = []
        should_break = False
        for i in xrange(loop):
            outputFile = os.path.join(cur_dir[i], "alignment.sto")
            tmphmmfile = os.path.join(cur_dir[i], "alignment.hmm")
            if hmmfile:
                # in this case copy it to the current dir
                shutil.copyfile(hmmfile, tmphmmfile)
            else:
                buildline = hmmbuild + \
                    " --amino %s %s" % (tmphmmfile, inputFile)
                executeCMD(buildline, 'hmmbuild')
            # will continue if not exception is found
            alignline = hmmalign + \
                " -o %s %s %s" % (outputFile, tmphmmfile, ungapedInputFile)
            executeCMD(alignline, 'hmmalign')
            # finding quality
            try:
                quality.append(accessQuality(outputFile, minqual, strategie))
            except IOError:
                raise IOError(
                    "File '%s' not found. Either alignment failed, or you are using a wrong hmm version with HMMER" % outputFile)

            inputFile = remove_gap_only_columns(outputFile, 'stockholm')
            if i > 0 and improve_is_stagned(outlist[-1], inputFile, len(outlist) * 1.0 / loop):
                should_break = True
                logging.debug(
                    "Stopping hmm loop : no alignment improvement after %d/%d iteration" % (len(outlist), loop))
            outlist.append(inputFile)
            if should_break:
                break

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
                os.remove(ungapedInputFile)
        return alignment

    @classmethod
    def clean_stop(clc, protseqs, dnaseqs, stop='*'):
        """Clean Stop from the sequence"""
        # this remove all the stop codon at the ends
        for gene, pseq in protseqs.items():
            seq_len_per_spec = {}
            for prot in pseq:
                if prot.seq.endswith(stop):
                    prot.seq = prot.seq[0:-1]
                seq_len_per_spec[prot.id] = len(prot.seq)

            for dseq in dnaseqs[gene]:
                if seq_len_per_spec.get(dseq.id, 0) * 3 == len(dseq.seq) - 3:
                    dseq.seq = dseq.seq[0:-3]
        return protseqs, dnaseqs

    def _prot_is_aligned(self):
        """Check if protein sequence are aligned"""
        alignment = {}
        try:
            for gene in self.genes:
                alignment[gene] = MultipleSeqAlignment(self.sequences[gene])
        except ValueError:
            return False, None

        return True, alignment

    @classmethod
    def translate(clc, core_inst, gcode=1, codemap={}):
        """Translate nucleotide sequences into protein sequence given a genetic code"""
        codontable = CodonTable.unambiguous_dna_by_id[1]
        try:
            codontable = CodonTable.unambiguous_dna_by_id[abs(gcode)]
        except:
            logging.warn("Wrong genetic code, resetting it to 1")

        if isinstance(core_inst, basestring):
            core_inst = CoreFile(core_inst, alphabet=generic_nucleotide)

        core_prot_inst = {}

        def stop_checker(x):
            return '*' if x[0].upper() in x[1].stop_codons else 'X'

        for gene, seqs in core_inst.items():
            translated = []
            for s in seqs:
                seq_code = codontable
                try:
                    has_code = codemap.get(s.id, None)
                    if has_code:
                        seq_code = CodonTable.unambiguous_dna_by_id[
                            abs(has_code)]
                except KeyError:
                    seq_code = codontable
                if len(s) % 3 != 0:
                    raise ValueError(
                        "Frame-shifting detected in %s : [%s], current version does not supported it." % (s.id, gene))
                else:
                    aa_list = [seq_code.forward_table.get(
                        s.seq[i:i + 3].upper(), stop_checker((s.seq[i:i + 3], seq_code))) for i in xrange(0, len(s), 3)]
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
        self.consensus = consensus.get_alt_cons()
        self.reacodons = makehash(1, Counter)
        self.usedcodons = makehash(1, Counter)
        self.mixtecodon = makehash(1, Counter)
        # find codon usage in the sequence
        self.subsmat = settings.SUBMAT
        self.cons_for_lik = settings.USE_CONSENSUS_FOR_LIKELIHOOD

        self.specs_amino_count = makehash(1, Counter)
        self.positions = positions
        self.genelimit = genelimit
        self.t_genelimit = len(genelimit)
        self.codons_distribution = makehash(2, set)
        self.total_codons_distribution = makehash(2, set)
        self.codon_map_in_genome = makehash(1, Counter)

        self._spec_codon_usage()

    def _spec_codon_usage(self):
        """Check codon usage in each species"""
        # Check this function to enable stop codon
        # reassignment later
        for i, aa in self.consensus.items():
            for spec in self.codon_alignment.keys():
                spec_codon = self.codon_alignment[
                    spec].seq._data[i * 3:(i + 1) * 3]
                spec_aa = self.dct.forward_table.get(spec_codon, None)
                cur_pos = self.positions[i]
                if isinstance(aa, list):
                    aa = sorted(aa, key=lambda x: x[-1])  # .pop()
                    aa_b, aa_s = aa.pop()
                    aa = aa_b + "".join([x[0] for x in aa if x[1] == aa_s])
                if(spec_aa):
                    if spec_aa in aa:
                        aa = spec_aa
                    else:
                        aa = aa[0]
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


class SequenceSet(object):
    """Representation of the current genome set and their sequences"""

    def __init__(self, seqloader, phylotree, table_num, codemap={}):
        # using sequence set suppose that the alignment and the protein
        # concatenation was already done
        self.codontable = CodonTable.unambiguous_dna_by_id[abs(table_num)]
        self.codemap = codemap
        if codemap:
            logging.debug("Species genetic code passed as input: ")
            logging.debug(
                "\n" + "\n".join(["\t%s\t%s" % (k, v) for k, v in codemap.items()]))
        self.prot_dict, self.dna_dict, self.gene_limits = seqloader.concat()
        self.prot_align = MultipleSeqAlignment(self.prot_dict.values())
        self.common_spec_per_gene = seqloader.common_spec_per_gene
        self.phylotree = phylotree
        self.common_genome = []
        self.is_filtered = False
        self.restrict_to_common()
        # self.compute_current_mat()
        self.codon_align()

    def fusion(self, *seqsetlist, **kwargs):
        # we assume that the same tree is used everywhere
        # not testing
        missing = kwargs.get('missing', '-')
        treenw = kwargs.get('tree', None)
        curlen = self.prot_align.get_alignment_length()
        filt_len = len(self.filt_position)
        self.filt_prot_align = SeqIO.to_dict(self.filt_prot_align)
        self.codon_alignment = SeqIO.to_dict(self.codon_alignment)
        self.fcodon_alignment = SeqIO.to_dict(self.fcodon_alignment)
        cod_alpha = self.codon_alignment.values()[0].seq.alphabet
        for o_seqset in seqsetlist:
            o_len = o_seqset.prot_align.get_alignment_length()
            o_filt_len = len(o_seqset.filt_position)
            self.common_genome = set.union(
                self.common_genome, o_seqset.common_genome)
            o_dna = SeqIO.to_dict(o_seqset.codon_alignment)
            o_filtprot = SeqIO.to_dict(o_seqset.filt_prot_align)
            o_filtdna = SeqIO.to_dict(o_seqset.fcodon_alignment)

            for spec in self.common_genome:
                # protdict
                cur_seq = self.prot_dict.get(spec, SeqRecord(
                    Seq(missing * curlen, generic_protein), id=spec, name=spec))
                o_seq = o_seqset.prot_dict.get(spec, SeqRecord(
                    Seq(missing * o_len, generic_protein), id=spec, name=spec))
                gseq = cur_seq.seq + o_seq.seq
                self.prot_dict[spec] = SeqRecord(gseq, id=spec, name=spec)

                # filt alignment
                cur_seq = self.filt_prot_align.get(spec, SeqRecord(
                    Seq(missing * filt_len, generic_protein), id=spec, name=spec))
                o_seq = o_filtprot.get(spec, SeqRecord(
                    Seq(missing * o_filt_len, generic_protein), id=spec, name=spec))
                fseq = cur_seq.seq + o_seq.seq
                self.filt_prot_align[spec] = SeqRecord(
                    fseq, id=spec, name=spec)

                # codon alignment
                cur_seq = self.codon_alignment.get(spec, SeqRecord(
                    CodonSeq(missing * curlen * 3, alphabet=cod_alpha), id=spec, name=spec))
                o_seq = o_dna.get(spec, SeqRecord(
                    CodonSeq(missing * o_len * 3, alphabet=cod_alpha), id=spec, name=spec))
                codseq = CodonSeq(cur_seq.seq._data +
                                  o_seq.seq._data, alphabet=cod_alpha)
                self.codon_alignment[spec] = SeqRecord(
                    codseq, id=spec, name=spec)

                # filt codon alignment
                cur_seq = self.fcodon_alignment.get(spec, SeqRecord(
                    CodonSeq(missing * filt_len * 3, alphabet=cod_alpha), id=spec, name=spec))
                o_seq = o_filtdna.get(spec, SeqRecord(
                    CodonSeq(missing * o_filt_len * 3, alphabet=cod_alpha), id=spec, name=spec))
                codseq = CodonSeq(cur_seq.seq._data +
                                  o_seq.seq._data, alphabet=cod_alpha)
                self.fcodon_alignment[spec] = SeqRecord(
                    codseq, id=spec, name=spec)
                # we can ignore dna_dict
                # self.dna_dict.update(o_seqset.dna_dict)

            self.position = np.append(
                self.position, (o_seqset.position + curlen))
            self.filt_position = np.append(
                self.filt_position, (o_seqset.filt_position + curlen))
            lastgene_pos = self.gene_limits[-1][-1]
            self.gene_limits.extend(
                [(x[0], x[1] + lastgene_pos, x[2] + lastgene_pos) for x in o_seqset.gene_limits])
            curlen += o_len
            del o_filtprot
            del o_filtdna
            del o_dna

        self.prot_align = MultipleSeqAlignment(
            self.prot_dict.values(), alphabet=generic_protein)
        self.filt_prot_align = MultipleSeqAlignment(
            self.filt_prot_align.values(), alphabet=generic_protein)
        self.codon_alignment = codonalign.CodonAlignment(
            self.codon_alignment.values(), alphabet=cod_alpha)
        self.fcodon_alignment = codonalign.CodonAlignment(
            self.fcodon_alignment.values(), alphabet=cod_alpha)
        # patch tree
        # because it could have been pruned
        if treenw:
            self.phylotree = Tree(treenw)
            self.phylotree.prune(self.common_genome)

        return self

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
                [1 for gene, speclist in self.common_spec_per_gene.items() if spec in speclist])
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

        # print [(x, len(self.dna_dict[x]),
        # len(self.prot_dict[x].seq.ungap('-'))*3,
        # len(self.dna_dict[x])==(len(self.prot_dict[x].seq.ungap('-'))+1)*3)
        if not c_genome:
            help1 = "1) Protein length do not match with dna and stop already checked! Look for mistranslation of Frame-shifting"
            help2 = "2) Wrong tree or sequences name not matching!"
            # print [x +" "+ str(len(self.dna_dict[x].seq.ungap('-'))) +"
            # "+str(len(self.prot_dict[x].seq.ungap('-')) * 3) for x in
            # common_genome]
            raise ValueError(
                'ID intersection for dna, prot and species tree is empty. Possible cause:\n%s\n%s\n' % (help1, help2))
        common_genome = c_genome
        if len(dna_set) != speclen or len(species_set) != speclen:
            logging.debug(
                'Non-uniform sequence in sequences and tree. The tree will be pruned.')

        # prune tree to sequence list
        self.phylotree.prune(common_genome)
        self.common_genome = common_genome
        logging.debug("List of common genome :")
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
            seq_code = self.codontable
            try:
                has_code = self.codemap.get(g, None)
                if has_code:
                    seq_code = CodonTable.unambiguous_dna_by_id[abs(has_code)]
            except KeyError:
                seq_code = self.codontable

            codon_rec, undef_c = self._get_codon_record(
                self.dna_dict[g], self.prot_dict[g], seq_code, alphabet)
            all_undef_codon[g] = undef_c
            codon_aln.append(codon_rec)

        self.codon_alignment = codonalign.CodonAlignment(
            codon_aln, alphabet=alphabet)

        if all_undef_codon:
            undef_codon = set().union(*all_undef_codon.values())
            logging.debug("Total position lost due to undef codon : %d" %
                          len(undef_codon))
            undef_codon = sorted(list(undef_codon))
            # R_aa_position = np.asarray(xrange(self.prot_align.get_alignment_length()))
            # R_aa_position = np.setxor1d(tt_filter_position, np.asarray(undef_codon))
            for gene, seqrec in self.prot_dict.items():
                seqrec.seq._data = seqrec.seq._data.replace('X', gap_char)
                # self.prot_dict[gene] = seqrec
            self.prot_align = MultipleSeqAlignment(
                self.prot_dict.values(), alphabet=alpha)

            # remove all the position with undef codon from the dna_dict
            for codseqrec in self.codon_alignment:
                k = codseqrec.id
                self.dna_dict[k] = SeqRecord(
                    codseqrec.seq.toSeq().ungap(gap_char), id=k, name=k)

        # self.codon_alignment = codonalign.build(self.prot_align, self.dna_dict, codon_table=self.codontable)

    def _get_codon_record(self, dnarec, protrec, codontable, alphabet, gap_char='-'):
        """Get a codon seq record from dna and prot seqrecord"""
        nuc_seq = dnarec.seq.ungap(gap_char)
        codon_seq = ""
        aa_num = 0
        max_error = 100
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
        return SeqRecord(CodonSeq(codon_seq, alphabet, enable_undef=True), id=dnarec.id, name=dnarec.name), x_undecoded

    def filter_sequences(self, id_thresh, gap_thresh, ic_thresh, rmcnst, nofilter=False):
        """Filter prot/codons alignment"""
        self.filter_prot_alignment(id_thresh, gap_thresh, ic_thresh, rmcnst)
        if nofilter:
            self.fcodon_alignment = self.codon_alignment
        else:
            r = self.filter_codon_alignment()
            self.fcodon_alignment = next(r)
        self.is_filtered = True

    def filter_prot_alignment(self, id_thresh=None, gap_thresh=None, ic_thresh=None, rmcnst=True):
        """Filter protein alignment"""
        current_alignment = self.prot_align
        tt_filter_position = np.asarray(
            xrange(current_alignment.get_alignment_length()))
        self.position = np.asarray(
            xrange(current_alignment.get_alignment_length()))

        logging.debug("Alignment length  %d" %
                      current_alignment.get_alignment_length())

        if(gap_thresh):
            self._gap_alignment, _gap_filtered_position = self.clean_alignment(
                current_alignment, threshold=(abs(gap_thresh) <= 1 or 0.01) * abs(gap_thresh))
            current_alignment = self._gap_alignment
            logging.debug("Alignment length after filtering by gaps percent : %d" %
                          current_alignment.get_alignment_length())
            tt_filter_position = tt_filter_position[_gap_filtered_position]

        if(ic_thresh):
            align_info = AlignInfo.SummaryInfo(current_alignment)
            align_info.information_content()
            ic_vector = align_info.ic_vector
            # little hack for biopython < 1.69 and > 1.65
            if isinstance(ic_vector, dict):
                ic_vector = np.zeros(len(align_info.ic_vector))
                for (ic_i, ic_v) in align_info.ic_vector.items():
                    ic_vector[ic_i] = ic_v
            max_val = max(ic_vector) * \
                ((abs(ic_thresh) <= 1 or 0.01) * abs(ic_thresh))

            ic_pos = (np.asarray(ic_vector) >= max_val).nonzero()
            # ic_pos here is a tuple, so we should access the first element
            self._ic_alignment = self.filter_align_position(
                current_alignment, ic_pos[0])
            current_alignment = self._ic_alignment
            tt_filter_position = tt_filter_position[ic_pos]
            logging.debug("Alignment length after filtering with information content: %d" %
                          current_alignment.get_alignment_length())

        if(id_thresh):
            self._id_alignment, _id_filtered_position = self.filter_alignment(
                current_alignment, remove_identity=rmcnst, threshold=(abs(id_thresh) <= 1 or 0.01) * abs(id_thresh))
            current_alignment = self._id_alignment
            tt_filter_position = tt_filter_position[_id_filtered_position]
            logging.debug("Alignment length after filtering by sequence identity : %d" %
                          current_alignment.get_alignment_length())

        self.filt_prot_align = current_alignment
        self.filt_position = tt_filter_position

    def filter_codon_alignment(self, codon_alignment=None, ind_array=None, get_dict=False, alphabet=default_codon_alphabet):
        """Return the codon alignment from a list of in array"""
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
        if self.is_filtered:
            return self.codon_alignment, self.fcodon_alignment
        else:
            return self.codon_alignment, None

    def save_align(self, alignment, outfile, format='fasta'):
        """save alignment in a file """
        AlignIO.write(alignment, open(outfile, 'w'), format)

    def write_data(self, ori_alignment=None, id_filtered=None, gap_filtered=None, ic_filtered=None, tree=None, nofilter=False):
        """Save file depending on the file use as input"""
        if not nofilter:
            if id_filtered:
                self.save_align(self._id_alignment, id_filtered)
            if gap_filtered:
                self.save_align(self._gap_alignment, gap_filtered)
            if ic_filtered:
                self.save_align(self._ic_alignment, ic_filtered)
        if ori_alignment:
            self.save_align(self.prot_align, ori_alignment)
            ori_core = ori_alignment.replace('fasta', 'core')
            core = CoreFile.split_alignment(self.prot_align, self.gene_limits)
            CoreFile.write_corefile(core, ori_core)
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
                    # do not know why this raise an error on python3
                    seq += seqrec[int(i)]
                seqrec.letter_annotations = _RestrictedDict(length=len(seq))
                seqrec.seq = seq.upper()

        else:
            edited_alignment = clc.copy_codon_alignment(alignment, codontable)
            index_array = sorted(index_array)
            t = time.time()
            for seqrec in edited_alignment:
                # print(seqrec), len(seqrec), len(index_array)
                codseq = Seq('', generic_nucleotide)
                for pos in index_array:
                    codseq += seqrec.seq[pos * 3:(pos + 1) * 3]
                seqrec.seq = CodonSeq.from_seq(codseq)
        return edited_alignment

    @classmethod
    def clean_alignment(clc, alignment=None, characs=['-'], threshold=0.5):
        """Remove position of alignment which contain character from characs"""
        align_array = np.array([list(rec) for rec in alignment], dtype=str)
        indel_array = np.where(~(np.mean(np.in1d(align_array,
                                                 characs).reshape(align_array.shape), axis=0) >= threshold))[0].tolist()

        return clc.filter_align_position(alignment, indel_array), indel_array

    @classmethod
    def filter_alignment(clc, alignment, threshold=0.5, remove_identity=False, ambiguous='X', alphabet=generic_protein):
        """Filter an alignment using threshold as the minimum aa identity per columns"""
        # Smart move : Make a consensus sequence with threshold
        # remove all ambiguous positions
        consensus = clc.differed_consensus(alignment, threshold, remove_identity,
                                           ambiguous, alphabet=generic_protein)
        cons_array = [i for i, c in enumerate(consensus) if c != ambiguous]

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
    def differed_consensus(clc, alignment, threshold, remove_identity=False,
                           ambiguous='X', alphabet=generic_protein):
        """ Adapted from Biopython gap_consensus
        """
        consensus = ''

        # find the length of the consensus we are creating
        con_len = alignment.get_alignment_length()
        max_size = len(alignment)

        # go through each seq item
        for n in range(con_len):
            # keep track of the counts of the different atoms we get
            # we suppose all seqs in the alignment have the same length
            atom_dict = Counter(alignment[:, n])
            # for the sake of filtering, we only need one major aa per pos
            max_atom, count_atom = atom_dict.most_common(1)[0]

            if remove_identity and count_atom == con_len:
                consensus += ambiguous
            elif (float(count_atom) / float(max_size)) >= threshold:
                consensus += max_atom
            else:
                consensus += ambiguous
        return Seq(consensus, alphabet)

    @classmethod
    def get_consensus(clc, alignment, threshold, ambiguous='X', alphabet=generic_protein):
        """return consensus, using a defined threshold"""
        return ConsensusKeeper(alignment, threshold, ambiguous='X', alphabet=generic_protein)

    @classmethod
    def get_aa_filtered_alignment(clc, consensus, aa):
        """Get the filtered alignment for an amino acid"""
        aa = aa.upper()
        #cons_array = [i for i, c in enumerate(consensus) if c.upper() == aa]
        # return cons_array
        return consensus.get_maj_aa_pos(aa)

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
        self.confd = getattr(settings, 'CONF', 0.05)
        self.suspected_species = defaultdict(dict)
        self.aa2aa_rea = defaultdict(dict)
        self.settings = settings
        self.interesting_case = []
        self.reassignment_mapper = makehash()

    @classmethod
    def update_reas(clc, codon, cible_aa, speclist, codon_align, codonpos, filt_pos, genelim):
        """Update the list of codons reassignment and position"""
        rea_position_keeper = defaultdict(dict)
        genelim = sorted(genelim, key=lambda x: x[2])

        def _pos_in_gene(pos, glim=genelim, ind=0):
            while pos > glim[ind][-1]:
                ind += 1
            gene, gpos = glim[ind][0], pos - glim[ind][1]
            return ind, gene, gpos

        for spec in speclist:
            spec_dt = []
            i = 0
            for ps in codonpos:
                if codon_align[spec].seq[ps * 3:(ps + 1) * 3] == codon:
                    i, gene, pos = _pos_in_gene(filt_pos[ps], ind=i)
                    spec_dt.append((gene, pos))
            rea_position_keeper[spec]["%s:%s" % (codon, cible_aa)] = spec_dt
        return rea_position_keeper

    def export_position(self, rea_position_keeper, outfile):
        """Export reassignment position to a file"""
        with open(outfile, 'w') as OUT:
            json.dump(rea_position_keeper, OUT, indent=4)

    def compute_sequence_identity(self, matCalc=None):
        """Compute a distance matrix from the alignment"""
        if not matCalc:
            matCalc = DistanceCalculator(self.settings.MATRIX)
        if not self.settings.USE_GLOBAL:
            self.paired_distance = matCalc.get_distance(
                self.seqset.filt_prot_align)
        else:
            self.paired_distance = matCalc.get_distance(
                self.seqset.prot_align)

    def get_genomes(self, use_similarity=1):
        """ Get suspected genomes """
        matCalc = DistanceCalculator(self.settings.MATRIX)
        self.compute_sequence_identity(matCalc)
        # logging.debug("Distance matrix : ")
        # logging.debug(self.filtered_paired_distance)
        self.filtered_consensus = self.seqset.get_consensus(
            self.seqset.filt_prot_align, self.settings.AA_MAJORITY_THRESH)
        self.global_consensus = self.seqset.get_consensus(
            self.seqset.prot_align, self.settings.AA_MAJORITY_THRESH)

        number_seq = self.seqset.get_total_genomes()
        self.seq_names = list(self.seqset.common_genome)
        self.sim_json = defaultdict(list)

        self.seqset.aa_filt_prot_align = {}
        self.aa_paired_distance = {}
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
        aa2suspect_dist = makehash(1, list)
        for (j, i) in itertools.combinations(xrange(number_seq), r=2):
            paired = abs(
                use_similarity - self.paired_distance[self.seq_names[i], self.seq_names[j]])
            self.sim_json[self.seq_names[i]].append(
                {"distance": paired, "species": self.seq_names[j]})
            for aa in self.aa_paired_distance.keys():
                aapaired = abs(
                    use_similarity - self.aa_paired_distance[aa][self.seq_names[i], self.seq_names[j]])
                self.aa_sim_json[aa_letters_1to3[aa]].append(
                    {'distance': paired, 'aafiltered': aapaired, "species": "%s||%s" % (self.seq_names[i], self.seq_names[j])})
                if (paired > aapaired and use_similarity) or (paired < aapaired and not use_similarity):
                    aa2suspect[aa][self.seq_names[i]] = aa2suspect[
                        aa].get(self.seq_names[i], 0) + 1
                    aa2suspect[aa][self.seq_names[j]] = aa2suspect[
                        aa].get(self.seq_names[j], 0) + 1

                aa2suspect_dist[aa][self.seq_names[
                    i]].append([paired, aapaired])
                aa2suspect_dist[aa][self.seq_names[
                    j]].append([paired, aapaired])

        if self.settings.MODE == 'count':
            self.get_suspect_by_count(aa2suspect, number_seq)
        elif self.settings.MODE in ['wilcoxon', 'mannwhitney', 'ttest']:
            self.get_suspect_by_stat(
                aa2suspect_dist, number_seq, use_similarity, test=self.settings.MODE, confd=self.confd)
        else:
            self.get_suspect_by_clustering(aa2suspect_dist, number_seq)
        # logging.debug('The list of suspected species per amino acid is:')
        # logging.debug(self.suspected_species)

    def get_suspect_by_count(self, aa2suspect, seq_num):
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

        logging.debug('Clustering chosen, labels for each species :')
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
        self.aa_count_per_spec = defaultdict(Counter)
        f_prot_align = SeqIO.to_dict(self.seqset.filt_prot_align)
        self.total_aa_count = Counter()
        for spec in f_prot_align:
            self.aa_count_per_spec[spec] = Counter(f_prot_align[spec])
            self.total_aa_count.update(self.aa_count_per_spec[spec])
        return self.aa_count_per_spec

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

    def possible_aa_reassignation(self, speclist=None):
        """Find to which amino acid the reassignation are probable"""
        # TODO :  CHANGE THIS FUNCTION TO A MORE REALISTIC ONE
        aa_count_spec = self.get_aa_count_in_alignment()
        aa_count_cons = Counter(self.filtered_consensus.consensus)
        if not speclist:
            speclist = self.seqset.common_genome
        else:
            speclist = set.intersection(self.seqset.common_genome, speclist)
        for aa, suspect in self.suspected_species.items():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa]]
            suspected_list = sorted(suspect.keys())
            for spec in speclist:
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
                    prob = spec_aa_counter[cur_aa] / tot if tot > 0 else 0
                    if prob > 0 and cur_aa != '-' and cur_aa != aa and cur_aa not in self.settings.EXCLUDE_AA_FROM:

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

    def get_expected_prob_per_species(self, genome, aa1, aa2, use_cost=False, use_align=False, l=0.8):
        """Get expected prob for binomial test"""
        ncount = 0
        if not use_align:
            aa1F = self.aa_count_per_spec[genome][aa1]
            aa2F = self.aa_count_per_spec[genome][aa2]
            ncount = sum(self.aa_count_per_spec[genome].values())
        else:
            aa1F = self.total_aa_count[aa1]
            aa2F = self.total_aa_count[aa2]
            ncount = sum(self.total_aa_count.values())

        exp_prob = aa1F * aa2F * 2.0 / (ncount**2)
        if not use_cost:
            return exp_prob
        else:
            try:
                cost = self.settings.SUBMAT[(aa1, aa2)]
            except:
                cost = self.settings.SUBMAT[(aa2, aa1)]
            return exp_prob * np.exp(cost * l)

    def get_codon_usage(self):
        """Get Codon usage from species"""
        spec_data = {}
        codons_align = SeqIO.to_dict(self.seqset.codon_alignment)
        for (spec, codalign) in codons_align.items():
            cseq = str(codalign.seq)
            spec_data[spec] = Counter([cseq[x:x + 3]
                                       for x in range(0, len(cseq), 3)])
        try:
            del spec_data['---']
        except:
            pass
        return spec_data

    def save_json(self):
        """Save result into a json file"""
        with open(os.path.join(self.settings.OUTDIR, "reassignment.json"), "w") as outfile1:
            json.dump(self.reassignment_mapper, outfile1, indent=4)
        with open(os.path.join(self.settings.OUTDIR, "similarity.json"), "w") as outfile2:
            json.dump(self.sim_json, outfile2, indent=4)

    def save_predictions(self, preds, outfile):
        """Save all prediction in the root folder"""
        def _valid_state(spec, valid):
            state = "None"
            if valid:
                if spec in valid[0] and valid[1]:
                    state = "Both"
                elif spec in valid[0]:
                    state = "Clad"
                elif valid[1]:
                    state = "Algn"
            return state

        # reorganize predictions by codon, aa
        predict = defaultdict(list)
        max_spec_len = 0
        for (xlab, y, ypred, valid) in preds:
            xlab = xlab[y == 1, :]
            prob = ypred[y == 1, -1]
            for (i, line) in enumerate(xlab):
                reaformat = "%s (%s, %s)" % (line[1], line[2], line[3])
                max_spec_len = max(max_spec_len, len(line[0]))
                predict[reaformat].append(
                    [line[0], "%.3f" % prob[i], _valid_state(line[0], valid.get(line[1], None))])

        with open(outfile, 'w') as OUT:
            gcode = self.seqset.codontable.id
            gcode_descr = self.seqset.codontable.names[
                0] if self.seqset.codontable.names else ""
            OUT.write("#Reference Genetic Code : %d (%s)\n" %
                      (gcode, gcode_descr))
            OUT.write("#Exception: \n")
            for rea, total_rea in predict.items():
                OUT.write("%s\n" % rea)
                for spec_with_rea in total_rea:
                    OUT.write("\t%s\t%s\n" % (spec_with_rea[0].ljust(
                        max_spec_len), "\t".join(spec_with_rea[1:])))

    def save_all(self, predictions, rjson, savecodon=False, nofilter=False):
        """Save everything"""
        if savecodon:
            self.reassignment_mapper['codons'] = self.get_codon_usage()
        self.reassignment_mapper['aa'] = rjson
        self.save_json()
        if self.settings.SAVE_ALIGN:
            id_filtfile = os.path.join(
                self.settings.OUTDIR, "filt_alignment.fasta")
            ic_filtfile = os.path.join(self.settings.OUTDIR, "ic_filt.fasta")
            gap_filtfile = os.path.join(self.settings.OUTDIR, "gap_filt.fasta")
            ori_al = os.path.join(self.settings.OUTDIR, "ori_alignment.fasta")
            newick = os.path.join(self.settings.OUTDIR, "tree.nwk")
            self.seqset.write_data(ori_alignment=ori_al,
                                   id_filtered=id_filtfile, gap_filtered=gap_filtfile, ic_filtered=ic_filtfile, tree=newick, nofilter=nofilter)
        self.save_predictions(predictions, os.path.join(
            self.settings.OUTDIR, "predictions.txt"))

    def run_analysis(self, codon_align, fcodon_align, aasubset=None):
        """ Run the filtering analysis of the current dataset in sequenceset"""

        cur_aa2aa_rea = self.aa2aa_rea
        # we are going to limit to only the reassignment passed in input
        if aasubset:
            cur_aa2aa_rea = defaultdict(dict)
            for (_aa, _ab) in aasubset:
                if _ab in self.aa2aa_rea.get(_aa, []):
                    cur_aa2aa_rea[_aa].update({_ab: self.aa2aa_rea[_aa][_ab]})

        for aa1, aarea in cur_aa2aa_rea.items():
            aa_alignment = self.seqset.aa_filt_prot_align[aa_letters_1to3[aa1]]
            gcodon_rea = CodonReaData((aa1, aarea), self.seqset.prot_align, self.global_consensus, codon_align,
                                      self.seqset.codontable, self.seqset.position, self.seqset.gene_limits, self.settings)
            fcodon_rea = CodonReaData((aa1, aarea), self.seqset.filt_prot_align, self.filtered_consensus, fcodon_align,
                                      self.seqset.codontable, self.seqset.filt_position, self.seqset.gene_limits, self.settings)
            for aa2, species in aarea.items():
                # logging.debug("%s to %s" % (aa2, aa1))
                t = self.seqset.phylotree.copy("newick")
                fitch = SingleNaiveRec(t, species, aa_letters_1to3[aa2], aa_letters_1to3[
                    aa1], self.seqset.codontable, (gcodon_rea, fcodon_rea))
                slist = fitch.get_species_list(
                    self.settings.LIMIT_TO_SUSPECTED_SPECIES)
                alldata = {}
                for genome in self.seqset.common_genome:
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
                    g_rep_keeper = self.global_consensus.get_alt_cons()
                    for pos in xrange(self.global_consensus.allen):
                        pos_best_aa = sorted(
                            g_rep_keeper[pos], key=lambda x: x[-1])  # .pop()
                        aa_b, aa_s = pos_best_aa.pop()
                        pos_best_aa = [aa_b] + [x[0]
                                                for x in pos_best_aa if x[1] == aa_s]
                        if aa1 in pos_best_aa and rec[pos] == aa2:
                            leaf.count += 1

                    # filtered codon infos
                    freacodon = fcodon_rea.get_reacodons(genome, aa2)
                    fusedcodon = fcodon_rea.get_usedcodons(genome, aa2)
                    fmixtecodon = fcodon_rea.get_mixtecodons(genome, aa2)
                    # global codon infos
                    greacodon = gcodon_rea.get_reacodons(genome, aa2)
                    gusedcodon = gcodon_rea.get_usedcodons(genome, aa2)
                    gmixtecodon = gcodon_rea.get_mixtecodons(genome, aa2)

                    # settings parameters
                    reacodon, usedcodon, mcodon = freacodon, fusedcodon, fmixtecodon
                    if self.settings.USE_GLOBAL:
                        reacodon, usedcodon, mcodon = greacodon, gusedcodon, gmixtecodon

                    eprob = None
                    if len(reacodon.values()) == 1 or len(usedcodon.values()) == 1:
                        eprob = self.get_expected_prob_per_species(genome, aa2, aa1,
                                                                   use_cost=True, use_align=True)
                        # print("%s | %s to %s : %f"%(genome, aa2, aa1, eprob))
                        assert (
                            eprob - 1 < 0), "Strange value for eprob %f" % eprob
                    fisher_passed, pval = independance_test(
                        reacodon, usedcodon, genome, confd=self.confd, expct_prob=eprob, force_chi2=self.settings.FORCED_CHI2)
                    # print "%s\t%s\t%s ==> %s" %(aa2, aa1, genome,
                    # str(column_comparision) )
                    # if('lost' in leaf.features and fisher_passed and
                    # leaf.count > self.settings.COUNT_THRESHOLD
                    if('lost' in leaf.features and fisher_passed):
                        leaf.lost = False

                    gdata = {}
                    if self.settings.CODON_COUNT_THRESHOLD > 1:
                        if self.settings.USE_GLOBAL:
                            greacodon = Counter({_cod: _count for _cod, _count in greacodon.iteritems(
                            ) if _count > self.settings.CODON_COUNT_THRESHOLD})
                        else:
                            freacodon = Counter({_cod: _count for _cod, _count in freacodon.iteritems(
                            ) if _count > self.settings.CODON_COUNT_THRESHOLD})
                    g_rea_dist = gcodon_rea.get_rea_aa_codon_distribution(
                        genome, aa2)
                    g_total_rea_dist = gcodon_rea.get_total_rea_aa_codon_distribution(
                        genome, aa2)
                    gdata['global'] = {'rea_codon': greacodon, 'used_codon': gusedcodon, 'mixte_codon': gmixtecodon,
                                       'count': leaf.count, 'rea_distribution': g_rea_dist, 'total_rea_distribution': g_total_rea_dist}

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

                if(fitch.is_valid(self.settings.COUNT_THRESHOLD) or self.settings.SHOW_ALL):
                    logging.debug(
                        "- Case with predictions: %s to %s" % (aa2, aa1))
                    yield (fitch, alldata, slist)
            # free objects
            del aa_alignment
            del gcodon_rea
            del fcodon_rea


def executeCMD(cmd, prog):
    """Execute a command line in the shell"""
    logging.debug("The following will be executed : \n%s\n" % cmd)
    p = subprocess.Popen(
        cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    if err:
        logging.debug(err)
    if out:
        pass
        # logging.debug("%s : \n----------------- %s" % (prog, out))
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
            right_branch = node.children[1]
            left_branch_leaf = left_branch.get_leaf_names()
            right_branch_leaf = right_branch.get_leaf_names()
            seqnames = [min(map(lambda x: seq_order.index(x), left_branch_leaf)) +
                        1, min(map(lambda x: seq_order.index(x), right_branch_leaf)) + 1, ]
            # seqnames = [seq_order.index(left_branch_leaf.name)+1, seq_order.index(right_branch_leaf.name)+1]
            # seqnames = [left_branch_leaf.name, right_branch_leaf.name]
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


def improve_is_stagned(alignfile1, alignfile2, prob=1):
    """Check if alignment are identical"""
    al1 = AlignIO.read(alignfile1, "fasta")
    al1.sort()
    al2 = AlignIO.read(alignfile2, "fasta")
    al2.sort()
    identical = (al1.format("fasta") == al2.format("fasta"))
    rand = np.random.uniform()
    return identical and (rand <= prob)


def remove_gap_only_columns(alignfile, curformat):
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


def independance_test(rea, ori, genome, confd=0.05, expct_prob=0.5, force_chi2=False):
    """Perform a Fisher's Exact test"""
    codon_list = set(rea.keys())
    codon_list.update(ori.keys())
    codon_list = [x for x in codon_list if not (
        rea.get(x, 0) == 0 and ori.get(x, 0) == 0)]
    nelmt = len(codon_list)

    # cannot possibily be reassigned when not present
    if not rea:
        return False, 2

    # line 1 is rea, line 2 is observed
    # in this case, we can perform a fisher test
    if nelmt > 1:
        obs = np.zeros((nelmt, 2))
        for i in xrange(nelmt):
            obs[i, 0] = rea.get(codon_list[i], 0)
            obs[i, 1] = ori.get(codon_list[i], 0)
        try:
            # check mat, to select chi2 or fisher_test
            # faut just croiser les doigts :)
            can_chi2 = chi2_is_possible(obs)
            if force_chi2 and (can_chi2[0] or can_chi2[1]):
                c, pval, dof, t = ss.chi2_contingency(obs)
            else:
                if np.sum(obs) < 1000:
                    pval = fisher_exact(obs, midP=True, attempt=1)
                elif np.any(obs < 5) and np.sum(obs) <= 1500:
                    pval = fisher_exact(obs, midP=True, attempt=2)
                else:
                    c, pval, dof, t = ss.chi2_contingency(obs)
                # fallback to chi2 test if fisher is impossible
        except Exception as e:
            logging.debug(
                "**warning: %s (%s, %s) trying chi2 instead of FISHER EXACT" % (genome, ori, rea))
            c, pval, dof, t = ss.chi2_contingency(
                obs)  # nothing to prevent this if error
        return pval <= confd, max(pval, SYSMINVAL)

    # strangely, codon is used only in rea column
    # complete reassignment ??
    # to avoid returning 0, we return the smallest positive int
    elif len(rea.values()) > 0 and len(ori.values()) == 0:
        return True, SYSMINVAL
    # codon is used in both column
    elif len(rea.values()) > 0 and len(ori.values()) > 0:
        # fpval = abs(rea.values()[0] - ori.values()[0]) / \
            # (rea.values()[0] + ori.values()[0])
        # return rea.values()[0] >= ori.values()[0], fpval / tot_size
        n = rea.values()[0] + ori.values()[0]  # ignoring mixcodon ??
        pval = onevalbinomtest(rea.values()[0], n, expct_prob)
        return pval <= confd, max(pval, SYSMINVAL)
    # In this case, the codon is neither in the rea column nor in the
    # second column, strange result
    else:
        return False, 2


def onevalbinomtest(obs, n, prob):
    """Return a binomial test given success and prob"""
    return ss.binom_test(obs, n, prob)


def chi2_is_possible(obs, nonstrict=True):
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
    if count:
        return nonstrict, (count * 1.0 / size) < (1 - limit)
    return True, True


@timeit
def execute_alignment(cmdline, inp, out):
    """Construct alignment command line and execute it """
    prog = 'mafft'
    if 'muscle' in cmdline:
        prog = 'muscle'
        cmdline += " -in %s -out %s" % (inp, out)
    elif 'mafft' in cmdline:
        cmdline += " %s > %s" % (inp, out)
    else:
        raise ValueError(
            "Cannot execute %s. Programme not expected. You can provide your own alignment instead.")

    executeCMD(cmdline, prog)


def compute_SP_per_col(al1, al2, columns, nspec, scoring_matrix):
    """Compute a SP score per column"""
    def scoring_function(aa1, aa2, scoring_matrix, gap_cost=-1):
        score = 0
        if aa1 == aa2 == '-':
            score = 0
        elif scoring_matrix == 'identity':
            score = (aa1 == aa2) * 1
        elif '-' in [aa1, aa2]:
            score = gap_cost
        else:
            # controversial decision
            # give -1 to gap event
            try:
                score = scoring_matrix.get(
                    (aa1, aa2), scoring_matrix.get((aa2, aa1)))
            except:
                score = gap_cost
        return score

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


def compute_ic_content(alignment):
    if isinstance(alignment, MultipleSeqAlignment):
        align_info = AlignInfo.SummaryInfo(alignment)
    else:
        align_info = AlignInfo.SummaryInfo(
            MultipleSeqAlignment(alignment.values()))
    align_info.information_content()
    ic_vector = align_info.ic_vector
    # biopython wrong version hack
    if isinstance(ic_vector, dict):
        ic_vector = np.zeros(len(align_info.ic_vector))
        for (ic_i, ic_v) in align_info.ic_vector.items():
            ic_vector[ic_i] = ic_v
    return list(ic_vector)


def extract_codaln(codon_alignment, alignment, not_wanted_species):
    new_cod_alignment = []
    new_alignment = []
    for spec, nucseq in codon_alignment.items():
        if spec not in not_wanted_species:
            new_cod_alignment.append(nucseq)
            new_alignment.append(alignment[spec])
    return new_cod_alignment, new_alignment


def translate_codaln(codon_alignment, codontable, changes=None):
    if changes is None:
        changes = {}
    # return dict((spec, "".join(map(lambda x: codontable[x] \
    # if changes.get(spec, (None,))[0] != x else changes[spec][1],
    # [prot[i:i + 3] for i in xrange(0, len(prot), 3)]))) for spec, prot in codon_alignment.items())
    translated_al = {}
    position = set([])
    for spec, nucseq in codon_alignment.items():
        translated = ""
        nucseq_len = int(len(nucseq.seq) / 3)
        for i in xrange(0, nucseq_len):
            cod = nucseq[i * 3:(i + 1) * 3].seq.tostring()
            if changes.get(spec, (None,))[0] != cod:
                # use X for amino acid when stop codon is found
                translated += codontable.get(cod, 'X')
            else:
                # in this case, there were reassignment
                translated += changes[spec][1]
                position.add(i)
        translated_al[spec] = SeqRecord(
            Seq(translated, generic_protein), id=spec, name=spec)
    return translated_al, sorted(list(position))


def check_gain(codon, cible_aa, speclist, tree, codontable, codon_alignment,
               scoring_method="identity", alignment=None, ic_cont=None, method="wilcoxon", spec_filter=True):
    """Check if there is an actuall gain in global sequence quality after applying reassignment"""
    translate = translate_codaln
    extract_alignment = extract_codaln
    if not isinstance(codon_alignment, dict):
        codon_alignment = SeqIO.to_dict(codon_alignment)

    if not alignment:
        alignment, _ = translate(codon_alignment, codontable)

    if spec_filter:
        SingleNaiveRec._dollo(tree)
        fcod_aln, f_aln = extract_alignment(
            codon_alignment, alignment, speclist)
        fake_rea = set([])
        for spec in speclist:
            # not just spec but
            # either we get all sister here
            # or all node under the same predicted reassignment
            # which could take too long
            cur_par = (tree & spec).up
            while cur_par.up is not None and cur_par.reassigned != {1}:
                cur_par = cur_par.up
            spec_sis = [x.name for x in cur_par if x.name in speclist]
            cur_recs = []
            cur_recs_al = []
            codchange = {}
            for x in spec_sis:
                cur_recs.append(codon_alignment[x])
                cur_recs_al.append(alignment[x])
                codchange[x] = (codon, cible_aa)
            corf_aln, pos = translate(SeqIO.to_dict(fcod_aln + cur_recs),
                                      codontable, codchange)
            f_aln_s = MultipleSeqAlignment(f_aln + cur_recs_al)
            # spec_ic = compute_ic_content(f_aln_s)
            # spec_cor_ic = compute_ic_content(corf_aln)
            pval, cor_al_sp, al_sp = check_align_upgrade(
                corf_aln, SeqIO.to_dict(f_aln_s), scoring_method, method, pos)

            # simple validation test here
            validated = sum(cor_al_sp) > sum(al_sp)
            logging.debug("Validation : %s | %s : (%s to %s) | (%f --> %f)" %
                          (validated, spec, codon, cible_aa, sum(cor_al_sp), sum(al_sp)))

            if not validated:
                fake_rea.add(spec)

        speclist = list(set(speclist) - fake_rea)
        if not speclist:
            return None, None, None, None, None, []

    changes = {}
    for spec in speclist:
        changes[spec] = (codon, cible_aa)

    cor_alignment, position = translate(codon_alignment, codontable, changes)
    score_improve, cor_al_sp, al_sp = check_align_upgrade(
        cor_alignment, alignment, scoring_method, method, position)

    cor_ic_cont = compute_ic_content(cor_alignment)
    if not ic_cont:
        ic_cont = compute_ic_content(alignment)

    return score_improve, (al_sp, cor_al_sp), (ic_cont, cor_ic_cont), (alignment, cor_alignment), position, speclist


def get_rea_genome_list(xlabel, y):
    """Return a dict containing each directory
    and if whether or it has a reassignation"""
    rea_genome = {}
    all_pos = xlabel[y == 1, 0]
    for g in all_pos:
        rea_genome[g] = True
    return rea_genome


def get_codon_set_and_genome(y, xlabel, trueval=None):
    """Get for each codon the list of reassigned genome"""
    all_pos = xlabel[:, 0:2]
    if trueval:
        all_pos = xlabel[y == trueval, 0:2]
    result = defaultdict(list)
    for (g, cod) in all_pos:
        result[cod].append(g)
    return result


def get_suffix(x, add_label):
    """Get suffix """
    return x if add_label else ""


def get_report(fitchtree, gdata, codon_align, prediction, genelimit, filt_position, settings, output="", pie_size=45):
    """ Render tree to a pdf"""
    OUTDIR = purge_directory(os.path.join(settings.OUTDIR, fitchtree.ori_aa +
                                          "_to_" + fitchtree.dest_aa))

    c_rea = get_rea_genome_list(prediction[1], prediction[3])
    if not output:
        output = os.path.join(OUTDIR, "Codon_data." + settings.IMAGE_FORMAT)
    if(GRAPHICAL_ACCESS):
        lprovider = extendfn(layout_func, show_n=settings.ADD_NUMBER_PIE, c_rea=c_rea, fitchtree=fitchtree, mixte_codons=settings.SHOW_MIXTE_CODONS,
                             global_codon=settings.SHOW_GLOBAL_CODON_DATA, gdata=gdata, get_suffix=partial(get_suffix, add_label=settings.ADD_LABEL_TO_LEAF), pie_size=pie_size)
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

        ts.layout_fn = lprovider
        # header declaration
        h1 = TextFace(fitchtree.dest_aa, fsize=10, fgcolor="#aa0000")
        h2 = TextFace(fitchtree.ori_aa, fsize=10, fgcolor="#aa0000")
        h3 = TextFace(fitchtree.ori_aa + " O.", fsize=10, fgcolor="#aa0000")
        h1g = TextFace(fitchtree.dest_aa, fsize=10, fgcolor="#aa0000")
        h2g = TextFace(fitchtree.ori_aa, fsize=10, fgcolor="#aa0000")
        h3g = TextFace(fitchtree.ori_aa + " O.", fsize=10, fgcolor="#aa0000")
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
            ts.aligned_header.add_face(h1g, column=next_column + 1)
            ts.aligned_header.add_face(h2g, column=next_column + 2)
            next_column += 3
            if(settings.SHOW_MIXTE_CODONS):
                ts.aligned_header.add_face(h3g, column=next_column)
                next_column += 1

        f_pval_h = TextFace("FE. pval", fsize=10, fgcolor="#aa0000")
        f_pval_h.vt_align = 1
        f_pval_h.hz_align = 1
        ts.aligned_header.add_face(f_pval_h, column=next_column + 1)

        ts.title.add_face(TextFace(fitchtree.ori_aa + " --> " +
                                   fitchtree.dest_aa, fsize=14), column=0)
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

        fitchtree.tree.render(output, dpi=700, tree_style=ts)

        rkp, data_var, codvalid = None, None, {}
        if settings.VALIDATION:
            # now get sequence retranslation improvement
            table = fitchtree.dct.forward_table
            # add the gap to codon_table
            table['---'] = '-'
            table['...'] = '.'
            data_present, data_var, rkp, codvalid = codon_adjust_improve(
                fitchtree, codon_align, table, prediction, genelimit, filt_position, settings, outdir=OUTDIR)

        # get report output
        rep_out = os.path.join(OUTDIR, "Report_" +
                               fitchtree.ori_aa + "_to_" + fitchtree.dest_aa)
        pdf_format_data(fitchtree.ori_aa1, fitchtree.dest_aa1, gdata, prediction,
                        codvalid, 'filtered', rep_out + ".pdf", settings.VALIDATION)

        return data_var, OUTDIR, rkp, codvalid


def format_tree(codon, alignment, SP_score, ic_contents, pos=[], tree=None, cible=None,
                limits=(None, None, None), dtype="aa", codontable={}, codon_col={}, colors=None):
    """Format the rendering of tree data for alignment"""
    t = tree.copy('newick')

    gpos, limiter, start_holder = limits

    def _get_genes(seq, jump=1):
        """ Return a dict of seq from the sequence"""
        gseq = {}
        ppos = 0
        if limiter is None:
            gseq['All'] = "".join([seq[i] for i in pos])
        else:
            for (l, name) in limiter:
                gseq[name] = "".join([seq[i * jump:jump * (i + 1)]
                                      for i in pos[ppos:l]])
                ppos = l
        return gseq

    for node in t:
        node_seq = alignment[node.name].seq.tostring()
        if pos:
            node_seq = _get_genes(node_seq, 2 * (dtype == "codon") + 1)
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
            faces.add_face_to_node(
                AttrFace('name', fsize=14, fgcolor=colors.get(node.name, 'black')), node, 0, position="aligned")
            if hasattr(node, "sequence"):
                ind = 0
                for (k, name) in limiter:
                    seq = node.sequence[name]
                    seqface = faces.SequenceFace(seq, fsize=13)
                    if dtype != 'aa':
                        seqface = SequenceFace(seq, cible, seqtype=dtype, fsize=13,
                                               codontable=codontable, spec_codon_col=codon_col)
                    faces.add_face_to_node(
                        seqface, node, 2 + ind, aligned=True)
                    ind += 2

    ts.layout_fn = layout

    if dtype == 'aa':

        def gfunc(x):
            return '*' if x == 1 else ' '

        if pos:

            ts.title.add_face(TextFace('(%s) - SP score : %.0f | IC = %.2f' % (codon, sum(SP_score), sum(ic_contents)),
                                       fsize=14, fgcolor='red'), 0)
            ts.aligned_header.add_face(
                faces.RectFace(14, 14, 'white', 'white'), 1)

            ts.aligned_foot.add_face(
                faces.RectFace(14, 14, 'white', 'white'), 1)
            start = 0
            ind = 3
            for (l, name) in limiter:
                listdata = pos[start:l]
                ic_content = np.asarray(ic_contents)[listdata]
                # sp_score = np.asarray(SP_score)[pos]
                footer_seq = SequenceFace(
                    "".join([gfunc(st) for st in start_holder[start:l]]), None, dtype, fsize=13)
                start = l
                ic_plot = faces.SequencePlotFace(ic_content, hlines=[(int(min(ic_content) - 0.5))],
                                                 hlines_col=['white'], ylim=(int(min(ic_content) - 0.5),
                                                                             int(max(ic_content) + 1)), fsize=10, col_width=14,
                                                 header="IC", kind='bar')
                ts.aligned_header.add_face(ic_plot, ind)
                ts.aligned_foot.add_face(List90Face(listdata, fsize=10), ind)
                ts.aligned_foot.add_face(footer_seq, ind)
                ts.aligned_foot.add_face(TextFace(name), ind)
                ind += 2

    else:
        for (cod, col) in codon_col.items():
            ts.legend.add_face(faces.RectFace(50, 25, col, col), column=0)
            ts.legend.add_face(TextFace("  %s " % cod, fsize=8), column=1)
        ts.legend.add_face(faces.RectFace(50, 25, 'black', 'black'), column=0)
        ts.legend.add_face(TextFace("  %s\'s codons " %
                                    (cible), fsize=8), column=1)
        ts.legend_position = 3
        ind = 1
        lprev = 0
        for (l, name) in limiter:
            ts.aligned_foot.add_face(faces.RectFace(
                14 * (l - lprev) * 3, 14, '#ffffff', 'white'), ind + 1)
            ts.aligned_foot.add_face(TextFace(' ', fsize=14), ind)
            ts.aligned_foot.add_face(TextFace(name, fsize=14), ind + 1)
            lprev = l
            ind += 2

    t.dist = 0
    return t, ts


def identify_position_with_codon(fcodal, codon, spec_to_check):
    """Get all positions where a codon is used"""
    positions = []
    try:
        al_len = int(len(fcodal.values()[0]) / 3)
    except:
        al_len = fcodal.get_aln_length()
    for i in xrange(al_len):
        for spec in spec_to_check:
            if str(fcodal[spec].seq[i * 3:(i + 1) * 3]) == codon:
                positions.append(i)
                break
    return positions


def violin_plot(vals, output, score, codon, cible, imformat="pdf"):
    """Return violin plot of SP score """
    keys = vals.keys()
    data = [vals[k] for k in keys]
    pos = [y + 1 for y in range(len(data))]
    title = 'p-values({} --> {}) : {:.2e}'.format(codon,
                                                  cible, max(score, SYSMINVAL))
    output = output + "." + imformat

    if not SEABORN:
        figure = plt.figure()
        axes = figure.add_subplot(1, 1, 1)
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
        plt.setp(axes, xticks=pos, xticklabels=keys)
        axes.set_title(title)
        axes.set_ylabel('Alignment score/col')
        figure.savefig(output, format=imformat)
    else:
        ax = sns.violinplot(data=data)
        ax.set_xticklabels(keys)
        ax.set_title(title)
        ax.set_ylabel('Alignment score/col')
        fig = ax.get_figure()
        fig.savefig(output, format=imformat)
        fig.clf()

    return output, (codon, cible, score)


def codon_adjust_improve(fitchtree, codon_align, codontable, prediction, genelimit, filt_position, settings, outdir=""):
    """Get a representation of the improvement after translation"""

    cible_aa = fitchtree.dest_aa1
    X_data, X_labels, pred_prob, pred = prediction
    true_codon_set = get_codon_set_and_genome(pred, X_labels, 1)

    sc_meth = settings.MATRIX
    method = settings.MODE if settings.MODE in [
        'wilcoxon', 'mannwhitney', 'ttest'] else 'wilcoxon'
    outputs = []
    violinout = None
    data_var = {}
    ori_al = None
    ic = None
    tree = fitchtree.tree.copy("newick-extended")
    rea_pos_keeper = defaultdict(dict)
    codvalid = {}
    # some local function
    limit_finder_fn = extendfn(
        _limit_finder, gene_limit=genelimit, proche=settings.STARTDIST)
    format_tree_fn = extendfn(format_tree, tree=tree, cible=cible_aa)

    for codon in true_codon_set.keys():
        speclist = true_codon_set[codon]
        if len(speclist) > 0:
            # pos = identify_position_with_codon(codon_align, codon, speclist)
            # check_gain is called only on filtered alignment
            # maybe it's a better idea to call it on the original alignmnt
            score_improve, alsp, alic, als, pos, speclist = check_gain(codon, cible_aa, speclist, tree, codontable,
                                                                       codon_align, scoring_method=sc_meth,
                                                                       alignment=ori_al, ic_cont=ic, method=method)
            if not speclist:
                # this mean that all predictions were probably fake
                continue

            codon_pos = [filt_position[i] for i in pos]
            limits = limit_finder_fn(codon_pos)
            if settings.COMPUTE_POS:
                # only compute this if asked
                rea_pos = ReaGenomeFinder.update_reas(
                    codon, cible_aa, speclist, codon_align, pos, filt_position, genelimit)
                for cuspec, readt in rea_pos.items():
                    for k in readt.keys():
                        rea_pos_keeper[cuspec][k] = readt[k]
            sp, cor_sp = alsp
            ic, cor_ic = alic
            ori_al, new_al = als
            violinout, viout = violin_plot({'Original': sp, 'Corrected': cor_sp},
                                           os.path.join(
                                               outdir, "%s_violin" % codon),
                                           score_improve, codon, fitchtree.dest_aa, imformat=settings.IMAGE_FORMAT)
            tmpvalid = dict((x, 'crimson') for x in speclist)

            if settings.SAVE_ALIGN:
                ori_t, ori_ts = format_tree_fn(
                    codon, ori_al, sp, ic, pos, limits=limits, colors=tmpvalid)
                rea_t, rea_ts = format_tree_fn(
                    codon, new_al, cor_sp, cor_ic, pos, limits=limits, colors=tmpvalid)
                cod_t, cod_ts = format_tree_fn(codon, codon_align, None, None, pos, limits=limits,
                                               dtype="codon", codontable=codontable, codon_col=fitchtree.colors, colors=tmpvalid)

                cod_ts.title.add_face(TextFace(
                    "Prediction validation for " + codon + " to " + fitchtree.dest_aa, fsize=14), column=0)

                cod_out = os.path.join(outdir, "%s_codons.%s" %
                                       (codon, settings.IMAGE_FORMAT))
                ori_out = os.path.join(outdir, "%s_ori.%s" %
                                       (codon, settings.IMAGE_FORMAT))
                rea_out = os.path.join(outdir, "%s_rea.%s" %
                                       (codon, settings.IMAGE_FORMAT))

                ori_t.render(ori_out, dpi=700, tree_style=ori_ts)
                rea_t.render(rea_out, dpi=700, tree_style=rea_ts)
                cod_t.render(cod_out, dpi=700, tree_style=cod_ts)

            logging.debug('{} --> {} : {:.2e}'.format(*viout))

            codvalid[codon] = (tmpvalid, viout[-1] < settings.CONF)

            data_var[codon] = score_improve

            outputs.append(True)

    return len(outputs) > 0, data_var, rea_pos_keeper, codvalid


def pdf_format_data(ori_aa, dest_aa, gdata, prediction, codvalid, dtype, output="", validate=True):
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
        genome_list = codon_set[codon]
        for g in genome_list:
            pos = pred_dict[g][codon]
            dt_by_att = {
                'Rea. count': gdata[g][dtype]['rea_codon'].get(codon, 0),
                'Use. count': gdata[g][dtype]['used_codon'].get(codon, 0),
                'Tot. count': gdata[g]['codons'][dtype].get(codon, 0),
                'Rea. Genes': len(gdata[g][dtype]['rea_distribution'].get(codon, [])),
                'Tot. Genes': len(gdata[g][dtype]['total_rea_distribution'].get(codon, [])),
                'Telford': gdata[g]['score'][dtype].get(codon, 0),
                # position 1 is true prob
                'Reassigned': pred[pos] == 1,
                # position 1 is true prob
                'Prob': pred_prob[pos, 1],
            }
            if validate:
                if codvalid.get(codon, None):
                    dt_by_att['Clad. Val'] = bool(
                        codvalid[codon][0].get(g, ''))
                    dt_by_att['Trans. Val'] = codvalid[codon][1]
                else:
                    dt_by_att['Clad. Val'] = False
                    dt_by_att['Trans. Val'] = False
            pd_data[g][codon] = dt_by_att

    frames = []
    glist = []
    for (g, dt) in pd_data.items():
        glist.append(g)
        frames.append(pd.DataFrame.from_dict(dt, orient='index'))
    fdata = pd.concat(frames, keys=glist)
    html_data = fdata.to_html()
    data_var = {'title': 'Codon reassignment prediction in all genome',
                'report_table': html_data,
                }
    return export_from_html(data_var, output)


def print_data_to_txt(outputfile, header, X, X_label, Y, Y_prob, codon_data, cible, suptext=None, valid={}):
    """Print result in a text file"""
    out = Output(file=outputfile)
    out.write("### Random Forest prediction\n")
    out.write("\t".join(["genome", "codon",
                         "ori_aa", "rea_aa"] + list([x.replace(' ', '_') for x in header]) + ["prediction", "probability"] +
                        (["clad_valid", "trans_valid"] if valid else [])))

    total_elm = len(Y)
    for i in xrange(total_elm):
        end_data = [str(Y[i]), str(Y_prob[i][-1])]
        if valid:
            try:
                gtmp_val = str(valid[X_label[i, 1]][-1] * 1)
                tmp_val = str(
                    bool(valid[X_label[i, 1]][0].get(X_label[i, 0], '')) * 1)
            except:
                tmp_val = '0'
                gtmp_val = '0'
            end_data.extend([tmp_val, gtmp_val])
        out.write("\n" + "\t".join(list(X_label[i]) + [str(x)
                                                       for x in X[i]] + end_data))

    if codon_data:
        out.write("\n\n### Alignment improvement:\n")
        for (codon, score) in codon_data.items():
            out.write("{}\t{}\t{:.2e}\n".format(codon, cible, score))

    tmp = Y_prob[Y > 0, 1]
    if len(tmp) > 0:
        prerange = (min(tmp), max(tmp))
        write_d = "\n### Prediction Prob. range for Positive :[ %.3f - %.3f ]\n" % prerange
        out.write("\n{}\n".format(write_d))

    if suptext and isinstance(suptext, basestring):
        out.write("\n\n{}\n".format(suptext))
    out.close()


def makehash(depth=None, type=None):
    """Utility method to make a multilevel dict"""
    if (depth, type) == (None, None):
        return defaultdict(makehash)
    elif depth == 0:
        return defaultdict(type)
    else:
        return defaultdict(partial(makehash, depth - 1, type))


def parse_spec_to_code(filename):
    spec_to_code = {}
    with open(filename) as GIN:
        for line in GIN:
            line = line.strip()
            if not line.startswith('#'):
                spec, code = line.split()
                spec_to_code[spec] = int(code)
    return spec_to_code


def extendfn(func, *args, **keywords):
    def newfunc(*fargs, **fkeywords):
        newkeywords = keywords.copy()
        newkeywords.update(fkeywords)
        fargs += args
        return func(*fargs, **newkeywords)
    newfunc.func = func
    newfunc.args = args
    newfunc.keywords = keywords
    return newfunc


def _limit_finder(codon_pos, gene_limit, proche):
    """ Return the gene for each position of codon pos"""
    gene_pos = []
    gene_break = []
    gll = len(gene_limit)
    start = 0
    close_to_start = []
    for i in range(len(codon_pos)):
        brokn = False
        entering = start
        while codon_pos[i] >= gene_limit[start][2] and start < gll:
            start += 1
            brokn = True
        if brokn and i > 0:
            gene_break.append((i, gene_limit[entering][0]))
        if codon_pos[i] - gene_limit[start][1] <= proche:
            close_to_start.append(1)
        else:
            close_to_start.append(0)
        gene_pos.append(gene_limit[start][0])

    gene_break.append((len(codon_pos), gene_limit[start][0]))
    # logging.debug("Current codon pos: %s"%" ".join([str(x) for x in codon_pos]))
    # logging.debug("Gene break : %s"%" ".join([str(x) for x in gene_break]))
    # logging.debug("Gene pos : %s"%" ".join([str(x) for x in gene_pos]))

    return gene_pos, gene_break, close_to_start


def layout_func(node, show_n, c_rea, fitchtree, mixte_codons, global_codon, gdata, pie_size, get_suffix):
    """Layout function for tree drawing"""
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

        if(mixte_codons):
            spec_mixtecodon_f = gdata[node.name][
                'filtered']['mixte_codon']
            faces.add_face_to_node(PPieChartFace(spec_mixtecodon_f.values(), pie_size, pie_size, show_label=show_n,
                                                 colors=[fitchtree.colors[k] for k in spec_mixtecodon_f.keys()]),
                                   node, column=3, position="aligned")
            next_column = 4

        if(global_codon):
            spec_codonrea_g = gdata[node.name]['global']['rea_codon']
            spec_codonused_g = gdata[node.name]['global']['used_codon']

            # add separator
            faces.add_face_to_node(LineFace(
                pie_size, pie_size, None), node, column=next_column, position="aligned")

            faces.add_face_to_node(PPieChartFace(spec_codonrea_g.values(), pie_size, pie_size, show_label=show_n,
                                                 colors=[fitchtree.colors[k] for k in spec_codonrea_g.keys()]),
                                   node, column=next_column + 1, position="aligned")

            faces.add_face_to_node(PPieChartFace(spec_codonused_g.values(), pie_size, pie_size, show_label=show_n,
                                                 colors=[fitchtree.colors[k] for k in spec_codonused_g.keys()]),
                                   node, column=next_column + 2, position="aligned")

            next_column += 3
            if(mixte_codons):
                spec_mixtecodon_g = gdata[node.name][
                    'global']['mixte_codon']
                faces.add_face_to_node(PPieChartFace(spec_mixtecodon_g.values(), pie_size, pie_size, show_label=show_n,
                                                     colors=[fitchtree.colors[k] for k in spec_mixtecodon_g.keys()]),
                                       node, column=next_column, position="aligned")
                next_column += 1

        faces.add_face_to_node(LineFace(
            pie_size, pie_size, None), node, column=next_column, position="aligned")

        faces.add_face_to_node(TextFace("{:.2e}".format(gdata[node.name]['lost']['pval']), fsize=10, fgcolor="#000"),
                               node, column=next_column + 1, position="aligned")
