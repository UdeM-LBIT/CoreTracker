# This module is currently not implemented and
# is not working. In the future, an algorithm to check
# for RNA editing will be implemented here

from Bio.Seq import Seq
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Align import AlignInfo, MultipleSeqAlignment
import itertools
from collections import defaultdict as ddict


class CodonGraph(dict):

    def __init__(self, subsrea, table):
        self.subsrea = subsrea
        self.codontable = table
        self.codon_graph = ddict(list)
        self.create_graph()

    def create_graph(self):
        for codon in codontable.forward_dict.keys():
            codtrans = CodonTransition(codon, self.subsrea, self.codontable)
            for trans in codtrans.iter_possible_transition():
                self.codon_graph[codon].append(trans)
            del codtrans

    def __getitem__(self, codon):
        return self.codon_graph[codon]


class CodonTransition(object):

    def __init__(self, codon, subsrea, table):
        self.codon = list(codon)
        self.subsrea = subsrea
        self.table = table

    def iter_possible_transition(self):
        expc_trans = [i for i, x in enumerate(
            self.codon) if self.subsrea.is_edited(x)]
        pstates = [0, 1]
        for state in itertools.product([0, 1], repeat=len(expc_trans)):
            state = map(lambda x: state[x] *
                        expc_trans[x], range(len(expc_trans)))
            codon = self.codon[:]
            changed = False
            for pos in state:
                if pos:
                    changed = True
                    codon[pos] = self.subsrea.get_edited[codon[pos]]
            codon = "".join(codon)
            expected_aa = self.table.get(codon, None)
            if changed and expected_aa:
                # it's possible that the return is a stop codon
                yield (codon, expected_aa)


class SubsReaType(object):

    def __init__(self, inputnuc, outputnuc):
        self.inputnuc = inputnuc.upper()
        self.outputnuc = outputnuc.upper()

    def is_edited(self, nuc):
        return self.inputnuc == nuc

    def get_edited(self, nuc):
        if nuc.upper() == self.inputnuc:
            return self.outputnuc
        return nuc

    def edit_sequence(self, sequence, position=[]):
        if not position:
            seq = [self.get_edited(nuc) for nuc in sequence]
        else:
            seq = []
            for (i, nuc) in enumerate(sequence):
                if i in position:
                    seq.append(self.get_edited(nuc))
                else:
                    seq.append(nuc)
        return "".join(seq)

    def __eq__(self, other):
        return (self.inputnuc, self.outputnuc) == (other.inputnuc, other.outputnuc)

    def __hash__(self):
        return hash(self.inputnuc, self.outputnuc)


class RNAediting (object):
    """Attempt to check edition with a simple algorithm"""

    def __init__(self, codonalign, aalign, gcode=1, reatype=('C', 'U'), radius=None):
        # radius expected value is 9
        self.codontable = CodonTable.unambiguous_dna_by_id[gcode]
        self.codonalign = codonalign
        self.aalign = self._to_mutable(aalign)
        self.dict_align = SeqIO.to_dict(codonalign)
        self.alen = codonalign.get_aln_length()
        self.glist = self.dict_align.keys()
        self.radius = radius
        self._compute_info_content()
        self.reatype = SubsReaType(reatype[0], reatype[1])
        self.possible_edited_pos = {}
        self.codongraph = CodonGraph(self.subsrea, self.codontable)
        self.id_map = self._id_to_index()

    def _to_mutable(self, alignment):
        # getting a copy of the alignment
        aalign = alignment[:, :]
        # now convert seq to mutable
        for seqrec in self.aalign:
            seqrec.seq = seqrec.seq.tomutable()
        return aalign

    def _id_to_index(self):
        mapping = {}
        for id, seq in enumerate(self.aalign):
            mapping[seq.id] = id
        return mapping

    def get_genome_editing(self, genomes=[], thresh=0):
        plausible_edited_pos = ddict(list)
        accgenome = self.glist
        if genomes:
            accgenome = set(genomes).intersection(self.glist)
        for geno in accgenome:
            # here we should get position for the genome where we have a U
            codonseq = self.dict_align[geno]
            for pos in range(self.alen):
                codon = codonseq.get_codon(pos)
                transitions = self.codongraph[codon]
                for trans in transitions:
                    diff = self.delta_column_ic(geno, pos, trans)
                    if diff > thresh:
                        plausible_edited_pos[geno].append((pos, trans, diff))

    def delta_colum_ic(self, geno, pos, trans):
        current_ic = self.ic_cont[pos]
        col_align = self.aalign[:, pos:pos + 1]
        col_align[self.id_map[geno]].seq[0] = trans[-1]

        align_info = AlignInfo.SummaryInfo(col_align)
        new_ic = align_info.information_content(pseudo_count=1)
        return current_ic - new_ic

    def _compute_info_content(self):
        align_info = AlignInfo.SummaryInfo(self.aalign)
        align_info.information_content(pseudo_count=1)
        self.ic_cont = align_info.ic_vector

    def check_editing():
        pass

    def _cu_subs_checking(self, codon1, amino2):
        pass
