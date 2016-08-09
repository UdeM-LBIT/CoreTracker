from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

class RNAediting (object):
    """Attempt to check edition with a simple algorithm"""
    def __init__(mtalignment, genome, dt, radius,gcode=1,):
        self.codontable = CodonTable.unambiguous_dna_by_id[gcode]

    def check_editing():
        pass

    def _cu_subs_checking(self, codon1, amino2):
        pass
