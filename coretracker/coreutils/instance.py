from ..classifier import Classifier, MODELPATH
from ..classifier.models import ModelType
from letterconfig import aa_letters_1to3
from . import ReaGenomeFinder
import yaml

try:
    from joblib import dump, load
except ImportError:
    from sklearn.externals.joblib import dump, load


class AADumper(yaml.Dumper):

    def increase_indent(self, flow=False, indentless=False):
        return super(AADumper, self).increase_indent(flow, False)


class RunningInstance():

    def __init__(self, reafinder, args):
        self.rfinder = reafinder
        self.args = args

    @classmethod
    def from_seqsets(clc, seqsets, speclist, settings, args):
        if len(seqsets) == 0:
            raise ValueError("Expect at least one dataset")
        else:
            seq_init = seqsets[0]
            seq_init = seq_init.fusion(*seqsets[1:], tree=args.tree)
            reafinder = ReaGenomeFinder(seq_init, settings)
            del seqsets[:]
            reafinder.get_genomes()
            reafinder.possible_aa_reassignation(speclist=speclist)
            reafinder.set_rea_mapper()
            return clc(reafinder, args)

    def get_model(self, etiquette):
        clf = Classifier.load_from_file(
            MODELPATH % self.rfinder.settings.MODEL_TYPE)
        model = ModelType(self.rfinder.settings.MODEL_TYPE, etiquette)
        if clf is None or not clf.trained:
            raise ValueError("Classifier not found or not trained!")
        return clf, model

    def save_instance(self, outfile, protocol=-1):
        """Dump a rea object into an outfile"""
        rea_aa_dict = dict((aa_letters_1to3[k], map(
            lambda vv: aa_letters_1to3[vv], v.keys())) for k, v in self.rfinder.aa2aa_rea.items())
        with open(outfile + ".aa", 'w') as AA:
            yaml.dump(rea_aa_dict, AA,
                      default_flow_style=False, Dumper=AADumper)
        dump(self, outfile + ".pgz", compress=True, protocol=protocol)

    @classmethod
    def load_instance(clc, outfile):
        """Load rea object"""
        return load(outfile)
