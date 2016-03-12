import parameters
import Bio.SubsMat.MatrixInfo as MatrixInfo

class Settings():
    """Contains global settings for the current run of CoReTracker"""

    def __init__(self):
        pass

    def set(self, **kwargs):
        """Set settings values from a dict"""

        EXCLUDE_AA = kwargs.get('EXCLUDE_AA', parameters.EXCLUDE_AA)
        # list of aa to exclude
        self.EXCLUDE_AA_FROM = kwargs.get(
            'EXCLUDE_AA_FROM', parameters.EXCLUDE_AA_FROM)
        # list of accepted aa
        self.AA_LETTERS = "".join(
            [aa for aa in "ACDEFGHIKLMNPQRSTVWY" if aa not in EXCLUDE_AA])
        # Set output directory
        self.OUTDIR = kwargs.get('OUTDIR', parameters.OUTDIR)
        self.AA_MAJORITY_THRESH = kwargs.get(
            # Set min frequency per column to use an aa as the most predominant
            'AA_MAJORITY_THRESH', parameters.AA_MAJORITY_THRESH)
        # whether or not analysis should be restricted to suspected species
        self.LIMIT_TO_SUSPECTED_SPECIES = kwargs.get(
            'LIMIT_TO_SUSPECTED_SPECIES', parameters.LIMIT_TO_SUSPECTED_SPECIES)
        # old method of finding suspected species by count
        self.FREQUENCY_THRESHOLD = kwargs.get(
            'FREQUENCY_THRESHOLD', parameters.FREQUENCY_THRESHOLD)
        # minimum aa substitution to consider reassignment as validation
        # this is use to filter out false positive
        self.COUNT_THRESHOLD = kwargs.get(
            'COUNT_THRESHOLD', parameters.COUNT_THRESHOLD)
        # default genetic code to use
        self.GENETIC_CODE = kwargs.get(
            'GENETIC_CODE', parameters.GENETIC_CODE)
        # whether or not codon in mixte position should be shown in the report
        self.SHOW_MIXTE_CODONS = kwargs.get(
            'SHOW_MIXTE_CODONS', parameters.SHOW_MIXTE_CODONS)
        # whether or not filtered position should be shown
        self.SHOW_GLOBAL_CODON_DATA = kwargs.get(
            'SHOW_GLOBAL_CODON_DATA', parameters.SHOW_GLOBAL_CODON_DATA)
        # output format. Should be pdf for the moment
        self.IMAGE_FORMAT = kwargs.get(
            'IMAGE_FORMAT', parameters.IMAGE_FORMAT)
        # Add a suffix to each leaf if colors is not available
        self.ADD_LABEL_TO_LEAF = kwargs.get(
            'ADD_LABEL_TO_LEAF', parameters.ADD_LABEL_TO_LEAF)
        # Add number inside the piechart
        self.ADD_NUMBER_PIE = kwargs.get(
            'ADD_NUMBER_PIE', parameters.ADD_NUMBER_PIE)
        # whether or not consensus should be used for the likelihood
        self.USE_CONSENSUS_FOR_LIKELIHOOD = kwargs.get(
            'USE_CONSENSUS_FOR_LIKELIHOOD', parameters.USE_CONSENSUS_FOR_LIKELIHOOD)
        # if global alignment should be used to find the suspected list
        self.USE_GLOBAL = kwargs.get(
            'USE_GLOBAL', parameters.USE_GLOBAL)
        # choose algorithm for computing the suspected species
        self.mode = kwargs.get('MODE', parameters.MODE)
        # choose matrix type to use, default is blosum62
        self.method = kwargs.get('MATRIX', parameters.MATRIX)
        # The following are the binaries setting for HMMER package
        self.hmmbuild = kwargs.get('hmmbuild', 'hmmbuild')
        self.eslalimask = kwargs.get('eslalimask', 'esl-alimask')
        self.eslalimanip = kwargs.get('eslalimanip', 'esl-alimanip')
        self.hmmalign = kwargs.get('hmmalign', 'hmmalign')
        # alpha to use , default is 0.05
        self.conf = kwargs.get('CONF', parameters.CONF)
        self.SHOW_ALL = False
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
