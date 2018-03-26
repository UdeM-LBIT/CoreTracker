import parameters
import Bio.SubsMat.MatrixInfo as MatrixInfo

AVAILABLE_MAT = MatrixInfo.available_matrices + ['identity']


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
        self.AA_MAJORITY_THRESH = kwargs.get(
            # Set min frequency per column to use an aa as the most predominant
            'AA_MAJORITY_THRESH', parameters.AA_MAJORITY_THRESH)
        # whether or not analysis should be restricted to suspected species
        self.LIMIT_TO_SUSPECTED_SPECIES = False
        try:
            self.LIMIT_TO_SUSPECTED_SPECIES = kwargs.get(
                'LIMIT_TO_SUSPECTED_SPECIES', parameters.LIMIT_TO_SUSPECTED_SPECIES)
        except:
            pass
        # old method of finding suspected species by count
        self.FREQUENCY_THRESHOLD = kwargs.get(
            'FREQUENCY_THRESHOLD', parameters.FREQUENCY_THRESHOLD)
        # minimum aa substitution to consider reassignment as validation
        # this is use to filter out false positive
        self.COUNT_THRESHOLD = kwargs.get(
            'COUNT_THRESHOLD', parameters.COUNT_THRESHOLD)
        # minimum codon presence in reassigned position to make prediction
        self.CODON_COUNT_THRESHOLD = kwargs.get(
            'CODON_COUNT_THRESHOLD', parameters.CODON_COUNT_THRESHOLD)
        # default genetic code to use
        self.GENETIC_CODE = kwargs.get(
            'GENETIC_CODE', parameters.GENETIC_CODE)
        # whether or not codon in mixte position should be shown in the report
        self.SHOW_MIXTE_CODONS = kwargs.get(
            'SHOW_MIXTE_CODONS', parameters.SHOW_MIXTE_CODONS)
        # whether or not filtered position should be shown
        self.SHOW_GLOBAL_CODON_DATA = kwargs.get(
            'SHOW_GLOBAL_CODON_DATA', parameters.SHOW_GLOBAL_CODON_DATA)
        # Add a suffix to each leaf if colors is not available
        self.ADD_LABEL_TO_LEAF = kwargs.get(
            'ADD_LABEL_TO_LEAF', parameters.ADD_LABEL_TO_LEAF)
        # Add number inside the piechart
        self.ADD_NUMBER_PIE = kwargs.get(
            'ADD_NUMBER_PIE', parameters.ADD_NUMBER_PIE)
        # whether or not consensus should be used for the likelihood
        self.USE_CONSENSUS_FOR_LIKELIHOOD = kwargs.get(
            'USE_CONSENSUS_FOR_LIKELIHOOD', parameters.USE_CONSENSUS_FOR_LIKELIHOOD)
        # save filtered alignment
        self.SAVE_ALIGN = kwargs.get('SAVE_ALIGN', parameters.SAVE_ALIGN)
        # do not output results for no predicted reassignment
        self.SKIP_EMPTY = kwargs.get('SKIP_EMPTY', parameters.SKIP_EMPTY)
        # if global alignment should be used to find the suspected list
        self.USE_GLOBAL = kwargs.get(
            'USE_GLOBAL', parameters.USE_GLOBAL)
        # hidden parameter for debug purpose
        self.MODEL_TYPE = kwargs.get(
            'MODEL_TYPE', parameters.MODEL_TYPE)
        # hmm loop
        self.HMMLOOP = kwargs.get(
            'HMMLOOP', parameters.HMMLOOP)
        # choose algorithm for computing the suspected species
        self.MODE = kwargs.get('MODE', parameters.MODE)

        # choose matrix to compute substitution, default is blosum62
        self.MATRIX = kwargs.get('MATRIX', parameters.MATRIX)
        if self.MATRIX not in AVAILABLE_MAT:
            self.MATRIX = "blosum62"
        try:
            self.SUBMAT = getattr(MatrixInfo, self.MATRIX)
        except:
            self.SUBMAT = getattr(MatrixInfo, "blosum62")

        # alpha to use , default is 0.05
        self.CONF = kwargs.get('CONF', parameters.CONF) or 0.05
        self.STARTDIST = kwargs.get('STARTDIST', parameters.STARTDIST)
        self.SHOW_ALL = kwargs.get('SHOW_ALL', parameters.SHOW_ALL)
        self.FORCED_CHI2 = kwargs.get('FORCED_CHI2', parameters.FORCED_CHI2)
        # output format. Should be pdf for the moment
        self.IMAGE_FORMAT = "pdf"
        # The following are the binaries setting for HMMER package
        self.hmmbuild = kwargs.get('HMMBUILD', 'hmmbuild')
        # self.eslalimask = kwargs.get('eslalimask', 'esl-alimask')
        # self.eslalimanip = kwargs.get('eslalimanip', 'esl-alimanip')
        # this can be used next version for a better filtering
        # not need right now
        self.hmmalign = kwargs.get('HMMALIGN', 'hmmalign')
        # default values for supplemental data
        self.VALIDATION = True
        self.COMPUTE_POS = False
        self.SCALE = 1.0
        self.PROTFORMAT = kwargs.get('PROTFORMAT', parameters.PROTFORMAT)
        self.DNAFORMAT = kwargs.get('DNAFORMAT', parameters.DNAFORMAT)
        self.OUTDIR = ""

    def get_external_binaries(self):
        ext = [self.hmmbuild, self.hmmalign]
        return ext

    def fill(self, params):
        if isinstance(params, dict):
            self.set(**params)
        else:
            self.__dict__.update(params.__dict__)

    def update_params(self, **kwargs):
        for k, v in kwargs.items():
            if k == 'MATRIX' and v in AVAILABLE_MAT:
                self.__dict__[k] = v
                if v != "identity":
                    self.__dict__['MATRIX'] = getattr(MatrixInfo, v)
            else:
                self.__dict__[k] = v
