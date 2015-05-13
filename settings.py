# this is used to enable a light parallelism in some part of the code
# change Process_enabled to your convenience
# but remember : more process won't necessary make the code run faster !!
try:
	import psutil
	PROCESS_ENABLED = psutil.cpu_count()
except ImportError, e:
	PROCESS_ENABLED = 2

# set up a temporary directory
TMP = "tmp/"

# Mafft output file 
MAFFT_OUTPUT = TMP+"alignment.fasta"

# Provide a multiple alignment sequence to compute the Accepted Replacement Matrix.
# This will be provided as additionnal input to compute everything related to substitution matrix
# and ic content
ARM_SEQUENCE = None

# Set this to decide whether or not, coretracker should perform a new alignement at each
# node of the tree (when we calculate the substitution matrix)
REALIGN_AT_EACH_NODE = False

# Use this to enable or disable the use of expected frequence (computed from the ARM_SEQUENCE)
# when the ic_information per site is computed
USE_EXPECTED_FREQ_FOR_IC = False

# amino acid to exclude from the plot
# this will speed up a lot the data generation
EXCLUDE_AA = ""

# value to use for ic information threshold (in percent (%), the true threshold will then be determined
# by : max(IC_INFO_VECTOR)*(IC_INFO_THRESHOLD/ 100.0))
IC_INFO_THRESHOLD = 0.5

# Threshold require to determine if an aa is majoritary in a column
AA_MAJORITY_THRESH = 0.5

# skip mafft alignment, use this to gain time for debug purpose
SKIPMAFFT = True

# skip substitution matrix computing
SKIPSUBMATRIX = True

# limit substitution count and analyses to suspected species only 
LIMIT_TO_SUSPECTED_SPECIES = False

# Counter threshold, using a random threshold is a bad idea
# Find the probability of global>filtered happening randomly and use that
# to compare 
COUNTER_THRESHOLD = 0.1

# Substitution must appear at least COUNT_THRE time in the global alignment to be considered
COUNT_THRE = 2