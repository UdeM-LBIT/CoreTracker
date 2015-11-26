# this is used to enable a light parallelism in some part of the code
# change Process_enabled to your convenience
# but remember : more process won't necessary make the code run faster !!
try:
	import psutil
	PROCESS_ENABLED = psutil.cpu_count()
except ImportError, e:
	PROCESS_ENABLED = 2

# Output directory name
OUTDIR = "tmp/"

# Save fitered file at each step

# amino acid to exclude from the plot
# this will speed up a lot the data generation
# If you exclude an amino acid, there won't be a reassignation to that aa
EXCLUDE_AA = "DEFGHPRSVWY"

# Threshold require to determine if an aa is majoritary in a column
AA_MAJORITY_THRESH = 0.5

# skip mafft alignment, use this to gain time for debug purpose
SKIP_ALIGNMENT = True

# this is a temporary mode to compute suspected species
MODE = 0

# limit substitution count and analyses to suspected species only 
LIMIT_TO_SUSPECTED_SPECIES = False

# Dump result into a file for the web app
JSON_DUMP = True

# Minimum frequency of global>filtered required for each specie 
# Using a random threshold is a bad idea
# It's better to Find the probability of global>filtered happening randomly and use that
# to compare 
FREQUENCY_THRESHOLD = 0.4

# Substitution must appear at least COUNT_THRESHOLD time in the global alignment to be considered
COUNT_THRESHOLD = 3

# Genetic code 
# standard : 1
# vertebrate mitochondrial : 2
# yeast mitochondrial : 3
# mold mitochondrial : 4
# yeast mitochondrial CUN : Leucine : -3
# Bacterial and plant plastid : 11

#GENETIC_CODE = -3
GENETIC_CODE = 1 # is in utils

# Show Filtered codon data
SHOW_FILTERED_CODON_DATA = True # is in utils

#IMAGE FORMAT
# Accepted format : svg, png and pdf
IMAGE_FORMAT = "png"  # is in utils

#set this to true if you cannot output color
ADD_LABEL_TO_LEAF = False # is in utils