# Output directory name
OUTDIR = "tmp/"

# Save fitered file at each step

# amino acid to exclude from the plot
# this will speed up a lot the data generation
# If you exclude an amino acid, there won't be a reassignation to that aa
EXCLUDE_AA = "ACEFHPQVWY"

EXCLUDE_AA_FROM = "ACEFHPQVWY"

# Threshold require to determine if an aa is majoritary in a column
AA_MAJORITY_THRESH = 0.5

# skip mafft alignment, use this to gain time for debug purpose
SKIP_ALIGNMENT = True

# this is a temporary mode to compute suspected species
# possible values : count, wilcoxon, any other
MODE = 'count'

# Distance matrice:
MATRIX = 'identity'

# limit substitution count and analyses to suspected species only 
LIMIT_TO_SUSPECTED_SPECIES = False
 
# Codon to amino acid likelihood. Use consensus or all alignement
USE_CONSENSUS_FOR_LIKELIHOOD = False


# Dump result into a file for the web app
JSON_DUMP = True

# Minimum frequency of global>filtered required for each specie 
# Using a random threshold is a bad idea
# It's better to Find the probability of global>filtered happening randomly and use that
# to compare 
FREQUENCY_THRESHOLD = 0.4

# Substitution must appear at least COUNT_THRESHOLD time in the global alignment to be considered
COUNT_THRESHOLD = 3

SHOW_MIXTE_CODONS = True
# Genetic code 
# standard : 1
# vertebrate mitochondrial : 2
# yeast mitochondrial : 3
# mold mitochondrial : 4
# yeast mitochondrial CUN : Leucine : -3
# Bacterial and plant plastid : 11

#GENETIC_CODE = -3
GENETIC_CODE = 4 # is in utils

# Show Filtered codon data
SHOW_FILTERED_CODON_DATA = True # is in utils

#IMAGE FORMAT
# Accepted format : svg, png and pdf
IMAGE_FORMAT = "pdf"  # is in utils

#set this to true if you cannot output color
ADD_LABEL_TO_LEAF = False # is in utils