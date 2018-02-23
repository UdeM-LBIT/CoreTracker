# amino acid to exclude from the plot
# this will speed up a lot the data generation
# If you exclude an amino acid, there won't be a reassignation to that aa

EXCLUDE_AA = ""

# If you exclude a aa, there won't be a reassignation from that aa
EXCLUDE_AA_FROM = ""

# Threshold require to determine if an aa is majoritary in a column
AA_MAJORITY_THRESH = 0.5

# Mode to compute suspected species
# possible values : count, wilcoxon, mannwhitney, kmean, ttest
# kmean does not work well, use wilcoxon
MODE = 'wilcoxon'

# Distance matrice :
# possible values :  identity or any of the biopython available matrices
# see https://web.archive.org/web/19991014010917/http://www.embl-heidelberg.de/~vogt/matrices/mlist1.html
MATRIX = 'blosum62'

# Do not attempt to filter results by limiting
# analysis to more likely reassignment
SHOW_ALL = False

# Limit prediction to suspected species for each reassignment
LIMIT_TO_SUSPECTED_SPECIES = False

# Codon to amino acid likelihood. Use consensus or complete
# alignement for Telford score
USE_CONSENSUS_FOR_LIKELIHOOD = False

# Save filtered alignment
SAVE_ALIGN = True

# skip empty results
SKIP_EMPTY = True

# Minimum frequency of global>filtered required for each specie
# Using a random threshold is a bad idea
# It's better to Find the probability of global>filtered happening randomly and use that
# to compare
FREQUENCY_THRESHOLD = 0.4

# Substitution must appear at least COUNT_THRESHOLD time in the global
# alignment to be considered
COUNT_THRESHOLD = 2

# Codon must appear at least CODON_COUNT_THRESHOLD time in the filtered/global alignment
# to be considered for codon reassignment prediction
CODON_COUNT_THRESHOLD = 1

# Show codon in non conserved position for its amino acid
SHOW_MIXTE_CODONS = True

# Genetic code
# standard : 1
# vertebrate mitochondrial : 2
# yeast mitochondrial : 3
# mold mitochondrial : 4
# Bacterial and plant plastid : 11
GENETIC_CODE = 4

# Show Filtered codon data
SHOW_GLOBAL_CODON_DATA = False  # is in utils

# Whether or not Label should be added to leaves
ADD_LABEL_TO_LEAF = False

# Add the total count in each pie chart in the
# codon detail figure
ADD_NUMBER_PIE = False

# Use global alignment to compute distance matrice
# instead of filtered Alignment
# not recommended
USE_GLOBAL = False

# alpha : rejection of the null hypothesis
CONF = 0.05

# codon start range
STARTDIST = 20

# Learning model to use for prediction
MODEL_TYPE = '3'

# Number of iteration for the hmm
HMMLOOP = 10

# input sequence format
PROTFORMAT = 'core'
DNAFORMAT = 'core'

#

