# amino acid to exclude from the plot
# this will speed up a lot the data generation
# If you exclude an amino acid, there won't be a reassignation to that aa

EXCLUDE_AA = "ALCDEFHKPQSVWY"

# If you exclude a aa, there won't be a reassignation from that aa
EXCLUDE_AA_FROM = "ALCDEFHKPQSVWY"

# Threshold require to determine if an aa is majoritary in a column
AA_MAJORITY_THRESH = 0.5

# Mode to compute suspected species
# possible values : count, wilcoxon, mannwhitney, kmean, ttest
# kmean does not work well, use wilcoxon
MODE = 'wilcoxon'

# Distance matrice:
MATRIX = 'identity'

# limit substitution count and analyses to suspected species only
LIMIT_TO_SUSPECTED_SPECIES = False

# Codon to amino acid likelihood. Use consensus or complete
# alignement for Telford score
USE_CONSENSUS_FOR_LIKELIHOOD = False

# Dump all result into a json file
JSON_DUMP = False

# Minimum frequency of global>filtered required for each specie
# Using a random threshold is a bad idea
# It's better to Find the probability of global>filtered happening randomly and use that
# to compare
FREQUENCY_THRESHOLD = 0.4

# Substitution must appear at least COUNT_THRESHOLD time in the global
# alignment to be considered
COUNT_THRESHOLD = 3

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

# IMAGE FORMAT for all plot
# Accepted format : png and pdf and svg
# The report will always be in pdf format
IMAGE_FORMAT = "pdf"

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


# alpha : rejection of the null hypothesis
CONF = 0.05
