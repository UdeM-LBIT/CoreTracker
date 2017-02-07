#!/bin/bash

echo "Downloading latest pfam release ..."
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
echo "\n.. Downloaded !"
echo "\n\nDecompressing .gz..."
gzip --verbose --decompress Pfam-A.hmm.gz
echo "\n\n Done !"
echo "\nRunning hmmpress on Pfam-A.hmm..."
hmmpress Pfam-A.hmm
echo "\n\n Done !"