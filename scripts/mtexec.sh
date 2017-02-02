echo '====> METAZOA\n'
nice -15  python coretracker.py -t data/input/meta_extend/species_tree.nw -p data/input/meta_extend/prot_aligned.core -n data/input/meta_extend/nuc.core --gapfilter 0.4 --iccontent 0.3  --idfilter 0.5  --norefine --wdir data/outnew/meta_extend2 --debug --params param.yml --parallel 6  --expos

