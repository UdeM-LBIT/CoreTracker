echo '====> METAZOA\n'
nice -18  coretracker -t data/input/meta_extend/species_tree.nwk -p data/input/meta_extend/prot_aligned.core -n data/input/meta_extend/nuc.core --gapfilter 0.4 --iccontent 0.3  --idfilter 0.5  --norefine --wdir data/outnew/meta_extend --debug --params param.yml --parallel  4 

echo '====> PLANT\n'
nice -18 coretracker -t data/input/plant/plant_species_tree.nwk -p data/input/plant/all_plants_prot_ali.fas  -n data/input/plant/out/nuc.core --gapfilter 0.4 --iccontent 0.3  --idfilter 0.5  --norefine --wdir data/outnew/plant --params param.yml  --parallel  4 --debug

echo '====> YEAST\n'
nice -18 coretracker -t data/input/yeast/specietree.nw -p data/input/yeast/all_yeast_prot.fas -n data/input/yeast/all_yeast_nuc.fas --gapfilter 0.4 --iccontent 0.3  --idfilter 0.5  --norefine --wdir data/outnew/yeast  --params param.yml  --parallel  4 --debug