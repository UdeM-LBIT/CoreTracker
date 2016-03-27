#!/usr/bin/env python

from ete3 import Tree
from coretracker.coreutils import CoreFile
import collections
import argparse

# argument parser
parser = argparse.ArgumentParser(
    description='Fusion, merge 2 dataset for CoreTracker')

parser.add_argument('--out', dest="output", default="merged", help="Merged output file")
parser.add_argument('--trees', nargs=2, dest='trees', help="list of trees")
parser.add_argument('--nucs', nargs=2, dest='nucs', help="list of nuc files")

args = parser.parse_args()
nuc =  CoreFile(args.nucs[0], alphabet="nuc").get_sequences()
nuc1 =  CoreFile(args.nucs[1], alphabet="nuc").get_sequences()

for (k,v) in nuc1.items():
    cnuc = nuc.get(k, None)
    if cnuc:
        nuc[k].extend(v)
    else:
        nuc[k] = v

CoreFile.write_corefile(nuc, args.output+".core")

t = Tree()
t1 = Tree(args.trees[0])
t2 = Tree(args.trees[1])
t.add_child(t1)
t.add_child(t2)
t.write(features=[], outfile=args.output+".nw")
