#!/usr/bin/env python

from coretracker.coreutils import AncestralRecon as AR
from ete3 import Tree
import yaml
import json
import argparse
from collections import defaultdict
from Bio.Data import CodonTable

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

def load_input(infile, isjson=False):
    readict = None
    with open(infile, 'r') as IN:
        if isjson:
            readict = json.load()
        else:
            readict = yaml.load(IN, Loader=Loader)
    if not readict:
        raise ValueError("Could not load codon reassignment, empty dict")
    else:
        return readict

def load_csvs(infiles, thresh=0.5):

    #   0         1         2       3      -2           -1
    # 'genome' 'codon'  'ori_aa' 'rea_aa'  'predict'   'proba'
    readict = defaultdict(dict)
    for f in infiles:
        with open(f.strip(), 'r') as IN:
            read_until = False
            for line in IN:
                line = line.strip()
                if line.startswith('### Random'):
                    read_until = True
                elif line.startswith('#'):
                    read_until = False
                if read_until and line:
                    dt = line.strip().split()
                    genome, codon = dt[:2]
                    rea_aa =  dt[3]
                    predict = int(dt[-2]) if thresh==0.5 else float(dt[-1]) > thresh
                    if predict:
                        readict[genome][rea_aa].append(codon)
    if not readict:
        raise ValueError("Could not load codon reassignment, empty dict")
    else:
        return readict

def perform_anc_recon(readict, tree, artype="fitch", gtable={}):
    """Perform ancestral reconstruction"""
    anc = AR.FitchBased(tree, readict)
    if artype == 'dollo':
        anc = AR.DolloParsimony(tree, readict)
    else:
        raise NotImplementedError("%s is not available yet"%artype)
    stmat, stmap, codlist =  anc.make_codonrea_matrice(tree, readict, gtable )
    stmat = anc.label_internal(stmat, codlist, stmap, alphmap)
    new_rea_dict =  anc.mat_to_dict(codlist, gtable)
    return new_rea_dict


def plot_recon_tree(tree, new_rea_dict, outputfile, show_head_text=True):
    aalist = set(sum([x.keys() for x in new_rea_dict.values()],[]))
    legend_w =  40
    legend_h = 15

    def layout(node):
        if node.is_leaf():
            faces.add_face_to_node(Faces.ReaRectFace(aalist, d[node.name], margin_left=20, ncodons=6), node, column=1, position='aligned')
        else:
            faces.add_face_to_node(Faces.ReaRectFace(aalist, d[node.name], margin_left=10, ncodons=6), node, column=1, position='branch-right')
    ts = TreeStyle()
    ts.layout_fn = layout

    if show_head_text:
        headtext = Faces.List90Face(sorted(aalist), fsize=13, ftype="Arial", rotation=0, col_w=34)
        headtext.margin_left = 13
        ts.aligned_header.add_face(headtext, column=1)


    for aa in sorted(aalist):
        f =  faces.TextFace("  " + aa + ":  ", fsize=14)
        r = RectFace(w, h, "#000000" , Faces._aabgcolors[aa] ,label=aa)
        f.margin_top = 8
        r.margin_top = 8
        ts.legend.add_face(f, column=0)
        ts.legend.add_face(r, column=1)

    ts.legend_position =  4
    tree.render(outputfile, tree_style=ts, dpi=600)


if __name__ == '__main__':

    # argument parser
    parser = argparse.ArgumentParser(
        description='AnsRec, Ancestral codon reassignment ploting')

    parser.add_argument(
        '--algorithm', '-a', choices=('fitch', 'dollo', 'MML'), default='fitch', dest="algo", help="Ancestral reconstruction algorithm")

    parser.add_argument('--json', action='store_true', help="Use json instead of yaml file for '-i' option")

    input_type = parser.add_mutually_exclusive_group(required=True)
    input_type.add_argument('--input', '-i', dest='input',
                            help="Reassignment file. Either a yaml or a json file")
    input_type.add_argument('--cinput', '-c', dest='cinput',  nargs='+',
                            help="CoreTracker csv input")

    parser.add_argument('--tree', '-t', dest='tree', help="Input tree to draw on")
    parser.add_argument('--gcode', default=4, type=int, dest='gcode', help="Base genetic code to build on.")
    parser.add_argument('--out', '-o', dest="outfile", default="outfile.svg", help="Output file name, ad image extension (svg, pdf or png).")
    parser.add_argument('--json', action='store_true', help="Use json instead of yaml file for '-i' option")
    parser.add_argument('--show_head_text', action='store_true', help="Wheter we should show text (AA) on header container or not")
    args = parser.parse_args()

    tree =  Tree(args.tree)
    if args.cinput:
	readict = load_csvs(args.cinput)
    else:
        readict =  load_input(args.input, args.json)

    table = CodonTable.unambiguous_dna_by_id[abs(gcode)]
    new_rea_dict = perform_anc_recon(tree, readict, artype=args.algo, gtable=table.forward_table)
    plot_recon_tree(tree, new_rea_dict, args.outfile, args.show_head_text)
