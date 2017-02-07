#!/usr/bin/env python

import os
import sys
import re
import json
from coretracker.classifier import *
from coretracker.classifier import MODELPATH
from coretracker.coreutils import utils
import numpy as np

dirprefix = '/usr/local/www/CoreTracker/data/out/'
gcode = 4
cltype = 'new'
fulldir = os.listdir(dirprefix)
etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac",
             "N. rea", "N. used", "Cod. count", "Sub. count",
             "G. len", "codon_lik", "N. mixte", "id"]  # , 'total_aa']
selected_feats = [2, 3, 4, 5, 6, 7, 8, 9, 11]
if cltype == 'new':
    selected_feats = [0, 2, 3, 4, 5, 6, 7, 8, 9, 11]
    selected_et = [etiquette[i] for i in selected_feats]

aa_letters_1to3 = {

    'A': 'Ala', 'C': 'Cys', 'D': 'Asp',
    'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met',
    'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp',
    'Y': 'Tyr',
}

aa_letters_3to1 = dict((x[1], x[0]) for x in aa_letters_1to3.items())

aamap = lambda x: (aa_letters_3to1[x[0]], aa_letters_3to1[x[1]])


def listrea(realist):
    for el in realist:
        if os.path.isdir(os.path.join(dirprefix, tottype, el)):
            yield (aamap(el.split('_to_'))) + (el,)


def bottom_read_til(tfile, stop='### A'):
    fillist = []
    with open(tfile) as IN:
        for line in reversed(IN.readlines()):
            fillist.insert(0, line.rstrip())
            if line.startswith(stop):
                break
    return fillist

clf = Classifier.load_from_file(MODELPATH % cltype)

for tottype in fulldir:
    reafile = os.path.join(dirprefix, tottype, "reassignment.json")
    X, Xlab, _ = read_from_json(reafile, None, use_global=False)
    Xlab = np.asarray(Xlab)
    X, _ = getDataFromFeatures(X, etiquette, feats=selected_feats)
    pred_prob = clf.predict_proba(X)
    pred = clf.predict(X)
    readetected = os.listdir(os.path.join(dirprefix, tottype))
    for (aa1, aa2, eldir) in listrea(readetected):
        outdir = os.path.join(dirprefix, tottype, eldir)
        print(outdir)
        ind = np.bitwise_and(Xlab[:, 2] == aa1, Xlab[:, 3] == aa2)
        X_c = X[ind, :]
        Xlab_c = Xlab[ind, :]
        pred_prob_c = pred_prob[ind, 1]
        pred_c = pred[ind]
        tmp = pred_prob_c[pred_c == 1]
        prerange = (0, 1)
        if len(tmp) > 0:
            prerange = (min(tmp), max(tmp))
        suptext = "\n".join(bottom_read_til(
            os.path.join(outdir, eldir + "_data.txt")))
        suptext += "\n### Prediction Prob. range for Positive :[ %.3f - %.3f ]\n" % prerange
        utils.print_data_to_txt(os.path.join(outdir, aa_letters_1to3[aa1] + "_to_" + aa_letters_1to3[aa2] + "_corr_data.txt"),
                                selected_et, X_c, Xlab_c, pred_c, pred_prob_c, None, aa2, suptext)
