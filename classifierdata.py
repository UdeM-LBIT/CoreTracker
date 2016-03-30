
import os, sys, re, json
from coretracker.classifier import *
from coretracker.classifier import MODELPATH
from coretracker.classifier.classifier import print_data, get_sensibility_and_precision, split_zeros_pos
from Bio.Data import CodonTable
from functools import partial
import numpy as np
from sklearn.utils import shuffle
from collections import defaultdict

metazoareafile = "/usr/local/www/CoreTracker/out/metajak/reassignment.json"
yeastreafile = "/usr/local/www/CoreTracker/out/yeast/reassignment.json"
datafile = "/usr/local/www/CoreTracker/input/classdata.txt.csv"

def train_load(Xlab, n, genetic_code):
    codontable = CodonTable.unambiguous_dna_by_id[genetic_code]
    train_list = np.zeros(n)
    Y = np.zeros(n)
    speclist = []
    datadict = makehash(depth=3, type=int)
    with open(datafile, 'r') as dfile:
        for line in dfile:
            data =  line.strip().split(',')
            genome = data[0]
            speclist.append(re.compile(genome))
            dlen = len(data)
            i = 1
            while i<len(data) and data[i]:
                codon, aa =  data[i:i+2]
                i+=2
                codon = codon.upper().replace('U', 'T')
                if len(codon) == 3 and codontable.forward_table[codon] != aa.upper():
                    datadict[genome][codon][aa] = 1
    for i in xrange(n):
        g, codon, ori, rea = Xlab[i]
        gkey = [x for x in speclist if x.match(g)]
        if gkey:
            gkey = gkey[0].pattern
            success = 0
            try:
                success =  datadict[gkey][codon].get(rea,0)
            except:
                pass
            train_list[i] = 1
            if success>0:
                Y[i] = 1
    return train_list==1, Y

def makehash(depth, type):
    if depth == 0:
        return defaultdict(type)
    else:
        return defaultdict(partial(makehash, depth-1, type))

def parse_input_file(file, genetic_code):
    train_list = makehash(depth=3, type=int)
    rest_of_list = makehash(depth=2, type=list)
    codontable = CodonTable.unambiguous_dna_by_id[genetic_code]
    with open(file, 'r') as FILE:
        for line in FILE:
            if not line.startswith('#'):
                sl = line.rstrip().split("\t")
                genome = ""
                codon = ""
                ori_aa = ""
                rea_aa = ""
                status = 0;
                comment = "";
                color = "";
                try:
                    genome = sl[0]
                    codon = sl[1];
                    ori_aa = sl[2]
                    rea_aa = sl[3]
                    status = sl[4]
                    color = sl[5]
                    comment = sl[6]
                except:
                    pass
                codon = codon.upper().replace('U', 'T')
                if len(codon) == 3 and codontable.forward_table[codon] != rea_aa.upper():

                    if comment == 'VN':
                        train_list[genome][codon][rea_aa] = 0

                    elif comment == 'VP':# and rea_aa.upper() != 'M':
                        train_list[genome][codon][rea_aa] = 1

                    else:
                        # add anything else to this dict
                        rest_of_list[genome][codon][ori_aa].append([rea_aa, status, color, comment])

    return train_list, rest_of_list

def test_and_print(Xtest, Ytest, Xlab, clf):
    # get result information
    Ypred = clf.predict(Xtest)
    c.get_stat(Xtest, Ytest)
    get_sensibility_and_precision(Ypred, Ytest, Xlab, Xtest)

def continuous_classification(c, X_train, Y_train, Xlab_train, X_test, Y_test ):
    true_data = X_train[meta_Y_train==1, :]
    false_data = X_train[meta_Y_train==0, :]
    false_datalabel = Xlab_train[Y_train==0, :]
    false_label = Y_train[Y_train==0]

    ssize = true_data.shape[0]
    tmp_false = true_data
    tmp_false_lab = Y_train[Y_train==1]
    i = 1
    for tmpx, _, tmpy in split_zeros_pos(false_datalabel, false_data, false_label, split_size=ssize):
        if tmp_false is None:
            tmp_false = tmpx
            tmp_false_lab = tmpy
        else:
            tmp_false = np.vstack((tmp_false, tmpx))
            tmp_false_lab = np.concatenate((tmp_false_lab, tmpy))

        cur_x, cur_y = shuffle(tmp_false, tmp_false_lab, random_state=12345)
        c.train(cur_x, cur_y)
        print("%d, Accuracy = %f"%(i, c.get_score(X_test, Y_test)))
        i += 1

if __name__ == '__main__':

    genetic_code = 1
    use_pca = False

    etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac", "N. rea", "N. used", "Cod. count", "Sub. count", "G. len", "codon_lik", "N. mixte" ,"id"]
    selected_feats = [0,2,3,4,5,6,7,8,9,11]

    meta_X, meta_Xlab, _ = read_from_json(metazoareafile, None, use_global=False)
    n = meta_X.shape[0]
    meta_train_list, meta_Y = train_load(meta_Xlab, n, genetic_code)
    meta_X_train = meta_X[meta_train_list]
    meta_Y_train = meta_Y[meta_train_list]
    meta_Xlab_train = np.asarray(meta_Xlab)[meta_train_list, :]
    meta_X_train, _ = getDataFromFeatures(meta_X_train, etiquette, feats=selected_feats)


    yeast_X, yeast_Xlab, _ = read_from_json(yeastreafile, None, use_global=False)
    n = yeast_X.shape[0]
    yeast_train_list, yeast_Y = train_load(yeast_Xlab, n, genetic_code)
    yeast_X_train = yeast_X[yeast_train_list]
    yeast_Y_train = yeast_Y[yeast_train_list]
    yeast_Xlab_train = np.asarray(yeast_Xlab)[yeast_train_list, :]
    yeast_X_train, etiquette = getDataFromFeatures(yeast_X_train, etiquette, feats=selected_feats)
    #print_data(yeast_X_train, yeast_Xlab_train, yeast_Y_train, etiquette)

    # initialize classifier
    parameters = {'random_state':12345, 'max_features':None, 'oob_score':True}#, 'class_weight': 'balanced'}
    c = Classifier('rf', parameters, False)

    # concat data for test and train
    Xtrain = np.vstack((meta_X_train, yeast_X_train))
    Ytrain = np.hstack((meta_Y_train, yeast_Y_train))
    #c.cross_validation(Xtrain, Ytrain, tsize=0.3)

    #continuous_classification(c, meta_X_train, meta_Y_train, meta_Xlab_train, yeast_X_train, yeast_Y_train)
    c.train(meta_X_train, meta_Y_train)

    # train and predict
    print('====> yeast')
    test_and_print(yeast_X_train, yeast_Y_train, yeast_Xlab_train, c)

    print('\n\n====> metazoa')
    test_and_print(meta_X_train, meta_Y_train, meta_Xlab_train, c)

    c.feature_importance(outfile="importance.png", features_list=etiquette)
    c.save_model(MODELPATH%'new')
    #X_pca_train = get_features(meta_X_train,  meta_Y_train, ncomp=6)
    #draw_pca_data(X_pca_train, meta_Xlab_train, meta_Y_train, outfile="PCA.png")
