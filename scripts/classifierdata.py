
import os, sys, re, json
from coretracker.classifier import *
from coretracker.classifier import MODELPATH
from coretracker.classifier.classifier import print_data, get_sensibility_and_precision, split_zeros_pos, get_aa_cross_val
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from Bio.Data import CodonTable
from functools import partial
import numpy as np
from sklearn.utils import shuffle
import seaborn as sns
sns.set()
from collections import defaultdict, Counter
import itertools
from matplotlib import pyplot as plt
from sklearn.cross_validation import cross_val_score
from sklearn.metrics import f1_score, roc_auc_score, accuracy_score

from unbalanced_dataset import UnderSampler

metazoareafile = "/usr/local/www/CoreTracker/data/out/metazoa/reassignment.json"
yeastreafile = "/usr/local/www/CoreTracker/data/out/yeast/reassignment.json"
datafile = "/usr/local/www/CoreTracker/data/input/classdata.txt.csv"

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

def test_and_print(Xtest, Ytest, Xlab, clf, seuil=0.5):
    # get result information
    print("Threshold used : %f"%seuil)
    Ypred = clf.predict_proba(Xtest)
    Ypred = (Ypred[:, -1]>seuil).astype(int)
    clf.get_stat(Xtest, Ytest)
    get_sensibility_and_precision(Ypred, Ytest, Xlab, Xtest)

def onhotencode(X, totsize=64):
    last_row = (X[:, -1]).astype(int)
    rlen = len(last_row)
    enc = np.zeros((rlen, totsize))
    ind = np.ravel_multi_index([np.arange(rlen), last_row], (rlen, totsize))
    enc.flat[ind] = 1
    return np.column_stack((X[:,:-1], enc))
    return X

def undersampling(x, y, ratio=15, xlab=None, xtest=None, ytest=None, ylab=None, lab=None, prefix='', stype="under"):
    # 'Random under-sampling'

    parameters = {'random_state':12345, 'max_features':None, 'oob_score':True}#, 'class_weight': 'balanced'}
    c = Classifier('rf', parameters, False)
    xtrain = x[:, :-1]
    # undersampler
    US = UnderSampler(ratio=ratio, verbose=False ,random_state=12345)
    usx, usy = US.fit_transform(x, y)
    if not xlab:
        xlab = usx[:, -1]
        usx = usx[:, :-1]
    seuil = threshold_tuning(usx, usy, xlab, title=prefix+'undersample')
    s2 = threshold_tuning(usx, usy, xlab, metric=roc_auc_score, title=prefix+'undersample_roc')
    s2 = threshold_tuning(usx, usy, xlab, metric=accuracy_score, title=prefix+'undersample_acc')
    tree = estimator_tree_tuning(usx, usy, pvalidator=xlab, title=prefix+'OverSampler_trees')
    if xtest is not None:
        c.train(usx, usy)
        print('====> undersample')
        print("------------Result on yeast")
        test_and_print(xtest, ytest, ylab, c,seuil)
        print("\n------------Result on metazoa")
        test_and_print(xtrain, y, lab, c,seuil)
    return c

def cutoff_predict(clf, X, cutoff):
    return (clf.predict_proba(X)[:,1]>cutoff).astype(int)

def custom_f1(cutoff, metric=f1_score):
    def f1_cutoff(clf, X, y):
        ypred = cutoff_predict(clf, X, cutoff)
        return metric(y,ypred)
    return f1_cutoff


def threshold_tuning(x, y, pvalidator, title='None', metric=f1_score):
    parameters = {'random_state':12345, 'max_features':None, 'oob_score':True}#, 'class_weight': 'balanced'}
    xaxis = []
    yaxis = []
    mean_score = []
    thresh_cutter = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    print(title)
    if not isinstance(pvalidator, int):
        pvalidator =  cross_val_data(x, y, pvalidator)
    for cutoff in thresh_cutter:
        clf = RandomForestClassifier(n_estimators=800, n_jobs=-1, max_leaf_nodes=1000, **parameters)
        validated = cross_val_score(clf, x, y, cv=pvalidator, scoring=custom_f1(cutoff, metric))
        xaxis.append(cutoff)
        yaxis.append(validated)
        mean_score.append(np.mean(validated))
    visualize_plot(xaxis, yaxis, 'threshold', 'f1_score (%f,%f)'%x.shape, title)
    seuil = np.argmax(mean_score)
    return xaxis[seuil]


def estimator_tree_tuning(x, y, minv=100, maxv=1500, step=100, pvalidator=10, scoret='f1', title="None"):
    # initialize classifier
    parameters = {'random_state':12345, 'max_features':None, 'oob_score':True}#, 'class_weight': 'balanced'}
    xaxis = []
    yaxis = []
    percentage = 1- (np.count_nonzero(y))*1.0/len(y)
    if not isinstance(pvalidator, int):
        pvalidator =  cross_val_data(x, y, pvalidator)
    for val in range(minv, maxv, step):
        clf = RandomForestClassifier(n_estimators=val, n_jobs=-1, max_leaf_nodes=1000, **parameters)
        validated = cross_val_score(clf, x, y, cv=pvalidator, scoring=scoret)
        xaxis.append(val)
        yaxis.append(validated)

    visualize_plot(xaxis, yaxis, 'number of trees', 'score', title, percentage)

def visualize_plot(xaxis, yaxis, xlab, ylab, title, hline=None):
    sns.boxplot(data=yaxis)
    if hline:
        plt.axhline(y=hline, ls='--')
    plt.xlabel(xlab)
    plt.xticks(range(len(xaxis)), [str(a) for a in xaxis])
    plt.ylabel(ylab)
    plt.title('Classification with %s'%title)
    plt.savefig('data/images/'+title+'.png')
    plt.clf()


def cross_val_data(x, y, xlab):
    aacrossval = []
    AAlist= ['K', 'S', 'M', 'G', 'N']
    for aa in AAlist:
        aacrossval.append(get_aa_cross_val(xlab, x, y, aa, tsize=0.2))
    return aacrossval


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

def oneHotFeatImport(cl, outfile="importance", features_list=[], shift_pos=9):
    """Show each feature importance"""
    if (cl.method in ['rf', 'etc']):
        importances = cl.clf.feature_importances_
        sum_imp = np.sum(importances[shift_pos:])
        importances = np.hstack((importances[:shift_pos], sum_imp))
        print(np.sum(importances))
        if len(features_list) > 0 and len(features_list) != len(importances):
            raise ValueError("Number of features does not fit!")

        indices = np.argsort(importances)[::-1]
        n_feats = len(features_list)
        np.savetxt(outfile+".txt", np.array([ np.hstack((tree.feature_importances_[:shift_pos], np.sum(tree.feature_importances_[shift_pos:]))) \
                                for tree in cl.clf.estimators_]), delimiter=',', fmt='%1.3e')
        std = np.std(
            [ np.hstack((tree.feature_importances_[:shift_pos], np.sum(tree.feature_importances_[shift_pos:]))) \
                                    for tree in cl.clf.estimators_], axis=0)
        plt.figure()
        plt.title("Feature importances")
        plt.bar(range(n_feats), importances[
                indices], width=0.5, color="b", yerr=std[indices], align="center")
        if len(features_list) > 0:
            features_list = np.asarray(features_list)[indices]
            plt.xticks(range(n_feats), features_list, rotation='vertical')
        plt.xlim([-1, n_feats])
        plt.margins(0.2)

        plt.subplots_adjust(bottom=0.15)
        plt.savefig(outfile+".svg", bbox_inches='tight')
    else:
        raise NotImplementedError(
            "Not supported for classifier other than Ensembl Tree")


if __name__ == '__main__':

    if False:
        genetic_code = 4
        use_pca = False

        etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac", "N. rea", "N. used", "Cod. count", "Sub. count", "G. len", "codon_lik", "N. mixte" ,"id"]
        selected_feats = [0,2,3,4,5,6,7,8,9,11]

        meta_X, meta_Xlab, _ = read_from_json(metazoareafile, None, use_global=False)
        n = meta_X.shape[0]
        meta_train_list, meta_Y = train_load(meta_Xlab, n, genetic_code)
        meta_X_train = meta_X[meta_train_list]
        meta_Y_train = meta_Y[meta_train_list]
        meta_Xlab_train = np.asarray(meta_Xlab)[meta_train_list, :]

        meta_X_train, et = getDataFromFeatures(meta_X_train, etiquette, feats=selected_feats)

        #aacrossval = cross_val_data(meta_X_train, meta_Y_train, meta_Xlab_train)

        meta_X_train_1hot = onhotencode(meta_X_train)
        X =  np.column_stack((meta_X_train, meta_Xlab_train[:,-1]))
        X_1hot =  np.column_stack((meta_X_train_1hot, meta_Xlab_train[:,-1]))

        # YEAST DATA SET LOADING
        yeast_X, yeast_Xlab, _ = read_from_json(yeastreafile, None, use_global=False)
        n = yeast_X.shape[0]
        yeast_train_list, yeast_Y = train_load(yeast_Xlab, n, genetic_code)
        yeast_X_train = yeast_X[yeast_train_list]
        yeast_Y_train = yeast_Y[yeast_train_list]
        yeast_Xlab_train = np.asarray(yeast_Xlab)[yeast_train_list, :]
        yeast_X_train, etiquette = getDataFromFeatures(yeast_X_train, etiquette, feats=selected_feats)
        yeast_X_train_1hot = onhotencode(yeast_X_train)
        #undersampling(X, meta_Y_train, ratio=10, lab=meta_Xlab_train, xtest=yeast_X_train, ytest=yeast_Y_train, ylab=yeast_Xlab_train, prefix="Ori_")
        print '\n\n---------------------------using one hot encoding ------------------------ \n\n'
        classifier = undersampling(X_1hot, meta_Y_train, ratio=10, lab=meta_Xlab_train, xtest=yeast_X_train_1hot, ytest=yeast_Y_train, ylab=yeast_Xlab_train, prefix="OneHot_")
        classifier.save_model(MODELPATH%'3')


    etiquette = ["Fitch", "suspected", "Fisher pval", "Gene frac.", "N. rea", "N. used", "Cod. count", "Sub. count", "G. len", "Telford", "N. mixte" ,"Codon ID"]
    selected_feats = [0,2,3,4,5,6,7,8,9,11]
    et = [etiquette[x] for x in selected_feats]
    print et
    clf =  Classifier.load_from_file(MODELPATH%'3')
    oneHotFeatImport(clf, outfile="importance", features_list=et, shift_pos=9)
