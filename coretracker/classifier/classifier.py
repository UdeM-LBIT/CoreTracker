from __future__ import division

import itertools
import json
import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from Bio.Data import CodonTable
from scipy import misc
from sklearn import cross_validation, feature_selection, preprocessing, svm
from sklearn.cross_validation import train_test_split
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.externals import joblib
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LassoCV, LogisticRegression
from sklearn.metrics import (average_precision_score, brier_score_loss,
                             f1_score, precision_recall_curve, precision_score,
                             recall_score, roc_auc_score, roc_curve)
from sklearn.naive_bayes import GaussianNB
from sklearn.pipeline import FeatureUnion, Pipeline
from sklearn.utils import shuffle

codon_identifier = dict(("".join(v), k)
                        for k, v in enumerate(itertools.product('ATGC', repeat=3)))


class Classifier(object):
    """This is a classifier for codon reassignment"""

    def __init__(self, method, classifier_spec={}, scale=False, n_estimators=1000):
        if method == 'rf':
            self.clf = RandomForestClassifier(
                n_estimators=n_estimators, n_jobs=-1, max_leaf_nodes=1000, **classifier_spec)
        elif method == 'svc':
            self.clf = svm.SVC(probability=True, **classifier_spec)
        elif method == 'etc':
            self.clf = ExtraTreesClassifier(
                n_estimators=n_estimators, **classifier_spec)
        elif method == 'gnb':
            self.clf = GaussianNB()
        else:
            raise NotImplementedError(
                "The method you chose (%s) is not implemented" % method)
        self.method = method
        self.trained = False
        self.scale = scale

    @classmethod
    def load_from_file(clc, loadfile):
        """Load model from a file"""
        try:
            clf = joblib.load(loadfile)
            return clf
        except IOError:
            print('Problem with file %s, can not open it' % loadfile)
        except Exception as e:
            raise e
        return None

    def save_model(self, outfile):
        """Save model to a file"""
        joblib.dump(self, outfile)

    def train(self, X=None, Y=None):
        """Train the model"""
        if self.scale:
            X = preprocessing.scale(X)
        self.clf.fit(X, Y)
        self.X = X
        self.y = Y
        self.trained = True

    @classmethod
    def from_classifier(clc, clfier):
        newclf = clc(clfier.method, {}, clfier.scale)
        newclf.__dict__.update(clfier.__dict__)
        return newclf

    def get_score(self, X, Y):
        """Return score for classification on X"""
        if self.scale:
            X = preprocessing.scale(X)
        return self.clf.score(X, Y)

    def predict(self, X):
        """Predict values for X"""
        if not self.trained:
            raise ValueError("Classifier is not trained")
        if self.scale:
            X = preprocessing.scale(X)
        return self.clf.predict(X)

    def predict_proba(self, X):
        """Return probability for each class prediction"""
        return self.clf.predict_proba(X)

    def feature_importance(self, outfile="importance.png", features_list=[]):
        """Show each feature importance"""
        if (self.method in ['rf', 'etc']):
            importances = self.clf.feature_importances_
            if len(features_list) > 0 and len(features_list) != len(importances):
                raise ValueError("Number of features does not fit!")

            indices = np.argsort(importances)[::-1]
            n_feats = len(features_list)
            np.savetxt(outfile + ".txt", np.array([tree.feature_importances_
                                                   for tree in self.clf.estimators_]), delimiter=',', fmt='%1.3e')
            std = np.std(
                [tree.feature_importances_ for tree in self.clf.estimators_], axis=0)
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
            plt.savefig(outfile, bbox_inches='tight')
        else:
            raise NotImplementedError(
                "Not supported for classifier other than Ensembl Tree")

    def cross_validation(self, X, Y, X_test=None, Y_test=None, tsize=0.3):
        """Cross validation on X and Y, using a sub sample"""
        if X_test is None or Y_test is None:
            X_train, X_test, Y_train, Y_test = train_test_split(
                X, Y, test_size=tsize)
        else:
            X_train = X
            Y_train = Y
        self.train(X_train, Y_train)
        Y_predicted = self.predict(X_test)
        self.get_stat(X_test, Y_test)
        return Y_predicted, self.get_score(X_test, Y_test)

    def plot_precision_recall(self, X_test, y_test, infos="", outfile="precision_recall.png"):
        """plot precicion-recall curve"""
        if self.trained:
            try:
                y_score = self.clf.decision_function(X_test)
            except:
                y_score = self.clf.predict_proba(X_test)[:, 1]
            precision, recall, _ = precision_recall_curve(y_test, y_score)
            average_precision = average_precision_score(
                y_test, y_score, average="micro")
            # Plot Precision-Recall curve for each class
            plt.clf()
            plt.plot(recall, precision,
                     label='Average Precision-recall curve (area = {0:0.2f})'
                     ''.format(average_precision))
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.title('Precision-Recall curve for %s (%s)' %
                      (self.method, infos))
            plt.legend(loc="lower right")
            plt.savefig(outfile)
        else:
            raise ValueError("Classifier is not trained")

    def get_stat(self, X_test, y_test):
        """Print list of score for the current classifier"""
        y_pred = self.predict(X_test)
        if hasattr(self.clf, "predict_proba"):
            prob_pos = self.clf.predict_proba(X_test)[:, 1]
        else:  # use decision function
            prob_pos = self.clf.decision_function(X_test)
            prob_pos = (prob_pos - prob_pos.min()) / \
                (prob_pos.max() - prob_pos.min())

        clf_score = brier_score_loss(y_test, prob_pos)
        print("%s:" % self.method)
        print("\tBrier: %1.3f" % (clf_score))
        print("\tPrecision: %1.3f" % precision_score(y_test, y_pred))
        print("\tRecall: %1.3f" % recall_score(y_test, y_pred))
        print("\tF1: %1.3f" % f1_score(y_test, y_pred))
        print("\tROC AUC score: %1.3f\n" % roc_auc_score(y_test, prob_pos))


def read_from_json(data, labels=None, use_global=True, use_pvalue=True):
    """Parse X array from data"""
    if isinstance(data, basestring):
        with open(data) as jfile:
            data = json.load(jfile)

    if labels and isinstance(labels, basestring):
        with open(labels) as jfile2:
            labels = json.load(jfile2)
    # matrice format
    # global
    min_value = np.finfo(np.float).min
    dtype = 'global'
    if not use_global:
        dtype = 'filtered'

    fisher_type = 'pass'
    if use_pvalue:
        fisher_type = 'pval'
    X = []
    X_label = []
    Y = []
    # each entry format :
    # [fitch, suspected, gene_frac, rea_frac, used_frac, subs_count, codon_lik_for_rea_aa]

    for aa2, val in data['aa'].items():
        for aa1, glist in val.items():
            for genome, gdata in glist.items():
                type_check = gdata[dtype]
                codon_total = gdata['codons'][dtype]
                fitch = gdata['fitch']
                suspected = gdata['suspected'] < 0.05
                rea_codon = type_check['rea_codon']
                mixte_codon = type_check['mixte_codon']
                used_codon = type_check['used_codon']
                # gene_in_genome = data['genes'][genome]
                was_lost = gdata['lost'][fisher_type]
                total_aa = np.sum(codon_total.values())
                # mixte_codon = type_check['mixte_codon']
                subs_count = type_check['count']
                for codon in codon_total.keys():
                    gene_count = 0
                    total_gene_count = 0
                    try:
                        gene_count = len(type_check[
                            'rea_distribution'].get(codon, []))
                        total_gene_count = len(type_check[
                            'total_rea_distribution'].get(codon, []))
                    except:
                        gene_count = type_check[
                            'rea_distribution'].get(codon, 0)
                        total_gene_count = type_check[
                            'total_rea_distribution'].get(codon, 0)

                    codon_count = codon_total[codon]
                    try:
                        codon_lik = gdata['score'][
                            dtype].get(codon, min_value)
                        if codon_lik == np.inf:
                            codon_lik = min_value
                    except:
                        codon_lik = min_value

                    # si il y a mutation, frequence d'utilisation du codon
                    rea_frac = rea_codon.get(
                        codon, 0) * (1.0 / codon_count) if codon_count > 0 else 0
                    # frequence d'utilisation du codon dans les positions
                    # ou l'acide amine est predo
                    used_frac = used_codon.get(
                        codon, 0) * (1.0 / codon_count) if codon_count > 0 else 0
                    mixte_frac = mixte_codon.get(
                        codon, 0) * (1.0 / codon_count) if codon_count > 0 else 0
                    gene_frac = gene_count * \
                        (1.0 / total_gene_count) if total_gene_count > 0 else 0
                    codon_id = codon_identifier[codon.replace('U', 'T')]
                    genome_len = data["genome"][dtype][genome]
                    # only add codon that are reassigned else it does not
                    # make sense, right?

                    if rea_codon.get(codon, 0) > 0:
                        entry = [fitch, suspected, was_lost, gene_frac, rea_frac, used_frac,
                                 codon_count * 1.0 / total_aa, subs_count * 1.0 / total_aa, genome_len,
                                 codon_lik, mixte_frac, codon_id]
                        X.append(entry)
                        X_label.append([genome, codon, aa2, aa1])
                        codon_mapper = None
                        if labels:
                            try:
                                codon_mapper = labels[genome][codon]
                            except:
                                pass
                            if codon_mapper is None:
                                Y.append(-1)
                            else:
                                Y.append(codon_mapper.get(aa1, -1))
                                # keep only what is well defined
                                # anything else will have a class of -1 to
                                # mean unsure
    if labels:
        assert (len(X) == len(
            Y)), "We should not have different length for data and for input"
    return np.array(X), np.asarray(X_label), np.array(Y)


def get_labels_from_csvfile(csvfile, genetic_code):
    """Get labels from a csvfile. Won't be used anymore"""
    def makehash():
        return defaultdict(makehash)
    codontable = CodonTable.unambiguous_dna_by_id[genetic_code]
    labels = makehash()
    with open(csvfile) as labfile:
        for line in labfile:
            line = line.strip()
            if line and not line.startswith('#'):
                genome, codon, reassignment = line.split()
                codon = codon.upper().replace('U', 'T')
                if len(codon) == 3 and codontable.forward_table[codon] != reassignment.upper():
                    labels[genome][codon] = reassignment
    return labels


def get_2D_distinct(Xdata, Xlabel, y, etiquette, outfile="2Dcompare.png", features=[]):
    """Project features in 2D to check data separation"""
    color = np.empty(y.shape[0], dtype="S7")
    color[y == 0] = '#dddddd'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'G'))] = '#00ff00'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'N'))] = '#000080'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'M'))] = '#FF55A3'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'S'))] = '#FFD800'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'K'))] = '#00ffd8'
    indices = np.argsort(y)
    Xdata = Xdata[indices, :]
    y = y[indices]
    ncomp = len(features)
    if ncomp == 0 or ncomp > len(etiquette):
        ncomp = len(etiquette)
        features = range(len(etiquette))
    else:
        Xdata = Xdata[:, features]

    total_size = int(misc.comb(ncomp, 2))

    i = int(np.floor(np.sqrt(total_size)))
    j = int(np.ceil(total_size / i))

    plt.close('all')
    f, axarr = plt.subplots(i, j)
    for xax in xrange(ncomp):
        for yax in xrange(xax + 1, ncomp):
            total_size -= 1
            i, j = np.unravel_index(total_size, axarr.shape)
            axarr[i, j].scatter(Xdata[:, xax], Xdata[:, yax], c=color)
            axarr[i, j].set_xlabel(etiquette[features[xax]], fontsize=6)
            axarr[i, j].set_ylabel(etiquette[features[yax]], fontsize=6)
    for axe in np.ravel(axarr):
        axe.tick_params(axis='both', which='both', bottom='off', labelbottom='off',
                        labeltop='off', top='off', right='off', left='off', labelleft='off')
    plt.tight_layout()
    plt.savefig(outfile)


def get_features(Xdata, y=None, ncomp=2, kbest=0):
    """Feature selection using PCA or Kbest variance selection"""
    if ncomp > 0 and kbest > 0:
        pca = PCA(n_components=ncomp)
        selection = SelectKBest(f_classif, k=(
            int(kbest) if int(kbest) < Xdata.shape[1] else 'all'))
        combined_features = FeatureUnion(
            [("pca", pca), ("univ_select", selection)])
        X_features = combined_features.fit_transform(Xdata, y)

    elif ncomp > 0:
        pca = PCA(n_components=ncomp)
        X_features = pca.fit_transform(Xdata, y)

    elif kbest > 0:
        selection = SelectKBest(k=int(kbest) if int(
            kbest) < Xdata.shape[1] else 'all')
        X_features = selection.fit_transform(Xdata, y)

    return X_features


def draw_pca_data(X_features, Xlabel, y, outfile="PCA.png"):
    """ Draw pca data and save in a file"""
    color = np.empty(y.shape[0], dtype="S7")
    color[y == 0] = '#dddddd'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'G'))] = '#00ff00'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'N'))] = '#000080'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'M'))] = '#FF55A3'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'S'))] = '#FFD800'
    color[np.logical_and((y == 1), (Xlabel[:, 3] == 'K'))] = '#00ffd8'

    ncomp = X_features.shape[1]
    total_size = int(misc.comb(ncomp, 2))
    i = int(np.floor(np.sqrt(total_size)))
    j = int(np.ceil(total_size / i))

    plt.close('all')
    f, axarr = plt.subplots(i, j)
    if total_size > 1:
        for xax in xrange(ncomp):
            for yax in xrange(xax + 1, ncomp):
                total_size -= 1
                i, j = np.unravel_index(total_size, axarr.shape)
                axarr[i, j].scatter(X_features[:, xax],
                                    X_features[:, yax], c=color)
                axarr[i, j].set_title('%d vs %d' % (xax, yax), fontsize=6)

        for ax in axarr.flatten():
            ax.tick_params(axis='both', which='both', bottom='off', labelbottom='off',
                           labeltop='off', top='off', right='off', left='off', labelleft='off')
    else:
        axarr.scatter(X_features[:, 0],
                      X_features[:, 1], c=color)
        axarr.set_title('%d vs %d' % (1, 2))
        axarr.tick_params(axis='both', which='both', bottom='off', labelbottom='off',
                          labeltop='off', top='off', right='off', left='off', labelleft='off')
    plt.tight_layout()
    plt.savefig(outfile)


def print_data(X, X_label, Y, etiquette=None):
    """Print data"""
    if etiquette is None:
        etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac", "N. rea",
                     "N. used", "Cod. count", "Sub. count", "G. len", "codon_lik", "N. mixte", "id"]
    etiquette = list(etiquette)

    print("\n" + "\t".join(["genome", "codon",
                            "ori_aa", "rea_aa"] + etiquette))
    for i in xrange(len(X_label)):
        if Y[i] == 1:
            print("\t".join(list(X_label[i]) + [str(x) for x in X[i]]))


def getDataFromFeatures(Xdata, etiquette, feats=[]):
    """Extract Data based on list of features"""
    if len(feats) == 0:
        return Xdata, etiquette
    else:
        return Xdata[:, feats], np.asarray(etiquette)[feats]


def get_sensibility_and_precision(pred_y, true_y, X_labels=None, X=None, log=True):
    """Get sensibility and precision after classification"""
    nel = len(true_y)
    assert nel == len(pred_y), 'Vector should be the same size\n'
    true_pos, true_neg, false_pos, false_neg = 0.0, 0.0, 0.0, 0.0
    false_neg_list, false_pos_list = [], []
    for i in xrange(len(pred_y)):
        if pred_y[i] == 0 and true_y[i] == 1:
            false_neg += 1
            false_neg_list.append(i)
        elif pred_y[i] == 1 and true_y[i] == 0:
            false_pos += 1
            false_pos_list.append(i)
        elif pred_y[i] == 1 and true_y[i] == 1:
            true_pos += 1
        elif pred_y[i] == 0 and true_y[i] == 0:
            true_neg += 1

    print("Test size is: %d\nTrue Positive is: %d\nTrue negative is: \
          %d\nFalse positive is: %d\nFalse negative is:%d" % (
        nel, true_pos, true_neg, false_pos, false_neg))
    print('-------------------------------------------')
    print("Sensibility is %f" % (true_pos / (true_pos + false_neg)
                                 if (true_pos + false_neg) > 0 else 1))
    print("Specificity is %f" % (true_neg / (true_neg + false_pos)
                                 if (true_neg + false_pos) > 0 else 1))
    print("Accuracy is %f" % ((true_neg + true_pos) / nel))
    print("Precision is %f\n\n" %
          (true_pos / (true_pos + false_pos) if (true_pos + false_pos) > 0 else 1))
    if log:
        if X_labels is not None and X is not None:
            if len(false_neg_list) > 0:
                print("List of false negatives")
                for i in false_neg_list:
                    print("\t".join(X_labels[i]))
                    print("\t".join([str(x) for x in X[i]]))

            if len(false_pos_list) > 0:
                print("\nList of False positives")
                for i in false_pos_list:
                    print("\t".join(X_labels[i]))
                    print("\t".join([str(x) for x in X[i]]))


def split_zeros_pos(L, X, Y, split_size=300):
    """Split data into training and test set"""
    zero_pos = np.where(Y == 0)[0]
    np.random.shuffle(zero_pos)
    tlen = len(zero_pos)
    nsplit = np.floor(tlen / split_size)
    if tlen <= split_size:
        yield (X[zero_pos], L[zero_pos], Y[zero_pos])
    else:
        delim = tlen - tlen % split_size
        last_chunck = zero_pos[delim:]
        chuncks = np.split(zero_pos[:delim], nsplit)
        chuncks.append(last_chunck)
        for ch in chuncks:
            yield (X[ch], L[ch], Y[ch])


def get_aa_cross_val(L, X, Y, AA, tsize=None, rstate=-1):
    """Get test data from dataset"""
    test_position = []
    aa_y = np.zeros(Y.shape)
    for i in xrange(len(Y)):
        if L[i][-1] == AA:
            aa_y[i] = 1
            test_position.append(i)

    if tsize:
        t_len = int(tsize * len(Y))
        # positions that are 0 without being the one for AA
        zero_pos = np.where(np.logical_and(Y == 0, aa_y == 0))[0]
        clen = t_len - len(test_position)
        if clen > 0:
            random_zero_pos = np.random.choice(zero_pos, clen, replace=False)
            test_position.extend(random_zero_pos)

    test_position = np.random.permutation(test_position)
    mask = np.ones(Y.shape, dtype=bool)
    mask[test_position] = False
    train_position = np.array(range(len(mask)))[mask]

    if rstate > 0:
        return shuffle(train_position, random_state=rstate), shuffle(test_position, random_state=rstate)
    # in this case, suppose we want only the train and test index
    else:
        return train_position, test_position
