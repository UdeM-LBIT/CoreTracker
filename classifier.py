from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn import feature_selection
from sklearn import svm
from sklearn.externals import joblib
from Bio.Data import CodonTable
from sklearn import preprocessing
from sklearn import cross_validation
from sklearn.utils import shuffle
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt

from sklearn.cross_validation import train_test_split

class Classifier:
    """This is a classifier for codon rassignment"""
    def __init__(self, method, classifier_spec, scale=False):
        if method == 'rf':
            self.clf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, max_leaf_nodes=100 ,**classifier_spec)
        elif method == 'svc':
            self.clf = svm.SVC(**classifier_spec)
        elif method == 'etc':
            self.clf = ExtraTreesClassifier(n_estimators=1000, **classifier_spec)

        self.trained = False
        self.scale = scale


    def load_from_file(self, loadfile):
        try:
            self.clf = joblib.load(loadfile)
        except IOError:
            print('Problem with file %s, can not open it'%loadfile)
        except Exception as e:
            print(e)

    def save_model(self, outfile):
        """Save model to a file"""
        joblib.dump(self.clf, outfile)


    def train(self, X=None, Y=None, debug=20):
        """Train the model"""
        if self.scale:
            X = preprocessing.scale(X)
        self.clf.fit(X, Y)
        self.trained = True

    def feature_importance(self, outfile="importance.png", features_list=None):
        """Show each feature importance"""
        importances = self.clf.feature_importances_

        if features_list is not None and len(features_list) != len(importances):
            raise ValueError("Number of features does not fit!")

        indices = np.argsort(importances)[::-1]
        n_feats =  len(features_list)
        std = np.std([tree.feature_importances_ for tree in self.clf.estimators_], axis=0)

        plt.figure()
        plt.title("Feature importances")
        plt.bar(range(n_feats), importances[indices], width=0.5, color="b", yerr=std[indices], align="center")
        if features_list:
            features_list = np.asarray(features_list)[indices]
            plt.xticks(range(n_feats), features_list, rotation='vertical')
        plt.xlim([-1, n_feats])
        plt.margins(0.2)

        plt.subplots_adjust(bottom=0.15)
        plt.savefig(outfile, bbox_inches='tight')


    @classmethod
    def get_test_from_dataset(clc, L, X, Y, AA, tsize=0.2, rstate=0):
        test_position = []
        for i in xrange(len(Y)):
            if L[i][-1] == AA:
                test_position.append(i)
        
        if tsize:
            t_len = int(tsize*len(Y))
            zero_pos = np.where(Y==0)[0]

            test_len = len(test_position)
            random_zero_pos = np.random.choice(zero_pos, t_len, replace=False)
            text_position.extend(random_zero_pos)
 
        test_position = np.random.permutation(test_position)
        mask = np.ones(Y.shape, dtype=bool)
        mask[test_position] = False
        
        return shuffle(X[mask], Y[mask], random_state=rstate), shuffle(X[~mask], Y[~mask], np.array(range(len(mask)))[~mask], random_state=rstate)

    
    def cross_validation(self, X, Y, X_test=None, Y_test=None, tsize=0.3):
        """Cross validation on X and Y, using a sub sample"""
        if X_test is None or Y_test is None:
            X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=tsize)
        else :
            X_train = X
            Y_train = Y
        self.train(X_train, Y_train)
        Y_predicted = self.predict(X_test)
        return Y_predicted, self.get_score(X_test, Y_test)


    def get_sensibility_and_precision(self, pred_y, true_y, X_labels=None, X=None, log=True):
        nel = len(true_y)
        assert nel==len(pred_y), 'Vector should be the same size\n'
        
        true_pos, true_neg, false_pos, false_neg = 0.0, 0.0, 0.0, 0.0
        false_neg_list, false_pos_list = [], []
        for i in xrange(len(pred_y)):
            if pred_y[i] == 0 and true_y[i]==1:
                false_neg += 1
                false_neg_list.append(i)
            elif pred_y[i] == 1 and true_y[i]==0:
                false_pos += 1
                false_pos_list.append(i)
            elif pred_y[i] == 1 and true_y[i]==1:
                true_pos += 1
            elif pred_y[i] == 0 and true_y[i]==0:
                true_neg += 1

        if log:
            print "Test size is: %d\nTrue Positive is: %d\nTrue negative is: %d\nFalse positive is: %d\nFalse negative is:%d"%(nel, true_pos, true_neg, false_pos, false_neg)
            print '-------------------------------------------'
            print "Sensibility is %f"%(true_pos/(true_pos+false_neg))
            print "Specificity is %f"%(true_neg/(true_neg+false_pos))
            print "Accuracy is %f"%((true_neg+true_pos)/nel)
            print "Precision is %f\n\n"%(true_pos/(true_pos+false_pos))
            if X_labels is not None and X is not None:
                print "List of false negative"
                for i in false_neg_list:
                    print "\t".join(X_labels[i]+X[i])
                print "List of False positive"
                for i in false_pos_list:
                    print "\t".join(X_labels[i]+X[i])


    def get_score(self, X, Y):
        if self.scale:
            X = preprocessing.scale(X)
        return self.clf.score(X, Y)

    def predict(self, X):
        if not self.trained:
            raise ValueError("Classifier is not trained")
        if self.scale:
            X = preprocessing.scale(X)
        return self.clf.predict(X)

    def read_from_json(self, data, labels, use_global=True, use_pvalue=False):
        """Parse X array from data"""
        if isinstance(data, basestring):
            with open(data) as jfile:
                data =  json.load(jfile)

        if isinstance(labels, basestring):
            with open(labels) as jfile2:
                labels =  json.load(jfile2)
        # matrice format
        # global
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
        # [fitch, suspected, gene_frac, rea_frac, used_frac, subs_count]
        #[fitch, suspected, gene_frac, rea_frac, used_frac, subs_count, codon_lik_for_rea_aa]
        for aa2, val in data['aa'].items():
            for aa1, glist in val.items():
                for genome, gdata in glist.items():
                    type_check = gdata[dtype]
                    codon_total = gdata['codons'][dtype]
                    fitch =  gdata['fitch']
                    suspected =  gdata['suspected']
                    rea_codon =  type_check['rea_codon']
                    used_codon = type_check['used_codon']
                    gene_in_genome = data['genes'][genome]
                    was_lost =  gdata['lost'][fisher_type]
                    #mixte_codon = type_check['mixte_codon']
                    subs_count = type_check['count']
                    for codon in codon_total.keys():
                        gene_count = type_check['rea_distribution'].get(codon, 0)
                        codon_count = codon_total[codon]
                        codon_lik = gdata['score'][dtype].get(codon, 0)
                        # si il y a mutation, frequence d'utilisation du codon
                        rea_frac =  rea_codon.get(codon, 0)*1.0/codon_count if codon_count> 0 else 0
                        # frequence d'utilisation du codon dans les positions ou l'acide amine est predo
                        used_frac = used_codon.get(codon, 0)*1.0/codon_count if codon_count> 0 else 0
                        gene_frac =  gene_count*1.0 / gene_in_genome if gene_in_genome > 0 else 0

                        entry = [fitch, suspected, was_lost, gene_frac, rea_frac, used_frac, subs_count, codon_lik]
                        X.append(entry)
                        X_label.append([genome, codon, aa2, aa1])
                        codon_mapper = labels.get(genome, None)
                        if codon_mapper is None:
                            Y.append(0)
                        else:
                            codon_mapper = codon_mapper.get(codon, None)
                            if codon_mapper and codon_mapper == aa1:
                                Y.append(1)
                            else:
                                Y.append(0)

        assert (len(X) == len(Y)), "We should not have different length for data and for label"
        return np.array(X), X_label, np.array(Y)


def get_labels_from_csvfile(csvfile, genetic_code):
    
    def makehash():
        return defaultdict(makehash)
    
    codontable = CodonTable.unambiguous_dna_by_id[genetic_code]
    labels = makehash()
    with open(csvfile) as labfile:
        for line in labfile:
            line = line.strip()
            if line and not line.startswith('#'):
                genome, codon, reassignment = line.split()
                codon = codon.upper().replace('U','T')
                if len(codon) == 3 and codontable.forward_table[codon]!= reassignment.upper():
                    labels[genome][codon] = reassignment
    return labels


if __name__ == '__main__':
    import sys
    assert len(sys.argv)>=3
    jsonfile =  sys.argv[1]
    labels = None
    
    genetic_code = 1

    if len(sys.argv) == 4:
        genetic_code = int(sys.argv[3])
    labels = get_labels_from_csvfile(sys.argv[2], genetic_code)
    
    c = Classifier('rf', {}, False)
    X, Y, X_labels = c.read_from_json(jsonfile, labels, use_global=True):

    (X_train, Y_train), (X_test, Y_test, mask) = c.get_test_from_dataset(X_labels, X, Y, 'N', tsize=0.3)
    
    pred_y, score = c.cross_validation(X_train, Y_train, X_test, Y_test)
    assert  np.array_equal(X[mask], X_test), "Difference between mask and test dataset"
    Xlab = np.asarray(X_labels)
    c.get_sensibility_and_precision(pred_y, Y_test, Xlab[mask], X_test)

    flist = [ "fitch", "suspected", "was_lost", "gene_frac", "rea_frac", "used_frac", "subs_count", "codon_lik"]
    c.feature_importance(features_list=flist)
    print("\t".join(["genome", "codon", "ori_aa", "rea_aa"] + flist))
    for i in xrange(len(X_labels)):
        if Y[i] == 1:
            print("\t".join(X_labels[i]+[str(x) for x in X[i]]))
