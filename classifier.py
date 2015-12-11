from sklearn.ensemble import RandomForestClassifier
from sklearn import feature_selection
from sklearn import svm
from sklearn.externals import joblib
from Bio.Data import CodonTable
from sklearn import preprocessing
from sklearn import cross_validation
from collections import defaultdict
import json
import numpy as np


from sklearn.cross_validation import train_test_split

class Classifier:
    """This is a classifier for codon rassignment"""
    def __init__(self, method, classifier_spec, scale=False):
        if method == 'rf':
            self.clf = RandomForestClassifier(n_estimators=100, n_jobs=-1, **classifier_spec)
        elif method =='svc':
            self.clf = svm.SVC(**classifier_spec)
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
        if self.scale:
            X = preprocessing.scale(X)
        self.clf.fit(X, Y)
        self.trained = True

    def cross_validation(self, X, Y, tsize=0.3):
        """Cross validation on X and Y, using a sub sample"""
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=tsize)
        self.train(X_train, y_train)
        return self.get_score(X_test, y_test)

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

    def read_from_json(self, data, labels, use_global=True):
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
                    codon_total = gdata['codons']
                    fitch =  gdata['fitch']
                    suspected =  gdata['suspected']
                    rea_codon =  type_check['rea_codon']
                    used_codon = type_check['used_codon']
                    gene_in_genome = data['genes'][genome]
                    was_lost =  gdata['lost']
                    #mixte_codon = type_check['mixte_codon']
                    subs_count = type_check['count']
                    for codon in codon_total.keys():
                        gene_count = type_check['rea_distribution'].get(codon, 0)
                        codon_count = codon_total[codon]
                        codon_lik = gdata['score'].get(codon, 0)
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
    print "Score on test set, using a Random Forest without tuning : "
    c = Classifier('rf', {}, True)
    X, X_labels, Y = c.read_from_json(sys.argv[1], labels)
    p = c.cross_validation(X, Y)
    print(p)

    print "Score on test set, using a SVC classifer with RBF without tuning : "
    c2 = Classifier('svc', {}, True)
    p = c2.cross_validation(X, Y)
    print(p)
    #model = svm.SVC()
    #print cross_validation.cross_val_score(model, X, Y, scoring='precision')
    
    #print("\t".join(["genome", "codon", "ori_aa", "rea_aa", "fitch", "suspected", "was_lost", "gene_frac", "rea_frac", "used_frac", "subs_count", "codon_lik"]))
    #for line in X:
    #    print("\t".join([str(x) for x in line]))


