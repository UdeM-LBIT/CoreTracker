
import os, sys
from classifier import *
from Bio.Data import CodonTable
from functools import partial

yeastfile = "/usr/local/www/CoreTracker/yeast.csv" 
yeastreafile = "/usr/local/www/CoreTracker/tmp/yeast/reassignment.json"

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


def test_on_unlabeled_data(classifier, X, X_label, unlabeled_data, etiquette=[], features=[]):
    """ Test classifier on data that was not labeled """
    if etiquette and features:
        X, etiquette = getDataFromFeatures(X_train, etiquette, feats=features)
    
    predicted_y = classifier.predict(X)
    for i in predicted_y:
        genome, codon, ori_aa, rea_aa = X_label[i]


def get_yeast_label(file, genetic_code):
    codontable = CodonTable.unambiguous_dna_by_id[genetic_code]
    train_list = makehash(depth=3, type=int)
    with open(file) as INPUT:
        for line in INPUT:
            line = line.strip()
            if line and not line.startswith('#'):
                genome, codon,  ori_aa, rea_aa, status = line.split("\t")
                status = int(status)
                if len(codon) == 3 and codontable.forward_table[codon] != rea_aa.upper():
                    train_list[genome][codon][rea_aa] = status

    return train_list


def get_train_data(X, X_labels, Y, testaa='M'):
    X_train = X[Y>=0]
    Y_train = Y[Y>=0]
    X_lab_train = np.asarray(X_labels)
    X_lab_train = X_lab_train[Y>=0]
    X_lab_train.shape
    test_ind = (X_lab_train[:, -1] == testaa)

    return X_train, X_lab_train, Y_train, test_ind


if __name__ == '__main__':

    assert len(sys.argv)>=3
    
    jsonfile =  sys.argv[1] 
    genetic_code = 1
    use_pca = False
    
    if len(sys.argv) == 4:
        genetic_code = int(sys.argv[3])

    if len(sys.argv) == 5:
        use_pca = bool(sys.argv[4])

    labels, not_labeled = parse_input_file(sys.argv[2], genetic_code)

    etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac", "N. rea", "N. used", "Cod. count", "Sub. count", "G. len", "codon_lik", "N. mixte" ,"id"]
    
    selected_feats = [2,3,4,5, 6, 7,8,9,11]
    X, X_labels, Y = read_from_json(jsonfile, labels, use_global=False)

    X_train, X_lab_train, Y_train, test_ind = get_train_data(X, X_labels, Y)

    X_train, etiquette = getDataFromFeatures(X_train, etiquette, feats=selected_feats)

    #print_data(X_train, X_lab_train, Y_train, etiquette)
    #sys.exit()
    
    # initialize classifier 
    parameters = {'random_state':12345, 'max_features':None, 'oob_score':False}
    c = Classifier('rf', parameters, False)
    
    # get train pos and test posi
    train_pos = np.where(test_ind==0)[0]
    test_pos = np.where(test_ind==1)[0]
    
    # get test and 
    x_test = X_train[test_pos]
    y_test = Y_train[test_pos]
    xlab_test = X_lab_train[test_pos]

    # train and predict    
    c.train(X_train[train_pos], Y_train[train_pos])
    y_pred = c.predict(x_test)
 
    # get result information
    c.get_stat(x_test, y_test)
    get_sensibility_and_precision(y_pred, y_test, xlab_test, x_test)

    c.feature_importance(outfile="importance.png", features_list=etiquette)
    get_2D_distinct(X_train, X_lab_train, Y_train, etiquette)
    X_pca_train = get_features(X_train,  Y_train, ncomp=6)
    draw_pca_data(X_pca_train, X_lab_train, Y_train, outfile="PCA.png")
    
    sys.exit()
    

    # get yeast data
    yeast_label = get_yeast_label(yeastfile, genetic_code)
    X_yeast, X_labels_yeast, Y_yeast = read_from_json(yeastreafile, yeast_label, use_global=False)
    testlab = json.load(open(yeastreafile))
    genome_list = set(testlab['genes'].keys())
    yeast_genome_list = set(yeast_label.keys())
    X_yeast_test, Xlab_yeast_test, Y_yeast_test, test_ind = get_train_data(X_yeast, X_labels_yeast, Y_yeast)


    etiquette = ["fitch", "suspected", "Fisher pval", "Gene frac", "N. rea", "N. used", "Cod. count", "Sub. count", "G. len", "codon_lik", "N. mixte" ,"id"] 

    X_yeast_test, etiquette = getDataFromFeatures(X_yeast_test, etiquette, feats=selected_feats)

    y_pred = c.predict(X_yeast_test)
 
    # get result information
    print("Result on yeast dataset prediction")
    get_sensibility_and_precision(y_pred, Y_yeast_test, Xlab_yeast_test, X_yeast_test)





    #print("\n"+"\t".join(["genome", "codon", "ori_aa", "rea_aa"]))
    #for i in xrange(len(X_labels)):
        #print("\t".join(list(X_labels[i])+[str(Y[i])]))

    # Xlab = np.asarray(X_labels)
    # X, etiquette = getDataFromFeatures(X, etiquette, feats=[1,2,3,4,6,9] )
    
    # (X_train, Y_train, train_mask), (X_test, Y_test, mask) = c.get_test_from_dataset(Xlab, X, Y, 'G', tsize=0.25)
    
    # y_pred, score = c.cross_validation(X_train, Y_train, X_test, Y_test)
    # c.get_sensibility_and_precision(y_pred, Y_test, Xlab[mask], X_test)
    # print_data(c, X, Xlab, Y,etiquette=etiquette)


    # zero_iter = c.split_zeros_pos(Xlab[train_mask], X_train, Y_train, split_size=100)
    # true_x = X_train[np.ma.make_mask(Y_train)]
    # true_x_len = true_x.shape[0]
    # true_x_lab = Xlab[train_mask]
    # true_x_lab = true_x_lab[np.ma.make_mask(Y_train)]
    # x_split, xlab_split, y_split = next(zero_iter)
    
    # x_split, xlab_split, y_split = shuffle(np.vstack((x_split,true_x)), np.vstack((xlab_split,true_x_lab)), np.hstack((y_split, np.ones(true_x_len))), random_state=0)

    # get_2D_distinct(x_split, xlab_split, y_split, etiquette, features=[])

    
    # warm_start = True
    # estimator_leap = 250 
    # estimator_start = 250

    # if warm_start:
    #     parameters = {'random_state':12345, 'max_features':None, 'oob_score':False, 'warm_start':warm_start}
    #     c_warm = Classifier('rf', parameters, False, n_estimators=estimator_start)
    #     zero_iter = c_warm.split_zeros_pos(Xlab[train_mask], X_train, Y_train, split_size=100)
    #     true_x = X_train[np.ma.make_mask(Y_train)]
    #     true_x_len = true_x.shape[0]
        
    #     x_split, xlab_split, y_split = next(zero_iter)
    #     x_split_train, y_split_train = shuffle(np.vstack((x_split,true_x)), np.hstack((y_split, np.ones(true_x_len))), random_state=0)
    #     c_warm.train(x_split_train, y_split_train)
    #     y_pred = c_warm.predict(X_test)
       
    #     for x_split, xlab_split, y_split in zero_iter:
    #         x_split_train, y_split_train = shuffle(np.vstack((x_split,true_x)), np.hstack((y_split, np.ones(true_x_len))), random_state=0)
    #         estimator_start += estimator_leap
    #         c_warm.clf.set_params(n_estimators=estimator_start)
    #         c_warm.train(x_split_train, y_split_train)        
    #         print "\n EST : %d -----------------------------------------------\n\n" % estimator_start

    #     y_pred = c_warm.predict(X_test)
    #     c_warm.get_sensibility_and_precision(y_pred, Y_test, Xlab[mask], X_test)