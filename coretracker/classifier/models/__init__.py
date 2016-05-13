import numpy as np

class ModelType(object):
    """Representation of each classification model"""
    mod1 = [2,3,4,5,6,7,8,9,11]
    mod2 = [0,2,3,4,5,6,7,8,9,11]
    mod3 = [0,2,3,4,5,6,7,8,9,11]
    oneencoded_models  = ['3']
    default_sfeat = {'1':mod1, '2':mod2, '3':mod3}
    def __init__(self, m, etiquette, axisselector=None, sfeat=[], encode=False):
        if m not in  self.default_sfeat.keys():
            raise ValueError('Selected model do not exist')
        self.model = str(m)
        self.selector = axisselector
        if not self.selector:
            self.selector = self.get_data_from_feature
        self.etiquette = etiquette
        if not sfeat:
            try:
                self.sfeat = self.default_sfeat[self.model]
            except:
                self.sfeat = range(len(etiquette))
        else:
            self.sfeat = sfeat
        self.encode = encode
        if str(self.model) in self.oneencoded_models:
            self.encode = True

    def format_data(self, data):
        """Format data into training and printing data"""
        training_data, selected_et = self.selector(data, self.etiquette, self.sfeat)
        if self.encode and 'id' in self.printing_data:
            training_data =  self._idonehotencode(training_data)
        printing_data = data[:, self.sfeat]
        return training_data, printing_data, selected_et

    def _idonehotencode(self, data, totsize=64):
        """On hot encoding of the last column"""
        last_row = (data[:, -1]).astype(int)
        rlen = len(last_row)
        enc = np.zeros((rlen, totsize))
        ind = np.ravel_multi_index([np.arange(rlen), last_row], (rlen, totsize))
        enc.flat[ind] = 1
        return np.column_stack((data[:,:-1], enc))
        return data

    @staticmethod
    def get_data_from_feature(data, etiquette, feats=[]):
        """Extract Data based on list of features"""
        if len(feats) == 0:
            return data, etiquette
        else:
            return data[:, feats], np.asarray(etiquette)[feats]
