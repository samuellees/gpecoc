import os.path
import os
import time
import numpy as np



feature_method = ["fclassify", "svmrfe", "forest", "bsswss"]

"""
# interface function of feature selection
"""
def select_features(filename, X, y, feature_num):
    fsel = FeatureSel(feature_num)
    #timer = time.strftime('%Y-%m-%d.%H:%M:%S',time.localtime(time.time()))
    #cache = filename + "_features_" + timer + ".txt"
    cache = filename + "_features_" + ".txt"
    sel_features = fsel.select(X, y, cache)
    
    feature_F1 = fsel.process_select_F1(X, y)
    feature_F2 = fsel.process_select_F2(X, y)
    feature_F3 = fsel.process_select_F3(X, y)
    feature_F4 = fsel.process_select_F4(X, y)
    
    return sel_features, feature_F1, feature_F2, feature_F3, feature_F4


class FeatureSel():
    """Fuse the genes selected by various feature selection technique
    """

    def __init__(self, f_num):
        self.sel_features = {}
        self.cache_changed = False
        self.f_num = f_num #

    def select(self, data, labels, cache_path=None):
        """

        Arguments:
        - `data`: expression of the data
        - `labels`: labels (class) of the data
        - `cache`: the path for result cache.
        """
        if cache_path and os.access(cache_path, os.F_OK):
            self._load_cache(cache_path)
            # If some new feature selection methods is involved.
            # we have to use them  and save the new selected features.
            # methods with cache loaded will skip the process
            self.process_select(data, labels)
            if self.cache_changed:
                #print "save1"
                self._save_cache(cache_path)
        else:
            self.process_select(data, labels)
            if cache_path:
                #print "save2"
                self._save_cache(cache_path)
        return self.fused

    def _save_cache(self, cache_path):
        self.makedir_for_cache(cache_path)
        f = file(cache_path, 'w')
        for (k, v) in self.sel_features.iteritems():
            f.write(str(k))
            # seperator between methods and selected features
            f.write("::")
            f.write(str(v.tolist()))
            f.write("\n")

    def _load_cache(self, cache_path):
        f = file(cache_path, 'r')
        for line in f:
            (k, v) = line.split("::")
            self.sel_features[k] = np.asarray(eval(v))

    def makedir_for_cache(self, file_path):
        dirname = os.path.dirname(file_path)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)

    def process_select(self, data, labels):
        """Do the real work of feature selection, with various techniques.

        Arguments:
        - `data`: expression of the data
        - `labels`: labels (class) of the data
        """
        self.process_select_fclassify(data, labels)
        #print "1"
        #self.process_select_irelief(data, labels)
        self.process_select_svmrfe(data, labels)
        #print "2"
        self.process_select_forest(data, labels)
        #print "3"
        self.process_select_bsswss(data, labels)
        return self._fuse_features()

    def _fuse_features(self):
        """Select a set of features from candidates.

        """
        self.fused = []
        for v in self.sel_features.itervalues():
            self.fused += v.tolist()
        self.fused = np.unique(self.fused)
        return self.fused

    def process_select_svmrfe(self, data, labels):
        if 'svm' in self.sel_features:
            return

        self.cache_changed = True

        from sklearn.svm import SVC
        from sklearn.feature_selection import RFE

        svc = SVC(kernel="linear", C=1)
        rfe = RFE(estimator=svc, n_features_to_select=100, step=2)
        rfe.fit(data, labels)

        #        selected = np.arange(features.shape[1])[rfe.support_]
        selected = np.argsort(rfe.ranking_)
        selected = selected[:self.f_num]
        self.sel_features['svm'] = selected
        return selected

    def process_select_fclassify(self, data, labels):
        if 'ftest' in self.sel_features:
            return

        self.cache_changed = True

        from sklearn.feature_selection import SelectKBest, f_classif
        selector = SelectKBest(f_classif, self.f_num)
        selector.fit(data, labels)

        selected = np.argsort(selector.scores_)[::-1]
        selected = selected[:self.f_num]
        self.sel_features["ftest"] = selected
        return selected

    def process_select_forest(self, data, labels):
        if 'forest' in self.sel_features:
            return

        self.cache_changed = True

        from sklearn.ensemble import RandomForestClassifier
        forest = RandomForestClassifier(random_state=0)
        forest.fit(data, labels)

        selected = np.argsort(forest.feature_importances_)[::-1]
        selected = selected[:self.f_num]
        self.sel_features["forest"] = selected
        return selected

    def process_select_irelief(self, data, labels):
        if 'irelief' in self.sel_features:
            return

        self.cache_changed = True

        from mlpy.irelief import IRelief
        selector = IRelief(sigma=2)
        selector.learn(data, labels)

        selected = np.argsort(selector.weights())[::-1]
        selected = selected[:self.f_num]
        self.sel_features["irelief"] = selected
        return selected
        
    def process_select_bsswss(self, data, labels):
        if 'bsswss' in self.sel_features:
            return

        self.cache_changed = True

        def bss_wss_value(f, labels):
            names = sorted(set(labels))
            wss, bss = np.array([]), np.array([])
            for name in names:
                f_k = f[labels == name]
                f_m = f_k.mean()
                d_m = (f_m - f.mean()) ** 2
                d_z = (f_k - f_m) ** 2
                bss = np.append(bss, d_m.sum())
                wss = np.append(wss, d_z.sum())
            z, m = bss.sum(), wss.sum()
            bsswss = z / m if m > 0 else 0
            return bsswss

        i = 0
        x, y = [], []
        for f in data.transpose():
            x.append(i)
            y.append(bss_wss_value(f, labels))
        selected = np.argsort(y)[::-1]
        selected = selected[:self.f_num]
        
        self.sel_features["bsswss"] = selected
        return selected   
        
        
        
    def process_select_F1(self, data, labels):
        if 'ftest' in self.sel_features:
            F1 = (self.sel_features['ftest']).tolist()
        #print F1
        return F1  
    def process_select_F2(self, data, labels):
        if 'svm' in self.sel_features:
            F2 = (self.sel_features['svm']).tolist()
        #print F2
        return F2  
    def process_select_F3(self, data, labels):
        if 'forest' in self.sel_features:
            F3 = (self.sel_features['forest']).tolist()
        #print F3
        return F3  
    def process_select_F4(self, data, labels):
        if 'bsswss' in self.sel_features:
            F4 = (self.sel_features['bsswss']).tolist()
        #print F4
        return F4  