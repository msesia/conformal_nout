import numpy as np
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import train_test_split

class DataSet:
    def __init__(self, X, Y):
        idx_in = np.where(Y==0)[0]        
        idx_out = np.where(Y==1)[0]        
        self.X_in  = X[idx_in] 
        self.Y_in  = Y[idx_in] 
        self.X_out = X[idx_out]
        self.Y_out = Y[idx_out]
        

    def sample(self, n_train, n_cal, n_test, prop_out, random_state=2024):
        # Sample the reference data
        n_ref = n_train+n_cal
        X_ref, X_in_hout = train_test_split(self.X_in, train_size=n_ref, random_state=random_state)
        # Split the reference data into training and calibration sets
        X_train, X_cal = train_test_split(X_ref, train_size=n_train, random_state=101+random_state)
        # Sample the test data
        n_test_out = np.round(n_test * prop_out).astype(int)
        n_test_in = n_test - n_test_out

        if n_test_in > 0:            
            X_test_in, _ = train_test_split(X_in_hout, train_size=n_test_in, random_state=102+random_state)
            Y_test_in = np.zeros((len(X_test_in),))
        else:
            X_test_in = np.ones((0, X_ref.shape[1]))
            Y_test_in = np.zeros((0, ))

        if n_test_out > 0:
            X_test_out, _ = train_test_split(self.X_out, train_size=n_test_out, random_state=103+random_state)
            Y_test_out = np.ones((len(X_test_out),))
        else:
            X_test_out = np.ones((0, X_ref.shape[1]))
            Y_test_out = np.ones((0, ))

        X_test = np.concatenate([X_test_in, X_test_out], axis=0)
        Y_test = np.concatenate([Y_test_in, Y_test_out], axis=0)
        idx_shuffle = np.random.choice(n_test, size=n_test, replace=False) 
        X_test = X_test[idx_shuffle]
        Y_test = Y_test[idx_shuffle]
        
        return X_train, X_cal, X_test, Y_test.astype(int)


def fdr_filter_bh(pvals, alpha):
    reject, _, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    return reject

def fdr_filter_storey_bh(pvals, alpha, lambda_par=0.5):
    pi_hat = (1.0 + np.sum(pvals>lambda_par)) / (len(pvals)*(1.0 - lambda_par))
    alpha_eff = alpha/pi_hat
    reject, _, _, _ = multipletests(pvals, alpha=alpha_eff, method='fdr_bh')
    return reject    

def eval_discoveries(discovered, is_outlier):
    idx_in = np.where(is_outlier==0)[0]
    idx_out = np.where(is_outlier==1)[0]
    fp = np.sum(discovered[idx_in])
    tp = np.sum(discovered[idx_out])
    fdp = fp/np.maximum(1,fp+tp)
    power = tp/np.maximum(1, np.sum(is_outlier))
    return fdp, power
