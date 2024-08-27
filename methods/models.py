import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification
from sklearn.base import BaseEstimator
from scipy import optimize
import pdb

class DataModel:
    def __init__(self, p, amplitude, random_state=None):
        if random_state is not None:
            np.random.seed(random_state)
        self.p = p
        self.a = amplitude

    def _sample_clean(self, n):
        return None

    def _sample_outlier(self, n):
        return None

    def sample(self, n, purity, offset=None, random_state=None):
        if random_state is not None:
            np.random.seed(random_state)
        p = self.p
        purity = np.clip(purity, 0, 1)
        n_clean = np.round(n * purity).astype(int)
        n_outlier = n - n_clean
        is_outlier = np.zeros((n,))
        if n_clean > 0:
            X_clean = self._sample_clean(n_clean)
            Y_clean = np.zeros((n_clean,))
        else:
            X_clean = np.zeros((0,p))
            Y_clean = np.zeros((0,))
        if n_outlier > 0:
            X_outlier = self._sample_outlier(n_outlier)
            Y_outlier = np.ones((n_outlier,))
        else:
            X_outlier = np.zeros((0,p))
            Y_outlier = np.ones((0,))
        X = np.concatenate([X_clean, X_outlier], axis=0)
        Y = np.concatenate([Y_clean, Y_outlier], axis=0)
        idx_shuffle = np.random.choice(n, size=n, replace=False)

        return X[idx_shuffle], Y[idx_shuffle].astype(int)



class GaussianMixture(DataModel):
    def _sample_clean(self, n):
        p = self.p
        X = np.random.randn(n, p)
        return X

    def _sample_outlier(self, n):
        p = self.p
        X = self.a + np.random.randn(n, p)
        return X


class ConcentricCircles(DataModel):
    def _sample_clean(self, n):
        p = self.p
        X = np.random.randn(n, p)
        return X

    def _sample_outlier(self, n):
        p = self.p
        X = np.sqrt(self.a)*np.random.randn(n, p)
        return X


class AdversarialModel(DataModel):
    def __init__(self, p, amplitude, occ_model: BaseEstimator, random_state=None):
        super().__init__(p, amplitude, random_state=random_state)
        self.occ_model = occ_model
        assert(self.a>=1)

    def _sample_clean(self, n):
        p = self.p
        X = np.random.randn(n, p)
        return X

    def _sample_outlier(self, n):
        p = self.p
        n_train = 1000
        n_candidates = int(np.ceil(self.a * n))

        # Generate clean training data
        X_train = self._sample_clean(n_train)

        # Train the OCC model on clean data
        self.occ_model.fit(X_train)

        # Generate candidate outliers by adding small perturbations
        candidate_outliers = self._sample_clean(n_candidates)

        # Compute decision function scores for candidate outliers
        candidate_scores = self.occ_model.decision_function(candidate_outliers)

        # Select the adversarial outliers closest to the decision boundary
        adversarial_indices = np.argsort(np.abs(candidate_scores))[:n]
        X = candidate_outliers[adversarial_indices]

        # Shuffle the order of the outliers
        np.random.shuffle(X)

        return X


class ConcentricCircles2(DataModel):
    def _sample_clean(self, n):
        p = self.p
        X = np.random.randn(n, p)
        return X

    def _sample_outlier(self, n):
        p = self.p
        rescale = np.ones((1,p))
        rescale[0,0:int(p/2)] = np.sqrt(self.a)
        X = rescale*np.random.randn(n, p)
        return X


class ConcentricCirclesMixture(DataModel):
    def __init__(self, p, amplitude, random_state=None):
        super().__init__(p, amplitude, random_state=random_state)
        self.Z = np.random.uniform(low=-3, high=3, size=(p,p))

    def _sample_clean(self, n):
        p = self.p
        X = np.random.randn(n, p)
        cluster_idx = np.random.choice(self.Z.shape[0], n, replace=True)
        X = X + self.Z[cluster_idx,]
        return X

    def _sample_outlier(self, n):
        p = self.p
        X = np.sqrt(self.a) * np.random.randn(n, p)
        cluster_idx = np.random.choice(self.Z.shape[0], n, replace=True)
        X = X + self.Z[cluster_idx,]
        return X


class BinomialModel(DataModel):
    def __init__(self, p, amplitude, random_state=None):
        super().__init__(p, amplitude, random_state=random_state)
        self.beta_Z = np.sqrt(amplitude)*np.random.normal(size=(p,2))

    def calculate_offset(self, purity):
        X = self.sample_X(1000)
        def foo(offset):
            Y = self.sample_Y(X, offset)
            return np.mean(Y) - (1.0-purity)
        offset = optimize.bisect(foo, -1000, 1000)
        return offset

    def sample_X(self, n):
        X = np.random.normal(0, 1, (n,self.p))
        X[:,0] = np.random.uniform(low=0, high=1, size=(n,))
        return X

    def compute_prob(self, X, offset):
        f = np.matmul(X,self.beta_Z)
        f[:,0] = f[:,0] - offset/2
        f[:,1] = f[:,1] + offset/2
        prob = np.exp(f)
        prob_y = prob / np.expand_dims(np.sum(prob,1),1)
        return prob_y

    def sample_Y(self, X, offset):
        prob_y = self.compute_prob(X, offset)
        g = np.array([np.random.multinomial(1,prob_y[i]) for i in range(X.shape[0])], dtype = float)
        classes_id = np.arange(2)
        y = np.array([np.dot(g[i],classes_id) for i in range(X.shape[0])], dtype = int)
        return y

    def sample(self, n, purity=1, offset=None, random_state=None):
        n_clean = np.round(n * purity).astype(int)
        n_out = n - n_clean

        if offset is None:
            offset = self.calculate_offset(purity)
            print("Purity: {:.3f}. Offet: {:.3f}.".format(purity, offset))
        X = self.sample_X(2*n)
        Y = self.sample_Y(X, offset)

        if n_clean>0:
            idx_in = np.random.choice(np.where(Y==0)[0], n_clean, replace=False)
            X_in = X[idx_in]
            Y_in = np.zeros((n_clean, ))
        else:
            X_in = np.zeros((0, X.shape[1]))
            Y_in = Y[idx_in]

        if n_out>0:
            idx_out = np.random.choice(np.where(Y==1)[0], n_out, replace=False)
            X_out = X[idx_out]
            Y_out = Y[idx_out]
        else:
            X_out = np.zeros((0, X.shape[1]))
            Y_out = np.zeros((0, ))

        X = np.concatenate([X_in, X_out], axis=0)
        Y = np.concatenate([Y_in, Y_out], axis=0)
        idx_shuffle = np.random.choice(n, size=n, replace=False)

        return X[idx_shuffle], Y[idx_shuffle].astype(int)
