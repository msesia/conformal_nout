import numpy as np
import pandas as pd
from sklearn.datasets import load_digits, fetch_covtype, fetch_openml
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from scipy.io import arff
import mat73
import re
import pickle
from scipy.io import arff, loadmat
import pdb
from joblib import Memory
memory = Memory('./tmp')
fetch_openml_cached = memory.cache(fetch_openml)

def unpickle(file):
    """load the cifar-10 data"""

    with open(file, 'rb') as fo:
        data = pickle.load(fo, encoding='bytes')
    return data


def load_cifar_10_data(data_dir, negatives=False):
    """
    Return train_data, train_filenames, train_labels, test_data, test_filenames, test_labels
    """

    meta_data_dict = unpickle(data_dir + "/batches.meta")
    cifar_label_names = meta_data_dict[b'label_names']
    cifar_label_names = np.array(cifar_label_names)

    # training data
    cifar_train_data = None
    cifar_train_filenames = []
    cifar_train_labels = []

    # cifar_train_data_dict
    # 'batch_label': 'training batch 5 of 5'
    # 'data': ndarray
    # 'filenames': list
    # 'labels': list

    for i in range(1, 6):
        cifar_train_data_dict = unpickle(data_dir + "/data_batch_{}".format(i))
        if i == 1:
            cifar_train_data = cifar_train_data_dict[b'data']
        else:
            cifar_train_data = np.vstack((cifar_train_data, cifar_train_data_dict[b'data']))
        cifar_train_filenames += cifar_train_data_dict[b'filenames']
        cifar_train_labels += cifar_train_data_dict[b'labels']

    #cifar_train_data = cifar_train_data.reshape((len(cifar_train_data), 3, 32, 32))
    #if negatives:
    #    cifar_train_data = cifar_train_data.transpose(0, 2, 3, 1).astype(np.float32)
    #else:
    #    cifar_train_data = np.rollaxis(cifar_train_data, 1, 4)
    cifar_train_filenames = np.array(cifar_train_filenames)
    cifar_train_labels = np.array(cifar_train_labels)

    # test data
    # cifar_test_data_dict
    # 'batch_label': 'testing batch 1 of 1'
    # 'data': ndarray
    # 'filenames': list
    # 'labels': list

    cifar_test_data_dict = unpickle(data_dir + "/test_batch")
    cifar_test_data = cifar_test_data_dict[b'data']
    cifar_test_filenames = cifar_test_data_dict[b'filenames']
    cifar_test_labels = cifar_test_data_dict[b'labels']

    cifar_test_data = cifar_test_data.reshape((len(cifar_test_data), 3, 32, 32))
    if negatives:
        cifar_test_data = cifar_test_data.transpose(0, 2, 3, 1).astype(np.float32)
    else:
        cifar_test_data = np.rollaxis(cifar_test_data, 1, 4)
    cifar_test_filenames = np.array(cifar_test_filenames)
    cifar_test_labels = np.array(cifar_test_labels)

    return cifar_train_data, cifar_train_filenames, cifar_train_labels, \
        cifar_test_data, cifar_test_filenames, cifar_test_labels

def load_real_data(base_path, data_name):
        # Load the data
        if data_name=="images_flowers":
            data_raw = pd.pandas.read_csv(base_path + data_name + ".csv", sep=",", header=None)
            Y = np.array(data_raw.iloc[:,0])
            X = np.array(data_raw.iloc[:,1:])
            labels_inlier = ["roses"]
            labels_outlier_train = []
            labels_outlier_test = ["sunflowers","dandelion","daisy","tulips"]

        elif data_name=="images_animals":
            data_raw = pd.pandas.read_csv(base_path + data_name + ".csv", sep=",", header=None)
            Y = np.array(data_raw.iloc[:,0])
            X = np.array(data_raw.iloc[:,1:])
            labels_inlier = ["hamster", "guinea pig"]
            labels_outlier_train = []
            labels_outlier_test = ["lynx","wolf","coyote","cheetah","jaguer","chimpanzee","orangutan","cat"]

        elif data_name=="images_cars":
            data_raw = pd.pandas.read_csv(base_path + data_name + ".csv", sep=",", header=None)
            Y = np.array(data_raw.iloc[:,0])
            X = np.array(data_raw.iloc[:,1:])
            labels_inlier = ["car"]
            labels_outlier_train = []
            labels_outlier_test = ["fruit", "dog", "motorbike", "person", "cat", "flower", "airplane"]

        elif data_name in ["mammography", "shuttle", "cover", "pendigits"]:
            mat = loadmat(base_path + data_name + ".mat")
            X = mat['X']
            Y = mat['y']
            labels_inlier = [0]
            labels_outlier_train = []
            labels_outlier_test = [1]

        elif data_name == "aloi":
            data_raw = arff.loadarff(base_path + data_name + ".arff")
            data_raw = pd.DataFrame(data_raw[0])
            Y = np.array(data_raw.iloc[:,-2])
            X = np.array(data_raw.iloc[:,0:-2])
            labels_inlier = [b'no']
            labels_outlier_train = []
            labels_outlier_test = [b'yes']

        elif data_name=="creditcard":
            data_raw = pd.pandas.read_csv(base_path + data_name + ".csv", sep=",", header=None)
            Y = np.array(data_raw.iloc[:,-1])
            X = np.array(data_raw.iloc[:,0:-1])
            labels_inlier = ['normal']
            labels_outlier_train = []
            labels_outlier_test = ['fraud-1','fraud-0']

        elif data_name=="annthyroid":
            mat = loadmat(base_path + "annthyroid.mat")
            X = mat['X']
            Y = mat['y'].flatten()
            labels_inlier = [0]
            labels_outlier_train = []
            labels_outlier_test = [1]

        is_inlier = np.array([y in labels_inlier for y in Y]).astype(int)
        is_outlier = 1-is_inlier
        is_outlier_train = np.array([y in labels_outlier_train for y in Y]).astype(int)
        is_outlier_test = np.array([y in labels_outlier_test for y in Y]).astype(int)

        print("Loaded data set with {:d} samples: {:d} inliers and {:d} outliers, of which {:d} are available for training."\
              .format(len(Y), np.sum(is_inlier), np.sum(is_outlier), np.sum(is_outlier_train)))

        return X, is_outlier
