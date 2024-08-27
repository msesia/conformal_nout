import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import pdb

# List of classifiers
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC


import matplotlib.pyplot as plt
import seaborn as sns

import os, sys
sys.path.append("../methods")

from models import GaussianMixture, ConcentricCircles, ConcentricCirclesMixture, BinomialModel, AdversarialModel
import conformal
from my_utils import DataSet, fdr_filter_bh, fdr_filter_storey_bh, eval_discoveries
from utils_data import load_real_data

#########################
# Experiment parameters #
#########################

if True: # Input parameters
    # Parse input arguments
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    model_num = 1
    if len(sys.argv) != 14:
        print("Error: incorrect number of parameters.")
        quit()

    setup = int(sys.argv[1])
    data_name = sys.argv[2]
    n_train = int(sys.argv[3])
    n_cal = int(sys.argv[4])
    n_test = int(sys.argv[5])
    p = int(sys.argv[6])
    a = float(sys.argv[7])
    prop_out = float(sys.argv[8])
    classifier = sys.argv[9]
    tune_size = float(sys.argv[10])
    alpha = float(sys.argv[11])
    selection_method = sys.argv[12]
    random_state = int(sys.argv[13])

else: # Default parameters
    setup = 1
    data_name = "circles"
    n_train = 1000
    n_cal = 1000
    n_test = 1000
    p = 100
    a = 0.9
    prop_out = 0.5
    classifier = "auto"
    tune_size = 0.25
    alpha = 0.1
    selection_method = "none"
    random_state = 2022


# Fixed experiment parameters
allow_sign_flip = True
if setup in [5]:
    num_repetitions = 100
else:
    num_repetitions = 10
#n_perm = 200
n_perm = -1 # 200
B_perm = 1000

tables_t2_path = "../tables/table_t2.csv"
tables_t3_path = "../tables/table_t3.csv"
tables_t4_path = "../tables/table_t4.csv"
tables_fisher_path = "../tables/table_fisher.csv"

def load_table(path, alpha):
    if os.path.exists(path):
        table_raw = pd.read_csv(path)
        column_name = "a"+str(alpha)
        if column_name in table_raw:
            table = pd.DataFrame({'m':table_raw["m"], 'n':table_raw["n"], 'critical':table_raw[column_name].astype(float)})
            print("Loaded table of critical values for alpha={:.2f} from {:s}.".format(alpha, path))
        else:
            table = None
            print("Could not find critical values for alpha={:.2f} in {:s}.".format(alpha, path))
    else:
        table = None
        print("Could not find table file {:s}.".format(path))
    return table

table_t2 = load_table(tables_t2_path, alpha)
table_t3 = load_table(tables_t3_path, alpha)
table_t4 = load_table(tables_t4_path, alpha)
table_fisher = load_table(tables_fisher_path, alpha)

# Candidate values of K for WMW local test
#method_candidates = ["wmw-K2", "wmw-K3", "wmw-K4", "fisher", "simes", "storey_simes"]
method_candidates = ["wmw-K2", "lmp", "fisher", "simes", "storey_simes"]

#########################
# Data-generating model #
#########################

if data_name=="lhco":
    data_path_local = "/media/msesia/Samsung1/data/physics/"
    data_path_cluster = "/project/sesia_1123/lhco/"
    if os.path.isdir(data_path_local):
        data_path = data_path_local
    elif os.path.isdir(data_path_cluster):
        data_path = data_path_cluster

    # Load everything into memory
    df_raw = pd.read_hdf(data_path + "events_anomalydetection_v2.features.h5")
    X_raw = np.array(df_raw.iloc[:,0:14])
    Y_raw = np.array(df_raw.iloc[:,14]).astype(int)

else:

    data_path_local = "/media/msesia/Samsung1/data/"
    data_path_local_2 = "../../data/"
    data_path_cluster = "/project/sesia_1123/data/outliers/"
    if os.path.isdir(data_path_local):
        data_path = data_path_local
    elif os.path.isdir(data_path_local_2):
        data_path = data_path_local_2
    elif os.path.isdir(data_path_cluster):
        data_path = data_path_cluster
    else:
        print("Error! Data directory not found.")
        exit(-1)

    if data_name in ["mammography", "creditcard", "shuttle", "cover", "pendigits", "aloi"]:
        X_raw, Y_raw = load_real_data(data_path, data_name)

    elif data_name=="moons":
        n_tot = 10 * (n_train + n_cal + n_test)
        from sklearn.datasets import make_moons
        X_raw, Y_raw = make_moons(n_samples=n_tot, shuffle=True, noise=0.2, random_state=random_state)

    elif data_name=="circles-mixed":
        model = ConcentricCirclesMixture(p, a, random_state=random_state)
        n_tot = 10 * (n_train + n_cal + n_test)
        X_raw, Y_raw = model.sample(n_tot, 0.5, offset=0)

    elif data_name=="adversarial":
        occ_model_adversary = OneClassSVM(kernel='rbf', degree=3)
        model = AdversarialModel(p, a, occ_model=occ_model_adversary, random_state=random_state)
        n_tot = 10 * (n_train + n_cal + n_test)
        X_raw, Y_raw = model.sample(n_tot, 0.5, offset=0)

    elif data_name=="binomial":
        model = BinomialModel(p, a, random_state=random_state)
        n_tot = 10 * (n_train + n_cal + n_test)
        X_raw, Y_raw = model.sample(n_tot, 0.5, offset=0)

    elif data_name.startswith("mixture-"):
        data_name_str, prop_binomial = data_name.split("-")
        prop_binomial = float(prop_binomial)
        model_1 = BinomialModel(p, a, random_state=random_state)
        model_2 = ConcentricCirclesMixture(p, 0.7, random_state=random_state)
        n_tot = 10 * (n_train + n_cal + n_test)
        n_tot_1 = int(prop_binomial * n_tot)
        n_tot_2 = int((1-prop_binomial) * n_tot)
        if n_tot_1 > 0:
            X_raw_1, Y_raw_1 = model_1.sample(n_tot_1, 0.5, offset=0)
        if n_tot_2 > 0:
            X_raw_2, Y_raw_2 = model_2.sample(n_tot_2, 0.5, offset=0)

        if n_tot_1==0:
            X_raw_1 = np.zeros((0,X_raw_2.shape[1]))
            Y_raw_1 = np.zeros((0,))
        if n_tot_2==0:
            X_raw_2 = np.zeros((0,X_raw_1.shape[1]))
            Y_raw_2 = np.zeros((0,))

        X_raw = np.concatenate([X_raw_1, X_raw_2], axis=0)
        Y_raw = np.concatenate([Y_raw_1, Y_raw_2], axis=0)

    # elif data_name=="lhco":
    #     data_path_local = "/media/msesia/Samsung/data/physics/events_anomalydetection_v2.features.h5"
    #     data_path_cluster = "/project/sesia_1123/lhco/events_anomalydetection_v2.features.h5"
    #     if os.path.isdir(data_path_local):
    #         data_path = data_path_local
    #     elif os.path.isdir(data_path_cluster):
    #         data_path = data_path_cluster

    #     n_tot = 10 * (n_train + n_cal + n_test)
    #     from sklearn.datasets import make_moons
    #     X_raw, Y_raw = make_moons(n_samples=n_tot, shuffle=True, noise=0.2, random_state=2024)

    else:
        print("Error! Unknown data")
        exit(-1)

dataset = DataSet(X_raw, Y_raw)

###############
# Classifiers #
###############

if classifier in ["auto", "occ-auto"]:
    clf_occ_list = {'IF':IsolationForest(random_state=random_state),
                    'SVM':OneClassSVM(kernel='rbf', degree=3),
                    'LOF':LocalOutlierFactor(novelty=True)}
elif classifier == "occ-if":
    clf_occ_list = {'IF':IsolationForest(random_state=random_state)}
elif classifier == "occ-svm":
    clf_occ_list = {'SVM':OneClassSVM(kernel='rbf', degree=3)}
elif classifier == "occ-lof":
    clf_occ_list = {'LOF':LocalOutlierFactor(novelty=True)}
else:
    clf_occ_list = {}

if classifier in ["auto", "bc-auto"]:
    clf_bc_list = {'MLP':MLPClassifier(max_iter=1000, random_state=random_state, verbose=False),
                   'RF':RandomForestClassifier(random_state=random_state),
                   'ABC':AdaBoostClassifier(random_state=random_state)}
elif classifier == "bc-mlp":
    clf_bc_list = {'MLP':MLPClassifier(max_iter=1000, random_state=random_state, verbose=False)}
elif classifier == "bc-rf":
    clf_bc_list = {'RF':RandomForestClassifier(random_state=random_state)}
elif classifier == "bc-abc":
    clf_bc_list = {'ABC':AdaBoostClassifier(random_state=random_state)}
else:
    clf_bc_list = {}


###############
# Output file #
###############
outfile_prefix = "results/setup" + str(setup) + "/" +str(data_name) + "_n"+str(n_train) + "_"+str(n_cal) + "_"+str(n_test)
outfile_prefix += "_p" + str(p) + "_a" + str(a) + "_pt" + str(prop_out)
outfile_prefix += "_" + classifier + "_ts"+str(tune_size) + "_alpha" + str(alpha) + "_" + str(selection_method) + "_s" + str(random_state)
outfile = outfile_prefix + ".txt"
print("Output file: {:s}".format(outfile), end="\n")

# Header for results file
def add_header(df):
    df["Setup"] = setup
    df["Data"] = data_name
    df["n_train"] = n_train
    df["n_cal"] = n_cal
    df["n_test"] = n_test
    df["p"] = p
    df["Signal"] = a
    df["prop_out"] = prop_out
    df["Classifier"] = classifier
    df["Alpha"] = alpha
    df["tune_size"] = tune_size
    df["selection"] = selection_method
    df["Seed"] = random_state
    return df



###################
# Run experiments #
###################


def run_experiment(dataset, random_state):
    # Sample data
    X_in_train, X_in_cal, X_test, Y_test = dataset.sample(n_train, n_cal, n_test, prop_out, random_state=random_state)

    # Initialize result data frame
    results = pd.DataFrame({})

    scores_cal, scores_test, sel_classifier, sel_method = conformal.calculate_scores_auto(X_in_train, X_in_cal, X_test,
                                                                                          clf_occ_list=clf_occ_list, clf_bc_list=clf_bc_list,
                                                                                          method_candidates=method_candidates, tune_size=tune_size,
                                                                                          random_state=random_state, allow_sign_flip=allow_sign_flip,
                                                                                          alpha=alpha, selection_method=selection_method,
                                                                                          n_perm=n_perm, B=B_perm,
                                                                                          table_t2=table_t2, table_t3=table_t3, table_t4=table_t4, table_fisher=table_fisher)
    print("Selected classifier: {:s}".format(sel_classifier))
    print("Selected method: {:s}".format(sel_method))

    # Calculate the effective calibration set size
    # Note: some calibration samples may be used for model selection
    n_cal_eff = len(scores_cal)

    # Calculate the conformal p-values
    pvals = conformal.calculate_pvalues(scores_cal, scores_test)

    # Debug
    if False:
        sns.histplot(x=scores_test, hue=Y_test); plt.show()
        sns.histplot(x=pvals, hue=Y_test); plt.show()

    ####################
    # Subset selection #
    ####################
    selected = conformal.subset_selection(scores_test, selection_method)

    ###################################
    # Estimate the number of outliers #
    ###################################

    total_nout = np.sum(Y_test)
    if selected is None:
        selected_num = len(Y_test)
        selected_nout = total_nout
    else:
        selected_num = len(selected)
        selected_nout = np.sum(Y_test[selected])
    print("Total number of outliers: {:d}; selection size: {:d}; number of selected outliers: {:d}.".format(total_nout, selected_num, selected_nout))

    # Calculate confidence lower bound for the number of true outliers in closed testing procedure using Wilcoxon-Mann-Whitney local test
    # applied to conformal p-values.
    lb_wmw_2, pval_wmw_2 = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="wmw-K2", n_perm=n_perm, B=B_perm, selected=selected, table_t2=table_t2)

    lb_wmw_3, pval_wmw_3 = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="wmw-K3", n_perm=n_perm, B=B_perm,
                                                           selected=selected, table_t3=table_t3)
    lb_wmw_4, pval_wmw_4 = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="wmw-K4", n_perm=n_perm, B=B_perm,
                                                           selected=selected, table_t4=table_t4)
    lb_lmp, pval_lmp = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="lmp", n_perm=n_perm, B=B_perm,
                                                       selected=selected)
    lb_auto, pval_auto = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method=sel_method, n_perm=n_perm, B=B_perm,
                                                         selected=selected,
                                                         table_t2=table_t2, table_t3=table_t3, table_t4=table_t4, table_fisher=table_fisher)

    print("WMW-2 lower bound: {:d}.".format(lb_wmw_2))
    print("WMW-3 lower bound: {:d}.".format(lb_wmw_3))
    print("WMW-4 lower bound: {:d}.".format(lb_wmw_4))
    print("LMP lower bound: {:d}.".format(lb_lmp))
    print("Automatic lower bound: {:d}.".format(lb_auto))

    # Calculate confidence lower bound for the number of true outliers using other methods
    lb_simes, pval_simes = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="simes", selected=selected)
    lb_storey_simes, pval_storey_simes = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="storey_simes", selected=selected)
    lb_fisher, pval_fisher = conformal.estimate_num_outliers(scores_cal, scores_test, alpha, method="fisher", n_perm=n_perm, B=B_perm,
                                                             selected=selected, table_fisher=table_fisher)

    print("Simes lower bound: {:d}.".format(lb_simes))
    print("Storey-Simes lower bound: {:d}.".format(lb_storey_simes))
    print("Fisher lower bound: {:d}.".format(lb_fisher))

    #####################
    # Outlier discovery #
    #####################

    reject_bh = fdr_filter_bh(pvals, alpha=alpha)
    disc_bh = np.sum(reject_bh)
    fdp_bh, power_bh = eval_discoveries(reject_bh, Y_test)
    print("Number of discoveries (BH): {:d}. FDP: {:.3f}, POW: {:.3f}".format(disc_bh, fdp_bh, power_bh))

    reject_sbh = fdr_filter_storey_bh(pvals, alpha=alpha)
    disc_sbh = np.sum(reject_sbh)
    fdp_sbh, power_sbh = eval_discoveries(reject_sbh, Y_test)
    #print("Number of discoveries (SBH): {:d}. FDP: {:.3f}, POW:{:.3f}".format(disc_sbh, fdp_sbh, power_sbh))


    ####################
    # Assemble results #
    ####################

    res_new = pd.DataFrame({'n_cal_eff':[n_cal_eff],
                            'n_out':[total_nout], 'n_out_sel':[selected_nout], 'selected_num':[selected_num],
                            'lb_simes':[lb_simes], 'lb_storey_simes': [lb_storey_simes], 'lb_fisher':[lb_fisher],
                            'lb_auto':[lb_auto], 'lb_wmw_k2':[lb_wmw_2], 'lb_wmw_k3':[lb_wmw_3], 'lb_wmw_k4':[lb_wmw_4],
                            'lb_lmp':[lb_lmp],
                            'pval_simes':[pval_simes], 'pval_storey_simes': [pval_storey_simes], 'pval_fisher':[pval_fisher],
                            'pval_auto':[pval_auto], 'pval_wmw_k2':[pval_wmw_2], 'pval_wmw_k3':[pval_wmw_3], 'pval_wmw_k4':[pval_wmw_4],
                            'pval_lmp':[pval_lmp],
                            'disc_bh':[disc_bh], 'fdp_bh':[fdp_bh], 'power_bh':[power_bh],
                            'disc_sbh':[disc_sbh], 'fdp_sbh':[fdp_sbh], 'power_sbh':[power_sbh],
                            'selected_classifier':[sel_classifier], 'selected_method':[sel_method],
                            'seed':[random_state]})

    return res_new

# Initialize result data frame
results = pd.DataFrame({})

for r in range(num_repetitions):
    print("\nStarting repetition {:d} of {:d}...\n".format(r+1, num_repetitions))
    sys.stdout.flush()
    # Change random seed for this repetition
    random_state_new = 10*num_repetitions*random_state + r
    # Run experiment and collect results
    results_new = run_experiment(dataset, random_state_new)
    results_new = add_header(results_new)
    results_new["Repetition"] = r
    results = pd.concat([results, results_new])
    # Save results
    results.to_csv(outfile, index=False)
    print("\nResults written to {:s}\n".format(outfile))
    sys.stdout.flush()

print("\nAll experiments completed.\n")
sys.stdout.flush()

#pdb.set_trace()
