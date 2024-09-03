import numpy as np
from sklearn.model_selection import train_test_split

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri

r = robjects.r
numpy2ri.activate()
nout = importr("nout")

import pdb
import matplotlib.pyplot as plt

def calculate_scores_auto(X_in_train, X_in_cal, X_test, clf_occ_list={}, clf_bc_list={}, method_candidates=[], tune_size=0.25,
                          random_state=2024, alpha=0.1, selection_method=None, allow_sign_flip=True, n_perm=100, B=1000,
                          table_t2=None, table_t3=None, table_t4=None, table_fisher=None):

    assert not ( (len(clf_occ_list)==0)  and (len(clf_bc_list)==0))
    assert len(method_candidates)>0

    # Randomly split the calibration set
    X_in_cal_1, X_in_cal_2 = train_test_split(X_in_cal, train_size=tune_size, random_state=random_state)

    print("------------------------------------------------------------")
    print("Starting automatic model/method selection.")
    print("Training set size: {:d}\nCalibration_1 set size: {:d}".format(len(X_in_train), len(X_in_cal_1)))
    print("Calibration_2 set size: {:d}\nTest set size: {:d}".format(len(X_in_cal_2), len(X_test)))

    # Combine the test and calibration_2 data sets
    X_test_tmp = np.concatenate([X_in_cal_2, X_test], axis=0)

    # Initialize the list of lower bounds computed with different methods
    num_occ = len(clf_occ_list)
    num_bc = len(clf_bc_list)
    lb_occ = -np.ones((np.maximum(1,num_occ),len(method_candidates)))
    lb_bc = -np.ones((np.maximum(1,num_bc), len(method_candidates),))
    pvals_occ = np.zeros((np.maximum(1,num_occ),len(method_candidates)))
    pvals_bc = np.zeros((np.maximum(1,num_bc), len(method_candidates),))
    flipped_occ = np.zeros((np.maximum(1,num_occ),len(method_candidates)))
    flipped_bc = np.zeros((np.maximum(1,num_bc), len(method_candidates)))

    # Apply the one-class method, using calibration_1 and test_tmp
    if num_occ>0:
        print("Computing scores for {:d} OCC models and {:d} methods...".format(num_occ, len(method_candidates)))
        for occ_idx in range(num_occ):
            clf_occ_key = [*clf_occ_list][occ_idx]
            clf_occ = clf_occ_list[clf_occ_key]
            sc_c_tmp_occ, scores_test_tmp_occ = calculate_scores_occ(X_in_train, X_in_cal_1, X_test_tmp, clf_occ)
            for k_index in range(len(method_candidates)):
                method = method_candidates[k_index]
                # Check whether the signs should be flipped
                if allow_sign_flip:
                    sc_c_tmp_occ_flip, scores_test_tmp_occ_flip, flipped_occ[occ_idx, k_index] = flip_score_signs(sc_c_tmp_occ, scores_test_tmp_occ)
                else:
                    flipped_occ[occ_idx, k_index] = False

                # Perform subset selection (using flipped signs)
                selected_tmp = subset_selection(scores_test_tmp_occ_flip, selection_method)

                # Estimate the number of outliers
                lb_occ[occ_idx, k_index], pvals_occ[occ_idx, k_index] = estimate_num_outliers(sc_c_tmp_occ_flip, scores_test_tmp_occ_flip, alpha,
                                                                                              method=method, selected=selected_tmp, n_perm=n_perm, B=B,
                                                                                              table_t2=table_t2, table_t3=table_t3, table_t4=table_t4,
                                                                                              table_fisher=table_fisher)

    if num_bc>0:
        print("Computing scores for {:d} BC models and {:d} methods...".format(num_bc, len(method_candidates)))
        for bc_idx in range(num_bc):
            clf_bc_key = [*clf_bc_list][bc_idx]
            clf_bc = clf_bc_list[clf_bc_key]
            # Apply the binary method, using calibration_1 and test_tmp
            sc_c_tmp_bc, scores_test_tmp_bc = calculate_scores_binary(X_in_train, X_in_cal_1, X_test_tmp, clf_bc)

            # Evaluate the lower bounds using different methods
            for k_index in range(len(method_candidates)):
                method = method_candidates[k_index]

                # Check whether the signs should be flipped
                if allow_sign_flip:
                    sc_c_tmp_bc_flip, scores_test_tmp_bc_flip, flipped_bc[bc_idx, k_index] = flip_score_signs(sc_c_tmp_bc, scores_test_tmp_bc)
                else:
                    flipped_bc[bc_idx, k_index] = False

                # Perform subset selection (using flipped signs)
                selected_tmp = subset_selection(scores_test_tmp_bc, selection_method)

                # Estimate the number of outliers
                lb_bc[bc_idx, k_index], pvals_bc[bc_idx, k_index] = estimate_num_outliers(sc_c_tmp_bc_flip, scores_test_tmp_bc_flip, alpha,
                                                                                          method=method, selected=selected_tmp, n_perm=n_perm, B=B,
                                                                                          table_t2=table_t2, table_t3=table_t3, table_t4=table_t4,
                                                                                          table_fisher=table_fisher)

    # Use pvalues to break ties
    overall_occ = lb_occ - pvals_occ
    overall_bc = lb_bc - pvals_bc

    # Select the best type of classifier
    if np.max(overall_occ) >= np.max(overall_bc):
        idx_selection = np.unravel_index(np.argmax(overall_occ), np.array(overall_occ).shape)
        selected_method = method_candidates[idx_selection[1]]
        key_selected = [*clf_occ_list][idx_selection[0]]
        clf_occ_selected = clf_occ_list[key_selected]
        selected_classifier = "auto-occ-"+key_selected
        selected_flipped = flipped_occ[idx_selection]
        print("Selected OCC model {:s}. Best method: {:s}".format(selected_classifier, selected_method))
    else:
        idx_selection = np.unravel_index(np.argmax(overall_bc), np.array(overall_bc).shape)
        selected_method = method_candidates[idx_selection[1]]
        key_selected = [*clf_bc_list][idx_selection[0]]
        clf_bc_selected = clf_bc_list[key_selected]
        selected_classifier = "auto-bc-"+key_selected
        selected_flipped = flipped_bc[idx_selection]
        print("Selected BC model {:s}. Best method: {:s}".format(selected_classifier, selected_method))

    # Calculate final scores with selected classifier
    X_in_train_full = np.concatenate([X_in_train, X_in_cal_1], axis=0)
    print("Calculating conformity scores using selected model and method...")
    print("Training set size: {:d}\nCalibration set size: {:d}".format(len(X_in_train_full), len(X_in_cal_2)))
    print("Test set size: {:d}".format(len(X_test)))
    if selected_classifier.startswith("auto-occ"):
        scores_cal, scores_test = calculate_scores_occ(X_in_train_full, X_in_cal_2, X_test, clf_occ_selected)
    else:
        scores_cal, scores_test = calculate_scores_binary(X_in_train_full, X_in_cal_2, X_test, clf_bc_selected)

    # Flip the scores (if needed)
    if selected_flipped:
        print("Flipped the sign of the conformity scores!")
        scores_test = - scores_test
        scores_cal = - scores_cal

    print("------------------------------------------------------------\n")
    return scores_cal, scores_test, selected_classifier, selected_method



def calculate_scores_occ(X_in_train, X_in_cal, X_test, clf):
    print("Fitting a one-class classifier using {:d} inlier data points... ".format(len(X_in_train)), end='')
    clf.fit(X_in_train)
    print("Fitting completed.")
    scores_cal = -clf.score_samples(X_in_cal)
    scores_test = -clf.score_samples(X_test)
    scores_cal += np.random.uniform(low=0.0, high=1e-6, size=scores_cal.shape)
    scores_test += np.random.uniform(low=0.0, high=1e-6, size=scores_test.shape)
    return scores_cal, scores_test


def calculate_scores_binary(X_in_train, X_in_cal, X_test, clf):
    # Create mixed training set
    X_pu_train_mixed = np.concatenate([X_in_cal, X_test], axis=0)
    mixed_is_test = np.concatenate([np.zeros((len(X_in_cal),)), np.ones((len(X_test),))], axis=0)

    # Create semi-supervised data set
    X_pu_train = np.concatenate([X_in_train, X_pu_train_mixed], axis=0)
    Y_pu_train = np.concatenate([np.zeros((X_in_train.shape[0],)), np.ones((X_pu_train_mixed.shape[0],))], axis=0)

    print("Fitting a binary classifier using {:d} inliers and {:d} unlabeled data points... ".format(len(X_in_train), len(X_pu_train_mixed)), end='')
    clf.fit(X_pu_train, Y_pu_train)
    print("Fitting completed.")

    # Compute the conformity scores
    scores_mixed = clf.predict_proba(X_pu_train_mixed)[:,1]
    scores_mixed += np.random.uniform(low=0.0, high=1e-6, size=scores_mixed.shape)
    scores_cal = scores_mixed[mixed_is_test==0]
    scores_test = scores_mixed[mixed_is_test==1]
    return scores_cal, scores_test


def flip_score_signs(scores_cal, scores_test):
    pvals = calculate_pvalues(scores_cal, scores_test)
    if np.median(pvals) > 0.5:
        scores_test = - scores_test
        scores_cal = - scores_cal
        flipped = True
    else:
        flipped = False
    return scores_cal, scores_test, flipped


def calculate_pvalues(scores_cal, scores_test):
    n_test = len(scores_test)
    pvals = np.array([(1+np.sum(scores_test[i]<=scores_cal)) / (1+len(scores_cal)) for i in range(n_test)])
    return pvals


def estimate_num_outliers_adaptive(scores_cal, scores_test, alpha, method_list, tune_size=0.25, selected=None, random_state=2024,
                                   n_perm=100, B=1000, table_t2=None, table_t3=None, table_t4=None, table_fisher=None):
    # Randomly split the calibration set
    scores_cal_1, scores_cal_2 = train_test_split(scores_cal, train_size=tune_size, random_state=random_state)
    # Combine the test and calibration_2 data sets
    scores_test_tmp = np.concatenate([scores_cal_2, scores_test], axis=0)

    num_methods = len(method_list)
    lb_list = [None]*num_methods
    for i in range(num_methods):
        method = method_list[i]
        lb_list[i], _ = estimate_num_outliers(scores_cal_1, scores_test_tmp, alpha, method=method, selected=selected,
                                              n_perm=n_perm, B=B, table_t2=table_t2, table_t3=table_t3, table_t4=table_t4, table_fisher=table_fisher)

    # Pick the best method
    i_selected = np.argmax(lb_list)
    method_selected = method_list[i_selected]
    lb, _ = estimate_num_outliers(scores_cal_2, scores_test, alpha, method=method_selected, selected=selected, n_perm=n_perm, B=B,
                                  table_t2=table_t2, table_t3=table_t3, table_t4=table_t4, table_fisher=table_fisher)
    return lb


def subset_selection(scores_test, selection_method):
    n = len(scores_test)

    if (selection_method == "none") or (selection_method is None):
        selected = None
    else:
        selection_cutoff = float(selection_method.split("-")[1]) / 100
        top_n = int(np.ceil(n*selection_cutoff))
        selected = np.argsort(scores_test)[::-1][0:top_n]
    return selected


def estimate_num_outliers(scores_cal, scores_test, alpha, method="wmw", selected=None, n_perm=100, B=1000, table_t2=None, table_t3=None,
                          table_t4=None, table_fisher=None):
    m = len(scores_cal)

    if selected is None:
        selected = robjects.NULL
    else:
        selected = selected + 1

    if method in ["wmw", "wmw-K2", "wmw-K3", "wmw-K4", "lmp", "fisher"]:
        if (method=="wmw") or (method=="wmw-K2"):
            statistic = "T2"
            if table_t2 is None:
                critical_values = []
            else:
                critical_values = np.array(table_t2.loc[table_t2['m'] == m]['critical'])
        elif method=="wmw-K3":
            statistic = "T3"
            if table_t3 is None:
                critical_values = []
            else:
                critical_values = np.array(table_t3.loc[table_t3['m'] == m]['critical'])
        elif method=="wmw-K4":
            statistic = "T4"
            if table_t4 is None:
                critical_values = []
            else:
                critical_values = np.array(table_t4.loc[table_t4['m'] == m]['critical'])
        elif method=="fisher":
            statistic = "fisher"
            if len(scores_cal)>=200:
                n_perm = -1
            if table_fisher is None:
                critical_values = []
            else:
                critical_values = np.array(table_fisher.loc[table_fisher['m'] == m]['critical'])
        elif method=="lmp":
            statistic = "lmp"
            n_perm = -1
            critical_values = []

        if method == "wmw":
            tmp = nout.d_selection_higher(scores_cal, scores_test, S=selected, local_test="wmw", alpha=alpha, n_perm=n_perm, B=B, critical_values=critical_values)
        elif method in ["wmw-K2", "wmw-K3", "wmw-K4"]:
            k = int(method.replace("wmw-K", ""))
            tmp = nout.d_selection_higher(scores_cal, scores_test, S=selected, local_test="higher", k=k, alpha=alpha, n_perm=n_perm, B=B, critical_values=critical_values)
        elif method == "lmp":
            tmp = nout.d_selection_G(scores_cal, scores_test, S=selected, monotone=True, fit_method="betamix", prop_cal=0.5, alpha=alpha, n_perm=n_perm, B=B, B_MC=10^4)
        elif method == "fisher":
            tmp = nout.d_selection_fisher(scores_cal, scores_test, S=selected, n_perm = 0)
        else:
            print("Error: Unknown method!")
            exit(0)

        if selected is None:
            lb = tmp[0][0]
            pval = tmp[1][0]
        else:
            lb = tmp[0][0]
            pval = tmp[1][0]


    elif method=="simes":
        tmp = nout.d_selection_simes(scores_cal, scores_test, S=selected, alpha=alpha)
        lb = tmp[0][0]
        pval = tmp[1][0]

    elif method=="storey_simes":
        tmp = nout.d_selection_storey(scores_cal, scores_test, S=selected, alpha=alpha)
        lb = tmp[0][0]
        pval = tmp[1][0]

#    elif method=="bh":
#        lb = nout.d_benjhoch(S_Y=scores_test, S_X=scores_cal, alpha=alpha)[0]
#    elif method=="storey_bh":
#        lb = nout.d_StoreyBH(S_Y=scores_test, S_X=scores_cal, alpha=alpha)[0]
    else:
        print("Error! Unknown method {:s}.".format(method))
        exit(-1)

    return int(lb), float(pval)
