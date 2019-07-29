# This is a set of functions used by gpm_train, gpm_test or both.

import numpy as np
import glob
import os
from sklearn.metrics import (roc_auc_score, accuracy_score, mean_squared_error,
                             roc_curve, auc, precision_recall_curve,
                             average_precision_score)
from scipy import interp


def load_data(name_prefix, filedir, n_samples):
    '''
    This function is used to collect the .npy output from the homology
    search and return it as one big matrix [n_proteins x n_features].
    '''
    feature_files = glob.glob('{}*.npy'.format(os.path.join(filedir, name_prefix)))[:n_samples]
    print(len(feature_files))
    count = 0
    accessions = []
    feature_mat = np.zeros((len(feature_files), 4))
    for n, feature_file in enumerate(feature_files):
        protein_summary = np.load(feature_file)
        acc = feature_file.split('.')[-2].split('_')[-1]
        accessions.append(acc)
        feature_mat[n] = protein_summary
        count += 1
        if count % 200 == 0:
            print(count)
    notnan_idx = ~np.isnan(feature_mat[:, -1])

    return accessions, feature_mat[notnan_idx]


def plot_pr_cv(ytest, ypred, aps=[], recalls=[], precisions=[]):
    #TODO Docstring
    import matplotlib.pyplot as plt
    plt.figure(2)
    ax = plt.subplot(111)
    mean_precision = np.linspace(0, 1, 100)
    precision, recall, _ = precision_recall_curve(ytest, ypred)
    recalls.append(interp(mean_precision, recall, precision))
    average_precision = average_precision_score(ytest, ypred)
    aps.append(average_precision)
    precisions.append(interp(mean_precision, precision, recall))
    ax.step(recall[:-1], precision[:-1], alpha=0.2, where='post')
    return aps, recalls, precisions


def plot_roc_cv(ytest, ypred, tprs=[], fprs=[], aucs=[]):
    #TODO Docstring
    import matplotlib.pyplot as plt
    plt.figure(1)
    ax1 = plt.subplot(111)
    mean_fpr = np.linspace(0, 1, 100)
    fpr, tpr, thresholds = roc_curve(ytest, ypred)
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    ax1.plot(fpr, tpr, lw=1, alpha=0.3)

    print(fprs, tprs, aucs)
    return tprs, fprs, aucs


def finish_ap_plot(tprs, fprs, aucs):
    #TODO Docstring
    import matplotlib.pyplot as plt
    plt.figure(2)
    ax = plt.subplot(1, 1, 1)
    mean_tpr = np.linspace(0, 1, 100)
    mean_fpr = np.mean(fprs, axis=0)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.step(mean_tpr[:-1], mean_fpr[:-1], color='b',
            label=r'Mean PR (AP = %0.3f $\pm$ %0.3f)' % (mean_auc, std_auc),
            lw=2, alpha=0.8)
    std_fpr = np.std(fprs, axis=0)
    fprs_upper = np.minimum(mean_fpr + std_fpr, 1)
    fprs_lower = np.maximum(mean_fpr - std_fpr, 0)
    plt.fill_between(mean_tpr[:-1], fprs_lower[:-1], fprs_upper[:-1],
                     color='grey', alpha=0.2, label=r'$\pm$ 1 std. dev.')

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('Recovery')
    plt.ylabel('Precision')
    plt.legend(loc="lower left")
    plt.savefig('RPcurve.png')


def finish_roc_plot(tprs, fprs, aucs):
    #TODO Docstring
    import matplotlib.pyplot as plt
    plt.figure(1)
    ax1 = plt.subplot(1, 1, 1)
    ax1.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey', alpha=0.8)
    mean_fpr = np.linspace(0, 1, 100)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax1.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.3f $\pm$ %0.3f)' % (mean_auc, std_auc),
             lw=2, alpha=0.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.2,
                     label=r'$\pm$ 1 std. dev.')

    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig('ROCcurve.png')
    plt.close()


def preprocess_data(X, min_all, max_all):
    #TODO Docstring
    X = add_pseudocount(X)
    X = keep_desired_features(X)
    X = logit(X)
    X = normalise(X, min_all, max_all)

    return X


def add_pseudocount(X):
    #TODO Docstring
    if np.shape(X) == (4,):
        X = X.reshape(1, 4)
    if X[0, 0] == 0:
        X = np.array(X, dtype=float)
    X[: ,0] += (1 / X[:, 3])
    return X


def keep_desired_features(X):
    #TODO Docstring
    return X[:, [0, 1, 2]]


def normalise(X, min_all, max_all):
    '''Normalise all features to 0-1'''
    if np.shape(X) == (3,):
        X = X.reshape((1, 3))
    for i in range(3):
        X[:, i] = (X[:, i] - min_all[i]) / (max_all[i] - min_all[i])
    return X


def unnormalise(X, min_all, max_all):
    '''Unused'''
    # TODO Then why is it here?
    Xunnorm = (X * (max_all[2] - min_all[2])) + min_all[2]
    return Xunnorm


def logit(X):
    '''Apply logarithm to all features'''
    if np.shape(X) == (3,):
        X = X.reshape((1, 3))
    X = np.log(X)
    return X


def eval_result(p_predicted, y_true):
    '''
    Given y_true and y_predicted, this function
    returns the area under the curve, accuracy
    and mean squared error
    '''

    mse = mean_squared_error(y_true, p_predicted)
    auc = roc_auc_score(y_true, p_predicted)
    acc = accuracy_score(y_true, p_predicted.round())

    return auc, acc, mse
