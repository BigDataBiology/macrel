'''This is the main script to create fancy pictures and prepare protein data
for further analysis.'''

import argparse
import os
import sys
import numpy as np
import helperfunctions as helperfunctions
import gpm_predict as gpm_predict

import matplotlib as mpl
mpl.use('Agg')

MODEL = sys.argv[1]

def predict_spuriosity(querytype, trained_model_path, npy_files_path,
                      n_samples=None):
    '''
    Classify any data previously processed by SPURIO

    Input: querytype: The querytype chosen in the SPURIO run (str, e.g. 'sp_train').
                      This will lead to the right .npy files in /output/np
           n_samples: Number of samples to use. (int or 'all')

    Output: predictions_querytype.txt: summary of processed proteins and
                                       their spuriosity score
            data plot in matplotlib
    '''

    # Load data (will take around 1min per 5k datapoints.)
    names, X = helperfunctions.load_data(querytype,
                                         filedir=npy_files_path,
                                         n_samples=n_samples)
    print('Data loaded. Classifying...')

    # Load classifier
    clf_full = np.load(trained_model_path)
    clf, min_, max_ = clf_full[0], clf_full[1], clf_full[2]

    # Predict labels
    #TODO This calculation is unused. Is it supposed to be?
    y_pred = gpm_predict.gpm_predict(X, names, clf, min_, max_,
                  outlist = 'predictions_newer_{}'.format(querytype),
                  write_to_list = True)

if __name__ == "__main__":
    predict_spuriosity('spuriomod', MODEL, 'npdir')
