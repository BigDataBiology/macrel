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

dqt = "unknown_af"
npy = sys.argv[1]
nmb = sys.argv[2]
modele = sys.argv[3]
SUMM_DIR='spurio_summaries'

def run_prediction(DEFAULT_MODEL_PATH, fileX, accs, querytype):

     clf_all = np.load(DEFAULT_MODEL_PATH)
     clf, min_, max_ = clf_all[0], clf_all[1], clf_all[2]

     X = helperfunctions.load_data(fileX,
                                   filedir="npdir",
                                   n_samples=1)

         prob = gpm_predict.gpm_predict(X, accs, clf, min_, max_,
                                        outlist='predictions.txt',
                                        write_to_list=True)

     print('Protein {} (Query {}) done!'
           ' Spurious probability: {}'.format(accs[0],
           querytype,
           np.round(prob[0], 6)))

run_prediction(modele, npy, nmb, dqt)
