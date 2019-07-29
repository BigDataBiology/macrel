import numpy as np
import copy
import os

import helperfunctions as helperfunctions

def gpm_predict(X, accs, clf, min_, max_, outlist='predictions.txt',
                write_to_list=True):
    '''
    This function is used to estimate the spuriosity of a protein.
    It is called a) to estimate test error
             b) by spurio after each homology search
             c) on its own, to estimate a probabilities from np files.

    Input: - numpy array of true proteins, dimensions n_proteins x 4:
             Accession Name, Stops_per_msa, Homologous seqs, Sequence length
           - numpy array with classifier and normalisation values

    Output: - 1d array of predicted values
            - optional: Writes a line for each protein into an output txt file
    '''

    # If a single protein is to be predicted (e.g. during SPURIO main run),
    # reshape!
    if np.shape(X) == (4,):
        X = X.reshape((1, 4))

    Xorig = copy.deepcopy(X)

    # Preprocess data
    X = helperfunctions.preprocess_data(X, min_, max_)

    # Predict labels
    y_preds = clf.predict_proba(X)[:, 0]

    if write_to_list:
        # Sort results by spuriosity score
        idx_sort = np.argsort(y_preds)
        y_preds_sorted = y_preds[idx_sort]
        accs_sorted = np.array(accs)[idx_sort]
        Xorig_sorted = Xorig[idx_sort]

        # Save to txtfile
        if not os.path.isfile(outlist):
            with open(outlist, 'a') as outlist_:
                outlist_.write('UniProt Accession \t Spuriosity score \t STOPS per AA'
                               '\t Homologous Matches \t Sequence Length \n \n')
        for i in range(len(y_preds)):
            res_full = ('\t'.join([accs_sorted[i].ljust(10),
                                   '%.3f' % np.abs(y_preds_sorted[i]),
                                   ('%.5f' % np.abs(Xorig_sorted[i, 0])).ljust(6),
                                   str(int(Xorig_sorted[i, 1])),
                                   str(int(Xorig_sorted[i, 2])),
                                   ]))
            with open(outlist, 'a') as outlist_:
                print(res_full, file=outlist_)

    return y_preds
