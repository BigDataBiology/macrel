import pandas as pd
import numpy as np
import pickle
import gzip

def predict(model1, model2, data):
    model1 = pickle.load(gzip.open(model1, 'rb'))
    model2 = pickle.load(gzip.open(model2, 'rb'))

    features = data.iloc[:, 3:]

    # predict_proba will raise an Exception if passed empty arguments
    if len(features):
        amp_prob = model1.predict_proba(features).T[0]
        hemo_prob = model2.predict_proba(features).T[0]
    else:
        amp_prob = np.array([])
        hemo_prob = np.array([])
    is_amp = np.where(amp_prob > .5, model1.classes_[0], model1.classes_[1])
    is_hemo = np.where(hemo_prob > .5, model2.classes_[0], model2.classes_[1])

    final = pd.DataFrame({'Sequence': data['sequence'],
                'AMP_family':
                    (np.where(data.eval("acidicAA > basicAA"),
                         np.where(data['sequence'].map(lambda s: ('C' in s)), 'ADP', 'ALP'),
                         np.where(data['sequence'].map(lambda s: ('C' in s)), 'CDP', 'CLP'))),
                'is_AMP': is_amp,
                'AMP_probability' : amp_prob,
                'Hemolytic': is_hemo,
                'Hemolytic_probability': hemo_prob})
    rfinal = final.query('is_AMP == "AMP"').drop('is_AMP', axis=1)
    return rfinal

