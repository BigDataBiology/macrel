import pandas as pd
import numpy as np
import pickle
import gzip
import logging

def predict(model1, model2, data, keep_negatives=False):
    '''
    Run prediction models

    Parameters
    ----------
    model1: filepath
    model2: filepath
    data: pandas DataFrame
    keep_negatives: boolean, optional (default: False)
        whether to keep negative results

    Returns
    -------
    Table with prediction labels
    '''
    model1 = pickle.load(gzip.open(model1, 'rb'))
    model2 = pickle.load(gzip.open(model2, 'rb'))

    # limit should be 100, but let's give the user 10% leeway before we warn them
    if data.sequence.map(len).max() >= 110:
        logger = logging.getLogger('macrel')
        logger.warning('Warning: some input sequences are longer than 100 amino-acids.'
                ' Macrel models were developed and tested for short peptides (<100 amino-acids).'
                ' Applying them on longer ones will return a result, but these should be considered less reliable.')
    elif data.sequence.map(len).min() <= 8:
        logger = logging.getLogger('macrel')
        logger.warning('Warning: some input sequences are smaller than 10 amino-acids.'
                ' Macrel models were developed and tested for short peptides (>= 10 amino-acids).'
                ' Applying them on shorter ones will return all of them as non-amps.')

    namp_seq = data[data['group'] == 'NAMP']['sequence'].tolist()
    namp_header = data[data['group'] == 'NAMP'].index.tolist()
    data.dropna(inplace=True)
    features = data.iloc[:, 3:]
 
    # predict_proba will raise an Exception if passed empty arguments
    if len(features):
        amp_prob = model1.predict_proba(features).T[0]
        hemo_prob = model2.predict_proba(features).T[0]
    else:
        amp_prob = np.array([])
        hemo_prob = np.array([])
    
    is_amp = np.where(amp_prob > .5, model1.classes_[0], model1.classes_[1])
    is_amp = (is_amp == "AMP")
    is_hemo = np.where(hemo_prob > .5, model2.classes_[0], model2.classes_[1])
    
    df = pd.DataFrame()
    df['Sequence'] = namp_seq
    df['AMP_family'] = ['-' for x in namp_header]
    df['is_AMP'] = [0 for x in namp_header]
    df['AMP_probability'] = [0 for x in namp_header]
    df['Hemolytic'] = ['-' for x in namp_header]
    df['Hemolytic_probability'] = ['-' for x in namp_header]
    df.index = namp_header

    final = pd.DataFrame({'Sequence': data['sequence'],
                'AMP_family':
                    (np.where(data.eval("acidicAA > basicAA"),
                         np.where(data['sequence'].map(lambda s: ('C' in s)), 'ADP', 'ALP'),
                         np.where(data['sequence'].map(lambda s: ('C' in s)), 'CDP', 'CLP'))),
                'is_AMP': is_amp,
                'AMP_probability' : amp_prob,
                'Hemolytic': is_hemo,
                'Hemolytic_probability': hemo_prob})

    final = pd.concat([final, df])

    if not keep_negatives:
        final = final[final['is_AMP'] == 1].drop('is_AMP', axis=1)
        final.Hemolytic_probability = pd.to_numeric(final.Hemolytic_probability)
    elif keep_negatives:
        newhemoprob = [ float(f'{x:.3f}') if x != '-' else '-' for x in final.Hemolytic_probability ]
        final.Hemolytic_probability = newhemoprob
        newisamp = [ 1 if x == 1 else 0 for x in final.is_AMP  ]
        final.is_AMP = newisamp

    return final

