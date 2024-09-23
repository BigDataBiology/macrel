import pandas as pd
import numpy as np
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
    import onnxruntime as rt
    with gzip.open(model1, 'rb') as f:
        model1 = rt.InferenceSession(f.read(), providers=["CPUExecutionProvider"])

    with gzip.open(model2, 'rb') as f:
        model2 = rt.InferenceSession(f.read(), providers=["CPUExecutionProvider"])

    # limit should be 100, but let's give the user 10% leeway before we warn them
    if data.sequence.map(len).max() >= 110:
        logger = logging.getLogger('macrel')
        logger.warning('Warning: some input sequences are longer than 100 amino-acids.'
                ' Macrel models were developed and tested for short peptides (<100 amino-acids).'
                ' Applying them on longer ones will return a result, but these should be considered less reliable.')
    features = data.iloc[:, 2:]

    # predict_proba will raise an Exception if passed empty arguments
    if len(features):
        [amp_prob] = model1.run(['output_probability'], {'input_features': features.values.astype(np.float32)})
        [hemo_prob] = model2.run(['output_probability'], {'input_features': features.values.astype(np.float32)})

    else:
        amp_prob = np.array([])
        hemo_prob = np.array([])
    amp_prob = np.array([x['AMP'] for x in amp_prob])
    is_amp = (amp_prob > .5)

    hemo_prob = np.array([x['Hemo'] for x in hemo_prob])
    is_hemo = np.where(hemo_prob > .5, 'Hemo', 'NonHemo')

    final = pd.DataFrame({'Sequence': data['sequence'],
                'AMP_family':
                    (np.where(data.eval("acidicAA > basicAA"),
                         np.where(data['sequence'].map(lambda s: ('C' in s)), 'ADP', 'ALP'),
                         np.where(data['sequence'].map(lambda s: ('C' in s)), 'CDP', 'CLP'))),
                'is_AMP': is_amp,
                'AMP_probability' : amp_prob,
                'Hemolytic': is_hemo,
                'Hemolytic_probability': hemo_prob})
    if not keep_negatives:
        final = final.query('is_AMP').drop('is_AMP', axis=1)
    return final

