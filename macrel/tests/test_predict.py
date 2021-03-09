import numpy as np
from macrel import AMP_features, AMP_predict
from macrel.main import data_file

def test_predict():
    fs = AMP_features.features('tests/peptides/expep.faa.gz')
    fsp = AMP_predict.predict( data_file("models/AMP.pkl.gz"),
                        data_file("models/Hemo.pkl.gz"),
                        fs)
    fsn = AMP_predict.predict( data_file("models/AMP.pkl.gz"),
                        data_file("models/Hemo.pkl.gz"),
                        fs, keep_negatives=True)
    assert len(fsp) < len(fsn)
    assert not np.all(fsn.is_AMP)

def test_predict_lt50():
    fs = AMP_features.features('tests/peptides/expep.faa.gz')
    fsp = AMP_predict.predict( data_file("models/AMP_lt50.pkl.gz"),
                        data_file("models/Hemo_lt50.pkl.gz"),
                        fs)
    fsn = AMP_predict.predict( data_file("models/AMP_lt50.pkl.gz"),
                        data_file("models/Hemo_lt50.pkl.gz"),
                        fs, keep_negatives=True)
    assert len(fsp) < len(fsn)
    assert not np.all(fsn.is_AMP)
