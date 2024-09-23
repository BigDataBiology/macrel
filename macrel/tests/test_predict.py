import numpy as np
from macrel import AMP_features, AMP_predict
from macrel.main import data_file
from os import path

def test_predict():
    fs = AMP_features.fasta_features('tests/peptides/expep.faa.gz')
    fsp = AMP_predict.predict(data_file("models/AMP.onnx.gz"),
                        data_file("models/Hemo.onnx.gz"),
                        fs)
    fsn = AMP_predict.predict( data_file("models/AMP.onnx.gz"),
                        data_file("models/Hemo.onnx.gz"),
                        fs, keep_negatives=True)
    assert len(fsp) < len(fsn)
    assert not np.all(fsn.is_AMP)

def test_predict_very_short():
    '''Test the prediction of very short sequences (used to crash)'''
    fs = AMP_features.fasta_features(
            path.join(path.dirname(__file__),
                    'data',
                    'very_short.faa'))
    assert len(fs) == 2
    fsn = AMP_predict.predict(data_file("models/AMP.onnx.gz"),
                        data_file("models/Hemo.onnx.gz"),
                        fs, keep_negatives=True)
    assert np.any(fsn.is_AMP)
