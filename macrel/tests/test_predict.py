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


_REPO_ROOT = path.dirname(path.dirname(path.dirname(path.abspath(__file__))))


def _predict_example_peptides():
    fs = AMP_features.fasta_features(
            path.join(_REPO_ROOT, 'tests', 'peptides', 'expep.faa.gz'))
    return AMP_predict.predict(data_file("models/AMP.onnx.gz"),
                        data_file("models/Hemo.onnx.gz"),
                        fs, keep_negatives=True)


def test_predict_probabilities_valid():
    '''Output probabilities must be valid probabilities

    onnxruntime 1.26 returned values outside [0, 1] (e.g., negative), which
    silently turned every prediction into a negative.
    '''
    fsn = _predict_example_peptides()
    for col in ['AMP_probability', 'Hemolytic_probability']:
        assert np.all(fsn[col] >= 0.), col
        assert np.all(fsn[col] <= 1.), col
    assert fsn.is_AMP.any()
    assert not fsn.is_AMP.all()
    assert (fsn.Hemolytic == 'Hemo').any()
    assert (fsn.Hemolytic == 'NonHemo').any()


def test_predict_matches_reference():
    '''Predictions must match the checked-in reference output

    This catches numerical regressions in the model runtime (e.g., onnxruntime
    updates) without requiring the external tools used by run-tests.sh.
    '''
    import pandas as pd
    expected = pd.read_table(
            path.join(_REPO_ROOT, 'tests', 'peptides', 'expected.prediction'),
            comment='#',
            index_col='Access')
    fsn = _predict_example_peptides()
    fsp = fsn.query('is_AMP')

    assert set(fsp.index) == set(expected.index)
    fsp = fsp.loc[expected.index]
    assert np.all(fsp.Sequence == expected.Sequence)
    assert np.all(fsp.AMP_family == expected.AMP_family)
    assert np.all(fsp.Hemolytic == expected.Hemolytic)
    # Reference values are written with 3 decimal places
    assert np.allclose(fsp.AMP_probability, expected.AMP_probability, atol=1e-3, rtol=0)
    assert np.allclose(fsp.Hemolytic_probability, expected.Hemolytic_probability, atol=1e-3, rtol=0)
