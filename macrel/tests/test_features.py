import numpy as np
import pandas as pd
from macrel import macrel_features
from macrel.AMP_features import features

def test_macrel_features():
    assert np.allclose(
                macrel_features.compute_all('FLPVLAGLTPSIVPKLVCLLTKKC'),
                [2.914112328517451,
                 9.66778564453125,
                 154.16666666666666,
                 42.63333333333334,
                 -1.2358333333333336,
                 0.3929166666666667,
                 0.5199225805121164])
    assert np.allclose(
                macrel_features.compute_all('HGRHVTLKDIVLDLQPPDPVGLHCYEQLVDSSEDEVDEVDGQDSQPLKQHFQIVTC'),
                [-8.11413482164507,
                  4.14520263671875,
                  91.9642857142857,
                  39.5142857142857,
                  2.0424999999999995,
                  -0.04285714285714287,
                  0.4804582672677425],)

def test_features():
    computed = features('tests/peptides/expep.faa.gz')
    expected = pd.read_table('macrel/tests/data/features.tsv.gz', sep='\t', index_col=0)

    for c in computed.columns:
        if c in ['sequence', 'group']:
            assert np.all(computed[c] == expected[c])
        else:
            print(c)
            assert np.allclose(computed[c], expected[c])

