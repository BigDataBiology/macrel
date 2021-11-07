import numpy as np
import pandas as pd
from macrel import macrel_features
from macrel.AMP_features import fasta_features

def test_macrel_features():
    assert np.allclose(
                macrel_features.compute_all('FLPVLAGLTPSIVPKLVCLLTKKC'),
                    [0.2916666666666667,
                      0.5416666666666666,
                      0.4583333333333333,
                      0.041666666666666664,
                      0.75,
                      0.25,
                      0.125,
                      0.125,
                      0.0,
                      2.914112328517451,
                      9.66778564453125,
                      154.16666666666666,
                      42.63333333333334,
                      -1.2358333333333336,
                      0.3929166666666666,
                      0.5199225805121164,
                      4.166666666666666,
                      62.5,
                      12.5,
                      8.333333333333332,
                      4.166666666666666,
                      12.5])
    assert np.allclose(
                macrel_features.compute_all('HGRHVTLKDIVLDLQPPDPVGLHCYEQLVDSSEDEVDEVDGQDSQPLKQHFQIVTC'),
                    [0.17857142857142858,
                      0.5178571428571429,
                      0.26785714285714285,
                      0.10714285714285714,
                      0.4642857142857143,
                      0.5357142857142857,
                      0.3392857142857143,
                      0.125,
                      0.21428571428571427,
                      -8.11413482164507,
                      4.14520263671875,
                      91.9642857142857,
                      39.5142857142857,
                      2.0424999999999995,
                      -0.04285714285714285,
                      0.4804582672677426,
                      3.571428571428571,
                      5.357142857142857,
                      1.7857142857142856,
                      3.571428571428571,
                      26.785714285714285,
                      1.7857142857142856])

def test_features():
    computed = fasta_features('tests/peptides/expep.faa.gz')
    expected = pd.read_table('macrel/tests/data/features.tsv.gz', sep='\t', index_col=0)

    for c in computed.columns:
        if c in ['sequence', 'group']:
            assert np.all(computed[c] == expected[c])
        else:
            print(c)
            assert np.allclose(computed[c], expected[c])

