import numpy as np
from macrel import macrel_features


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
                macrel_features.compute_all('MHGRHVTLKDIVLDLQPPDPVGLHCYEQLVDSSEDEVDEVDGQDSQPLKQHFQIVTC'),
                [-8.11413482164507,
                  4.14520263671875,
                  91.9642857142857,
                  39.5142857142857,
                  2.0424999999999995,
                  -0.04285714285714287,
                  0.4804582672677425],)
    assert np.allclose(
                macrel_features.compute_all('HGRHVTLKDIVLDLQPPDPVGLHCYEQLVDSSEDEVDEVDGQDSQPLKQHFQIVTC'),
                [-8.11413482164507,
                  4.14520263671875,
                  91.9642857142857,
                  39.5142857142857,
                  2.0424999999999995,
                  -0.04285714285714287,
                  0.4804582672677425],)
