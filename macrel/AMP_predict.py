import pandas as pd
import numpy as np
import rpy2.robjects
from rpy2.robjects import numpy2ri, pandas2ri

def predict(model1, model2, data):
    pandas2ri.activate()
    numpy2ri.activate()

    r = rpy2.robjects.r
    r.require('randomForest')
    # When the model returns probability == 0.5, then the choice of class is
    # random, so the seed has an impact:
    r('set.seed(95014)')

    rpy2.robjects.globalenv['model1'] = model1
    rpy2.robjects.globalenv['model2'] = model2

    rpy2.robjects.globalenv['data'] = data

    r('model1 <- readRDS(model1)')
    r('model2 <- readRDS(model2)')

    is_amp = pd.Categorical(r('predict(model1, data)'))
    amp_prob = r('predict(model1, data, type="prob")[,1]')
    is_hemo = r('predict(model2, data)')
    hemo_prob = r('predict(model2, data, type="prob")[,1]')

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


def main(args):
    if len(args) < 5:
        import sys
        sys.stderr.write("This is an internal FACS script and is not meant to be used independently")
        sys.exit(1)

    ifile = args[1]
    model1 = args[2]
    model2 = args[3]
    ofile = args[4]

    data = pd.read_table(ifile, index_col=0)
    rfinal = predict(model1, model2, data)
    rfinal.to_csv(ofile, sep='\t', header=False)


if __name__ == '__main__':
    from sys import argv
    main(argv)


