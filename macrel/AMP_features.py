#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import pandas as pd
import numpy as np

from .fasta import fasta_iter
from .database import GROUPS_SA, GROUPS_HB
from .macrel_features import ctdd, amino_acid_composition, normalize_seq, compute_all


def features(ifile):
    seqs = []
    headers = []
    features = []
    for h, seq in fasta_iter(ifile):
        seq = normalize_seq(seq)
        seqs.append(seq)
        headers.append(h)
        features.append(
                np.concatenate((
                    amino_acid_composition(seq),
                    compute_all(seq),
                    ctdd(seq, GROUPS_SA + GROUPS_HB))))

    features = pd.DataFrame(features, index=headers, columns=[
            "tinyAA",
            "smallAA",
            "aliphaticAA",
            "aromaticAA",
            "nonpolarAA",
            "polarAA",
            "chargedAA",
            "basicAA",
            "acidicAA",
            "charge",
            "pI",
            "aindex",
            "instaindex",
            "boman",
            "hydrophobicity",
            "hmoment",
            "SA.Group1.residue0",
            "SA.Group2.residue0",
            "SA.Group3.residue0",
            "HB.Group1.residue0",
            "HB.Group2.residue0",
            "HB.Group3.residue0",
            ])
    features.insert(0, 'group', 'Unk')
    features.insert(0, 'sequence', seqs)
    return features

def main(args):
    if len(args) < 3:
        sys.stderr.write("This is an internal MACREL script and is not meant to be used independently")
        sys.exit(1)

    ifile = args[1]
    ofile = args[2]
    feat = features(ifile)
    feat.to_csv(ofile, sep='\t', index_label='Access')

if __name__ == '__main__':
    main(sys.argv)
