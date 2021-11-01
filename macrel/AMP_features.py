#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import pandas as pd
import numpy as np
from .feature_functions import *
from .fasta import fasta_iter
from .database import GROUPS_SA, GROUPS_HB

def features(ifile):
    groups = [set(g) for g in (GROUPS_SA+GROUPS_HB)]
    seqs = []
    headers = []
    encodings = []
    aaComp = []
    desc_features = []
    for h, seq in fasta_iter(ifile):
        if seq[-1] == '*':
            seq = seq[:-1]
        if seq[0] == 'M':
            seq = seq[1:]
        seqs.append(seq)
        headers.append(h)
        encodings.append(ctdd(seq, groups))
        aaComp.append(amino_acid_composition(seq))
        if len(seq) < 3:
            import logging
            logger = logging.getLogger('macrel')
            logger.warning("Warning: input sequence '{}' is shorter longer than 3 amino acids."
                " Macrel models were developed and tested for peptides with at least 10 amino acids and "
                " results on very short peptides these should be considered unreliable.".format(h))

            encodings[-1] *= 0
            aaComp[-1] *= 0
            # This is a major hack, but otherwise, the R code will fail
            if len(seq) == 1:
                seq += 'XX'
            elif len(seq) == 2:
                seq += 'X'
            seqs[-1] = seq
        ch = pep_charge(seq, 7)
        pI = isoelectric_point(seq)
        aIndex = aliphatic_index(seq)
        instaIndex = instability_index(seq)
        iboman = boman_index(seq)
        hydroph = hydrophobicity(seq)
        h_moment = hmoment(seq, angle=100, window=11)
        desc_features.append([ch, pI, aIndex, instaIndex, iboman, hydroph, h_moment])

    features = np.hstack([aaComp, desc_features, encodings])

    # This is arguably a Pandas bug (at least inconsistency), but
    # pd.DataFrame([], ...) works, while pd.DataFrame(np.array([]), ...) does
    # not:
    if len(features) == 0:
        features = []

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
    features = features(ifile)
    features.to_csv(ofile, sep='\t', index_label='Access')

if __name__ == '__main__':
    main(sys.argv)
