#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import gzip
import lzma
import pandas as pd
import numpy as np

from .fasta import fasta_iter
from .database import GROUPS_SA, GROUPS_HB
from .macrel_features import MacrelFeatures

def openfile(ifile):
    if ifile.endswith('.gz'):
        records = SeqIO.parse(gzip.open(ifile, 'rt', encoding='utf-8'), 'fasta')
    elif ifile.endswith('.xz'):
        records = SeqIO.parse(lzma.open(ifile, 'rt', encoding='utf-8'), 'fasta')
    else:
        records = SeqIO.parse(ifile, 'fasta')

    return records
    
def features(ifile):
    seqs = []
    headers = []
    encodings = []
    aaComp = []
    desc_features = [] 
    for h, seq in fasta_iter(ifile):
        features = MacrelFeatures(seq)
        seqs.append(features.checkseq()) 
        headers.append(h)
        encodings.append(features.ctdd(GROUPS_SA + GROUPS_HB))
        aaComp.append(features.amino_acid_composition())
        desc_features.append(features.compute_all())

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
    feat = features(ifile)
    feat.to_csv(ofile, sep='\t', index_label='Access')

if __name__ == '__main__':
    main(sys.argv)
