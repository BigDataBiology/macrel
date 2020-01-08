#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import pandas as pd
import numpy as np
import rpy2
import rpy2.robjects
from rpy2.robjects import numpy2ri
from .fasta import fasta_iter
numpy2ri.activate()
r = rpy2.robjects.r
r.library('Peptides')

GROUPS_SA = ['ALFCGIVW', 'RKQEND', 'MSPTHY'] #solventaccess
GROUPS_HB = ['ILVWAMGT', 'FYSQCN', 'PHKEDR'] # HEIJNE&BLOMBERG1979

def ctdd(sequence, groups):
    code = []
    for group in groups:
        for i, aa in enumerate(sequence):
            if aa in group:
                code.append((i + 1)/len(sequence) * 100)
                break
        else:
            code.append(0)
    return code


def features(ifile):
    groups = [set(g) for g in (GROUPS_SA+GROUPS_HB)]
    seqs = []
    headers = []
    encodings = []
    for h,seq in fasta_iter(ifile):
        if seq[-1] == '*':
            seq = seq[:-1]
        seqs.append(seq)
        headers.append(h)
        encodings.append(ctdd(seq, groups))

    # We can do this inside the loop so that we are not forced to pre-load all
    # the sequences into memory. However, it becomes much slower
    rpy2.robjects.globalenv['seq'] = seqs
    aaComp = r('aaComp(seq)')
    rfeatures = r('''
    ch <- charge(seq=seq, pH=7, pKscale="EMBOSS")
    pI <- pI(seq=seq, pKscale="EMBOSS")
    aIndex <- aIndex(seq=seq)
    instaIndex <- instaIndex(seq=seq)
    boman <- boman(seq=seq)
    hydrophobicity <- hydrophobicity(seq=seq, scale="Eisenberg")
    hmoment <- hmoment(seq=seq, angle=100, window=11)
    cbind(ch, pI, aIndex, instaIndex, boman, hydrophobicity, hmoment)
    ''')
    aaComp = np.array([np.array(v) for v in aaComp])
    aaComp = aaComp[:,:,1]

    features = np.hstack([aaComp, rfeatures, encodings])
    # The column names must match those in the saved model
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
            "SA.G1.residue0",
            "SA.G2.residue0",
            "SA.G3.residue0",
            "hb.Group.1.residue0",
            "hb.Group.2.residue0",
            "hb.Group.3.residue0",
            ])
    features.insert(0, 'group', 'Unk')
    features.insert(0, 'sequence', seqs)
    return features

def main(args):
    if len(args) < 3:
        sys.stderr.write("This is an internal FACS script and is not meant to be used independently")
        sys.exit(1)

    ifile = args[1]
    ofile = args[2]
    features = features(ifile)
    features.to_csv(ofile, sep='\t', index_label='access')

if __name__ == '__main__':
    main(sys.argv)
