#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import logging
import pandas as pd
import numpy as np

from .fasta import fasta_iter
from .macrel_features import normalize_seq, compute_all

# The 20 standard amino acids. Sequences containing any other character
# (ambiguity codes such as X/B/Z, selenocysteine U, pyrrolysine O, internal
# stop codons *, ...) are skipped: feature computation is not meaningful for
# them and some physicochemical scales are undefined.
CANONICAL_AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWY')


def fasta_features(ifile):
    '''Compute features for all sequences in a given FASTA file

    Empty sequences and sequences containing non-standard amino acids are
    skipped with a warning.
    '''
    logger = logging.getLogger('macrel')
    seqs = []
    headers = []
    features = []
    n_empty = 0
    n_noncanonical = 0
    for h, seq in fasta_iter(ifile):
        if not seq:
            n_empty += 1
            continue
        seq = normalize_seq(seq)
        if not seq:
            n_empty += 1
            continue
        if not set(seq) <= CANONICAL_AMINO_ACIDS:
            n_noncanonical += 1
            continue
        seqs.append(seq)
        headers.append(h)
        features.append(compute_all(seq))

    if n_empty:
        logger.warning(
            'Skipped %d empty sequence(s) while computing features.', n_empty)
    if n_noncanonical:
        logger.warning(
            'Skipped %d sequence(s) containing non-standard amino acids '
            '(only the 20 standard amino acids are supported).', n_noncanonical)

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
    feat = fasta_features(ifile)
    feat.to_csv(ofile, sep='\t', index_label='Access')

if __name__ == '__main__':
    main(sys.argv)
