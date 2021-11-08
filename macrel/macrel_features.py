'''
Functions were adapted from:

  -- ModlAMP(https://github.com/alexarnimueller/modlAMP)
  -- originally available in Macrel v.1 (https://github.com/celiosantosjr/macrel)

Some functions were adapted initially because differences in the scales 
then to avoid calling an entire module for 1 single function, we adapted
some others into here.

The new module for macrel_features, was implemented using some
functions from [modlAMP](https://github.com/alexarnimueller/modlAMP/)
which is under BSD3 license which is completely overlapped by Macrel
licensing under the MIT license.
'''

import math
import numpy as np
from .database import eisenberg, instability2, _aa_groups
from .database import pos_pks, neg_pks, boman_scale
from .database import CTDD_groups
from collections import Counter


def normalize_seq(seq):
    if seq[0] == 'M':
        seq = seq[1:]
    if seq[-1] == '*':
        seq = seq[:-1]
    return seq


def ctdd(seq, groups):
    code = []
    for group in groups:
        for i, aa in enumerate(seq):
            if aa in group:
                code.append((i + 1)/len(seq) * 100)
                break
        else:
            code.append(0)
    return code


def amino_acid_composition(seq, aa_content):
    return [
            sum(aa_content.get(aa, 0.0) for aa in g)/len(seq)
            for g in _aa_groups]

pos_pks10 = {k:10**pk for k,pk in pos_pks.items()}
neg_pks10 = {k:10**pk for k,pk in neg_pks.items()}

def pep_charge_aa(aa_content, ph):
    ph10 = 10**ph
    net_charge = 0.0

    for aa, pK10 in pos_pks10.items():
        c = aa_content.get(aa)
        if c:
            # original was
            # c_r = 10 ** (pK - ph)
            # But this requires one exponentiation per iteration, while this
            # approach is slightly faster
            c_r = pK10 / ph10
            partial_charge = c_r / (c_r + 1.0)
            net_charge += c * partial_charge

    for aa, pK10 in neg_pks10.items():
        c = aa_content.get(aa)
        if c:
            # See above
            c_r = ph10 / pK10
            partial_charge = c_r / (c_r + 1.0)
            net_charge -= c * partial_charge

    return net_charge


def isoelectric_point(aa_content, ph=7.0):
    charge = pep_charge_aa(aa_content, ph)

    if charge > 0.0:
        while charge > 0.0:
            ph += 1.0
            charge = pep_charge_aa(aa_content, ph)
        ph1 = ph - 1.0
        ph2 = ph
    else:
        while charge < 0.0:
            ph -= 1.0
            charge = pep_charge_aa(aa_content, ph)
        ph1 = ph
        ph2 = ph + 1.0

    while ph2 - ph1 > 0.0001 and charge != 0.0:
        ph = (ph1 + ph2) / 2.0
        charge = pep_charge_aa(aa_content, ph)
        if charge > 0.0:
            ph1 = ph
        else:
            ph2 = ph
    return ph



def instability_index(seq):
    stabindex = 0.0
    for i in range(len(seq) - 1):
        stabindex += instability2.get(seq[i:i+2], 0.0)

    return (10.0 / len(seq)) * stabindex


def hydrophobicity(seq, aa_content):
    return sum(
            v * eisenberg.get(k, 0.0) for k,v in aa_content.items()
            ) / len(seq)


def aliphatic_index(seq, aa_content):
    return 100.0*(
            aa_content.get('A', 0.0) +
            2.9 * aa_content.get('V', 0.0) +
            3.9 * (aa_content.get('I', 0.0) + aa_content.get('L', 0.0))
            ) / len(seq) # formula for calculating the AI (Ikai, 1980)


def boman_index(seq, aa_content):
    return sum(
            v * boman_scale.get(k, 0.0) for k,v in aa_content.items()
            ) / len(seq)


def hmoment(seq, angle = 100, window = 11):
    '''
    # http://emboss.bioinformatics.nl/cgi-bin/emboss/hmoment
    # SEQUENCE: FLPVLAGLTPSIVPKLVCLLTKKC
    # ALPHA-HELIX ANGLE=100 : 0.52
    # BETA-SHEET  ANGLE=160 : 0.271
    # 
    # ALPHA HELIX VALUE
    # hmoment(seq = "FLPVLAGLTPSIVPKLVCLLTKKC", angle = 100, window = 11)
    # 0.5199226
    # 
    # BETA SHEET VALUE
    # hmoment(seq = "FLPVLAGLTPSIVPKLVCLLTKKC", angle = 160, window = 11)
    # 0.2705906
    '''
    wdw = min(window, len(seq))  # if sequence is shorter than window, take the whole sequence instead
    mtrx = np.array([eisenberg[aa] for aa in seq], dtype=np.float64)  #[::-1]
    mwdw = np.array([mtrx[i:i + wdw] for i in range(len(mtrx) - wdw + 1)])
    rads = angle * (np.pi / 180) * np.arange(wdw)  # calculate actual moment (radial)

    vcos = np.dot(mwdw, np.cos(rads))
    vsin = np.dot(mwdw, np.sin(rads))
    # The code below is optimized to avoid copies
    vcos **= 2.
    vsin **= 2.
    moms = vsin
    moms += vcos
    return math.sqrt(moms.max()) / wdw

def compute_all(seq):
    aa_content = dict(Counter(seq))
    aa_content['Nterm'] = 1
    aa_content['Cterm'] = 1
    return np.array(
             amino_acid_composition(seq, aa_content) + [
               pep_charge_aa(aa_content, ph=7.0),
               isoelectric_point(aa_content, ph=7.0),
               aliphatic_index(seq, aa_content),
               instability_index(seq),
               boman_index(seq, aa_content),
               hydrophobicity(seq, aa_content),
               hmoment(seq, angle=100, window=11)] +
             ctdd(seq, CTDD_groups))

