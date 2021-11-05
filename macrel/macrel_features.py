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

import numpy as np
from .database import eisenberg, instability, _aa_groups
from .database import pos_pks, neg_pks, boman_scale
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
    return np.array(code)


def amino_acid_composition(seq):
    return np.array(
            [sum(map(g.__contains__, seq)) for g in _aa_groups],
            dtype=float)/len(seq)


def pep_charge(seq, ph=7.0):
    aa_content = dict(Counter(seq))
    aa_content['Nterm'] = 1
    aa_content['Cterm'] = 1
    return pep_charge_aa(aa_content, ph)

def pep_charge_aa(aa_content, ph):

    net_charge = 0.0

    for aa, pK in pos_pks.items():
        if aa_content.get(aa):
            c_r = 10 ** (pK - ph)
            partial_charge = c_r / (c_r + 1.0)
            net_charge += aa_content.get(aa) * partial_charge

    for aa, pK in neg_pks.items():
        if aa_content.get(aa):
            c_r = 10 ** (ph - pK)
            partial_charge = c_r / (c_r + 1.0)
            net_charge -= aa_content.get(aa) * partial_charge

    return net_charge


def isoelectric_point(seq, ph=7.0):
    aa_content = dict(Counter(seq))
    aa_content['Nterm'] = 1
    aa_content['Cterm'] = 1
    ph, ph1, ph2 = float(), float(), float()
    desc = []
    charge = pep_charge_aa(aa_content, ph)

    if charge > 0.0:
        ph1 = ph
        charge1 = charge
        while charge1 > 0.0:
            ph = ph1 + 1.0
            charge = pep_charge_aa(aa_content, ph)
            if charge > 0.0:
                ph1 = ph
                charge1 = charge
            else:
                ph2 = ph
                break
    else:
         ph2 = ph
         charge2 = charge
         while charge2 < 0.0:
             ph = ph2 - 1.0
             charge = pep_charge_aa(aa_content, ph)
             if charge < 0.0:
                ph2 = ph
                charge2 = charge
             else:
                ph1 = ph
                break

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
        stabindex += instability[seq[i]][seq[i+1]]

    return (10.0 / len(seq)) * stabindex


def hydrophobicity(seq):
    hydrophobicity = 0.0
    for aa in seq:
        if eisenberg.get(aa):
            hydrophobicity += eisenberg[aa] / len(seq)

    return hydrophobicity


def aliphatic_index(seq):
    aliph = 'AVIL'
    d = {aa: float(seq.count(aa)) / len(seq) for aa in aliph}  # count aa
    return 100*(d['A'] + 2.9 * d['V'] + 3.9 * (d['I'] + d['L']))  # formula for calculating the AI (Ikai, 1980)


def boman_index(seq):
    val = []
    for a in seq:
        val.append(boman_scale[a])

    return sum(val) / len(val)


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
    mtrx = [eisenberg[aa] for aa in seq]  #[::-1]
    mwdw = [mtrx[i:i + wdw] for i in range(len(mtrx) - wdw + 1)]
    mwdw = np.asarray(mwdw)
  
    rads = angle * (np.pi / 180) * np.asarray(range(wdw))  # calculate actual moment (radial)
    vcos = (mwdw * np.cos(rads)).sum(axis=1)
    vsin = (mwdw * np.sin(rads)).sum(axis=1)
    moms = np.sqrt(vsin ** 2 + vcos ** 2) / wdw
    
    return np.max(moms)

def compute_all(seq):
    return [pep_charge(seq, ph=7.0),
            isoelectric_point(seq, ph=7.0),
            aliphatic_index(seq),
            instability_index(seq),
            boman_index(seq),
            hydrophobicity(seq),
            hmoment(seq, angle=100, window=11)]

