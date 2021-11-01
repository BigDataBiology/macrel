'''
Functions were adapted from:

  -- ModlAMP(https://github.com/alexarnimueller/modlAMP)
  -- Peptides R Package (https://github.com/dosorio/Peptides), and
  -- originally available in Macrel v.1 (https://github.com/celiosantosjr/macrel)

Some functions were adapted initially because differences in the scales 
then to avoid calling an entire module for 1 single function, we adapted
some others into here.
'''

def amino_acid_composition(seq):
    from .database import _aa_groups
    import numpy as np
    
    # See groups above
    return np.array(
            [sum(map(g.__contains__, seq)) for g in _aa_groups],
            dtype=float)/len(seq)


def ctdd(sequence, groups):
    import numpy as np
    
    code = []
    for group in groups:
        for i, aa in enumerate(sequence):
            if aa in group:
                code.append((i + 1)/len(sequence) * 100)
                break
        else:
            code.append(0)
    return np.array(code)
    
    
def pep_charge(seq, ph=7.0):
    from collections import Counter
    from .database import pos_pks, neg_pks

    aa_content = dict(Counter(seq))
    aa_content['Nterm'] = 1
    aa_content['Cterm'] = 1

    net_charge = 0.0

    for aa, pK in pos_pks.items():
        c_r = 10 ** (pK - ph)
        partial_charge = c_r / (c_r + 1.0)
        if aa_content.get(aa):
            net_charge += aa_content.get(aa) * partial_charge

    for aa, pK in neg_pks.items():
        c_r = 10 ** (ph - pK)
        partial_charge = c_r / (c_r + 1.0)
        if aa_content.get(aa):
            net_charge -= aa_content.get(aa) * partial_charge

    return round(net_charge, 3)
    
def isoelectric_point(seq, ph=7.0):
    ph1, ph2 = float(), float()
    desc = []
    charge = pep_charge(seq, ph)

    if charge > 0.0:
        ph1 = ph
        charge1 = charge
        while charge1 > 0.0:
            ph = ph1 + 1.0
            charge = pep_charge(seq, ph)
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
             charge = pep_charge(seq, ph)
             if charge < 0.0:
                ph2 = ph
                charge2 = charge
             else:
                ph1 = ph
                break

    while ph2 - ph1 > 0.0001 and charge != 0.0:
        ph = (ph1 + ph2) / 2.0
        charge = pep_charge(seq, ph)
        if charge > 0.0:
            ph1 = ph
        else:
            ph2 = ph

    return ph


def instability_index(seq):
    from .database import instability as dimv

    stabindex = 0.0
    for i in range(len(seq) - 1):
        stabindex += dimv[seq[i]][seq[i+1]]

    return (10.0 / len(seq)) * stabindex


def hydrophobicity(seq):
    from .database import eisenberg

    hydrophobicity = 0.0

    for aa in seq:
        if eisenberg.get(aa):
            hydrophobicity += eisenberg[aa] / len(seq)

    return hydrophobicity


def aliphatic_index(seq):
    aliph = 'AVIL'
    d = {aa: (seq.count(aa) / len(seq)) for aa in aliph}  # count aa
    
    return d['A'] + 2.9 * d['V'] + 3.9 * (d['I'] + d['L'])  # formula for calculating the AI (Ikai, 1980)


def boman_index(seq):
    from .database import boman_scale

    val = []
    for a in seq:
        val.append(boman_scale[a])

    return sum(val) / len(val)


def hmoment(seq, angle = 100, window = 11):
  from .database import eisenberg
  import numpy as np

  wdw = min(window, len(seq))  # if sequence is shorter than window, take the whole sequence instead
  mtrx = [eisenberg[aa] for aa in seq]
  mwdw = [[sum(mtrx[i:i + wdw])] for i in range(len(mtrx) - wdw + 1)]
  mwdw = np.asarray(mwdw)
  
  rads = angle * (np.pi / 180) * np.asarray(range(wdw))  # calculate actual moment (radial)
  vcos = (mwdw * np.cos(rads)).sum(axis=1)
  vsin = (mwdw * np.sin(rads)).sum(axis=1)
  moms = np.sqrt(vsin ** 2 + vcos ** 2) / wdw
    
  return np.max(moms)
