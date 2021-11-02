import logging
import numpy as np
from .database import eisenberg, instability, _aa_groups
from .database import pos_pks, neg_pks, boman_scale
from collections import Counter

class MacrelFeatures:

    '''
    Functions were adapted from:

      -- ModlAMP(https://github.com/alexarnimueller/modlAMP)
      -- Peptides R Package (https://github.com/dosorio/Peptides), and
      -- originally available in Macrel v.1 (https://github.com/celiosantosjr/macrel)

    Some functions were adapted initially because differences in the scales 
    then to avoid calling an entire module for 1 single function, we adapted
    some others into here.
    '''
    
    def __init__(self, seq):
        self.seq = seq

    def checkseq(self):
        if self.seq[-1] == '*' and self.seq[0] == 'M':
            self.seq = self.seq[1:-1]
        elif self.seq[0] == 'M':
                self.seq = self.seq[1:]
        elif self.seq[-1] == '*':
            self.seq = self.seq[:-1]
        return self.seq
  

    def ctdd(self, groups):
        code = []
        for group in groups:
            for i, aa in enumerate(self.seq):
                if aa in group:
                    code.append((i + 1)/len(self.seq) * 100)
                    break
            else:
                code.append(0)
        return np.array(code)


    def amino_acid_composition(self):
        # See groups above
        return np.array(
                [sum(map(g.__contains__, self.seq)) for g in _aa_groups],
                dtype=float)/len(self.seq)

    
    def pep_charge(self, ph=7.0):
        aa_content = dict(Counter(self.seq))
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

        return net_charge
    
    
    def isoelectric_point(self, ph=7.0):
        ph, ph1, ph2 = float(), float(), float()
        desc = []
        charge = self.pep_charge(ph)

        if charge > 0.0:
            ph1 = ph
            charge1 = charge
            while charge1 > 0.0:
                ph = ph1 + 1.0
                charge = self.pep_charge(ph)
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
                 charge = self.pep_charge(ph)
                 if charge < 0.0:
                    ph2 = ph
                    charge2 = charge
                 else:
                    ph1 = ph
                    break
    
        while ph2 - ph1 > 0.0001 and charge != 0.0:
            ph = (ph1 + ph2) / 2.0
            charge = self.pep_charge(ph)
            if charge > 0.0:
                ph1 = ph
            else:
                ph2 = ph
    
        return ph


    def instability_index(self):
        stabindex = 0.0
        for i in range(len(self.seq) - 1):
            stabindex += instability[self.seq[i]][self.seq[i+1]]
    
        return (10.0 / len(self.seq)) * stabindex


    def hydrophobicity(self):
        hydrophobicity = 0.0
        for aa in self.seq:
            if eisenberg.get(aa):
                hydrophobicity += eisenberg[aa] / len(self.seq)
    
        return hydrophobicity
    

    def aliphatic_index(self):
        aliph = 'AVIL'
        d = {aa: float(self.seq.count(aa)) / len(self.seq) for aa in aliph}  # count aa
        
        return 100*(d['A'] + 2.9 * d['V'] + 3.9 * (d['I'] + d['L']))  # formula for calculating the AI (Ikai, 1980)


    def boman_index(self):
        val = []
        for a in self.seq:
            val.append(boman_scale[a])

        return sum(val) / len(val)


    def hmoment(self, angle = 100, window = 11):
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
        wdw = min(window, len(self.seq))  # if sequence is shorter than window, take the whole sequence instead
        mtrx = [eisenberg[aa] for aa in self.seq]  #[::-1]
        mwdw = [mtrx[i:i + wdw] for i in range(len(mtrx) - wdw + 1)]
        mwdw = np.asarray(mwdw)
      
        rads = angle * (np.pi / 180) * np.asarray(range(wdw))  # calculate actual moment (radial)
        vcos = (mwdw * np.cos(rads)).sum(axis=1)
        vsin = (mwdw * np.sin(rads)).sum(axis=1)
        moms = np.sqrt(vsin ** 2 + vcos ** 2) / wdw
        
        return np.max(moms)
    
    def compute_all(self):
        seq = self.checkseq()
        return [self.pep_charge(),
                self.isoelectric_point(),
                self.aliphatic_index(),
                self.instability_index(),
                self.boman_index(),
                self.hydrophobicity(),
                self.hmoment()]
                
