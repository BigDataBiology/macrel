import pandas as pd
import numpy as np
import rpy2
import rpy2.robjects
from rpy2.robjects import numpy2ri
numpy2ri.activate()
r = rpy2.robjects.r
r.library('Peptides')


def fasta_iter(fname):
    header = None
    chunks = []
    with open(fname) as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                header = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)

def main(args):
    ifile = args[1]
    ofile = args[2]

    seqs = []
    headers = []
    for h,seq in fasta_iter(ifile):
        seqs.append(seq)
        headers.append(h)

    # We can do this inside the loop so that we are not forced to pre-load all
    # the sequences into memory. However, the overhead is significant and it
    # becomes much slower
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

    features = np.hstack([aaComp, rfeatures])
    features = pd.DataFrame(features, index=headers, columns=["tinyAA",
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
            "hmoment"])
    features.insert(0, 'group', 'Unk')
    features.insert(0, 'sequence', seqs)
    features.to_csv(ofile, sep='\t')

if __name__ == '__main__':
    from sys import argv
    main(argv)
