#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import re
import math

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

def ctdd(sequence, groups):
    code = []
    for group in groups:
        for i, aa in enumerate(sequence):
            if aa in group:
                code.append(str( (i + 1)/len(sequence) * 100 ))
                break
        else:
            code.append("0")
    return code


def main(args):
    if len(args) < 6:
        sys.stderr.write("This is an internal FACS script and is not meant to be used independently")
        sys.exit(1)

    ifile = args[1]
    ofile = args[2]
    name = args[3]
    groups = args[4:]

    groupsStr = ''.join(groups)
    groupsStr = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', groupsStr)
    if len(groupsStr) != 20 or len(set(groupsStr)) != 20:
        print('\nERROR: The amino acid must be no-repeat in each groups and the sum is 20!\n')

    groups = [set(g) for g in groups]
    with open(ofile, 'wt') as out:
        for g in range(len(groups)):
            if g > 0:
                out.write('\t')
            out.write(name + str(g+1) + '.residue0')
        out.write('\n')
        for h,f in fasta_iter(ifile):
            encodings = ctdd(f, groups)
            out.write('\t'.join(encodings))
            out.write('\n')

if __name__ == '__main__':
    main(sys.argv)
