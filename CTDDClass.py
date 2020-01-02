#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import re
import math

USAGE = """
USAGE:
    python CTDDClass.py input.fasta output amino_acids_group_1 amino_acids_group_2 ... amino_acids_group_N

    input.fasta:                 the input protein sequence file in fasta format.    
    output:                      the encoding file.
    amino_acids_group_x          the amino acids groups.

EXAMPLE:
    python CTDDClass.py input.faa CTDDClass.tsv RKEDQN GASTPHY CLVIMFW
"""



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



def count(group, sequence):
    total = sum((a in group) for a in sequence)
    if total == 0:
        return ["0" for _ in range(5)]
    
    cutoffs = [1, math.floor(0.25 * total), math.floor(0.50 * total), math.floor(0.75 * total), total]
    cutoffs = [num if num >= 1 else 1 for num in cutoffs]

    code = []
    next_cut = 0
    count = 0
    for i, aa in enumerate(sequence):
        if aa in group:
            count += 1
            while count == cutoffs[next_cut]:
                code.append( str((i + 1)/len(sequence) * 100) )
                next_cut += 1
                if next_cut == len(cutoffs):
                    break
    return code


def ctdd(sequence, groups):
    code = []
    for group in groups:
        code.extend(count(group, sequence))
    return code


def main(args):
    if len(args) < 5:
        print(USAGE)
        sys.exit(1)

    ifile = args[1]
    ofile = args[2]
    groups = args[3:]

    groupsStr = ''.join(groups)
    groupsStr = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', groupsStr)
    if len(groupsStr) != 20 or len(set(groupsStr)) != 20:
        print('\nERROR: The amino acid must be no-repeat in each groups and the sum is 20!\n')

    groups = [set(g) for g in groups]
    with open(ofile, 'wt') as out:
        out.write('#')
        for g in range(len(groups)):
            for d in ['0', '25', '50', '75', '100']:
                out.write('\tGroup.' + str(g+1) + '.residue' + d)
        out.write('\n')
        for h,f in fasta_iter(ifile):
            encodings = ctdd(f, groups)
            out.write('\t'.join([h] + encodings))
            out.write('\n')

if __name__ == '__main__':
    main(sys.argv)
