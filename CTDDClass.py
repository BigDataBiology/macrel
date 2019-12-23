#!/usr/bin/env python
# _*_coding:utf-8_*_

import sys
import os
import re
import math
import collections

# the current python file location
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)

USAGE = """
USAGE:
    python CTDDClass.py input.fasta output amino_acids_group_1 amino_acids_group_2 ... amino_acids_group_N

    input.fasta:                 the input protein sequence file in fasta format.    
    output:                      the encoding file.
    amino_acids_group_x          the amino acids groups.

EXAMPLE:
    python CTDDClass.py example/test-protein.txt CTDDClass.tsv RKEDQN GASTPHY CLVIMFW
"""

# -------------------------------------------------------------------------------------------------------------------------------------


# merged from readFasta.py
def readFasta(file):
    if not os.path.exists(file):
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open(file) as f:
        records = f.read()

    if re.search('>', records) is None:
        print('The input file seems not in fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    fastas = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
        fastas.append([name, sequence])
    return fastas


# merged from saveCode.py
def savetsv(encodings, file='encoding.tsv'):
    lengthOfEncodings = len(encodings)
    with open(file, 'w') as f:
        if lengthOfEncodings == 0:
            f.write('Descriptor calculation failed.')
        else:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for encoding in encodings[1:]:
                for j in range(0, len(encoding) - 1):
                    f.write(str(encoding[j]) + '\t')
                f.write(str(encoding[len(encoding)-1]) + '\n')
    return None

# -------------------------------------------------------------------------------------------------------------------------------------


def count(group, sequence):
    number = 0
    lengthOfSequence = len(sequence)

    letterFrequency = collections.Counter(sequence)
    for letter, frequency in letterFrequency.items():
        if letter in group:
            number = number + frequency

    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [num if num >= 1 else 1 for num in cutoffNums]

    code = []
    for num in cutoffNums:
        count = 0
        for i in range(lengthOfSequence):
            if sequence[i] in group:
                count += 1
                if count == num:
                    code.append((i + 1) / lengthOfSequence * 100)
                    break
        if count == 0:
            code.append(0)
    return code


def ctdd(fastas, groups):
    encodings = []
    header = ['#']

    for g in range(len(groups)):
        for d in ['0', '25', '50', '75', '100']:
            header.append('Group.' + str(g+1) + '.residue' + d)
    encodings.append(header)

    for fasta in fastas:
        name, sequence = fasta[0], fasta[1]
        code = [name]
        for group in groups:
            code = code + count(group, sequence)
        encodings.append(code)

    return encodings


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print(USAGE)
        sys.exit(1)

    groups = sys.argv[3:]

    groupsStr = ''.join(groups)
    groupsStr = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', groupsStr)
    if len(groupsStr) != 20 or len(set(groupsStr)) != 20:
        print('\nERROR: The amino acid must be no-repeat in each groups and the sum is 20!\n')

    fastas = readFasta(sys.argv[1])
    if (0 == len(fastas)):
        print("The content of input file is empty!\n")
    encodings = ctdd(fastas, groups)
    savetsv(encodings, sys.argv[2])
