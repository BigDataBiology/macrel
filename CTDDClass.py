#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re, math
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

def Count(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence) * 100)
                    break
        if myCount == 0:
            code.append(0)
    return code
	
# -------------------------------------------------------------------------------------------------------------------------------------
	
# merged from readFasta.py
def readFasta(file):
    if os.path.exists(file) == False:
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open(file) as f:
        records = f.read()

    if re.search('>', records) == None:
        print('The input file seems not in fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    myFasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
    return myFasta

# merged from saveCode.py
def savetsv(encodings, file = 'encoding.tsv'):
    with open(file, 'w') as f:
        if encodings == 0:
            f.write('Descriptor calculation failed.')
        else:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    return None

# -------------------------------------------------------------------------------------------------------------------------------------

def CTDDClass(fastas, groups):
    encodings = []
    header = ['#']

    for g in range(len(groups)):
        for d in ['0', '25', '50', '75', '100']:
            header.append('Group.' + str(g+1) + '.residue' + d)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for group in groups:
            code = code + Count(group, sequence)
        encodings.append(code)

    return encodings


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print(USAGE)
        sys.exit(1)

    groups = sys.argv[3:]
    myStr = ''.join(groups)
    myStr = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', myStr)
    if len(myStr) != 20 or len(set(myStr)) != 20:
        print('\nERROR: The amino acid must be no-repeat in each groups and the sum is 20!\n')
    fastas = readFasta(sys.argv[1])
    encodings = CTDDClass(fastas, groups)
    savetsv(encodings, sys.argv[2])