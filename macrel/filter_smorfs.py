#!/usr/bin/env python

# TODO This function is copy & paste from AMP-features.py, but should
# eventually be merged
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

def filter_smorfs(ifile, ofile):
    with open(ofile, 'wt') as output:
        for h,seq in fasta_iter(ifile):
            if len(seq) > 100: continue
            output.write(">{}\n{}\n".format(h,seq))

