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
def main(args):
    if len(args) < 2:
        import sys
        sys.stderr.write("This is an internal FACS script and is not meant to be used independently")
        sys.exit(1)

    ifile = args[1]

    for h,seq in fasta_iter(ifile):
        if len(seq) > 100: continue
        print(">{}\n{}".format(h,seq))

if __name__ == '__main__':
    from sys import argv
    main(argv)
