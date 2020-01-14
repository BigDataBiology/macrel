from .fasta import fasta_iter

def filter_smorfs(ifile, ofile):
    with open(ofile, 'wt') as output:
        for h,seq in fasta_iter(ifile):
            if len(seq) > 100: continue
            output.write(">{}\n{}\n".format(h,seq))

