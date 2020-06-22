from .fasta import fasta_iter
from .utils import open_output

def filter_smorfs(ifile, ofile, uniq, full_headers=False):
    '''Remove larger ORFs, leaving only smORFs behind'''
    seen = set()
    with open_output(ofile, mode='wt') as output:
        for h,seq in fasta_iter(ifile, full_headers):
            if len(seq) > 100: continue
            if uniq:
                if seq in seen: continue
                h = 'smORF_{}'.format(len(seen))
                seen.add(seq)
            output.write(">{}\n{}\n".format(h,seq))

