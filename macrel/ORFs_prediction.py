def create_pyrodigal_orffinders():
    import pyrodigal
    # Prior to version 3.0, pyrodigal called the GeneFinder class
    # OrfFinder. This little bit of code allows us to use either
    if hasattr(pyrodigal, 'OrfFinder'):
        GeneFinder = pyrodigal.OrfFinder
    else:
        GeneFinder = pyrodigal.GeneFinder
    gorf = GeneFinder(closed=True,
                      min_gene=33,
                      max_overlap=0)
    morf_finder = GeneFinder(meta=True,
                      closed=True,
                      min_gene=33,
                      max_overlap=0)
    return gorf, morf_finder


def ppyrodigal_out(contig, ind, idx, pred):
    orfid = f'{contig}_{ind}'
    seconid = f'ID={idx}_{ind}'
    part = ''.join([str(int(pred.partial_begin)),
                    str(int(pred.partial_end))])
    part = f'partial={part}'
    st = f'start_type={pred.start_type}'
    motif = f'rbs_motif={pred.rbs_motif}'
    sp = f'rbs_spacer={pred.rbs_spacer}'
    gc = f'gc_cont={pred.gc_cont:.3f}'
    last = ';'.join([seconid, part, st, motif, sp, gc])
    header = f'>{orfid} # {pred.begin} # {pred.end} # {pred.strand} # {last}'
    return f'{header}\n{pred.translate()}\n'


def predict_genes(infile, ofile):
    import random
    import pandas as pd
    from .fasta import fasta_iter
    from .utils import open_output

    clen = []
    gorf, morf_finder = create_pyrodigal_orffinders()

    # predict genes
    with open_output(ofile) as odb:
        for idx, (h, s) in enumerate(fasta_iter(infile)):
            orfs, smorfs = [0, 0]
            if len(s) <= 100_000:
                # if contig length less than 100kbp then not suitable for training
                # predict genes using metagenome pretrained models
                for i, pred in enumerate(morf_finder.find_genes(s)):
                    t = ppyrodigal_out(h, i+1, idx+1, pred)
                    odb.write(t)
                    orfs += 1
                    if len(pred.translate()) <= 100: smorfs += 1
            else:
                # if contig length is above or 100kbp then suitable for training of
                # its own model, therefore proceed in a genome wise way
                # each train procedure updates the predictor automatically
                gorf.train(s)
                for i, pred in enumerate(gorf.find_genes(s)):
                    t = ppyrodigal_out(h, i+1, idx+1, pred)
                    odb.write(t)
                    orfs += 1
                    if len(pred.translate()) <= 100: smorfs += 1

            clen.append([h, len(s), orfs, smorfs])

    return pd.DataFrame(clen, columns=['contig', 'length', 'ORFs', 'smORFs'])

