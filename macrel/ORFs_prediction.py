def create_pyrodigal_orffinder():
    # generating orf_finder
    import pyrodigal
    gorf = pyrodigal.OrfFinder(closed=True,
                               min_gene=33,
                               max_overlap=0)       
    morf_finder = pyrodigal.OrfFinder(meta=True,
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
    gc = f'gc_cont={round(pred.gc_cont, 3)}'
    last = ';'.join([seconid, part, st, motif, sp, gc])
    header = f'>{orfid} # {pred.begin} # {pred.end} # {pred.strand} # {last}'
    return (header, pred.translate())


def predict_genes(infile, ofile):
    import random
    import pandas as pd
    from .fasta import fasta_iter
    random.seed(1991)
    predictions = []
    gorf, morf_finder = create_pyrodigal_orffinder()
    # predict genes
    for idx, (h, s) in enumerate(fasta_iter(infile)):
        if len(s) <= 100_000:
            # if contig length less than 100kbp then not suitable for training
            # predict genes using metagenome pretrained models
            for i, pred in enumerate(morf_finder.find_genes(s)):
                predictions.append(ppyrodigal_out(h, i+1, idx+1, pred))
        else:
            # if contig length is above or 100kbp then suitable for training of
            # its own model, therefore proceed in a genome wise way
            # each train procedure updates the predictor automatically
            gorf.train(s)
            for i, pred in enumerate(gorf.find_genes(s)):
                predictions.append(ppyrodigal_out(h, i+1, idx+1, pred))
    # saving sequences
    with open(ofile, 'w') as odb:
        for h, s in predictions:
            odb.write(f'{h}\n{s}\n')
    
