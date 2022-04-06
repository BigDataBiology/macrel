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
    

def ppyrodigal_out(contig, numb, pred):
    return [f'{contig}_{numb}', contig, pred.begin, pred.end,
            pred.strand, pred.confidence(),
            pred.partial_begin, pred.partial_end, 
            pred.gc_cont, pred.translation_table,
            pred.rbs_motif, pred.rbs_spacer,
            pred.start_type, pred.sequence(),
            pred.translate()]


def predict_genes(infile):
    import random
    import pandas as pd
    from .fasta import fasta_iter
    random.seed(1991)
    predictions = []
    gorf, morf_finder = create_pyrodigal_orffinder()
    # predict genes
    for h, s in fasta_iter(infile):
        if len(s) <= 100_000:
            # if contig length less than 100kbp then not suitable for training
            # predict genes using metagenome pretrained models
            for i, pred in enumerate(morf_finder.find_genes(s)):
                predictions.append(ppyrodigal_out(h, i, pred))
        else:
            # if contig length is above or 100kbp then suitable for training of
            # its own model, therefore proceed in a genome wise way
            # each train procedure updates the predictor automatically
            gorf.train(s)
            for i, pred in enumerate(gorf.find_genes(s)):
                predictions.append(ppyrodigal_out(h, i, pred))
    # converting to df
    df = pd.DataFrame(predictions,
                      columns=['artificial_name',
                               'contig', 'start',
                               'end', 'strand',
                               'confidence', 'partial_begin',
                               'partial_end', 'gc_content',
                               'translation_table', 'rbs_motif',
                               'rbs_spacer', 'start_type',
                               'sequence', 'peptide'])
    return df
    

def write_seqs(x):
    # X is a dataframe with two main columns id, peptide
    return f'>{x[0]}\n{x[1]}\n'
    
    
def retrieve_smorfs(df, cluster, smorfs_out, cluster_out):
    df = df[['artificial_name', 'peptide']]
    smorfs = df[df.peptide.str.len() <= 100]
    if cluster:
        smorfs = smorfs.groupby('peptide').agg(lambda x: ','.join(x))
        smorfs = smorfs.reset_index()
        smorfs['id'] = 'smORF_'+smorfs.index.astype('str')
        smorfs = smorfs.rename({'artificial_name': 'genes'}, axis=1)
        smorfs[['id', 'genes']].to_csv(cluster_out, 
                                       sep='\t',
                                       header=True,
                                       index=None)
    else:    
        smorfs.rename({'artificial_name': 'id'},
                      axis=1,
                      inplace=True)
    smorfs = smorfs[['id', 'peptide']]
    with open(smorfs_out, 'w') as ofile:
        lines = smorfs.apply(write_seqs, axis=1)
        ofile.writelines(lines)
        
