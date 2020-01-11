from macrel.fasta import fasta_iter
from macrel.AMP_features import features
normalized_fname = 'data/AMP_NAMP.training.faa'
with open(normalized_fname, 'wt') as output:
    for i, (_, seq) in enumerate(fasta_iter('data/M_model_train_AMP_sequence.fasta')):
        output.write(f">AMP_{i}\n{seq}\n")
    for i, (_, seq) in enumerate(fasta_iter('data/M_model_train_nonAMP_sequence.fasta')):
        output.write(f">NAMP_{i}\n{seq}\n")
fs = features(normalized_fname)
fs['group'] = fs.index.map(lambda ix: ix.split('_')[0])
fs.to_csv('data/AMP.train.tsv', sep='\t')

normalized_fname_validation = 'data/AMP_NAMP.validation.faa'
with open(normalized_fname_validation, 'wt') as output:
    for i, (_, seq) in enumerate(fasta_iter('data/Supp-S2_AMP.faa')):
        output.write(f">AMP_{i}\n{seq}\n")
    for i, (_, seq) in enumerate(fasta_iter('data/Supp-S2_NAMP.faa')):
        output.write(f">NAMP_{i}\n{seq}\n")
        
fs_v = features(normalized_fname_validation)
fs_v['group'] = fs_v.index.map(lambda ix: ix.split('_')[0])
fs_v.to_csv('data/AMP.validation.tsv', sep='\t')
