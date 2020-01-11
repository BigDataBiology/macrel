from macrel.fasta import fasta_iter
from macrel.AMP_features import features


normalized_fname = 'data/Hemo.training.normalized.faa'
with open(normalized_fname, 'wt') as output:
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.training.pos.faa')):
        output.write(f">Hemo_{i}\n{seq}\n")
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.training.neg.faa')):
        output.write(f">nHemo_{i}\n{seq}\n")
fs = features(normalized_fname)
fs['group'] = fs.index.map(lambda ix: ix.split('_')[0])
fs.to_csv('data/Hemo.train.tsv', sep='\t')

normalized_fname = 'data/Hemo.validation.normalized.faa'
with open(normalized_fname, 'wt') as output:
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.validation.pos.faa')):
        output.write(f">Hemo_{i}\n{seq}\n")
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.validation.neg.faa')):
        output.write(f">nHemo_{i}\n{seq}\n")
fs_v = features(normalized_fname)
fs_v['group'] = fs_v.index.map(lambda ix: ix.split('_')[0])
fs_v.to_csv('data/Hemo.validation.tsv', sep='\t')
