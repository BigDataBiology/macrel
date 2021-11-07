from macrel.fasta import fasta_iter
from macrel.AMP_features import fasta_features
from os import makedirs

makedirs('preproc/', exist_ok=True)

normalized_fname = 'preproc/Hemo.train.normalized.faa'
with open(normalized_fname, 'wt') as output:
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.training.pos.faa')):
        output.write(f">Hemo_{i}\n{seq}\n")
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.training.neg.faa')):
        output.write(f">NonHemo_{i}\n{seq}\n")
fs = fasta_features(normalized_fname)
fs['group'] = fs.index.map(lambda ix: ix.split('_')[0])
fs.to_csv('preproc/Hemo.train.tsv', sep='\t')

normalized_fname = 'preproc/Hemo.test.normalized.faa'
with open(normalized_fname, 'wt') as output:
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.validation.pos.faa')):
        output.write(f">Hemo_{i}\n{seq}\n")
    for i, (_, seq) in enumerate(fasta_iter('data/hemo.validation.neg.faa')):
        output.write(f">NonHemo_{i}\n{seq}\n")
fs_v = fasta_features(normalized_fname)
fs_v['group'] = fs_v.index.map(lambda ix: ix.split('_')[0])
fs_v.to_csv('preproc/Hemo.test.tsv', sep='\t')
