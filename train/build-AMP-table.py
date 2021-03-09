from macrel.fasta import fasta_iter
from macrel.AMP_features import features
from os import makedirs

makedirs('preproc/', exist_ok=True)
normalized_fname = 'preproc/AMP_NAMP.train.faa'
normalized_fname2 = 'preproc/AMP_NAMP_lt50.train.faa'

# The AmPEP data has duplicates! The same exact same sequences appear on both
# the positive and negative classes:
seen = set()
with open(normalized_fname, 'wt') as output:
    with open(normalized_fname2, 'wt') as output2:
        for i, (_, seq) in enumerate(fasta_iter('data/M_model_train_AMP_sequence.fasta')):
            if len(seq) <= 50:
                output2.write(f">AMP_{i}\n{seq}\n")
            output.write(f">AMP_{i}\n{seq}\n")
            seen.add(seq)          
        for i, (_, seq) in enumerate(fasta_iter('data/M_model_train_nonAMP_sequence.fasta')):
            if seq in seen: continue
            if len(seq) <= 50:
                output2.write(f">NAMP_{i}\n{seq}\n")
            output.write(f">NAMP_{i}\n{seq}\n")
            seen.add(seq)

fs = features(normalized_fname)
fs['group'] = fs.index.map(lambda ix: ix.split('_')[0])
fs.to_csv('preproc/AMP.train.tsv', sep='\t')

fs = features(normalized_fname2)
fs['group'] = fs.index.map(lambda ix: ix.split('_')[0])
fs.to_csv('preproc/AMP_lt50.train.tsv', sep='\t')

normalized_fname_test = 'preproc/AMP_NAMP.test.faa'
normalized_fname_test2 = 'preproc/AMP_NAMP_lt50.test.faa'
with open(normalized_fname_test, 'wt') as output:
    with open(normalized_fname_test2, 'wt') as output2:
        for i, (_, seq) in enumerate(fasta_iter('data/Supp-S2_AMP.faa')):
            output.write(f">AMP_{i}\n{seq}\n")
            if len(seq) <= 50:
                output2.write(f">AMP_{i}\n{seq}\n")
        for i, (_, seq) in enumerate(fasta_iter('data/Supp-S2_NAMP.faa')):
            output.write(f">NAMP_{i}\n{seq}\n")
            if len(seq) <= 50:
                output2.write(f">NAMP_{i}\n{seq}\n")

fs_t = features(normalized_fname_test)
fs_t['group'] = fs_t.index.map(lambda ix: ix.split('_')[0])
fs_t.to_csv('preproc/AMP.test.tsv', sep='\t')

fs_t = features(normalized_fname_test2)
fs_t['group'] = fs_t.index.map(lambda ix: ix.split('_')[0])
fs_t.to_csv('preproc/AMP_lt50.test.tsv', sep='\t')

normalized_fname_test = 'preproc/AMP_NAMP.train.bench.faa'
normalized_fname_test = 'preproc/AMP_NAMP_lt50.train.bench.faa'
with open(normalized_fname_test, 'wt') as output:
    with open(normalized_fname_test2, 'wt') as output2:
        for i, (_, seq) in enumerate(fasta_iter('data/Supp-S1_AMP.faa')):
            output.write(f">AMP_{i}\n{seq}\n")
            if len(seq) <= 50:
                output2.write(f">AMP_{i}\n{seq}\n")
        for i, (_, seq) in enumerate(fasta_iter('data/Supp-S1_NAMP.faa')):
            output.write(f">NAMP_{i}\n{seq}\n")
            if len(seq) <= 50:
                output2.write(f">NAMP_{i}\n{seq}\n")

fs_bench = features(normalized_fname_test)
fs_bench['group'] = fs_bench.index.map(lambda ix: ix.split('_')[0])
fs_bench.to_csv('preproc/AMP.train_bench.tsv', sep='\t')

fs_bench = features(normalized_fname_test2)
fs_bench['group'] = fs_bench.index.map(lambda ix: ix.split('_')[0])
fs_bench.to_csv('preproc/AMP_lt50.train_bench.tsv', sep='\t')
