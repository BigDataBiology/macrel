from macrel import ampsphere

seq1 = 'KRVKSFFKGYMRAIEINAALMYGYRPK'
seq2 = 'GRVIGKQGRIAKAIRVVMRAAAVRVDEKVLVEID'
def test_exact():
    r = ampsphere.get_ampsphere_exact_match(seq1, 'seq1')
    assert r.index[0] == 'seq1'
    assert r.iloc[0]['result'] == 'AMP10.000_002'
    r = ampsphere.get_ampsphere_exact_match(seq1 + 'HELO', 'seq2')
    assert r.index[0] == 'seq2'
    assert r.iloc[0]['result'] is None

def test_mmseqs():
    r = ampsphere.get_ampsphere_mmseqs_match(seq1, 'seq1')
    assert r.index[0] == 'seq1'
    assert r.iloc[0]['target_identifier'] == 'AMP10.000_002'

def test_hmmer():
    r = ampsphere.get_ampsphere_hmmer_match(seq2, 'seq2')
    assert r.index[0] == 'seq2'
