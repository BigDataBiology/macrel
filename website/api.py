import flask
from flask import request
import subprocess
import tempfile
import pandas as pd
from os import path
import gzip

app = flask.Flask('Macrel')

@app.route('/predict', methods=['POST'])
def predict():
    args = ['macrel']
    with tempfile.TemporaryDirectory() as tdir:
        if request.form.get('dataType') == 'peptides':
            args.append('peptides')
            ifname = path.join(tdir, 'pep.faa')
        else:
            args.append('contigs')
            ifname = path.join(tdir, 'contigs.fna')

        if request.form.get('textData'):
            with open(ifname, 'wt') as out:
                out.write(request.form.get('textData'))
        else:
            # The input filename is never used, also can be security
            # vulnerability, so use "input.fa":
            ifname = path.join(tdir, 'input.fa')
            request.files['file'].save(ifname)

        args.extend([
            '--fasta', ifname,
            '--output', tdir + '/macrel.out',
            ])
        subprocess.check_call(args)
        tfile = path.join(tdir, 'macrel.out', 'macrel.out.prediction.gz')
        data = pd.read_table(tfile)
        data.rename(inplace=True, columns={
            'Unnamed: 0': 'access',
            'group': 'amp_family',})

        r = flask.jsonify({
            'code': 1,
            'message': 'Macrel OK',
            'rawdata': gzip.open(tfile, 'rt').read(),
            'data': {
                'objects': [
                    data.iloc[i].to_dict()
                        for i in range(len(data))]
                },
            })
        r.headers.add('Access-Control-Allow-Origin', '*')
        r.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
        r.headers.add('Access-Control-Allow-Methods', 'GET,POST')
        return r

if __name__ == '__main__':
    app.run()
