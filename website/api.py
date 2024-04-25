import flask
from flask import request
import subprocess
import tempfile
import pandas as pd
from os import path
from macrel.macrel_version import __version__ as macrel_version
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
            '--keep-negatives',
            ])
        try:
            subprocess.check_call(args)
            tfile = path.join(tdir, 'macrel.out', 'macrel.out.prediction.gz')
            data = pd.read_table(tfile)
            data.rename(inplace=True, columns={
                'Unnamed: 0': 'access',
                'group': 'amp_family',})

            # is_AMP is a boolean, but is of type bool_ which is not JSON
            # serializable. Since it's redundant with AMP_probability, we can
            # drop it
            data.drop(['is_AMP'], axis=1, inplace=True)

            r = flask.jsonify({
                'code': 1,
                'message': 'Macrel OK',
                'rawdata': gzip.open(tfile, 'rt').read(),
                'macrel_version': macrel_version,
                'data': {
                    'objects': [
                        data.iloc[i].to_dict()
                            for i in range(len(data))]
                    },
                })
        except Exception as e:
            print(e)
            r = flask.jsonify({
                'code': 0,
                'message': 'Error running Macrel',
                })
        return r

if __name__ == '__main__':
    app.run()
