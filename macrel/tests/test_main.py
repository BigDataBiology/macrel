import sys
from os import path

from macrel import main


def test_get_smorfs_file_output(tmp_path, monkeypatch):
    '''`get-smorfs --file-output <file>` used to crash writing the README

    `args.output` is None for a plain --file-output, so the README write
    raised TypeError before the pipeline could run (see improvements #1).
    '''
    fasta = path.join(path.dirname(__file__), '..', '..',
                      'example_seqs', 'excontigs.fna.gz')
    out_file = tmp_path / 'smorfs.faa'
    monkeypatch.setattr(sys, 'argv',
            ['macrel', 'get-smorfs',
             '-f', fasta,
             '--file-output', str(out_file)])
    main.main(sys.argv)
    assert out_file.exists()
    assert out_file.read_text().startswith('>')
