import shutil
import tempfile
import requests
import pandas as pd
from os import path
import logging

from macrel.utils import open_output
from macrel.macrel_version import __version__

REQUESTS_HEADER = {
        'User-Agent': f'macrel/{__version__} (python-requests)'
        }

def get_cache_directory(args):
    '''Get cache directory'''
    from os import environ
    if args.cache_dir:
        cache_dir = args.cache_dir
    else:
        if 'XDG_CACHE_HOME' in environ:
            cache_dir = path.join(environ['XDG_CACHE_HOME'], 'macrel')
        else:
            cache_dir = path.join(path.expanduser("~"), '.cache', 'macrel')
    if not path.exists(cache_dir):
        from os import makedirs
        makedirs(cache_dir, exist_ok=True)
        with open_output(path.join(cache_dir, 'CACHEDIR.TAG'), 'wt') as f:
            f.write('Signature: 8a477f597d28d172789f06886806bc55\n')
            f.write('# This file is a cache directory tag created by macrel.\n')
            f.write('# For information about cache directory tags, see:\n')
            f.write('#\thttp://www.brynosaurus.com/cachedir/\n')
    return cache_dir

def maybe_download_ampsphere_mmseqs(args):
    import tarfile
    target = path.join(get_cache_directory(args), 'AMPSphere_latest.mmseqsdb')
    if path.exists(target):
        logging.debug(f'AMPSphere MMSeqs2 database already downloaded to {target}')
        if args.re_download_database:
            logging.debug(f'Forced redownload enabled, re-downloading AMPSphere MMSeqs2 database')
            shutil.rmtree(target)
        else:
            return target
    AMPSPHERE_MMSEQS2_URL = 'https://ampsphere-api.big-data-biology.org/v1/downloads/AMPSphere_latest.mmseqsdb.tar.xz'
    with tempfile.TemporaryDirectory() as tmpdir:
        tfile = path.join(tmpdir, 'AMPSphere_latest.mmseqsdb.tar.xz')
        r = requests.get(AMPSPHERE_MMSEQS2_URL, stream=True, headers=REQUESTS_HEADER)
        with open(tfile, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        logging.debug(f'Downloaded AMPSphere MMSeqs2 database to {tfile}')
        with tarfile.open(tfile) as tar:
            tar.extractall(tmpdir)
        logging.debug(f'Extracted AMPSphere MMSeqs2 database to {tmpdir}')
        logging.debug(f'Moving AMPSphere MMSeqs2 database to {target}')
        return shutil.move(path.join(tmpdir, 'mmseqs_db'), target)


def maybe_download_ampsphere_hmm(args):
    target_dir = path.join(get_cache_directory(args), 'hmm_db')
    target = path.join(target_dir, 'AMPSphere_latest.hmm')
    if path.exists(target):
        logging.debug(f'AMPSphere HMM database already downloaded to {target}')
        if args.re_download_database:
            logging.debug(f'Force redownload enabled, re-downloading AMPSphere HMM database')
            shutil.rmtree(target)
        else:
            return target
    HMM_URL = 'https://ampsphere-api.big-data-biology.org/v1/downloads/AMPSphere_latest.hmm'
    if not shutil.which('hmmpress'):
        logging.error('HMMER not found. Please install it first (you can use `conda install -c bioconda hmmer`)')
        sys.exit(1)
    with tempfile.TemporaryDirectory() as tmpdir:
        tfile = path.join(tmpdir, 'AMPSphere_latest.hmm')
        r = requests.get(HMM_URL, stream=True, headers=REQUESTS_HEADER)
        with open(tfile, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        logging.debug(f'Downloaded AMPSphere HMM database to {tfile}')
        _logged_subprocess_call(['hmmpress', tfile])
        logging.debug(f'Indexed AMPSphere HMM database')
        logging.debug(f'Moving AMPSphere HMM database to {target}')
        # move all files in tmpdir to target
        from os import listdir, makedirs
        makedirs(target_dir, exist_ok=True)
        for f in listdir(tmpdir):
            shutil.move(path.join(tmpdir, f), target_dir)
        return target

def get_ampsphere_hmmer_match_local(args, seqs):
    from macrel.fasta import fasta_iter
    hmm = maybe_download_ampsphere_hmm(args)
    if not hmm:
        logging.error('AMPSphere HMM database not found. Please download it first or provide the path to the database using the --cache-dir flag')
        sys.exit(1)
    with tempfile.TemporaryDirectory() as tmpdir:
        query_file = path.join(tmpdir, 'query.faa')
        output_file = path.join(tmpdir, 'output.tsv')
        with open(query_file, 'wt') as f:
            for (query_name, seq) in seqs:
                f.write(f'>{query_name}\n{seq}\n')
        _logged_subprocess_call(
            ['hmmsearch', '--tblout', output_file, hmm, query_file])
        return pd.read_csv(output_file, sep='\s+', comment='#', header=None,
                names=['target', 'target_accession', 'query_name',
                    'query_accession', 'evalue', 'score', 'bias',
                    'domain_score', 'domain_bias', 'exp', 'reg', 'clu',
                    'ov', 'env', 'dom', 'rep', 'inc', 'description'])\
            .set_index('query_name')

MMSEQS2_OUTPUT_FORMAT = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qseq,tseq'

def _logged_subprocess_call(cmd):
    import subprocess
    logging.debug(f'Running command: {" ".join(cmd)}')
    subprocess.check_call(cmd)

def get_ampsphere_mmseqs_match_local(args, seqs):
    mmseqs_bin = shutil.which('mmseqs')
    if not mmseqs_bin:
        logging.error('MMSeqs2 not found. Please install it first (you can use `conda install -c bioconda mmseqs2`)')
        sys.exit(1)
    logging.debug(f'Using MMSeqs2 binary at {mmseqs_bin}')

    mmseqs_db = maybe_download_ampsphere_mmseqs(args)
    if not mmseqs_db:
        logging.error('AMPSphere MMSeqs2 database not found. Please download it first or provide the path to the database using the --cache-dir flag')
        sys.exit(1)
    mmseqs_db = path.join(mmseqs_db, 'AMPSphere_latest.mmseqsdb')
    logging.info(f'Using AMPSphere MMSeqs2 database at {mmseqs_db}')
    with tempfile.TemporaryDirectory() as tmpdir:
        query_file = path.join(tmpdir, 'query.faa')
        output_file = path.join(tmpdir, 'output.tsv')
        with open(query_file, 'wt') as f:
            for (query_name, seq) in seqs:
                f.write(f'>{query_name}\n{seq}\n')
        _logged_subprocess_call(
            ['mmseqs', 'createdb', query_file, f'{query_file}.mmseqsdb'])
        _logged_subprocess_call(
            ['mmseqs', 'search',
                f'{query_file}.mmseqsdb', mmseqs_db,
             output_file + '.mmseqsdb',
             tmpdir + '/tmp'])

        _logged_subprocess_call(
            ['mmseqs', 'convertalis',
             f'{query_file}.mmseqsdb',
             mmseqs_db,
             output_file + '.mmseqsdb',
             output_file,
             '--format-output', MMSEQS2_OUTPUT_FORMAT])

        return pd.read_csv(output_file, sep='\t', header=None,
                    names=MMSEQS2_OUTPUT_FORMAT.split(','))\
            .set_index('query') \
            .sort_values('evalue').sort_index(kind='stable')

def maybe_download_ampsphere_faa(args):
    target = path.join(get_cache_directory(args), 'AMPSphere_v.2022-03.faa.gz')
    if path.exists(target):
        logging.debug(f'AMPSphere database already downloaded to {target}')
        if args.re_download_database:
            logging.debug(f'Force re-download enabled, re-downloading AMPSphere database')
        else:
            return target
    if args.no_download_database:
        return None
    URL = 'https://zenodo.org/records/6511404/files/AMPSphere_v.2022-03.faa.gz?download=1'
    r = requests.get(URL, stream=True, headers=REQUESTS_HEADER)
    with open_output(target, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return target

def get_ampsphere_exact_match_local(args, seqs):
    from macrel.fasta import fasta_iter
    faa = maybe_download_ampsphere_faa(args)
    if not faa:
        logging.error('AMPSphere database not found. Please download it first or provide the path to the database using the --cache-dir flag')
        sys.exit(1)
    seq2id = {}
    for h,seq in fasta_iter(faa):
        seq2id[seq] = h
    results = []
    for (query_name, seq) in seqs:
        # Walrus operator would be nice here, but it's Python 3.8+ (and we support 3.6)
        if seq in seq2id:
            hit = seq2id[seq]
            results.append((query_name, seq, hit))
        elif seq[0] == 'M' and (seq[1:] in seq2id):
            hit = seq2id[seq[1:]]
            results.append((query_name, seq, h))
        else:
            results.append((query_name, seq, 'No_Hit'))
    return pd.DataFrame(results, columns=['query_name', 'query', 'result'])\
            .set_index('query_name')


def get_ampsphere_exact_match(seq, query_name):
    '''Get exact match from AMPSphere API'''
    URL = f'https://ampsphere-api.big-data-biology.org/v1/search/sequence-match?query={seq}'
    response = requests.get(URL, headers=REQUESTS_HEADER)
    data = response.json()
    return pd.DataFrame.from_dict({query_name : data}, orient='index')

def get_ampsphere_mmseqs_match(seq, query_name):
    '''Get MMSeqs2 match from AMPSphere API'''
    query = f'>{query_name}\n{seq}'
    URL = f'https://ampsphere-api.big-data-biology.org/v1/search/mmseqs?query={query}'
    response = requests.get(URL, headers=REQUESTS_HEADER)
    data = response.json()
    return pd.DataFrame.from_dict(data)\
            .drop("alignment_strings", axis=1)\
            .set_index('query_identifier')

def get_ampsphere_hmmer_match(seq, query_name):
    '''Get HMMER match from AMPSphere API'''
    query = f'>{query_name}\n{seq}'
    URL = f'https://ampsphere-api.big-data-biology.org/v1/search/hmmer?query={query}'
    response = requests.get(URL, headers=REQUESTS_HEADER)
    data = response.json()
    if not data:
        return pd.DataFrame()
    return pd.DataFrame.from_dict(data)\
            .set_index('query_name')

