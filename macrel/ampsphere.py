import requests
import pandas as pd
from os import path
import logging

from macrel.utils import open_output

def get_cache_directory(args):
    '''Get cache directory'''
    from os import environ
    if args.cache_dir:
        cache_dir = args.cache_dir
    else:
        if 'XDG_CACHE_HOME' in environ:
            cache_dir = path.join(environ['XDG_CACHE_HOME'], 'ampsphere')
        else:
            cache_dir = path.join(path.expanduser("~"), '.cache', 'ampsphere')
    if not path.exists(cache_dir):
        from os import makedirs
        makedirs(cache_dir, exist_ok=True)
        with open_output(path.join(cache_dir, 'CACHEDIR.TAG'), 'wt') as f:
            f.write('Signature: 8a477f597d28d172789f06886806bc55\n')
            f.write('# This file is a cache directory tag created by macrel.\n')
            f.write('# For information about cache directory tags, see:\n')
            f.write('#\thttp://www.brynosaurus.com/cachedir/\n')
    return cache_dir

def maybe_download_ampsphere_faa(args):
    target = path.join(get_cache_directory(args), 'AMPSphere_v.2022-03.faa.gz')
    if path.exists(target):
        logging.debug(f'AMPSphere database already downloaded to {target}')
        if not args.no_download_database and args.force:
            logging.debug(f'Force download enabled, re-downloading AMPSphere database')
        else:
            return target
    if args.no_download_database:
        return None
    URL = 'https://zenodo.org/records/6511404/files/AMPSphere_v.2022-03.faa.gz?download=1'
    r = requests.get(URL, stream=True)
    with open_output(target, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return target

def get_ampsphere_exact_match_local(args, seqs):
    from macrel.fasta import fasta_iter
    faa = maybe_download_ampsphere_faa(args)
    if not faa:
        logging.error('AMPSphere database not found. Please download it first (use the --download-database flag) or provide the path to the database using the --cache-dir flag')
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
    response = requests.get(URL)
    data = response.json()
    return pd.DataFrame.from_dict({query_name : data}, orient='index')

def get_ampsphere_mmseqs_match(seq, query_name):
    '''Get MMSeqs2 match from AMPSphere API'''
    query = f'>{query_name}\n{seq}'
    URL = f'https://ampsphere-api.big-data-biology.org/v1/search/mmseqs?query={query}'
    response = requests.get(URL)
    data = response.json()
    return pd.DataFrame.from_dict(data)\
            .drop("alignment_strings", axis=1)\
            .set_index('query_identifier')

def get_ampsphere_hmmer_match(seq, query_name):
    '''Get HMMER match from AMPSphere API'''
    query = f'>{query_name}\n{seq}'
    URL = f'https://ampsphere-api.big-data-biology.org/v1/search/hmmer?query={query}'
    response = requests.get(URL)
    data = response.json()
    if not data:
        return pd.DataFrame()
    return pd.DataFrame.from_dict(data)\
            .set_index('query_name')

