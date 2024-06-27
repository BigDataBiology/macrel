import requests
import pandas as pd


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

