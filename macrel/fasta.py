def fasta_iter(fname):
    header = None
    chunks = []
    if fname.endswith('.gz'):
        import gzip
        op = gzip.open
    else:
        op = open
    with op(fname, 'rt') as f:
        for line in f:
            if line[0] == '>':
                if header is not None:
                    yield header,''.join(chunks)
                line = line[1:].strip()
                if not line:
                    header = ''
                else:
                    header = line.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)

