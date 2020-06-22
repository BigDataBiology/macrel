from contextlib import contextmanager

@contextmanager
def open_output(ofile, mode='wt'):
    if ofile == '-' or ofile == '/dev/stdout':
        with open('/dev/stdout', mode=mode) as out:
            yield out
    else:
        from atomicwrites import atomic_write
        with atomic_write(ofile, mode=mode, overwrite=True) as out:
            yield out
