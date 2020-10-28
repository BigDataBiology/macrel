from contextlib import contextmanager

@contextmanager
def open_output(ofile, mode='wt'):
    if ofile == '-' or ofile == '/dev/stdout':
        with open('/dev/stdout', mode=mode) as out:
            yield out
    else:
        try:
            from atomicwrites import atomic_write
        except ImportError:
            with open(ofile, mode=mode) as out:
                yield out
            return

        with atomic_write(ofile, mode=mode, overwrite=True) as out:
            yield out
