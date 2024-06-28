import argparse
import subprocess
import tempfile
import gzip
import logging
import os
from os import path, makedirs
import textwrap
from .utils import open_output

def error_exit(args, errmessage):
    import sys
    sys.stderr.write(errmessage+'\n')
    sys.exit(1)

def data_file(fname):
    return path.join(path.dirname(__file__),
                    'data',
                    fname)

def parse_args(args):
    from .macrel_version import __version__



    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='macrel v{}'.format(__version__), epilog=textwrap.dedent('''\
             Examples:
                 run Macrel on peptides:
                 macrel peptides --fasta example_seqs/expep.faa.gz --output out_peptides
                 
                 run Macrel on contigs:
                 macrel contigs --fasta example_seqs/excontigs.fna.gz --output out_contigs
                 
                 run Macrel on paired-end reads:
                 macrel reads -1 example_seqs/R1.fq.gz -2 example_seqs/R2.fq.gz --output out_metag --outtag example_metag
                 
                 run Macrel to get abundance profiles:
                 macrel abundance -1 example_seqs/R1.fq.gz --fasta example_seqs/ref.faa.gz --output out_abundance --outtag example_abundance

                 Query the AMPSphere database
                 macrel query-ampsphere --query-mode mmseqs --fasta peptides.faa
                 
                 For more information,please read the docs: https://macrel.readthedocs.io/en/latest/
             '''))

    parser.add_argument('command', nargs=1,
            help='Macrel command to execute (see documentation)')
    parser.add_argument('-t', '--threads', required=False, action='store',
            help='Number of threads to use',
            default='1', dest='threads')
    parser.add_argument('-o', '--output', required=False, default=None,
            help='path to the output directory', dest='output')
    parser.add_argument('--file-output', required=False, default=None,
            help='path to the output file', dest='output_file')
    parser.add_argument('--tag', required=False, default='macrel.out',
            help='Set output tag', dest='outtag')
    parser.add_argument('-f', '--fasta', required=False,
            help='path to the input FASTA file. This is used in both the peptides command (where the file is expected to contain short amino-acid sequences) and in the contigs command (where the file is expected to contain longer nucleotide contigs)', dest='fasta_file')
    parser.add_argument('-1', '--reads1', required=False, action='store',
            dest='reads1')
    parser.add_argument('-2', '--reads2', required=False, action='store',
            default=None, dest='reads2')
    parser.add_argument('--mem', required=False, action='store', default='0.75', )
    parser.add_argument('--cluster', required=False, action='store_true', default=False, dest='cluster',
            help='Whether to pre-cluster the smORFs (at 100%% identity) to avoid repeats')
    parser.add_argument('--force', required=False, action='store_true')
    parser.add_argument('--keep-fasta-headers', required=False, action='store_true', default=False, dest='keep_fasta_headers',
            help='Keep complete FASTA headers [get-smorfs command]')
    parser.add_argument('--tmpdir', required=False, default=None, dest='tmpdir', action='store',
            help='Temporary directory to use (default: $TMPDIR in the environment or /tmp)')
    parser.add_argument('--keep-negatives', required=False, default=False, dest='keep_negatives', action='store_true',
            help='Whether to keep non-AMPs in the output')
    parser.add_argument('--version', '-v', action='version',
                    version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('--log-file', required=False, default=None, dest='logfile',
            help='Path to the output logfile')
    parser.add_argument('--log-append', required=False, action='store_true',
            default=False, dest='logfile_append',
            help='If set, then the log file is appended to (default: overwrite existing file)')

    parser.add_argument('--query-mode', required=False, default='exact',
                        help='Query mode to use in the AMPSphere query (options: exact, mmseqs, hhm)', dest='query_mode')

    parser.add_argument('--local', required=False, default=False, action='store_true',
                        help='Use local AMPSphere database', dest='local')

    parser.add_argument('--re-download-database', required=False, default=False, action='store_true',
                        help='Download the AMPSphere database even if it already was downloaded before', dest='re_download_database')
    parser.add_argument('--no-download-database', required=False, default=False, action='store_true',
                        help='Do not download the AMPSphere database', dest='no_download_database')

    parser.add_argument('--cache-dir', required=False, default=None, dest='cache_dir',
                        help='Directory to use for caching AMPSphere data')

    return parser.parse_args()

def validate_args(args):
    '''Checks that args are consistent

    Exits with an error message if not the case
    '''
    from os import path
    if len(args.command) != 1:
        error_exit(args, 'Could not parse argument list (multiple commands given?)')
    [args.command] = args.command
    if args.command == 'get-examples':
        return
    if args.command == 'peptides':
        if not args.fasta_file:
            error_exit(args, "PEPTIDE File is necessary for 'peptides' command.")
    elif args.command == 'contigs':
        if not args.fasta_file:
            error_exit(args, "FASTA File is necessary for 'contigs' command.")
    elif args.command == 'reads':
        if not args.reads1:
            error_exit(args, "FQ file is necessary for 'reads' command.")
    elif args.command == 'abundance':
        if not args.reads1:
            error_exit(args, "FQ file is necessary for 'abundance' command.")
        if not args.fasta_file:
            error_exit(args, "Fasta file is necessary for 'abundance' command.")
    elif args.command == 'get-smorfs':
        if not args.fasta_file:
            error_exit(args, "Fasta file is necessary for 'get-smorfs' command.")
    elif args.command == 'query-ampsphere':
        if not args.fasta_file:
            error_exit(args, "Fasta file is necessary for 'query-ampsphere' command.")
        if args.query_mode not in ['exact', 'mmseqs', 'hmmer']:
            error_exit(args, f"Unknown query mode {args.query_mode} (options: exact, mmseqs, hmmer)")
    else:
        error_exit(args, "Unknown command {}".format(args.command))
    if not args.output and not args.output_file:
        error_exit(args, "Either --output or --file-output argument must be used")

    if args.keep_fasta_headers and args.command in ['peptides', 'abundance']:
        error_exit(args, '--keep-fasta-headers is not available for peptides and abundance commands')

    args.output_dir = args.output
    if args.output_dir:
        if not path.exists(args.output_dir):
            makedirs(args.output, exist_ok=True)
        elif args.force:
            import sys
            sys.stderr.write("Output folder already exists, but --force flag was used")
        else:
            error_exit(args, "Output folder [{}] already exists".format(args.output_dir))
    elif args.command != 'get-smorfs':
        error_exit(args, '--file-output is only possible for `get-smorfs` command')

    if args.logfile:
        if not path.exists(os.path.dirname(args.logfile)) and os.path.dirname(args.logfile) != '':
            makedirs(os.path.dirname(args.logfile), exist_ok=True)

def do_smorfs(args, tdir,logfile):
    from .ORFs_prediction import predict_genes
    from .filter_smorfs import filter_smorfs
    import sys
    
    if args.output_dir:
        all_peptide_file = path.join(args.output, args.outtag+'.all_orfs.faa')
        peptide_file = path.join(args.output, args.outtag+'.smorfs.faa')
    else:
        all_peptide_file = path.join(tdir, 'all_orfs.faa')
        peptide_file = (args.output_file if args.output_file != '-' else '/dev/stdout')

    # predict genes with pyrodigal
    clen = predict_genes(args.fasta_file, all_peptide_file)
    filter_smorfs(all_peptide_file, peptide_file, args.cluster, args.keep_fasta_headers)
    args.fasta_file = peptide_file
    return clen
    

def link_or_uncompress_fasta_file(orig, dest):
    '''
    If the input is compress, uncompress it. Otherwise, link it to `dest`
    '''
    if orig.endswith('.gz'):
        logging.debug('Uncompressing FASTA file ({})'.format(orig))
        with open(dest, 'wb') as ofile:
            with gzip.open(orig, 'rb') as ifile:
                while True:
                    chunk = ifile.read(32*1024)
                    if not chunk:
                        break
                    ofile.write(chunk)
    else:
        os.symlink(path.abspath(orig), dest)
    return dest

def do_abundance(args, tdir,logfile):
    do_read_trimming(args, tdir,logfile)
    # Even when paired reads, to do the abundance, 
    # only one pair is used to make the measures
    # comparable
    if args.reads2:
        trimmed_file = path.join(tdir, 'preproc.pair.1.fq.gz')
    else:
        trimmed_file = path.join(tdir, 'preproc.fq.gz')
    sam_file = path.join(tdir, 'paladin.out.sam')
    fasta_file = link_or_uncompress_fasta_file(
                        args.fasta_file,
                        path.join(tdir, 'paladin.faa'))

    subprocess.check_call([
        'paladin', 'index',

        # -r<#>  Reference type:
        #     1: Reference contains nucleotide sequences (requires corresponding .gff annotation)
        #     2: Reference contains nucleotide sequences (coding only, eg curated transcriptome)
        #     3: Reference contains protein sequences (UniProt or other source)
        #     4: Development tests
        '-r3',
        fasta_file],
        stdout=logfile)
    logging.debug('Mapping reads against references')
    with open(sam_file, 'wb') as sout:
        subprocess.check_call([
            'paladin', 'align',

            # -t INT        number of threads [1]
            '-t', str(args.threads),

            # -T INT        minimum score to output [15]
            '-T', '20',

            # -f INT        minimum ORF length accepted (as constant value) [250]
            '-f', '10',

            # -z INT[,...]  Genetic code used for translation (-z ? for full list) [1]
            '-z', '11',

            # -a            output all alignments for SE or unpaired PE
            '-a',

            #-V            output the reference FASTA header in the XR tag
            '-V',

            # -M            mark shorter split hits as secondary
            '-M',
            fasta_file,
            trimmed_file],
            stdout=sout)
    subprocess.check_call([
        'ngless',
        '--no-create-report',
        '--quiet',
        '-j', str(args.threads),
        data_file('scripts/count.ngl'),
        sam_file,
        path.join(args.output, args.outtag + '.abundance.txt')],
        stdout=logfile)

def do_read_trimming(args, tdir,logfile):
    ofile = path.join(tdir, 'preproc.fq.gz')
    if args.reads2:
        ngl_file = data_file('scripts/trim.pe.ngl')
        ngl_args = [args.reads1, args.reads2, ofile]
    else:
        ngl_file = data_file('scripts/trim.se.ngl')
        ngl_args = [args.reads1, ofile]
    subprocess.check_call([
        'ngless',
        '--no-create-report',
        '--quiet',
        '-j', str(args.threads),
        ngl_file,
        ] + ngl_args,
        stdout=logfile)

def do_assembly(args, tdir,logfile):
    if args.reads2:
        megahit_args = ['-1', path.join(tdir, 'preproc.pair.1.fq.gz'),
                        '-2', path.join(tdir, 'preproc.pair.2.fq.gz')]
    else:
        megahit_args = ['-r', path.join(tdir, 'preproc.fq.gz')]
    megahit_output = path.join(args.output, args.outtag + '.megahit_output')
    do_read_trimming(args, tdir,logfile)
    subprocess.check_call([
        'megahit',
        '--presets', 'meta-large',
        '-o', megahit_output, # output directory
        '-t', str(args.threads),
        '-m', str(args.mem),
        '--min-contig-len', '1000',
        ] + megahit_args,
        stdout=logfile)
    args.fasta_file = path.join(megahit_output, 'final.contigs.fa')

def do_predict(args, tdir):
    # These imports are slow, so we do them inside the functions
    from .AMP_features import fasta_features
    from .AMP_predict import predict
    import gzip
    fs = fasta_features(args.fasta_file)
    prediction = predict(
                         data_file("models/AMP.pkl.gz"),
                         data_file("models/Hemo.pkl.gz"),
                         fs,
                         args.keep_negatives)
    ofile = path.join(args.output, args.outtag + '.prediction.gz')
    with open_output(ofile, mode='wb') as raw_out:
        with gzip.open(raw_out, 'wt') as out:
            from .macrel_version import __version__
            out.write('# Prediction from macrel v{}\n'.format(__version__))
            prediction.to_csv(out, sep='\t', index_label='Access', float_format="%.3f")
    return prediction

def do_density(args, clen, prediction):
    tpred = prediction.reset_index()
    tpred['contig'] = tpred['index'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    tpred = tpred.query('AMP_probability > 0.5')
    tpred = tpred.groupby('contig').agg('size')
    tpred = tpred.reset_index()
    tpred = tpred.rename({0: 'AMPs'}, axis=1)
    clen = clen.merge(on='contig', right=tpred, how='outer').fillna(0)
    clen[clen.columns[1:]] = clen[clen.columns[1:]].astype(int)
    ofile = path.join(args.output, args.outtag + '.percontigs.gz')
    sample = clen.set_index('contig').sum(axis=0).tolist()
    sample_density = sample[-1] * 1e6 / sample[0]
    clen.sort_values('contig', inplace=True)
    with open_output(ofile, mode='wb') as raw_out:
        with gzip.open(raw_out, 'wt') as out:
            from .macrel_version import __version__
            out.write('# Prediction from macrel v{}\n'.format(__version__))
            out.write(f'# Macrel calculated for the sample a density of {sample_density:.3f} AMPs / Mbp.\n')
            clen.to_csv(out, sep='\t', index=False, float_format="%.3f")
    print(f'Macrel processed the sample and verified a density of {sample_density:.3f} AMPs / Mbp.')

def do_get_examples(args):
    from urllib.request import urlretrieve

    DATA_FILES = [
        'excontigs.fna.gz',
        'expep.faa.gz',
        'R1.fq.gz',
        'R2.fq.gz',
        'ref.faa.gz',
        ]
    BASEURL = 'https://github.com/BigDataBiology/macrel/raw/main/example_seqs/'
    if path.exists('example_seqs') and not args.force:
        error_exit(args, 'example_seqs/ directory already exists')
    makedirs('example_seqs', exist_ok=True)
    for f in DATA_FILES:
        print('Retrieving {}...'.format(f))
        urlretrieve(BASEURL + f, 'example_seqs/'+f)

def do_ampsphere_query(args):
    from . import ampsphere
    from .fasta import fasta_iter
    from time import sleep
    import pandas as pd
    if args.local:
        if args.query_mode == 'hmmer':
            error_exit(args, 'Local mode does not support HMMER queries')
        match_function = {
                'exact': ampsphere.get_ampsphere_exact_match_local,
                'mmseqs': ampsphere.get_ampsphere_mmseqs_match_local,
                }[args.query_mode]
        results = match_function(args, fasta_iter(args.fasta_file))
    else:
        match_function = {
            'exact': ampsphere.get_ampsphere_exact_match,
            'mmseqs': ampsphere.get_ampsphere_mmseqs_match,
            'hmmer': ampsphere.get_ampsphere_hmmer_match,
            }[args.query_mode]
        results = []
        logging.debug(f'Calling the AMPSphere API in {args.query_mode} mode')
        for h,seq in fasta_iter(args.fasta_file):
            results.append(match_function(seq, h))
            sleep(0.1)
            if len(results) == 20:
                import sys
                if sys.stdout.isatty():
                    print('Note that to avoid overloading the AMPSphere API, this script will pause a bit after every query')
                    if args.query_mode != 'hmmer':
                        print('You can use the --local flag to download and use a local version of the AMPSphere database')
        results = pd.concat(results)
        results.index.name = 'query_name'
        results.fillna({'result': 'No_Hit'}, inplace=True)
    ofile = (args.output_file if args.output_file != '-' else '/dev/stdout')
    if ofile is None:
        ofile = path.join(args.output,
                          f'{args.outtag}.ampsphere_{args.query_mode}.tsv.gz')
    with open_output(ofile, mode='wb') as raw_out:
        with gzip.open(raw_out, 'wt') as out:
            out.write(f'# AMPSphere query results (mode: {args.query_mode})\n')
            results.to_csv(out, sep='\t', index=True)

def main(args=None):
    if args is None:
        import sys
        args = sys.argv
    args = parse_args(args)
    validate_args(args)

    if args.logfile:
        logfile = open(args.logfile, ('a' if args.logfile_append else 'w'))
    else:
        logfile = None

    if args.command == 'get-examples':
        do_get_examples(args)
        return

    with tempfile.TemporaryDirectory(dir=args.tmpdir) as tdir:
        from .output import readme_output_abundance_mode,readme_output_contigs_mode,\
            readme_output_peptides_mode,readme_output_reads_mode, readme_output_ampsphere_mode

        creadme = {'reads': readme_output_reads_mode,
                   'contigs': readme_output_contigs_mode,
                   'get-smorfs': readme_output_contigs_mode,
                   'peptides': readme_output_peptides_mode,
                   'abundance': readme_output_abundance_mode,
                   'query-ampsphere': readme_output_ampsphere_mode,
                  }
        # print readme
        if args.output_file != '-':
            with open_output(os.path.join(args.output, 'README.md')) as ofile:
                ofile.write(creadme[args.command])
        
        # commands
        if args.command == 'reads':
            do_assembly(args, tdir,logfile)
        if args.command in ['reads', 'contigs', 'get-smorfs']:
            clen = do_smorfs(args, tdir,logfile)
        if args.command in ['reads', 'contigs', 'peptides']:
            prediction = do_predict(args, tdir)
        if args.command in ['reads', 'contigs']:
            if not args.cluster:
                do_density(args, clen, prediction)
        if args.command == 'abundance':
            do_abundance(args, tdir,logfile)
        if args.command == 'query-ampsphere':
            do_ampsphere_query(args)

if __name__ == '__main__':
    import sys
    main(sys.argv)
