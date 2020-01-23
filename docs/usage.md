# Manual

> Macrel uses a _subcommand interface_. You run `macrel COMMAND ...` with the
> COMMAND specifying which components of the pipeline you want to use.


## Subcommands

- `peptides`: to classify a fasta containing peptide sequences in a fasta file,
- `contigs`: to input pre-assembled contigs in a fasta file,
- `reads`: to input reads in fastQ format (Macrel accepts single- and paired-end reads),
- `abundance`: to measure abundance of a peptides fasta file using a given fastQ file,
- `get-smorfs`: to predict small genes from a contigs fasta file.

## Mandatory input flags

- `--out`: Folder destination to output [for all subcommands]
- `--file-output`: Can be used in the `get-smorfs` subcommand to output just
  the final file. If this is used, `--out` is no longer mandatory
- `--fasta`: Path to the input FASTA file. This is used in both the `peptides`
  subcommands (where the file is expected to contain short amino-acid
  sequences) and in the `contigs/get-smorfs` subcommands (where the file is
  expected to contain longer nucleotide contigs)
- `--reads1/--reads2`: Path to the input FastQ files [for `reads` subcommand]


## Optional flags

- `--threads` (`-t`): Specify the number of cpus used
- `--cluster`: Pre-cluster the smORFs (at 100%% identity) to avoid repeats
- `--tag`: If used, then all the output files will incorporate it.
- `--mem`: Used during assembly to specify the percent of RAM used (1-0)
- `--tmpdir`: Temporary directory to use (default: `$TMPDIR` in the environment or /tmp)
- `--force`: Potentially over-write existing output files
- `--keep-fasta-headers`: Keep complete FASTA headers [get-smorfs command]

Also, `macrel --version` will show the version of macrel installed.


### Examples

To run Macrel on peptides, use the `peptides` subcommand:

```bash
macrel peptides \
    --fasta example_seqs/expep.faa.gz \
    --output out_peptides \
    -t 4
```

In this case, we use `example_seqs/expep.faa.gz` as input sequence. This should
be an amino-acid FASTA file. The outputs will be written into a folder called
`out_peptides`, and Macrel will 4 threads.

To run Macrel on contigs, use the `contigs` subcommand:

```bash
macrel contigs \
    --fasta example_seqs/excontigs.fna.gz \
    --output out_contigs
```

In this example, we use the example file `excontigs.fna.gz` which is a FASTA
file with nucleotide sequences, writing the output to `out_contigs`.

To run Macrel on paired-end reads, use the `reads` subcommand:

```bash
macrel reads \
    -1 example_seqs/R1.fq.gz \
    -2 example_seqs/R2.fq.gz \
    --output out_metag \
    --outtag example_metag
```

The paired-end reads are given as paired files (here, `example_seqs/R1.fq.gz`
and `example_seqs/R2.fq.gz`). If you only have single-end reads, you can omit
the `-2` argument.

To run Macrel to get abundance profiles, you only need the short reads file
and a reference with peptide sequences. Use the `abundance` subcommand:


```bash
macrel abundance \
    -1 example_seqs/R1.fq.gz \
    --fasta example_seqs/ref.faa.gz \
    --output out_abundance \
    --outtag example_abundance
```
