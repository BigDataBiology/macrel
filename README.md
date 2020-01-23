# MACREL: Meta(genomic) AMP Classification and Retrieval

Pipeline to mine antimicrobial peptides (AMPs) from (meta)genomes.

[![Build Status](https://travis-ci.com/BigDataBiology/macrel.svg?branch=master)](https://travis-ci.com/BigDataBiology/macrel)
[![Documentation Status](https://readthedocs.org/projects/macrel/badge/?version=latest)](https://macrel.readthedocs.io/en/latest/?badge=latest)
[![license: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Install with Bioconda](https://anaconda.org/bioconda/macrel/badges/installer/conda.svg)](https://anaconda.org/bioconda/macrel)
[![Install with Bioconda](https://anaconda.org/bioconda/macrel/badges/downloads.svg)](https://anaconda.org/bioconda/macrel)

If you use this software in a publication please cite

>   FACS: antimicrobial peptide screening in genomes and metagenomes
>   Celio Dias Santos-Junior, Shaojun Pan, Xing-Ming Zhao, Luis Pedro Coelho
>   bioRxiv 2019.12.17.880385; doi:
>   [https://doi.org/10.1101/2019.12.17.880385](https://doi.org/10.1101/2019.12.17.880385)

(The preprint still uses the old name of the tool, _FACS_ and will be updated
soon).

**NOTES**: This is still a _work in progress_ and, while the results of the
tool should be correct, we continuously work on making Macrel easier to install
and use.

## License

GPLv3.

While Macrel as a whole is **GPL v3** licensed (to comply with it being used in
some of its dependencies, namely Peptides), the macrel-specific code is also
licensed under the **MIT** license.

## Install

The recommended method of installation is through
[bioconda](https://anaconda.org/bioconda/macrel):

```bash
conda install -c bioconda macrel
```

To install from source, [read the
docs](https://macrel.readthedocs.io/en/latest/install)

### Examples

> Macrel uses a _subcommand interface_. You run `macrel COMMAND ...` with the
> COMMAND specifying which components of the pipeline you want to use.

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
