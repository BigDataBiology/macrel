header = """If you find Macrel useful, please cite:

> Santos-Junior, C.D. et al. Macrel: antimicrobial peptide screening in
> genomes and metagenomes. The PeerJ 8:e10555
>  https://doi.org/10.7717/peerj.10555

For more information, please read [the macrel
documentation](https://macrel.readthedocs.io) and use the [AMPsphere mailing
list](https://groups.google.com/g/ampsphere-users) for questions."""

footer = f"""If you used the `--outtag` argument, then the above files will be named using
that tag, instead of `macrel.out`"""

prediction_table_doc = f"""- `macrel.out.prediction.gz`

Compressed tab-separated table with the following columns

1. `Access` Identififiers (arbitrarily assigned by Macrel)
2. `Sequence` Predicted amino acid sequence
3. `AMP_family`: AMP family classified into anionic/cationic and
   cysteine-containing/linear coded as: `A`/`C` and `D`/`L`, respectively,
   followed by a `P` (for peptide), _e.g._, a cationic cysteine-containing
   peptide would be `CDP`.
4. `AMP_probability`: Probability that the peptide is antimicrobial
5. `Hemolytic`: Classification into hemolytic (Hemo) or non-hemolytic (NonHemo)
6. `Hemolytic_probability`: Probability of hemolytic activity

The table contains a header (as comments marked with the `#` character at the
start of the line) identifying the version of macrel used to generate these
results.

Note that, by default, only peptides predicted to be AMPs are output. If the
`--keep-negatives` flag is used, however, all sequences will be present in the
table."""

predicted_faas_doc = """- `macrel.out.all_orfs.faa`

Fasta file containing all predicted genes from the contigs (including large
proteins).

- `macrel.out.smorfs.faa`

Fasta file containing the predicted small ORFs (length filtered ranging from 10
to 100 amino acids)."""

megahit_output_doc = """- `example_metag.megahit_output`

This folder contains the full outputs of running megahit for assembly."""

readme_output_abundance_mode = f"""{header}

## Outputs for `abundance` mode

- `macrel.out.out.abundance.txt`

This file contains a table with two columns:

1. `smORF_accession`: smORFs IDs (i.e., the name of the sequences in the reference file)
2. `MACREL`: Peptide abundance in number of (short) reads mapped

{footer}
"""

readme_output_contigs_mode = f"""{header}

## Outputs for `contigs` mode

{prediction_table_doc}
{predicted_faas_doc}

{footer}
"""


readme_output_peptides_mode = f"""{header}

## Outputs for `peptides` mode

{prediction_table_doc}

{footer}
"""

readme_output_reads_mode = f"""{header}

## Outputs for `reads` mode

{prediction_table_doc}
{predicted_faas_doc}
{megahit_output_doc}

{footer}
"""
