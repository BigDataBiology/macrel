header = """If you find Macrel useful, please cite:

> Santos-Junior, C.D. et al. Macrel: antimicrobial peptide screening in
> genomes and metagenomes. The PeerJ 8:e10555
> https://doi.org/10.7717/peerj.10555

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
table.
"""

predicted_faas_doc = """- `macrel.out.all_orfs.faa`

Fasta file containing all predicted genes from the contigs (including large
proteins).

- `macrel.out.smorfs.faa`

Fasta file containing the predicted small ORFs (length filtered ranging from 10
to 100 amino acids).

In the both files the original header structure from Prodigal indicates the AMP
gene origins, and is shown if the option `--keep-fasta-headers` is used. The 
header structure has a consistent notation, e.g.:

Header for the 2nd small protein in the contig k77_5:

>k77_5_2 # 1984 # 2010 # 1 # ID=2_2;partial=00;
start_type=ATG;rbs_motif=TAAAAAA;rbs_spacer=7bp;gc_cont=0.487

The next three fields in the header, delimited by "#" signs, are the: (1) leftmost coordinate,
(2) rightmost coordinate, and (3) the strand (1 for forward strand genes, -1 for reverse strand genes).

Following the coordinate information is a semicolon-delimited string with he following fields:

   - ID: unique identifier for each gene, consisting of the ordinal ID of the sequence and an ordinal ID of that
         gene within the sequence (separated by an underscore).
   
   - partial: "0" indicates the gene has a true boundary (a start or a stop), whereas a "1" indicates the gene is
              "unfinished" at that edge (i.e. a partial gene).
   
   - start_type: sequence of the start codon (usually ATG, GTG, TTG, or "Edge" if the gene has no start codon).
   
   - stop_type: sequence of the stop codon (usually TAA, TGA, TAG or "Edge" if the gene has no stop codon).
   
   - rbs_motif: RBS motif found by Prodigal (e.g. "AGGA" or "GGA", etc.)
   
   - rbs_spacer: number of bases between the start codon and the observed motif.
   
   - gc_cont: GC content of the gene sequence.
   
"""

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

