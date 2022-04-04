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

1. `Access` Identifiers (arbitrarily assigned by Macrel)
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

predicted_faas_doc = """- `macrel.out.all_orfs.fna`

Fasta file containing all predicted genes from the contigs (including large
proteins).

- `macrel.out.all_orfs.faa`

Fasta file containing the peptide sequences for all predicted genes
from the contigs (including large proteins).

- `macrel.out.smorfs.faa`

Fasta file containing the predicted small ORFs (length filtered ranging from 10
to 100 amino acids).

- `macrel.correspondence.tsv`

Table containing the cluster information as peptide id and corresponding genes
at 100% of identity at amino acids level.

- `macrel.genes_coordinates.tsv`

The gene prediction using Pyrodigal creates a complete set of information about
the gene itself and is given as a table:

   - artificial_name: unique identifier for each gene, consisting of the contig name followed by an ordinal ID
         of the sequence separated by an underscore.
         
   - contig: contig name as given in the original fasta
   
   - start/end: these fields give the leftmost coordinate and the rightmost coordinate of the predicted gene
   
   - strand: information about the sense (1 for forward strand genes, -1 for reverse strand genes).
   
   - confidence: confidence from 0 to 100 on the true gene prediction.
   
   - partial_begin/partial_end: False indicates the gene has a true boundary (a start or a stop), whereas a
                                True indicates the gene is "unfinished" at that edge (i.e. a partial gene).
   
   - gc_content: GC content of the gene sequence.
   
   - translation_table: NCBI translation table used in the prediction and translation
   
   - rbs_motif: RBS motif found by Prodigal (e.g. "AGGA" or "GGA", etc.)
   
   - rbs_spacer: number of bases between the start codon and the observed motif.
   
   - start_type: sequence of the start codon (usually ATG, GTG, TTG, or "Edge" if the gene has no start codon).

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

