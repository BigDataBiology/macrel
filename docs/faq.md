# Frequent Asked Questions - FAQ

### Does Macrel detect post-translationally modified peptides?

Currently, Macrel only detects mature peptides. We expect to develop and test
approaches to mask sequences such as leaders or signal peptides too in the next
versions of the program.

### Does Macrel detect peptides of any length?

Currently, Macrel detects peptides in the range of 10-100 residues. The
classifier will generally not return meaningful result for longer sequences.

### Can Macrel work well on isolate genomes or isolate genome reads?

Macrel is able to handle isolated genomes, as contigs, short reads (in this
case it will trim and assemble before starting the gene prediction) or directly
on predicted peptide sequences (obtained through some alternative pipeline).
