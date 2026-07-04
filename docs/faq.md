# Frequent Asked Questions - FAQ

### Does Macrel detect post-translationally modified peptides?

Currently, Macrel only detects mature peptides. Masking sequences such as
leaders or signal peptides is not yet supported.

### Does Macrel detect peptides of any length?

Macrel targets peptides in the range of 10-100 residues. The upper bound is
enforced: predicted small genes longer than 100 amino acids are filtered out,
and the classifier will generally not return a meaningful result for longer
sequences. The lower bound is a practical recommendation rather than a hard
filter — Macrel will not crash on very short sequences, but predictions for
peptides shorter than ~10 residues are unlikely to be meaningful.

### What happens to sequences with ambiguous or non-standard residues?

Empty sequences and sequences containing non-standard amino acid codes (such as
`X`, `B`, `Z`, `U`, `O`, or an internal `*`) are skipped during feature
computation. Macrel emits a warning for each skipped sequence and continues with
the remaining input rather than aborting.

### Can Macrel work well on isolate genomes or isolate genome reads?

Macrel is able to handle isolated genomes, as contigs, short reads (in this
case it will trim and assemble before starting the gene prediction) or directly
on predicted peptide sequences (obtained through some alternative pipeline).
