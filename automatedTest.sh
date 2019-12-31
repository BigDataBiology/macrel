#!/usr/bin/env bash

set -e

bash FACS.sh --mode p --fasta example_seqs/expep.faa.gz --outfolder out_peptides --outtag example -t 4
bash FACS.sh --mode c --fasta example_seqs/excontigs.fna.gz --outfolder out_contigs --outtag out_contigs -t 4
bash FACS.sh --mode r --fwd example_seqs/R1.fq.gz --rev example_seqs/R2.fq.gz --outfolder out_metag --outtag example_metag -t 4
bash FACS.sh --mode a --fwd example_seqs/R1.fq.gz --fasta example_seqs/ref.faa.gz --outfolder out_abundance --outtag example_abundance -t 4
bash FACS.sh --mode a --fwd example_seqs/R1.fq.gz --ref out_metag/example_metag.tsv.gz --outfolder out_abundance_pred --outtag example_abundance_pred -t 4

