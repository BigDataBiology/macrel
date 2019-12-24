Lib="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bash FACS.sh --mode p --fasta "$Lib"/example_seqs/expep.faa.gz --outfolder out_peptides --outtag example -t 4
bash FACS.sh --mode c --fasta "$Lib"/example_seqs/excontigs.fna.gz --outfolder out_contigs --outtag out_contigs -t 4
bash FACS.sh --mode r --fwd "$Lib"/example_seqs/R1.fq.gz --rev "$Lib"/example_seqs/R2.fq.gz --outfolder out_metag --outtag example_metag -t 4
bash FACS.sh --mode a --fwd "$Lib"/example_seqs/R1.fq.gz --fasta "$Lib"/example_seqs/ref.faa.gz --outfolder out_abundance --outtag example_abundance -t 4
#bash FACS.sh --mode a --fwd "$Lib"/example_seqs/R1.fq.gz --ref out_metag/example_metag.tsv.gz --outfolder out_abundance_pred --outtag example_abundance_pred -t 4