/home/celio/Desktop/Lib/FACS15.sh -m c --fasta rep_contigs.fasta.gz --outfolder /home/celio/Desktop/benchmark/analysis/genomes/ --outtag rep_contigs_prodigal -t 3 --block 100M --log repcontigs_prodigal.log --clust 0

/home/celio/Desktop/Lib/spurius/spurius.mod.sh -m ../Lib/spurius/clf1.npy -d ../Lib/spurius/dict_np.npy -DB ../Lib/spurius/db/ --FACS rep_contigs_prodigal.tsv.gz 


