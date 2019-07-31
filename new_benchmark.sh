FACS=/home/celio/Desktop/Lib/FACS15.sh
spurius=/home/celio/Desktop/Lib/spurius/spurius.mod.sh

#echo "Searching for AMP in the Representative Contigs"

#time $FACS15 -m c --fasta /home/celio/Desktop/benchmark/genomes/representatives.contigs.fasta.gz --outfolder /home/celio/Desktop/benchmark/analysis/genomes/ --outtag repcontigs -t 3 --block 100M --log repcontigs.log --clust 0 > /home/celio/Desktop/benchmark/analysis/genomes/bashlog.repcontigs.txt

#echo "Doing metagenomes lists"

#ls /home/celio/Desktop/benchmark/metagenomes/fastq-files/*_1.fastq.gz > list1
#ls /home/celio/Desktop/benchmark/metagenomes/fastq-files/*_2.fastq.gz > list2
#paste -d'\t' list1 list2 > readslist.txt; rm -rf list*

while read a b;
do
	echo "Searching for AMP in the metagenomes: $a -- $b"

	time $FACS -m r --fwd /home/celio/Desktop/benchmark/metagenomes/fastq-files/$a --rev /home/celio/Desktop/benchmark/metagenomes/fastq-files/$b --outfolder /home/celio/Desktop/benchmark/analysis/metagenomes/ --outtag ${a/_1.fastq.gz/} -t 3 --block 100M --log ${a/_1.fastq.gz/}.log --clust 0 > /home/celio/Desktop/benchmark/analysis/metagenomes/bashlog.${a/_1.fastq.gz/}.txt

	echo "Calculating spurious ORFs"

	$spurius -m ../Lib/spurius/clf1.npy -d ../Lib/spurius/dict_np.npy -DB ../Lib/spurius/db/ --FACS /home/celio/Desktop/benchmark/analysis/metagenomes/${a/_1.fastq.gz/}.tsv.gz
	mv spurious_out.txt /home/celio/Desktop/benchmark/analysis/metagenomes/${a/_1.fastq.gz/}.spuriuspred.txt

	echo "Calculating abundance of Representative contigs' AMPs using metagenomes: $a -- $b"

	time $FACS -m mr --fwd /home/celio/Desktop/benchmark/metagenomes/fastq-files/$a --rev /home/celio/Desktop/benchmark/metagenomes/fastq-files/$b --ref /home/celio/Desktop/benchmark/analysis/genomes/FACS15/rep_contigs_prodigal.tsv.gz --outfolder /home/celio/Desktop/benchmark/analysis/abundances/ --outtag ${a/_1.fastq.gz/} -t 3 --block 100M --log ${a/_1.fastq.gz/}.log > /home/celio/Desktop/benchmark/analysis/abundances/bashlog.${a/_1.fastq.gz/}.txt

done < readslist.txt
