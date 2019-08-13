while read a b;
do
	echo "$a -- $b"
	time /home/celio/Desktop/Lib/FACS14.sh -m r --fwd /home/celio/Desktop/benchmark/metagenomes/fastq-files/$a --rev /home/celio/Desktop/benchmark/metagenomes/fastq-files/$b --outfolder /home/celio/Desktop/benchmark/analysis/metagenomes/ --outtag ${a/_1.fastq.gz/} -t 3 --block 100M --log ${a/_1.fastq.gz/}.log --clust 0 > /home/celio/Desktop/benchmark/analysis/metagenomes/bashlog.${a/_1.fastq.gz/}.txt
done < readslist.txt
