#!/usr/bin/bash

FACS="/home/celio/Desktop/Lib/FACS11.sh"

echo "[Doing metagenomes lists]"
ls /home/celio/Desktop/benchmark/metagenomes/fastq-files/*_1.fastq.gz > list1
ls /home/celio/Desktop/benchmark/metagenomes/fastq-files/*_2.fastq.gz > list2
paste -d'\t' list1 list2 > input.file; rm -rf list*

sed -i 's/\/home\/celio\/Desktop\/benchmark\/metagenomes\/fastq-files\///g' input.file

echo "[Starting FACS]"
while read a b;
do
	echo "[Starting FACS with samples $a and $b]"
	time $FACS -m r --fwd /home/celio/Desktop/benchmark/metagenomes/fastq-files/"$a" --rev /home/celio/Desktop/benchmark/metagenomes/fastq-files/"$b" --outtag pe_"${a/_1.fastq.gz/}" --outfolder /home/celio/Desktop/benchmark/analysis/metagenomes/ -t 3 --block 100M --log pe_"${a/_1.fastq.gz/}".log > /home/celio/Desktop/benchmark/analysis/metagenomes/logbash.pe_"${a/_1.fastq.gz/}".txt
	echo "[Starting FACS_abundance with samples $a and $b]"
	time $FACS -m mr --ref /home/celio/Desktop/benchmark/analysis/metagenomes/pe_"${a/_1.fastq.gz/}".tsv.gz --fwd "$a" --rev "$b" --outtag "${a/_1.fastq.gz/}"_abmr --outfolder /home/celio/Desktop/benchmark/analysis/metagenomes -t 3 --block 100M > /home/celio/Desktop/benchmark/analysis/metagenomes/logbash.ab_"${a/_1.fastq.gz/}".txt
done < input.file

rm -rf input.file
