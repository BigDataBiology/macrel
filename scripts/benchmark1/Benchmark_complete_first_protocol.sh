mkdir analysis
mkdir analysis/genomes
mkdir analysis/metagenomes
mkdir analysis/comparison

echo "Searching for AMP in the Representative Contigs"

time /home/celio/Desktop/Lib/FACS14.sh -m c --fasta /home/celio/Desktop/benchmark/genomes/representatives.contigs.fasta.gz --outfolder /home/celio/Desktop/benchmark/analysis/genomes/ --outtag repcontigs -t 3 --block 100M --log repcontigs.log --clust 0 > /home/celio/Desktop/benchmark/analysis/genomes/bashlog.repcontigs.txt

echo "Doing metagenomes lists"

ls /home/celio/Desktop/benchmark/metagenomes/fastq-files/*_1.fastq.gz > list1
ls /home/celio/Desktop/benchmark/metagenomes/fastq-files/*_2.fastq.gz > list2
paste -d'\t' list1 list2 > readslist.txt; rm -rf list*

while read a b;
do
	echo "Searching for AMP in the metagenomes: $a -- $b"

	time /home/celio/Desktop/Lib/FACS14.sh -m r --fwd /home/celio/Desktop/benchmark/metagenomes/fastq-files/$a --rev /home/celio/Desktop/benchmark/metagenomes/fastq-files/$b --outfolder /home/celio/Desktop/benchmark/analysis/metagenomes/ --outtag ${a/_1.fastq.gz/} -t 3 --block 100M --log ${a/_1.fastq.gz/}.log --clust 0 > /home/celio/Desktop/benchmark/analysis/metagenomes/bashlog.${a/_1.fastq.gz/}.txt

	echo "Calculating abundance of Representative contigs' AMPs using metagenomes: $a -- $b"
	time /home/celio/Desktop/Lib/FACS14.sh -m mr --fwd /home/celio/Desktop/benchmark/metagenomes/fastq-files/$a --rev /home/celio/Desktop/benchmark/metagenomes/fastq-files/$b --ref /home/celio/Desktop/benchmark/analysis/genomes/representative_contigs.tsv.gz --outfolder /home/celio/Desktop/benchmark/analysis/abundances/ --outtag ${a/_1.fastq.gz/} -t 3 --block 100M --log ${a/_1.fastq.gz/}.log > /home/celio/Desktop/benchmark/analysis/abundances/bashlog.${a/_1.fastq.gz/}.txt

done < readslist.txt

echo "Doing comparisons"

zcat /home/celio/Desktop/benchmark/analysis/genomes/repcontigs.tsv.gz | awk '{print ">"$1"\n"$2}' > analysis/comparison/repcontigs.AMP.faa

for i in $(ls analysis/metagenomes/ | grep ".gz"); do zcat analysis/metagenomes/$i | awk '{print ">"$1"\n"$2}' > analysis/comparison/${i/.tsv.gz/.AMP.faa}; done

cd analysis/comparison/
~/diamond makedb --in repcontigs.AMP.faa --db REPCONTIGS
~/diamond makedb --in ADAM.faa --db ADAM
~/diamond makedb --in representativeproteins.faa --db REPprot

for i in from_*.AMP.faa; do ~/diamond blastp -d REPCONTIGS -o ${i/.AMP.faa/.m8} -q $i --sensitive; done
for i in from_*.AMP.faa; do ~/diamond blastp -d ADAM -o ${i/.AMP.faa/.m8.adam} -q $i --sensitive; done
for i in from_*.AMP.faa; do ~/diamond blastp -d REPprot -o ${i/.AMP.faa/REPprot.m8} -q $i; done

for i in $(ls from_*AMP.faa); do grep -c ">" $i >> file1; done
for i in $(ls *.m8); do  cut -f1 $i | sort | uniq | wc -l >> file2; done
for i in $(ls *.m8.adam); do  cut -f1 $i | sort | uniq | wc -l >> file3; done

paste -d'\t' file1 file2 > tmp; rm -rf file1 file2
paste -d'\t' tmp file3 > tmp2; rm -rf tmp file3
echo -e "AMPs\tDetect_on_REP\tDetect_on_ADAM" > header
cat header tmp2 > comparison; rm -rf header tmp2
