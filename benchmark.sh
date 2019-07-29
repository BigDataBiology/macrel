#!/usr/bin/env bash

set -e

BENCHMARK_DIR=$PWD/BENCHMARK

mkdir -p $BENCHMARK_DIR/analysis
mkdir -p $BENCHMARK_DIR/analysis/genomes
mkdir -p $BENCHMARK_DIR/analysis/metagenomes
mkdir -p $BENCHMARK_DIR/analysis/comparison

echo "Searching for AMP in the Representative Contigs"

time FACSv14.sh -m c \
    --fasta $BENCHMARK_DIR/genomes/representatives.contigs.fasta.gz \
    --outfolder $BENCHMARK_DIR/analysis/genomes/ \
    --outtag repcontigs \
    -t 3 \
    --block 100M \
    --log repcontigs.log \
    --clust 0 > $BENCHMARK_DIR/analysis/genomes/bashlog.repcontigs.txt

echo "Doing metagenomes lists"

ls $BENCHMARK_DIR/metagenomes/fastq-files/*_1.fastq.gz > list1
ls $BENCHMARK_DIR/metagenomes/fastq-files/*_2.fastq.gz > list2
paste -d'\t' list1 list2 > readslist.txt; rm -rf list*

while read a b;
do
	echo "Searching for AMP in the metagenomes: $a -- $b"

    base_a=$(basename $a)
	time FACSv14.sh \
        -m r \
        --fwd $a \
        --rev $b \
        --outfolder $BENCHMARK_DIR/analysis/metagenomes/ \
        --outtag ${base_a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${base_a/_1.fastq.gz/}.log \
        --clust 0 > $BENCHMARK_DIR/analysis/metagenomes/bashlog.${base_a/_1.fastq.gz/}.txt

	echo "Calculating abundance of Representative contigs' AMPs using metagenomes: $a -- $b"
	time FACSv14.sh -m mr \
        --fwd $a \
        --rev $b \
        --ref $BENCHMARK_DIR/analysis/genomes/representative_contigs.tsv.gz \
        --outfolder $BENCHMARK_DIR/analysis/abundances/ \
        --outtag ${base_a/_1.fastq.gz/} \
        -t 3 \
        --block 100M \
        --log ${base_a/_1.fastq.gz/}.log > $BENCHMARK_DIR/analysis/abundances/bashlog.${base_a/_1.fastq.gz/}.txt

done < readslist.txt

echo "Doing comparisons"

zcat $BENCHMARK_DIR/analysis/genomes/repcontigs.tsv.gz | awk '{print ">"$1"\n"$2}' > analysis/comparison/repcontigs.AMP.faa

for i in $(ls analysis/metagenomes/ | grep ".gz"); do zcat analysis/metagenomes/$i | awk '{print ">"$1"\n"$2}' > analysis/comparison/${i/.tsv.gz/.AMP.faa}; done

cd analysis/comparison/
diamond makedb --in repcontigs.AMP.faa --db REPCONTIGS
diamond makedb --in ADAM.faa --db ADAM
diamond makedb --in representativeproteins.faa --db REPprot

for i in from_*.AMP.faa; do diamond blastp -d REPCONTIGS -o ${i/.AMP.faa/.m8} -q $i --sensitive; done
for i in from_*.AMP.faa; do diamond blastp -d ADAM -o ${i/.AMP.faa/.m8.adam} -q $i --sensitive; done
for i in from_*.AMP.faa; do diamond blastp -d REPprot -o ${i/.AMP.faa/REPprot.m8} -q $i; done

for i in $(ls from_*AMP.faa); do grep -c ">" $i >> file1; done
for i in $(ls *.m8); do  cut -f1 $i | sort | uniq | wc -l >> file2; done
for i in $(ls *.m8.adam); do  cut -f1 $i | sort | uniq | wc -l >> file3; done

paste -d'\t' file1 file2 > tmp; rm -rf file1 file2
paste -d'\t' tmp file3 > tmp2; rm -rf tmp file3
echo -e "AMPs\tDetect_on_REP\tDetect_on_ADAM" > header
cat header tmp2 > comparison; rm -rf header tmp2


