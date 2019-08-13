#!/bin/bash
#zgrep ">" representatives.contigs.fasta.gz | sed 's/>//g' > headers_rep_genome$
###cat headers_rep_genomes.txt | awk -F "." '{print $2}' | sort | uniq > uniq_g$
# To generate uniq_genomes, I have retrieved the genomes from abund files.
#mkdir genomes_separated_files
#zgrep -A 1 -f ids.txt  representatives.contigs.fasta.gz > contigs.test.fna
for i in $(cat uniq_genomes);
do
        grep "$i" headers_rep_genomes.txt > tmp
        echo "Extracting sequence genome $i"
        grep -A 1 -w -f tmp contigs.test.fna > genomes_separated_files/$i.fna
        rm -rf tmp
done

