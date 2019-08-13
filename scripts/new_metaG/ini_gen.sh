#!/usr/bin/bash
## Generate ini files to produce metagenomes
for i in $(ls ../abund);
do
	echo "[general]
output_sample_name = ${i/.abund/}
insert_size = 30
insert_size_std = 1
short_read_length = 100
error_rate = 0.05" > ini_files/${i/.abund/.ini}
	awk '{print $1"\t"$2*1000}' ../abund/$i > $i
	while read a b;
	do
		echo -e "\n" >> ini_files/${i/.abund/.ini}
		echo "[/home/celio/Desktop/benchmark/metagenomes/genomes_separated_files/$a]" >> ini_files/${i/.abund/.ini}
		echo "coverage = $b" >> ini_files/${i/.abund/.ini}		
	done < $i
	rm -rf $i
done

for i in $(ls ../abund/); do awk '{print $1"\t"$2*1000}' ../abund/$i > ${i/.abund/.AbundanceFile.txt}; done
cd genomes_separated_files/
for i in $(ls *); do val=`grep -v ">" $i | wc | awk '{print $3-$1}'`; echo -e "$i\t$val\t1" >> ../sample_genomeInfo.txt; done
cd ../
