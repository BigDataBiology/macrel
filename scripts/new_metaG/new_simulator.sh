cp -r genomes_separated_files sep_gen
cd genomes_separated_files
for i in $(ls *.fna); do sed -i "s/>.*/>${i/.fna/}/g" $i; done
cd ../
for i in $(ls ab/);
do
	cp ab/$i genomes_separated_files/ab
	val=`awk '{print $1}' ab/$i | awk 'BEGIN { ORS = " " } { print }'`
	cd genomes_separated_files/
	iss generate -g $val --abundance_file ab -m HiSeq --n_reads 35000000 --gc_bias -o ../metagenomes_simulated/${i/.AbundanceFile.txt/}
	rm -rf ab
	cd ../
done
rm -rf genomes_separated_files
mv sep_gen genomes_separated_files

for i in $(ls ab/); do cp ab/$i genomes_separated_files/ab; val=`awk '{print $1}' ab/$i | awk 'BEGIN { ORS = " " } { print }'`; cd genomes_separated_files/; ; rm -rf ab; cd ../;  done
