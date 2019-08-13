for i in $(ls ../abund/); do awk '{print $1"\t"int($2*100000000)}' ../abund/$i | sort -k1,1 > ab/${i/.abund/.AbundanceFile.txt}; done # number of reads

for i in $(ls ab);
do
	awk '{print $1}' ab/$i > tmp # select list
	grep -f tmp genomes_genomeInfo.txt | sort -k1,1 > sample_genomeInfo.txt # select genomes sizes
	join ab/$i sample_genomeInfo.txt | awk '{print $1"\t"$2"\t"(150*$2)/$3}' > tmp
	while read a b c;
	do
		~/art_bin_MountRainier/art_illumina -ss HS25 -i sep_genomes/$a -p -l 150 -f $c -m 200 -s 10 -o ${a/.fna}
	done < tmp
	cat *1.fq | pigz --best > metagenomes_simulated/${i/.AbundanceFile.txt/_1.fastq.gz}
	cat *2.fq | pigz --best > metagenomes_simulated/${i/.AbundanceFile.txt/_2.fastq.gz}
	rm -rf tmp *.fq *.aln sample_genomeInfo.txt 
done


#val=`perl simulatorForSolexaReads.pl --outfolder=metagenomes_simulated --outprefix=${i/.AbundanceFile.txt/} --abundancefile=ab/$i --genomeInfo=sample_genomeInfo.txt --genomefolder=genomes_separated_files --readlength=75 --insertSize=180 --insertSD=0 --insertNumber=10000000 --qualityfiles="fakeQualFileNoError_75bps.txt"`
#echo $val > ${i/.AbundanceFile.txt/.sh}
#chmod +x ${i/.AbundanceFile.txt/.sh}
#./${i/.AbundanceFile.txt/.sh}
#rm -rf ${i/.AbundanceFile.txt/.sh} tmp sample_genomeInfo.txt
#
