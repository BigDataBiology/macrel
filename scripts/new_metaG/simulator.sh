for i in $(ls ab); do
awk '{print $1}' ab/$i > tmp
grep -f tmp genomes_genomeInfo.txt > sample_genomeInfo.txt 
val=`perl simulatorForSolexaReads.pl --outfolder=metagenomes_simulated --outprefix=${i/.AbundanceFile.txt/} --abundancefile=ab/$i --genomeInfo=sample_genomeInfo.txt --genomefolder=genomes_separated_files --readlength=75 --insertSize=180 --insertSD=0 --insertNumber=10000000 --qualityfiles="fakeQualFileNoError_75bps.txt"`
echo $val > ${i/.AbundanceFile.txt/.sh}
chmod +x ${i/.AbundanceFile.txt/.sh}
./${i/.AbundanceFile.txt/.sh}
rm -rf ${i/.AbundanceFile.txt/.sh} tmp sample_genomeInfo.txt
done
