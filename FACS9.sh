#!/usr/bin/env bash

############################################################################################################################
############################## FACS pipeline - Fast AMP Clustering and Screening        ####################################
############################################################################################################################
############################## Authors: Célio Dias Santos Júnior, Luis Pedro Coelho     ####################################
############################################################################################################################
############################## Institute > ISTBI - FUDAN University / Shanghai - China  ####################################
############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

##################################################### Variables ############################################################
Lib="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pigz="/usr/bin/pigz"
pandaseq="/home/celio/pandaseq/pandaseq"

outfolder="./"
block=100000000  #100Mb
j=`echo "$(nproc) * 9 / 10" | bc`
adapters="NexteraPE-PE.fa"
outtag="OUT"

show_help ()
{
	echo "
	############################################################################################################################
	############################## FACS pipeline - Fast AMP Clustering System               ####################################
	############################################################################################################################
	############################## Authors: Célio Dias Santos Júnior, Luis Pedro Coelho     ####################################
	############################################################################################################################
	############################## Institute > ISTBI - FUDAN University / Shanghai - China  ####################################
	############################################################################################################################

	Usage: FACS.sh [options]

	Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].

	Basic options:
	-h, --help	        Show this helpful help page
	--fwd                   Illumina sequencing file in Fastq format (R1), please leave it compressed and full address
	--rev		        Illumina sequencing file in Fastq format (R2), please leave it compressed and full address
	--outfolder		Folder where output will be generated [./]
	--outtag          	Tag used to name outputs [OUT]
	-t, --threads [N]	Number of threads [90% of avaiable threads]
	--block			Bucket size (take in mind it is measured in bits and also it determines the memory usage). [100MB]
	--adapters		Adapters currently available in Trimmomatic program, if not available you can create one accordingly
				to its manual. [NexteraPE-PE.fa]
	--log			Save results of FACS run to a log file in output folder.
	
"
}


while [[ $# -gt 0 ]]
do
	case $1 in
		-h|--help|-help|\?|-\?)
			show_help
			exit
		;;
		-t|-T|--threads|--Threads|--THREADS|--t|--T)
			j=${2}
		;;
		--fwd|--Fwd|--FWD|--F|--R1|--l)
			read_1=${2}
		;;
		--rev|--Rev|--REV|--R2|--r)
			read_2=${2}
		;;
		--outfolder|--Outfolder|--OutFolder|--OUTFOLDER)
			outfolder=${2}
		;;
		--outtag|--Outtag)
			outtag=${2}
		;;
		--Block|--block|--BLOCK)
			block=${2}
		;;
		--adapters|--Adapters|--ADAPTERS)
			adapters=${2}
		;;
		--log|--Log|--LOG)
			log=${2}
		;;
	esac
	shift
done

echo -e "[M :::: Here we specify your variables ]

Threads		$j
read_R1		$read_1
read_R2		$read_2
Adapters	$adapters
Folder		$outfolder
Tag		$outtag
Bucket		$block
Log		$outfolder/$log"

if [[ -n $read_1 || -n $read2 ]]
then
	echo ""
else
	echo "[M ::: Your reads files are not set properly, please review the command line]"
	show_help
	exit
fi

if [[ -n $outfolder ]]
	then
		if [ -d "$outfolder" ] 
			then
    				echo " " 
			else
				echo "Error: Directory $outfolder does not exists."
				show_help
				exit
		fi
	else
		show_help
		exit
fi

############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#####################################################    CMDs   ############################################################

export PATH=$PATH:$Lib

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################# 1. generating two files with mating pairs of reads #######################################
java -jar $Lib/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -phred33 -threads $j \
$read_1 \
$read_2 \
.${read_1/.fastq.gz/.paired.fastq.gz} \
.${read_1/.fastq.gz/.singles.fastq.gz} \
.${read_2/.fastq.gz/.paired.fastq.gz} \
.${read_2/.fastq.gz/.singles.fastq.gz} \
ILLUMINACLIP:$Lib/Trimmomatic-0.39/adapters/$adapters:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:80

if [ -s ".${read_2/.fastq.gz/.paired.fastq.gz}" ]
then
	$pigz -dc .${read_2/.fastq.gz/.paired.fastq.gz} > .${read_2/.fastq.gz/.paired.fastq}
	$pigz -dc .${read_1/.fastq.gz/.paired.fastq.gz} > .${read_1/.fastq.gz/.paired.fastq}
	rm -rf .${read_1/.fastq.gz/.singles.fastq.gz} .${read_2/.fastq.gz/.singles.fastq.gz} .${read_1/.fastq.gz/.paired.fastq.gz} .${read_2/.fastq.gz/.paired.fastq.gz}
else
	echo "[ WARNING ::: Your trimming procedures did not result into a true value ]"
	rm -rf .${read_1/.fastq.gz/.singles.fastq.gz} .${read_2/.fastq.gz/.singles.fastq.gz} .${read_1/.fastq.gz/.paired.fastq.gz} .${read_2/.fastq.gz/.paired.fastq.gz}
	exit
fi


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
########################################## 2. merging files ################################################################

echo "[ M ::: Merging paired end reads ]"
$pandaseq -A pear -f .${read_1/.fastq.gz/.paired.fastq} -r .${read_2/.fastq.gz/.paired.fastq} -T $j -w .read_.assembled.fa 2> .test
if [ -s ".read_.assembled.fa" ]
then
	rm -rf .read_.unassembled.forward.fastq .read_.unassembled.reverse.fastq .read_.dicarded.fastq .${read_1/.fastq.gz/.paired.fastq} .${read_2/.fastq.gz/.paired.fastq} .test
else
	echo "[ WARNING ::: Your merging did not result into a true value, the pandaseq message follows ]"
	cat .test
	rm -rf .read_.assembled.fa .read_.unassembled.forward.fastq .read_.unassembled.reverse.fastq .read_.dicarded.fastq .${read_1/.fastq.gz/.paired.fastq} .${read_2/.fastq.gz/.paired.fastq} .test
	exit
fi


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
############################################## 3. calling orfs #############################################################

echo "[ M ::: Calling ORFs ]"
$Lib/orfm -m 30 .read_.assembled.fa > .pep.faa
if [ -s ".pep.faa" ]
then
	rm -rf .read_.assembled.fa
else
	echo "[ WARNING ::: Your ORFs calling procedure did not result into a true value ]"
	rm -rf .read_.assembled.fa
	exit
fi

# Sorting big fasta
mkdir sorting_folder/
perl -ne '$i++,next if /^>/;print if $i' .pep.faa | parallel --pipe --block $block --recstart "\n" "cat > sorting_folder/small-chunk{#}"
rm -rf .pep.faa
for X in $(ls sorting_folder/small-chunk*); do sort -S 80% --parallel=$j < $X > $X.sorted; rm -rf $X; done
sort -S 80% --parallel=$j -T . -m sorting_folder/small-chunk* > .sorted-huge-file; rm -rf sorting_folder/
perl -i -n -e "print if /S/" .sorted-huge-file

# Eliminating duplicates
LC_ALL=C uniq -c .sorted-huge-file > .out; rm -rf .sorted-huge-file
mkdir $outtag.pep.faa.split/
awk 'FNR==NR{s+=$1;next;} {printf "%s\n%s\n",">pep_"FNR"|"$1*1000000/s,$2}' .out .out | parallel --pipe --block $block --recstart ">" "cat > $outtag.pep.faa.split/small-chunk{#}"
rm -rf .out

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
######################################## 5. calculating descriptors ########################################################

for i in $(ls $outtag.pep.faa.split/small-chunk*); do

echo "[ M :::: Counting distribution using SA scale -- $i ]"
python3 $Lib/CTDDClass.py $i .CTDDC-SA.tsv 'ALFCGIVW' 'RKQEND' 'MSPTHY' #solventaccess
awk '{print $2"\t"$7"\t"$12}' .CTDDC-SA.tsv > .tmp
sed -i '1,1d' .tmp
echo -e "SA.G1.residue0\tSA.G2.residue0\tSA.G3.residue0" > .header
cat .header .tmp > .CTDDC-SA.tsv; rm -rf .tmp .header

echo "[ M :::: Counting distribution using HB scale -- $i ]"
python3 $Lib/CTDDClass.py $i .CTDDC-hb.tsv 'ILVWAMGT' 'FYSQCN' 'PHKEDR' # HEIJNE&BLOMBERG1979
awk '{print $2"\t"$7"\t"$12}' .CTDDC-hb.tsv > .tmp
sed -i '1,1d' .tmp
echo -e "hb.Group.1.residue0\thb.Group.2.residue0\thb.Group.3.residue0" > .header
cat .header .tmp > .CTDDC-hb.tsv; rm -rf .tmp .header

echo "[ M :::: Calculating cheminformatics descriptors -- $i ]"
echo -e "header\tseq\tgroup" > .tmp
perl -pe 's/>(.*)/>\1\t/g; s/\n//g; s/>/\n/g' $i | sed '1,1d' | awk '{print $1"\t"$2"\t""Unk"}' >> .tmp
rm -rf $i
R --vanilla --slave --args .tmp .out.file < $Lib/features_04061950.R >/dev/null 2>/dev/null
rm -rf .tmp


echo "[ M :::: Formatting descriptors -- $i ]"
paste -d'\t' .out.file .CTDDC-SA.tsv .CTDDC-hb.tsv | sed 's/\"//g' > .tab.desc.tsv
rm -rf .out.file .CTDDC-SA.tsv .CTDDC-hb.tsv

echo "[ M :::: Predicting AMPs -- $i ]"
R --vanilla --slave --args .tab.desc.tsv $Lib/r22_largeTraining.rds $Lib/orfsvm_19desc.rds $i.fin < $Lib/AMP_HEMOscreening.R >/dev/null 2>/dev/null
rm -rf $Lib/__pycache__/ .tab.desc.tsv; done

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#################################### 6. Generating abundance profiles ######################################################
echo "[ M :::: Formatting results ]"
cat $outtag.pep.faa.split/* > .out.file
rm -rf $outtag.pep.faa.split/
sed -i 's/\"//g' .out.file
sed -i 's/|/\t/g' .out.file
AMP=`wc -l .out.file | awk '{print $1}'`

echo "[ M :::: Exporting results ]"
echo -e "access\tabundance\tsequence\tgroup\themolytic" > .header
cat .header .out.file | pigz --best > $outfolder/$outtag.tsv.gz; rm -rf .header

echo "[ M :::: Calculating statistics ]"
val=`awk -F '\t' '{a[$4] += $2} END{for (i in a) print i, a[i]}' .out.file`
pepval=`awk '{print $4}' .out.file | sort | uniq -c`
rm -rf .out.file

echo "[ M :::: Generating log ]"
## Preparing report
echo -e "######################################################## FACS report

[ M :::: Variables ]

Threads		$j
read_R1		$read_1
read_R2		$read_2
Adapters	$adapters
Folder		$outfolder
Tag		$outtag
Bucket		$block
Log		$outfolder/$log

========================================= Files were treated

*** Deduplicated ORFs were called

These ORFs were then screened for AMPs

*** $AMP peptides were called

These peptides were classified into 4 families which accounts:

-- Abundance in ppm:
$val

-- Abundance in peptides:
$pepval

** Legend: ADP - Anionic Dissulphide-bond forming peptide
	   ALP - Anionic Linear peptide
	   CDP - Cationic Dissulphide-bond forming peptide
	   CLP - Cationic Linear peptide
########################################################" > $outfolder/$log

############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
