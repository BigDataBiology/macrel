#!/usr/bin/bash

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
model="/home/celio/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod"
mgm="/home/celio/MetaGeneMark_linux_64/mgm/gmhmmp"
mmseqs="/home/celio/MMseqs2/build/bin/mmseqs"
megahit="/home/celio/MEGAHIT-1.2.3-beta-Linux-static/bin/megahit"
pandaseq="/usr/local/bin/pandaseq"
paladin="/home/celio/paladin/paladin"
samtools="/usr/bin/samtools"
eXpress="/home/celio/express-1.5.1-linux_x86_64/express"

# Default variables
outfolder="./"
block=100000000  #100Mb
j=$(echo "$(nproc) * 9 / 10" | bc)
outtag="FACS_OUT"
clust="0"
mem="0.75"

# Help message
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

	Usage: FACS.sh --mode c/r [options]

	Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].

	Basic options:

	-h, --help	        Show this help message

	-m			Mode of operation, type:
				\"c\" to work with contigs,
				\"r\" to work with paired-end reads, 
				\"mr\" to map reads against AMP output database and generate abundances table
	--fasta			Compressed (gzipped) contigs Fasta file
	--fwd                   Illumina sequencing file in Fastq format (R1), please leave it compressed and full address
	--rev		        Illumina sequencing file in Fastq format (R2), please leave it compressed and full address
	--ref                   Compressed (gzipped) reference fasta file with peptides to be mapped to
	--outfolder		Folder where output will be generated [Default: ./]
	--outtag          	Tag used to name outputs [Default: FACS_OUT]
	-t, --threads [N]	Number of threads [Default: 90% of avaiable threads]
	--block			Bucket size (take in mind it is measured in bits and also it determines the memory usage). [100MB]
	--log			Log file name. It allows FACS to save the run results to this log file in output folder.
	--clust			Insert this flag with 1 if you want to pre-cluster peptides with MMSeqs2 at
				95% of id and cov. [Defult: 0]
	--mem			Memory available to FACS ranging from 0 - 1. [Defult: 0.75]	
"
}

# Taking flags
while [[ $# -gt 0 ]]
do
	case $1 in
		-h|--help|-help|\?|-\?)
			show_help
			exit
		;;
		-m|-M|--mode|--Mode|--m|--M)
			mode=${2}
		;;
		-t|-T|--threads|--Threads|--THREADS|--t|--T)
			j=${2}
		;;
		--fasta|--Fasta|--FASTA|-fasta)
			fasta=${2}
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
		--log|--Log|--LOG)
			log=${2}
		;;
		-ref|--ref|--Ref|-Ref)
			Reference=${2}
		;;
		-clust|--Clust|--clust|-Clust)
			clust=${2}
		;;
		-mem|--Mem|--mem|-Mem)
			mem=${2}
		;;
	esac
	shift
done

sanity_check ()
# Checking input variables
{
if [[ -n $mode ]]
then
	if [[ $mode == "r" ]]
	then
		if [[ -n $read_1 ]]
		then
			if [[ -n $read_2 ]]
			then
				if [[ -s $read_1 ]] || [[ -s $read_2 ]]
				then 
					echo "[ M ::: FACS mode has been assigned as paired-end reads ]
[ M ::: FACS has found your reads files, starting work... ]"
					mode="pe"
					echo -e "[ M ::: Here we specify your variables ]

Mode		$mode
Threads		$j
read_R1		$read_1
read_R2		$read_2
Folder		$outfolder
Tag		$outtag
Bucket		$block
Log		$outfolder/$log"
				else
					echo "[ M ::: FACS mode has been assigned as Reads ]
[ W ::: ERR010 - FACS has not found your reads files. ]"
					show_help
					exit
				fi
			elif [[ -s $read_1 ]]
			then 
				echo "[ M ::: FACS mode has been assigned as single-end reads ]
[ M ::: FACS has found your reads files, starting work... ]"
				mode="se"
				echo -e "[ M ::: Here we specify your variables ]

Mode		$mode
Threads		$j
read_R1		$read_1
Folder		$outfolder
Tag		$outtag
Bucket		$block
Log		$outfolder/$log"
			else
				echo "[ M ::: FACS mode has been assigned as Reads ]
[ W ::: ERR010 - FACS has not found your reads files. ]"
				show_help
				exit
			fi
		else
			echo "[ M ::: FACS mode has been assigned as Reads ]
[ W ::: ERR010 - FACS has not found your reads files. ]"
			show_help
			exit
		fi
	elif	[[ $mode == "c" ]]
	then
		if [ -s "$fasta" ]
		then 
			echo "[ M ::: FACS mode has been assigned as Contigs ]
[ M ::: FACS has found your contigs file, starting work... ]"
			echo -e "[ M ::: Here we specify your variables ]

Mode			$mode
Threads			$j
Contigs			$fasta
Folder			$outfolder
Tag			$outtag
Bucket			$block
Log			$outfolder/$log"
		else
			echo "[ M ::: FACS mode has been assigned as Contigs ]
[ W ::: ERR011 - Your contigs file is not present, please review the command line ]"
			show_help
			exit
		fi
	elif [[ $mode == "mr" ]]
	then
		echo "[ M ::: FACS mode has been assigned as read mapper ]"
		if [ -s "$Reference" ]
		then
			if [[ -s $read_1 ]]
			then
				if [[ -s $read_2 ]]
				then
					mode="mpe"
					echo "[ M ::: FACS has found your reference data set, starting work... ]"
					echo -e "[ M ::: Here we specify your variables ]

** Mapper with paired-end reads
Mode			$mode
Threads			$j
Reference		$Reference
R1			$read_1
R2			$read_2
Folder			$outfolder
Tag			$outtag
Bucket			$block
Log			$outfolder/$log"
				else
					mode="mse"
					echo "[ M ::: FACS has found your reference data set, starting work... ]"
					echo -e "[ M ::: Here we specify your variables ]

** Mapper with single-end reads
Mode			$mode
Threads			$j
Reference		$Reference
R1			$read_1
Folder			$outfolder
Tag			$outtag
Bucket			$block
Log			$outfolder/$log"
				fi
			else
				echo "[ M ::: FACS mode has been assigned as mapper ]
[ W ::: ERR011 - FACS has not found your reads files ]"
				show_help
				exit
			fi
	fi
else
	echo "[ W ::: ERR010 - The user needs to specify a valid FACS mode, please review the command line]"
	show_help
	exit
fi
fi

if [[ -n $outfolder ]]
	then
		if [ -d "$outfolder" ] 
			then
    				echo "" 
			else
				echo "[ W ::: ERR017 - Directory $outfolder does not exist. ]"
				show_help
				exit
		fi
	else
		show_help
		exit
	fi
}

############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################################       ENV      ############################################################

export PATH=$PATH:$Lib

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################################    FUNCTIONS   ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################# 1. generating two files with matching pairs of reads #####################################
PEreads_trimming ()
{
echo "[ M ::: Trimming low quality bases ]"
echo "[ M ::: Sorting paired-end reads ]"
java -jar "$Lib"/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -phred33 -threads "$j" \
"$read_1" \
"$read_2" \
.read_1.paired.fastq.gz \
.read_1.singles.fastq.gz \
.read_2.paired.fastq.gz \
.read_2.singles.fastq.gz \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:80 >/dev/null 2>/dev/null

if [ -s ".read_2.paired.fastq.gz" ]
then
	rm -rf ."read_1.singles.fastq.gz" ."read_2.singles.fastq.gz"
else
	echo "[ W ::: ERR231 - Your trimming procedures did not result into a true value ]"
	rm -rf .read_1.singles.fastq.gz .read_2.singles.fastq.gz .read_1.paired.fastq.gz .read_2.paired.fastq.gz
	exit
fi
}

SEreads_trimming ()
{
echo "[ M ::: Trimming low quality bases ]"
java -jar "$Lib"/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads "$j" "$read_1" .read_1.paired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:80 >/dev/null 2>/dev/null

if [ -s ."read_1.paired.fastq.gz" ]
then
	touch ."read_1.paired.fastq.gz"
else
	echo "[ W ::: ERR231 - Your trimming procedures did not result into a true value ]"
	rm -rf .read_1.paired.fastq.gz
	exit
fi
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
########################################## 2. Generating peptides ##########################################################

ASSEMBLY ()
{
echo "[ M ::: Assembly using MEGAHIT ]"
if [[ $mode == "pe" ]]
then
	$megahit --presets meta-large -1 .read_1.paired.fastq.gz -2 .read_2.paired.fastq.gz --min-contig-len 1000 -o out -t "$j" -m "$mem"
elif [[ $mode == "se" ]]
then
	$megahit --presets meta-large -r .read_1.paired.fastq.gz --min-contig-len 1000 -o out -t "$j" -m "$mem"
else 
	echo "[ W ::: ERR222 - FACS followed by a weird way ]"
	exit
fi

if [[ -s out/final.contigs.fa ]]
then
	mv out/final.contigs.fa .callorfinput.fa
	rm -rf out/ .read_1.paired.fastq.gz .read_2.paired.fastq.gz 
else
	echo "[ W ::: ERR128 - Assembly returned ECC0 ]"
	rm -rf out/ .read_1.paired.fastq.gz .read_2.paired.fastq.gz 
	exit
fi
}

callorf ()
{
echo "[ M ::: Calling ORFs ]"
$mgm -A .out.faa -m "$model" -n -i 0 .callorfinput.fa
rm -rf .callorfinput.fa .callorfinput.fa.lst
sed -i 's/\t>.*//g' .out.faa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' .out.faa | awk -F "\t" '$2 ~ /^ *M/ || $2 ~ /^ *V/ || $2 ~ /^ *L/ || $2 ~ /^ *S/' | awk '{print $1"\n"$2}' | awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])<=100) {printf("%s\n%s\n",L[0],L[1]);}}' > .pep.faa

if [ -s ".pep.faa" ]
then
	rm -rf .out.faa
else
	echo "[ W ::: ERR122 - Your ORFs calling procedure did not result into a true value ]"
	rm -rf .pep.faa .out.faa
	exit
fi
}

memdev ()
{
echo "[ M ::: Sorting peptides list ]"
mkdir sorting_folder/
perl -ne '$i++,next if /^>/;print if $i' .pep.faa | parallel --pipe --block "$block" --recstart "\n" "cat > sorting_folder/small-chunk{#}"
rm -rf .pep.faa
for X in $(ls sorting_folder/small-chunk*); do sed -i '/x/d' "$X"; sort -S 80% --parallel="$j" < "$X" > "$X".sorted; rm -rf "$X"; done
sort -S 80% --parallel="$j" -T . -m sorting_folder/small-chunk* > .sorted-huge-file; rm -rf sorting_folder/
perl -i -n -e "print if /S/" .sorted-huge-file
if [[ -s .sorted-huge-file ]]
then
	echo "[ M ::: Eliminating duplicated sequences ]"
	LC_ALL=C uniq .sorted-huge-file | awk 'FNR==NR {print ">pep_"FNR"\n"$1}' > .pep.faa
	if [[ -s .pep.faa ]]
	then
		rm -rf .sorted-huge-file
	else
		echo "[ W ::: ERR566 - Deduplication has failed ]"
		exit
	fi
else
	echo "[ W ::: ERR510 - Memdev sorting stage has failed ]"
	exit
fi	
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
########################################## 3. Clustering peptides ##########################################################

clusteRC ()
{
if [ "$clust" == "0" ]
then
	echo "[ M ::: Skipping clustering using MMSeqs2 ]"

elif [ "$clust" == "1" ]
then 

	echo "[ M ::: Indexing database ]"
	$mmseqs createdb .pep.faa .DB > log.clus
	mv .pep.faa "$outfolder"/"$outtag".pep.faa
	$pigz --best "$outfolder"/"$outtag".pep.faa

	echo "[ M ::: Executing clustering ]"
	$mmseqs cluster --cov-mode 0 -c 0.95 --alignment-mode 3 --min-seq-id 0.95 .DB .clu tmp >> log.clus

	echo "[ M ::: Output file from ffindex ]"
	$mmseqs createseqfiledb .DB .clu .clu_seq >> log.clus
	$mmseqs result2flat .DB .DB .clu_seq .clu_seq.fasta >> log.clus

	echo "[ M ::: Generating TSV report from ffindex ]"
	$mmseqs createtsv .DB .DB .clu .clu.tsv >> log.clus

	echo "[ M ::: Extracting representative sequences ]"
	$mmseqs result2repseq .DB .clu .DB_clu_rep >> log.clus
	$mmseqs result2flat .DB .DB .DB_clu_rep .DB_clu_rep.fasta --use-fasta-header >> log.clus

	echo "[ M ::: Cleaning ]"

	if [ -s "$outfolder"/"$outtag".pep.faa.gz ]
	then
		mv .clu.tsv "$outfolder"/"$outtag".clusters.tsv
		mv .DB_clu_rep.fasta .pep.faa
		rm -rf .clu.* .clu_seq.* .DB.* .DB_clu_rep.dbtype .DB_clu_rep .DB_clu_rep.index .*.dbtype .*.index .*_h .*.lookup .DB_clu_rep .DB tmp/ log.clus
		$pigz --best "$outfolder"/"$outtag".clusters.tsv
	else
		echo "[ W ::: ERR125 - Clustering failed ]"
		rm -rf .clu.* .clu_seq.* .DB.* .DB_clu_rep.dbtype .DB_clu_rep .DB_clu_rep.index .*.dbtype .*.index .*_h .*.lookup .DB_clu_rep .DB tmp/
		cat log.clus
		rm -rf log.clus
		exit
	fi
else
	echo "[ W ::: ERR621 - Failed in interpret variable ]"
	exit
fi
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
######################################## 5. calculating descriptors ########################################################

descripter ()
{
echo "[ M ::: Reducing file sizes ]"
mkdir splitted/
parallel --pipe --block "$block" --recstart ">" "cat > splitted/small-chunk{#}" < .pep.faa
rm -rf .pep.faa

for i in splitted/small-chunk*; do

echo "[ M ::: Counting distribution using SA scale -- $i ]"
python3 "$Lib"/CTDDClass.py "$i" .CTDDC-SA.tsv 'ALFCGIVW' 'RKQEND' 'MSPTHY' #solventaccess
awk '{print $2"\t"$7"\t"$12}' .CTDDC-SA.tsv > .tmp
sed -i '1,1d' .tmp
echo -e "SA.G1.residue0\tSA.G2.residue0\tSA.G3.residue0" > .header
cat .header .tmp > .CTDDC-SA.tsv; rm -rf .tmp .header

echo "[ M ::: Counting distribution using HB scale -- $i ]"
python3 "$Lib"/CTDDClass.py "$i" .CTDDC-hb.tsv 'ILVWAMGT' 'FYSQCN' 'PHKEDR' # HEIJNE&BLOMBERG1979
awk '{print $2"\t"$7"\t"$12}' .CTDDC-hb.tsv > .tmp
sed -i '1,1d' .tmp
echo -e "hb.Group.1.residue0\thb.Group.2.residue0\thb.Group.3.residue0" > .header
cat .header .tmp > .CTDDC-hb.tsv; rm -rf .tmp .header

echo "[ M ::: Calculating cheminformatics descriptors -- $i ]"
echo -e "header\tseq\tgroup" > .tmp
sed '/>/d' "$i" > .seqs; grep '>' "$i" | sed 's/ .*//g' | sed 's/>//g' > .heade; paste -d'\t' .heade .seqs | awk '{print $1"\t"$2"\t""Unk"}' >> .tmp; rm -rf .heade .seqs
rm -rf "$i"
R --vanilla --slave --args .tmp .out.file < "$Lib"/features_04061950.R >/dev/null 2>/dev/null
rm -rf .tmp

echo "[ M ::: Formatting descriptors -- $i ]"
paste -d'\t' .out.file .CTDDC-SA.tsv .CTDDC-hb.tsv | sed 's/\"//g' > .tab.desc.tsv
rm -rf .out.file .CTDDC-SA.tsv .CTDDC-hb.tsv

echo "[ M ::: Predicting AMPs -- $i ]"
R --vanilla --slave --args .tab.desc.tsv "$Lib"/r22_largeTraining.rds "$Lib"/orfsvm_19desc.rds "$i".fin < "$Lib"/AMP_HEMOscreening02.R >/dev/null 2>/dev/null
rm -rf "$Lib"/__pycache__/
rm -rf .tab.desc.tsv; done
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#################################### 6. Generating abundance profiles ######################################################

cleanmessy ()
{
echo "[ M ::: Formatting results ]"
echo -e "access\tsequence\tgroup\themolytic" > .out.file
cat splitted/* >> .out.file
rm -rf splitted/
sed -i 's/\"//g' .out.file
AMP=$(sed '1,1d' .out.file | wc -l | awk '{print $1}')
echo "[ M ::: Calculating statistics ]"
pepval=$(awk '{print $3}' .out.file | sort | uniq -c | awk '{ print $2"\t"$1 }' | grep -v "group")
echo "[ M ::: Exporting results ]"
mv .out.file "$outfolder"/"$outtag".tsv
$pigz --best "$outfolder"/"$outtag".tsv
}

loggen()
{
echo "[ M ::: Generating log ]"

if [[ $mode == "r" ]]
then
	## Preparing report
	echo -e "######################################################## FACS report

[ M ::: Variables ]

Threads		$j
read_R1		$read_1
read_R2		$read_2
Folder		$outfolder
Tag		$outtag
Bucket		$block

========================================= Files were treated

*** Deduplicated ORFs were called

These ORFs were then screened for AMPs

*** $AMP peptides were called

These peptides were classified into 4 families which accounts:

-- Distribution of peptides:
$pepval

** Legend: ADP - Anionic Dissulphide-bond forming peptide
	   ALP - Anionic Linear peptide
	   CDP - Cationic Dissulphide-bond forming peptide
	   CLP - Cationic Linear peptide
########################################################" > "$outfolder"/.log
else
	## Preparing report
	echo -e "######################################################## FACS report

[ M ::: Variables ]

Threads		$j
Contigs		$fasta
Folder		$outfolder
Tag		$outtag
Bucket		$block

========================================= Files were treated

*** Deduplicated ORFs were called

These ORFs were then screened for AMPs

*** $AMP peptides were called

These peptides were classified into 4 families which accounts:

-- Distribution of peptides:
$pepval

** Legend: ADP - Anionic Dissulphide-bond forming peptide
	   ALP - Anionic Linear peptide
	   CDP - Cationic Dissulphide-bond forming peptide
	   CLP - Cationic Linear peptide
########################################################" > "$outfolder"/.log
fi

sed -i 's/\"//g' "$outfolder"/.log

if [[ -n $log ]]
then
	cat "$outfolder"/.log
	mv "$outfolder"/.log "$outfolder"/"$log"
else
	cat "$outfolder"/.log
	rm -rf "$outfolder"/.log
fi

}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###########################################      7. Abundances      ########################################################

mapping ()
{
echo "[ M ::: Indexing references ]"
$pigz -dc "$Reference" | awk '{ print ">"$1"\n"$2 }' | sed '1,2d' > .ref.fa
$paladin index -r 3 .ref.fa
if [ -s .ref.fa.amb ]
then
	if [[ $mode == "mse" ]]
	then
		touch .ref.fa
	elif [[ $mode == "mpe" ]]
	then
		echo "[ M ::: Merging paired end reads ]"
		$pigz -d .read_1.paired.fastq.gz
		$pigz -d .read_2.paired.fastq.gz
		$pandaseq -A pear -F -f .read_1.paired.fastq -r .read_2.paired.fastq -T "$j" -w .read_.assembled.fastq 2> .test
		if [ -s ".read_.assembled.fastq" ]
		then
			rm -rf .read_.unassembled.forward.fastq .read_.unassembled.reverse.fastq .read_.discarded.fastq .read_1.paired.fastq .read_2.paired.fastq .test
			mv .read_.assembled.fastq .read_1.paired.fastq
			$pigz --best .read_1.paired.fastq
		else
			echo "[ W ::: ERR301 - Your merging did not result into a true value, the pandaseq message follows ]"
			cat .test
			rm -rf .read_.assembled.fastq .read_.unassembled.forward.fastq .read_.unassembled.reverse.fastq .read_.discarded.fastq .read_1.paired.fastq .read_2.paired.fastq .test
			exit
		fi
	else
		echo "[ W ::: ERR33 - Please review command line // INTERNAL ERROR SANITARY PROCESS ]"
		exit
	fi
else
	rm -rf .ref.fa
	echo "[ W ::: ERR303 - Error in indexing ]"	
	exit
fi

echo "[ M ::: Mapping reads against references, be aware it can take a while ]"

echo "[ M ::: Starting the paladin ]"
$paladin align -t "$j" -T 20 -f 10 -z 11 -a -V -M .ref.fa .read_1.paired.fastq.gz | $samtools view -Sb > .m.bam

if [[ -s .m.bam ]]
then
	touch .m.bam
else
	echo "[ W ::: ERR052 - Mapping failed ]"
	.re* .m.bam
	exit
fi
}


ab_profiling ()
{
echo "[ M ::: Expressing results of abundance ]"
$eXpress --no-bias-correct -o "$outfolder"/ .ref.fa .m.bam
if [[ -s "$outfolder"/results.xprs ]] || [[ -s "$outfolder"/params.xprs ]] 
then
	mv "$outfolder"/results.xprs "$outfolder"/"$outtag".xprs
	mv "$outfolder"/params.xprs "$outfolder"/"$outtag".params.xprs
else
	echo "[ W ::: ERR054 - Abundance calling failed ]"
fi
rm -rf .ref.* .read_1.paired.fastq.* .m*
}

############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################################    CMDs   #################################################################

export PATH=$PATH:$Lib

sanity_check

if [[ $mode == "pe" ]]
then
	PEreads_trimming
	ASSEMBLY
	callorf
	memdev
	clusteRC
	descripter
	cleanmessy
	mode="r"
	loggen
elif [[ $mode == "se" ]]
then
	SEreads_trimming
	ASSEMBLY
	callorf
	memdev
	clusteRC
	descripter
	cleanmessy
	mode="r"
	loggen
elif [[ $mode == "c" ]]
then
	echo "[ M ::: Decompressing contigs ]"
	$pigz -dc "$fasta" > .callorfinput.fa
	callorf
	memdev
	clusteRC
	descripter
	cleanmessy
	loggen
elif [[ $mode == "mse" ]]
then
	SEreads_trimming
	mapping
	ab_profiling
elif [[ $mode == "mpe" ]]
then
	PEreads_trimming
	mapping
	ab_profiling
else
	echo "[ W ::: ERR010 - The user needs to specify a valid FACS mode, please review the command line]"
	show_help
	exit
fi
