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
	--fasta			Compressed (or not gzipped) contigs Fasta file
	--fwd                   Illumina sequencing file in Fastq format (R1), please leave it compressed and full address
	--rev		        Illumina sequencing file in Fastq format (R2), please leave it compressed and full address
	--ref                   Output of module \"c\" in its raw format [file type tsv and compressed ]
	--outfolder		Folder where output will be generated [Default: ./]
	--outtag          	Tag used to name outputs [Default: FACS_OUT]
	-t, --threads [N]	Number of threads [Default: 90% of avaiable threads]
	--block			Bucket size (take in mind it is measured in bits and also it determines the memory usage). [100MB]
	--log			Log file name. It allows FACS to save the run results to this log file in output folder.
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
		--fasta|-fasta|-fa)
			fasta=${2}
		;;
		--fwd|-fwd|--FWD|-r1|-R1|-l)
			read_1=${2}
		;;
		--rev|--Rev|-rev|-R2|-r2)
			read_2=${2}
		;;
		--outfolder|-outfolder|-of)
			outfolder=${2}
		;;
		--outtag|-tag)
			outtag=${2}
		;;
		--block|-b|-block)
			block=${2}
		;;
		--log|-log)
			log=${2}
		;;
		-ref|--ref|--Ref|-Ref)
			Reference=${2}
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
trimmomatic PE -phred33 -threads "$j" \
"$read_1" \
"$read_2" \
.read_1.paired.fastq.gz \
.read_1.singles.fastq.gz \
.read_2.paired.fastq.gz \
.read_2.singles.fastq.gz \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:75 >/dev/null 2>/dev/null

if [ -s ".read_2.paired.fastq.gz" ]
then
	rm -rf ".read_1.singles.fastq.gz" ".read_2.singles.fastq.gz"
else
	echo "[ W ::: ERR231 - Your trimming procedures did not result into a true value ]"
	rm -rf .read_1.singles.fastq.gz .read_2.singles.fastq.gz .read_1.paired.fastq.gz .read_2.paired.fastq.gz feat.R pred.R
	exit
fi
}

SEreads_trimming ()
{
echo "[ M ::: Trimming low quality bases ]"
trimmomatic SE -phred33 -threads "$j" "$read_1" .read_1.paired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 >/dev/null 2>/dev/null

if [ -s ".read_1.paired.fastq.gz" ]
then
	touch ."read_1.paired.fastq.gz"
else
	echo "[ W ::: ERR231 - Your trimming procedures did not result into a true value ]"
	rm -rf .read_1.paired.fastq.gz feat.R pred.R
	exit
fi
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
########################################## 2. Assembly of contigs ##########################################################

ASSEMBLY ()
{
echo "[ M ::: Assembly using MEGAHIT ]"
if [[ $mode == "pe" ]]
then
	megahit --presets meta-large -1 .read_1.paired.fastq.gz -2 .read_2.paired.fastq.gz -o out -t "$j" -m "$mem" --min-contig-len 1000
elif [[ $mode == "se" ]]
then
	megahit --presets meta-large -r .read_1.paired.fastq.gz -o out -t "$j" -m "$mem" --min-contig-len 1000
else 
	echo "[ W ::: ERR222 - FACS followed by a weird way ]"
	rm -rf feat.R pred.R out/
	exit
fi

if [[ -s out/final.contigs.fa ]]
then
	mv out/final.contigs.fa .callorfinput.fa
	rm -rf out/ .read_1.paired.fastq.gz .read_2.paired.fastq.gz 
else
	echo "[ W ::: ERR128 - Assembly returned ECC0 ]"
	rm -rf out/ .read_1.paired.fastq.gz .read_2.paired.fastq.gz feat.R pred.R
	exit
fi
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
########################################## 2. Generating peptides ##########################################################

callorf ()
{
echo "[ M ::: Calling ORFs ]"

rm -rf callorfs/
mkdir callorfs

cat .callorfinput.fa | parallel -j $j --block $block --recstart '>' --pipe $Lib/envs/FACS_env/bin/prodigal -c -m -n -p meta -f sco -a callorfs/{#}.pred.smORFs.fa >/dev/null 2>/dev/null

ls callorfs/*pred.smORFs.fa > t
if [ -s t ]
then
	rm -rf t
else
	echo "[ W ::: ERR910 - Error in producing predictions ]"
	rm -rf t callorfs/ .callorfinput.fa feat.R pred.R
	exit
fi

if [[ $mode == "c" ]]
then
	rm -rf .callorfinput.fa
else
	mv .callorfinput.fa $outfolder/$outtag.contigs.fna
	pigz --best $outfolder/$outtag.contigs.fna
fi

echo "[ M ::: Filtering smORFs ]" 

for i in callorfs/*;
do
	awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $i | awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])<=100) {printf("%s\n%s\n",L[0],L[1]);}}' > $i.filtered
	rm -rf $i
done

cat callorfs/* > .pep.faa.tmp
rm -rf callorfs/
rm -rf .callorfinput.fa

if [ -s ".pep.faa.tmp" ]
then

	echo "[ M ::: Performing reduction in sampling space ]"

	sed 's/ /_/g' .pep.faa.tmp | sed 's/;/|/g' | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | awk '{print $2"\t"$1}' > .tmp; rm -rf .pep.faa.tmp
	
	echo "[ M ::: Sorting peptides list ]"
	mkdir sorting_folder/
	
	cat .tmp | parallel --pipe  -j $j --block "$block" --recstart "\n" "cat > sorting_folder/small-chunk{#}" >/dev/null 2>/dev/null
	rm -rf .tmp

	for X in $(ls sorting_folder/small-chunk*); do sort -S 80% --parallel="$j" -k1,1 -T . < "$X" > "$X".sorted; rm -rf "$X"; done

	sort -S 80% --parallel="$j" -T . -k1,1 -m sorting_folder/small-chunk* > .sorted-huge-file; rm -rf sorting_folder/

	perl -n -e "print unless /^$/" .sorted-huge-file | awk -F'\t' -v OFS=';' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' | sed 's/;;/\t/g' | sed 's/>//g' | sed 's/\*//g' > .tmp2

	if [[ -s .tmp2 ]]
	then
		echo "[ M ::: Eliminated all duplicated sequences ]"
		rm -rf .sorted-huge-file
	else
		echo "[ W ::: ERR510 - Memdev sorting stage has failed ]"
		rm -rf .sorted-huge-file feat.R pred.R
		exit
	fi	

	echo "[ M ::: Outputting genes lists ]"

	awk '{print ">smORF_"NR"\t"$1"\t"$2}' .tmp2 | pigz --best > $outfolder/$outtag.ids.tsv.gz

	echo "[ M ::: Outputting fastas ]"

	awk '{print ">smORF_"NR"\n"$1}' .tmp2 > .pep.faa 
	rm -rf .tmp2
	
	echo "[ M ::: Collecting statistics ]"

	tail -2 .pep.faa | head -1 | sed 's/|/\t/g' | cut -f2 | sed 's/>smORF_//g' > .all.nmb
else
	echo "[ W ::: ERR122 - Your ORFs calling procedure did not result into a true value ]"
	rm -rf .pep.faa.tmp feat.R pred.R
	exit
fi


}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
######################################## 5. calculating descriptors ########################################################

descripter ()
{
echo "[ M ::: Reducing file sizes ]"
mkdir splitted/
parallel --pipe  -j $j --block "$block" --recstart ">" "cat > splitted/small-chunk{#}" < .pep.faa >/dev/null 2>/dev/null
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

	if [ -s .CTDDC-SA.tsv ]
	then
		if [ -s .CTDDC-hb.tsv ]
		then

			echo "[ M ::: Calculating cheminformatics descriptors -- $i ]"
			echo -e "header\tseq\tgroup" > .tmp
			sed '/>/d' "$i" > .seqs; grep '>' "$i" | sed 's/ .*//g' | sed 's/>//g' > .heade; paste -d'\t' .heade .seqs | awk '{print $1"\t"$2"\t""Unk"}' >> .tmp; rm -rf .heade .seqs
			rm -rf "$i"
			$Lib/envs/FACS_env/bin/R --vanilla --slave --args .tmp .out.file < feat.R #>/dev/null 2>/dev/null
			rm -rf .tmp
		
			echo "[ M ::: Formatting descriptors -- $i ]"
			if [ -s .out.file ]
			then
				paste -d'\t' .out.file .CTDDC-SA.tsv .CTDDC-hb.tsv | sed 's/\"//g' > .tab.desc.tsv
				rm -rf .out.file .CTDDC-SA.tsv .CTDDC-hb.tsv
		
			        echo "[ M ::: Predicting AMPs -- $i ]"
		        	$Lib/envs/FACS_env/bin/R --vanilla --slave --args .tab.desc.tsv "$Lib"/r22_largeTraining.rds "$Lib"/orfsvm_19desc.rds "$i".fin < pred.R #>/dev/null 2>/dev/null
		        	rm -rf "$Lib"/__pycache__/
				rm -rf .tab.desc.tsv
			else
				rm -rf .CTDDC-SA.tsv .CTDDC-hb.tsv $i
				echo "[ W ::: Skipping AMPs prediction for $i -- ERR 229 ]"
			fi

		else

	                rm -rf $i
        	        echo "[ W ::: Skipping AMPs prediction for $i -- ERR 230 ]"

		fi	

	else

		rm -rf $i
		echo "[ W ::: Skipping AMPs prediction for $i -- ERR 230 ]"

	fi
done
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#################################### 6. Generating abundance profiles ######################################################

cleanmessy ()
{
echo "[ M ::: Formatting results ]"
tot=$(cat .all.nmb)
rm -rf .all.nmb

if [ "$(ls -A splitted/)" ]; then
	echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > .out.file
	cat splitted/* >> .out.file
	rm -rf splitted/ 
	sed -i 's/\"//g' .out.file
	AMP=$(sed '1,1d' .out.file | wc -l | awk '{print $1}')
	echo "[ M ::: Calculating statistics ]"
	pepval=$(awk '{print $3"_"$5}' .out.file | sort | uniq -c | awk '{ print $2"\t"$1 }' | grep -v "AMP_family")
	echo "[ M ::: Exporting results ]"
	mv .out.file "$outfolder"/"$outtag".tsv
	pigz --best "$outfolder"/"$outtag".tsv
else
    	echo "[ W ::: ERR585 - We are sorry to inform. Your procedure did not return any AMP sequence. We are closing now ]"
        echo "[ W ::: We really tried. None of $tot called smORFs were assigned to a probability higher than 0.5 or survived until here ]"
	rm -rf splitted/ feat.R pred.R
	exit
fi
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

A total of $tot of smORFs were called.

These ORFs were then screened for AMPs.

*** $AMP peptides were called as potential AMPs

These peptides were classified into 4 families which accounts:

-- Distribution of peptides:

$pepval

** Legend: ADP     - Anionic Dissulphide-bond forming peptide
	   ALP     - Anionic Linear peptide
	   CDP     - Cationic Dissulphide-bond forming peptide
	   CLP     - Cationic Linear peptide
	   HEMO    - Hemolytic
	   NonHEMO - Non-Hemolytic
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

A total of $tot of smORFs were called.

These ORFs were then screened for AMPs.

*** $AMP peptides were called as potential AMPs

These peptides were classified into 4 families which accounts:

-- Distribution of peptides:

$pepval

** Legend: ADP     - Anionic Dissulphide-bond forming peptide
	   ALP     - Anionic Linear peptide
	   CDP     - Cationic Dissulphide-bond forming peptide
	   CLP     - Cationic Linear peptide
	   HEMO    - Hemolytic
	   NonHEMO - Non-Hemolytic
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
pigz -dc "$Reference" | awk '{ print ">"$1"\n"$2 }' | sed '1,2d' > .ref.fa
paladin index -r 3 .ref.fa

if [ -s .ref.fa.amb ]
then
	if [[ $mode == "mse" ]]
	then
		touch .ref.fa
	elif [[ $mode == "mpe" ]]
	then
		echo "[ M ::: Merging paired end reads ]"
		pigz -d .read_1.paired.fastq.gz
		pigz -d .read_2.paired.fastq.gz
		pandaseq -A pear -F -f .read_1.paired.fastq -r .read_2.paired.fastq -T "$j" -w .read_.assembled.fastq 2> .test
		if [ -s ".read_.assembled.fastq" ]
		then
			rm -rf .read_.unassembled.forward.fastq .read_.unassembled.reverse.fastq .read_.discarded.fastq .read_1.paired.fastq .read_2.paired.fastq .test
			mv .read_.assembled.fastq .read_1.paired.fastq
			pigz --best .read_1.paired.fastq
		else
			echo "[ W ::: ERR301 - Your merging did not result into a true value, the pandaseq message follows ]"
			cat .test
			rm -rf .read_.assembled.fastq .read_.unassembled.forward.fastq .read_.unassembled.reverse.fastq .read_.discarded.fastq .read_1.paired.fastq .read_2.paired.fastq .test feat.R pred.R
			exit
		fi
	else
		echo "[ W ::: ERR33 - Please review command line // INTERNAL ERROR SANITARY PROCESS ]"
		rm -rf feat.R pred.R
		exit
	fi
else
	rm -rf .ref.fa*
	echo "[ W ::: ERR303 - Error in indexing ]"	
	rm -rf feat.R pred.R
	exit
fi

echo "[ M ::: Mapping reads against references, be aware it can take a while ]"

echo "[ M ::: Starting the paladin ]"
paladin align -t "$j" -T 20 -f 10 -z 11 -a -V -M .ref.fa .read_1.paired.fastq.gz | samtools view -Sb > .m.bam

if [[ -s .m.bam ]]
then
	touch .m.bam
else
	echo "[ W ::: ERR052 - Mapping failed ]"
	rm -rf .re* .m.bam feat.R pred.R
	exit
fi
}


ab_profiling ()
{
echo "[ M ::: Expressing results of abundance ]"
eXpress --no-bias-correct -o "$outfolder"/ .ref.fa .m.bam
if [[ -s "$outfolder"/results.xprs ]] || [[ -s "$outfolder"/params.xprs ]] 
then
	mv "$outfolder"/results.xprs "$outfolder"/"$outtag".xprs
	mv "$outfolder"/params.xprs "$outfolder"/"$outtag".params.xprs
else
	echo "[ W ::: ERR054 - Abundance calling failed ]"
fi
rm -rf .ref.* .read_1.paired.fastq.* .m* feat.R pred.R
}

############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################################    CMDs   #################################################################

export PATH=$PATH:$Lib/envs/FACS_env/bin:$Lib/envs/FACS_env/conda-meta:$Lib/envs/FACS_env/etc:$Lib/envs/FACS_env/jmods:$Lib/envs/FACS_env/lib:$Lib/envs/FACS_env/libexec:$Lib/envs/FACS_env/mkspecs:$Lib/envs/FACS_env/plugins:$Lib/envs/FACS_env/resources:$Lib/envs/FACS_env/share:$Lib/envs/FACS_env/translations:$Lib/envs/FACS_env/x86_64-conda_cos6-linux-gnu:$Lib/envs/FACS_env/compiler_compat:$Lib/envs/FACS_env/conf:$Lib/envs/FACS_env/doc:$Lib/envs/FACS_env/include:$Lib/envs/FACS_env/legal:$Lib/envs/FACS_env/lib64:$Lib/envs/FACS_env/man:$Lib/envs/FACS_env/phrasebooks:$Lib/envs/FACS_env/qml:$Lib/envs/FACS_env/sbin:$Lib/envs/FACS_env/ssl:$Lib/envs/FACS_env/var:$Lib:~/miniconda3/pkgs/:$Lib/envs/FACS_env/lib/R/library/

sanity_check

sed "s|PEPPERIDY|$Lib/envs/FACS_env/lib/R/library/|g" $Lib/features_130819.R > feat.R
sed "s|PEPPERIDY|$Lib/envs/FACS_env/lib/R/library/|g" $Lib/Predict_130819.R > pred.R
chmod +x feat.R
chmod +x pred.R

date

if [[ $mode == "pe" ]]
then
	PEreads_trimming
	ASSEMBLY
	callorf
	descripter
	cleanmessy
	mode="r"
	loggen
elif [[ $mode == "se" ]]
then
	SEreads_trimming
	ASSEMBLY
	callorf
	descripter
	cleanmessy
	mode="r"
	loggen
elif [[ $mode == "c" ]]
then
	if [[ "$fasta" =~ \.gz$ ]];
	then
		echo "[ M ::: Decompressing contigs ]"
		pigz -dc "$fasta" > .callorfinput.fa
	else
		ln -s $fasta .callorfinput.fa
	fi
	callorf
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
	rm -rf feat.R pred.R
	exit
fi

rm -rf feat.R pred.R

date
