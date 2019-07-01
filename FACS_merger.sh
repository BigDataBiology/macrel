#!/usr/bin/env bash

############################################################################################################################
############################## FACS pipeline - FACS results merger                      ####################################
############################################################################################################################
############################## Authors: Célio Dias Santos Júnior, Luis Pedro Coelho     ####################################
############################################################################################################################
############################## Institute > ISTBI - FUDAN University / Shanghai - China  ####################################
############################################################################################################################
############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
##################################################### Variables ############################################################
pigz="/usr/bin/pigz"
output="out.gz"
j=`echo "$(nproc) * 9 / 10" | bc` # Use 90% of threads
ReferenceFolder="./"

show_help ()
{
	echo "
	############################################################################################################################
	############################## FACS pipeline - FACS results merger                      ####################################
	############################################################################################################################
	############################## Authors: Célio Dias Santos Júnior, Luis Pedro Coelho     ####################################
	############################################################################################################################
	############################## Institute > ISTBI - FUDAN University / Shanghai - China  ####################################
	############################################################################################################################

	Usage: FACS_merger.sh [options]
	
	Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].
	
	N=Integer (eg. 100)
	
	Tutorial options:
		-h, --help	Show this help message
	
	Basic options:
		--output [file]	File where output is sent, it needs to be gzipped (ending em .gz)
		-t, --threads [N]	Number of threads [90% of avaiable threads]
		--reference [folder]	Folder where your reference files are located, if none current folder will be used [Reference_seqs]
"
}

while [[ $# -gt 0 ]]
do
	case $1 in
		-h|--help|-help|\?|-\?)
			show_help
			exit
		;;
		--output|--Output|--OUTPUT|-output|-Output|-OUTPUT|-O|-o|-Out|-OUT|-out|--out|--Out|--OUT|--o|--O|--address|--Address|--ADDRESS)
			if [[ -n "${2}" ]]
			then
				output=${2}
			fi
			shift
		;;
		-t|-T|--threads|--Threads|--THREADS|--t|--T)
			if [[ -n "${2}" ]]
			then
				j=${2}
			fi
			shift
		;;
		--reference|--Reference|--REFERENCE|--references|--References|--REFERENCES|--ref|--Ref|--REF|--refs|--Refs|--REFS|-ref|-Ref|-REF|-reference|-references)
			if [[ -d "${2}" ]]
			then
				ReferenceFolder=$(cd "$(dirname "${2}")" && pwd)/$(basename "${2}")
			fi
			shift
		;;
	esac
	shift
done

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#####################################################    CMDs   ############################################################

echo "[M ::: Starting in the reference folder]"
cd $ReferenceFolder

echo "[M ::: Arranging headers]"
ls *.tsv.gz | awk -vRS="\n" -vORS="\t" '1' | sed 's/\.tsv\.gz//g' | awk '{print "Sequence""\t"$0}' > header

if [[ -s header ]]
then
	echo ".."
else
	echo "[M ::: Your references do not return any positive value, please redesignate the references folder]"
	exit
fi

echo "[M ::: Intermediate files]"
for i in $(ls *.tsv.gz); do $pigz -dc $i | awk -F'\t' '{print $3"\t"$2}' | sed '1,1d' | sort -k1,1 -S 80% --parallel=$j -T . > ${i/.tsv.gz/.tmp}; $pigz -dc $i | awk -F'\t' '{print $3}' | sed '1,1d' | sort -k1,1 -S 80% --parallel=$j -T . > ${i/.tsv.gz/.tmp2}; done

echo "[M ::: Sorting main list]"
sort -S 80% --parallel=$j -T . -m *.tmp2 | uniq > sorted-huge-file; rm -rf *.tmp2

echo "[M ::: Starting main loop]"
for i in $(ls *tmp); do join -1 1 -2 1 -a1 -e0 -o'0,2.2' sorted-huge-file $i | sed 's/ /\t/g' | awk '{print $2}' > ${i/.tmp/.tmp2}; mv ${i/.tmp/.tmp2} $i; done

echo "[M ::: Arranging files]"
paste sorted-huge-file *tmp > final.t; rm -rf *.tmp sorted-huge-file

echo "[M ::: Formatting]"
cat header final.t | pigz --best > $output

echo "[M ::: Cleaning]"
rm -rf header final.t
