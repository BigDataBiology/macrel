#!/bin/bash
##################################################################################
#### This is spurius module                                        ###############
#### It is based in  https://f1000research.com/articles/7-261      ###############
#### Authors: Celio Dias Santos Junior, Luis Pedro Coelho          ###############
#### Affiliations: Fudan University                                ###############
#### It predicts which proteins are novel, spurius and non spurius ###############
##################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##################################################################################
# Stating hard codes
Lib="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
db1="$Lib/db/dbfilter2_part1/dbfilter2_part1"
db2="$Lib/db/dbfilter2_part2/dbfilter2_part2"
db3="$Lib/db/fullsource_filter.fa"
tblastn="/usr/bin/tblastn"
formatter="/usr/bin/blast_formatter"
merger="$Lib/blastXMLmerge.py"
order_bed="$Lib/bedtools_prep.py"
bedtools="/usr/bin/bedtools"
calc1="$Lib/findorfs.py"
pred="$Lib/spu_pred.py"
model="$Lib/clf1.npy"
dictio="$Lib/dict_np.npy"
outtag="spuriomod"

# Help message
show_help ()
{
	echo "
	############################################################################################################################
	############################## FACS pipeline - SPURIUS mode                             ####################################
	############################################################################################################################
	############################## Authors: Célio Dias Santos Júnior, Luis Pedro Coelho     ####################################
	############################################################################################################################
	############################## Institute > ISTBI - FUDAN University / Shanghai - China  ####################################
	############################################################################################################################

	Usage: ./spurius.mod.sh --facs <FACS FILE> -DB /path/to/db/

	Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].

	Basic options:

	-h, --help	        Show this help message

	-m			AI model of operation - Addresss
	-d                      Dictionary - Address
	-DB			Database folder
	--FACS			Compressed (or not gzipped) proteins \".tsv.gz\" FACS output
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
		-m|-M|--model|--Model|--m|--M)
			model=${2}
		;;
		--facs|--FACS|-FACS|-facs)
			fasta=${2}
		;;
		--d|--D|-d|-D)
			dictio=${2}
		;;
		--db|--DB|-db|-DB)
			db=${2}
		;;
	esac
	shift
done

sanity_check ()
# Checking input variables
{
if [[ -n $fasta ]]
then
	echo "[ M ::: SPURIUS starting ]"
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

convertFACS ()
{
pigz -dc $fasta | awk '{print ">"$1"\n"$2}' | sed '1,2d' > .input
}

splitter ()
{
mkdir splitted
cat .input | parallel -N 100 --recstart '>' --pipe "cat > splitted/small-chunk{#}"
rm -rf .input
}

searching ()
{
mkdir searches/
for i in splitted/*;
do 
	echo "[ M ::: Processing Bucket -- $i ]"
	# Running qlen
	awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $i > tmp; mv tmp $i
	awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $i | perl -n -e "print unless /^$/" | awk '{print NR"\t"$1"\t"$2}' > .namelist
	# Running db1
	echo "[ M ::: Blast with DB1 ]"
	parallel --recstart '>' -N 1 --pipe $tblastn -query - -db $db/dbfilter2_part1/dbfilter2_part1 -out splitted/{#}.part1.asn -outfmt 11 -max_target_seqs 10000 -evalue 1 -num_threads 1 -dbsize 5448272552 < $i
	# Running db2
	echo "[ M ::: Blast with DB2 ]"
	parallel --recstart '>' -N 1 --pipe $tblastn -query - -db $db/dbfilter2_part2/dbfilter2_part2 -out splitted/{#}.part2.asn -outfmt 11 -max_target_seqs 10000 -evalue 1 -num_threads 1 -dbsize 5448272552 -db_gencode 4 < $i
	echo "[ M ::: Starting merging of files ]"
	while read number access qlen;
	do
		$formatter -archive splitted/$number.part1.asn -out searches/$access.part1.xml -outfmt 5
		$formatter -archive splitted/$number.part2.asn -out searches/$access.part2.xml -outfmt 5
		rm -rf splitted/$number.part*
		if [ -s searches/$access.part1.xml ]
		then
			if [ -s searches/$access.part2.xml ]
			then
				python3 $merger searches/$access.xml searches/$access.part1.xml searches/$access.part2.xml
				rm -rf searches/$access.part1.xml searches/$access.part2.xml
			else
				mv searches/$access.part1.xml searches/$access.xml
				rm -rf searches/$access.part1.xml searches/$access.part2.xml
			fi
		elif [ -s searches/$access.part2.xml ]
		then
			mv searches/$access.part2.xml searches/$access.xml
			rm -rf searches/$access.part1.xml searches/$access.part2.xml
		else
			rm -rf searches/$access.part1.xml searches/$access.part2.xml
			rem=$(echo -e "$access\t$qlen"); grep -v "$rem" .namelist > tmp; mv tmp .namelist
			rem="0"
		fi
		echo "[ M ::: Producing Bedfiles // sequence -- $access : $number ]"
		python3 $order_bed searches/$access.xml searches/$access.bedrequest $qlen
		echo "[ M ::: Getting fastas ]"
		$bedtools getfasta -s -fi $db//fullsource_filter.fa -bed searches/$access.bedrequest -fo searches/$access.fasta -name 
	done < .namelist
	touch reflist
	cat .namelist reflist > tmp; mv tmp reflist; rm -rf .namelist
	rm -rf $i
done
awk '{print NR"\t"$2"\t"$3}' reflist > tmp; mv tmp reflist
rm -rf splitted/
}

calc1 ()
{
echo "[ M ::: Producing NPY files ::: ]"
mkdir npdir
while read nmb access qlen;
do
	echo "[ M ::: Processing sequence -- $access : $nmb ]"
	python3 $calc1 $access searches/$access.fasta $qlen $nmb $outtag 1 $dictio >/dev/null 2>/dev/null
done < reflist
rm -rf searches/
}

pred ()
{
echo "[ M ::: Predicting spurius coefficients ]"
python3 $pred $model >/dev/null 2>/dev/null
mv predictions_newer_spuriomod spurious_out.txt
rm -rf npdir/
sed -i 's/UniProt //g' spurious_out.txt
}

############################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################################    CMDs   #################################################################

export PATH=$PATH:$Lib

sanity_check
convertFACS
splitter
searching
calc1
pred
rm -rf reflist
