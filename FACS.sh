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

if [[ -d $Lib/envs/FACS_env/bin ]]
then
    export PATH=$Lib/envs/FACS_env/bin:$PATH
fi



# Default variables
outfolder_DEFAULT="FACS_output"
outfolder=$outfolder_DEFAULT

block_DEFAULT="100M"
block=$block_DEFAULT

j=$(nproc)

outtag_DEFAULT="FACS_OUT"
outtag=$outtag_DEFAULT

mem_DEFAULT="0.75"
mem=$mem_DEFAULT

tp=$(mktemp --tmpdir --directory FACS.XXXXXXX)
cls="1"
ep="0"

log_DEFAULT="log.txt"
log=log_DEFAULT

BUG_URL="https://github.com/FACS-Antimicrobial-Peptides-Prospection/FACS/issues"

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

    Usage: FACS.sh --mode c/r/p/a [options]

    Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].

    Basic options:

    -h, --help          Show this help message

    -m          Mode of operation, one of:
                \"c\" to work with contigs,
                \"p\" to predict AMPs directly from a peptides FASTA file,
                \"r\" to work with paired-end reads,
                \"a\" to map reads against AMP output database and generate abundances table

    --fasta             Contigs or peptides fasta file (possibly compressed)
    --fwd               Short-read sequencing file in FASTQ format (forward reads file)
    --rev               Short-read sequencing file in FASTQ format (reverse reads file)
    --ref               Output of module \"c\" in its raw format [file type TSV and compressed ]
    --outfolder         Folder where output will be generated [$outfolder_DEFAULT]
    --outtag            Tag used to name outputs [Default: $outtag_DEFAULT]
    -t, --threads [N]   Number of threads [Default: number of processors]
    --block             Bucket size (in Bytes, but M/G/T/m/g/t postfixes are accepted). [$block_DEFAULT]
    --log               Log file name. FACS will save the run results to this log file in output folder [$log_DEFAULT].
    --mem               Memory available to FACS ranging from 0 - 1. [$mem_DEFAULT]
    --tmp               Temporary folder
    --cls               Cluster peptides: yes (1) or no (0). [Default: 1 - yes]
    --ep                Extra profilling (solubility, proteases susceptibility and antigenicity): yes (1) or no (0). [Default: 0 - no ]
"
}

# Taking flags
while [[ $# -gt 0 ]]
do
    case $1 in
        -h|--help|-help|\?|-\?)
            show_help
            exit 0
        ;;
        -m|-M|--mode|--Mode|--m|--M)
            mode=${2}
        ;;
        -tmp|-TMP|--tmp|--TMP|--tp)
            rmdir $tp
            tp=${2}
        ;;
        -t|-T|--threads|--Threads|--THREADS|--t|--T)
            j=${2}
        ;;
        --fasta|-fasta|-fa)
            fasta=$(readlink -f ${2})
            if [[ $? != 0 ]]; then
                echo "[ Could not find fasta file (${2}) ]"
                exit 1
            fi
        ;;
        --fwd|-fwd|--FWD|-r1|-R1|-l)
            read_1=$(readlink -f ${2})
            if [[ $? != 0 ]]; then
                echo "[ Could not find forward read file (${2}) ]"
                exit 1
            fi
        ;;
        --rev|--Rev|-rev|-R2|-r2)
            read_2=$(readlink -f ${2})
            if [[ $? != 0 ]]; then
                echo "[ Could not find reverse read file (${2}) ]"
                exit 1
            fi
        ;;
        --outfolder|-outfolder|-of)
            outfolder=$(readlink -f ${2})
            if [[ $? != 0 ]]; then
                echo "[ Could not find output folder (${2}) ]"
                exit 1
            fi
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
            reference=$(readlink -e ${2})
            if [[ $? != 0 ]]; then
                echo "[ Could not find reference file (${2}) ]"
                exit 1
            fi
        ;;
        -mem|--Mem|--mem|-Mem)
            mem=${2}
        ;;
        -cls|--cls)
            cls=${2}
        ;;
        -ep|--ep)
            ep=${2}
        ;;
    esac
    shift
done

if [[ -n $log ]]; then
    log=/dev/null
else
    log=$(readlink -f "$outfolder/$log")
fi

sanity_check ()
# Checking input variables
{

# TODO: Remove duplication in code below
if [[ -z $mode ]]; then
    >&2 echo "[ ERR010 - No FACS mode specified. Please review the command line]"
    exit 1
elif [[ $mode == "r" ]]; then
    if [[ -z $read_1 ]]; then
        >&2 echo "[ ERR010 - Reads file not specified "]
        exit 1
    fi
    if [[ -n $read_2 ]]; then
        if [[ -s $read_1 ]] || [[ -s $read_2 ]]
        then
            echo "[ FACS has found your reads files ($read_1/$read_2), starting work... ]"
            mode="pe"
            echo -e "[ Configuration ]
Mode        $mode
Threads     $j
read_R1     $read_1
read_R2     $read_2
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
        else
            >&2 echo "[ ERR010 - FACS has not found your reads files. ]"
            exit 1
        fi
    elif [[ -s $read_1 ]]
    then
        echo "[ FACS mode has been assigned as single-end reads ]
[ FACS has found your reads files, starting work... ]"
        mode="se"
        echo -e "[ Here we specify your variables ]

Mode        $mode
Threads     $j
read_R1     $read_1
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
    else
        >&2 echo "[ ERR010 - FACS has not found your reads file ($read_1). ]"
        exit 1
    fi
elif [[ $mode == "p" ]]; then
    mode="pep"

    if [ -s "$fasta" ]
    then
        echo "[ FACS has found your peptides fasta file, starting work... ]"
        echo -e "[ Configuration ]

Mode        p
Threads     $j
Contigs     $fasta
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
    else
        >&2 echo "[ ERR011 - Peptides file ($fasta) is not present ]"
        exit 1
    fi
elif [[ $mode == "c" ]]; then
    if [[ -s "$fasta" ]]; then
        echo "[ FACS has found your contigs file, starting work... ]"
        echo -e "[ Here we specify your variables ]
Mode        $mode
Threads     $j
Contigs     $fasta
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
    else
        >&2 echo "[ W ::: ERR011 - Contigs file ($fasta) not found ]"
        exit 1
    fi
elif [[ $mode == "a" ]]; then
    echo "[ M ::: FACS mode has been assigned as read mapper ]"
    if [ -s "$reference" ]
    then
        refmode="tsv"
        if [[ -s "$read_1" ]] && [[ -s "$read_2" ]]
        then
            mode="mpe"
            echo "[ FACS has found your reference data set, starting work... ]"
            echo -e "[ Here we specify your variables ]

** Mapper with paired-end reads
Mode        $mode
Threads     $j
Reference   $reference
R1          $read_1
R2          $read_2
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
        elif [[ -s "$read_1" ]]
        then
            mode="mse"
            echo "[ FACS has found your reference data set, starting work... ]"
            echo -e "[ Here we specify your variables ]

** Mapper with single-end reads
Mode        $mode
Threads     $j
Reference   $reference
R1          $read_1
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
        else
            >&2 echo " [ W ::: ERR011 - FACS has not found your reads files ]"
            exit 1
        fi
    elif [ -s "$fasta" ]
    then
        refmode="fasta"
        if [[ -s "$read_1" ]] && [[ -s "$read_2" ]]
        then
            mode="mpe"
            echo "[ FACS has found your reference data set, starting work... ]"
            echo -e "[ Here we specify your variables ]

** Mapper with paired-end reads
Mode        $mode
Threads     $j
Reference   $fasta
R1          $read_1
R2          $read_2
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
        elif [[ -s "$read_1" ]]
        then
            mode="mse"
            echo "[ FACS has found your reference data set, starting work... ]"
            echo -e "[ Here we specify your variables ]

** Mapper with single-end reads
Mode        $mode
Threads     $j
Reference   $fasta
R1          $read_1
Folder      $outfolder
Tag         $outtag
Bucket      $block
Log         $log"
        else
            >&2 echo "[ ERR011 - FACS has not found your reads files ]"
            exit 1
        fi
    else
        >&2 echo "[ ERR011.1 - FACS has not found your reference file ($reference) ]"
        exit 1
    fi
fi

if [[ ! -e $tp ]]; then
    echo "[ Temporary folder ($tp) does not exist. Creating it... ]"
    mkdir -p $tp
    cd $tp
elif [[ ! -d $tp ]]; then
    echo "[ Temporary folder ($tp) already exists but is not a directory ]" 1>&2
    exit 1
fi
cd $tp
export TMPDIR=$tp

if [[ -z $outfolder ]]; then
    >&2 echo "[ No output folder specified ]"
    exit 1
elif [[ ! -d "$outfolder" ]]; then
    echo "[ Directory $outfolder does not exist.  Creating it... ]" >> "$log"
    mkdir -p $outfolder || (>&2 echo "[ Error creating output folder ($outfolder) ]" ; exit 1)
fi

}

clean_temp ()
{
    cd ..
    rm -rf $tp
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################################    FUNCTIONS   ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
################################# 1. generating two files with matching pairs of reads #####################################

read_trimming ()
{
echo "[ M ::: Trimming low quality bases ]"
if [[ $1 == "paired" ]]
then
    ngl=trim.pe.ngl
    files="$read_1 $read_2"
else
    ngl=trim.se.ngl
    files="$read_1"
fi

ngless --no-create-report --quiet -j $j --temporary-directory $tp $Lib/$ngl $files
if [[ $? != 0 ]]
then
    >&2 echo "[ W ::: ERR231 - Read trimming failed ]"
    cd ../
    rm -rf $tp
    exit 1
fi
}

assemble ()
{
echo "[ M ::: Assembly using MEGAHIT ]"
if [[ $mode == "pe" ]]
then
    read_trimming "paired"
    megahit_args="-1 preproc.pair.1.fq.gz -2 preproc.pair.2.fq.gz"
elif [[ $mode == "se" ]]
then
    read_trimming "single"
    megahit_args="-r preproc.pair.1.fq.gz"
else
    >&2 echo "[ ERR222 - Internal Error: should never happen. Please report a bug at $BUG_URL ]"
    cd ../
    rm -rf $tp
    exit 1
fi

megahit --presets meta-large $megahit_args -o out -t "$j" -m "$mem" --min-contig-len 1000

if [[ $? != 0 ]]
then
    >&2 echo "[ Megahit failed "]
    cd ../
    rm -rf $tp
    exit 1
fi

if [[ -s out/final.contigs.fa ]]
then
    mv out/final.contigs.fa contigs.fna
    rm -rf out/ preproc.pair.1.fq.gz preproc.pair.2.fq.gz
else
    >&2 echo "[ W ::: ERR128 - Assembly failed ]"
    cd ../; rm -rf $tp/
    exit 1
fi
}

predict_smorfs ()
{

if [[ "$mode" != "pep" ]]
then
    echo "[ M ::: Predicting ORFs ]"

    rm -rf callorfs/
    mkdir callorfs

    parallel -j $j --block $block --recstart '>' --pipe \
        prodigal_sm -c -m -n -p meta -f sco -a "callorfs/{#}.pred.smORFs.fa" \
        <contigs.fna >/dev/null
    if [[ $? != 0 ]]
    then
        >&2 echo "[ ERR910 - Error in prodigal_sm ]"
        clean_temp
        exit 1
    fi

    if [[ $mode == "c" ]]
    then
        rm -rf contigs.fna
    else
        mv contigs.fna "$outfolder/$outtag.contigs.fna"
        pigz --best "$outfolder/$outtag.contigs.fna" &
        waitfor=$!
    fi

    echo "[ M ::: Filtering smORFs ]"

    for i in callorfs/*;
    do
        python "$Lib/filter_smorfs.py" "$i" >> pep.faa.tmp
        rm -rf $i
    done

    rm -rf callorfs/

    if [[ -n $waitfor ]]; then
        wait $waitfor
        unset waitfor
    fi
else
    python "$Lib/filter_smorfs.py" contigs.fna > pep.faa.tmp
    rm -rf contigs.fna
fi

if [[ ! -s "pep.faa.tmp" ]]; then
    >&2 echo "[ ERR122 - ORFs calling procedure failed ]"
    clean_temp
    exit 1
fi
if [[ $cls == 1 ]]; then
    echo "[ Clustering peptides ]"

    mkdir sorting_folder/

    sed 's/ /_/g' pep.faa.tmp | sed 's/;/|/g' \
        | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' \
        | awk '{print $2"\t"$1}' \
        | parallel --pipe  -j $j --block "$block" --recstart "\n" "cat > sorting_folder/small-chunk{#}" >/dev/null
    rm -rf pep.faa.tmp

    for X in $(ls sorting_folder/small-chunk*); do sort -S 80% --parallel="$j" -k1,1 -T . < "$X" > "$X".sorted; rm -rf "$X"; done

    sort -S 80% --parallel="$j" -T . -k1,1 -m sorting_folder/small-chunk* > .sorted-huge-file; rm -rf sorting_folder/

    perl -n -e "print unless /^$/" .sorted-huge-file \
        | awk -F'\t' -v OFS=';' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' \
        | sed 's/;;/\t/g' \
        | sed 's/>//g' \
        | sed 's/\*//g' \
        > .tmp2

    if [[ -s .tmp2 ]]
    then
        echo "[ M ::: Eliminated all duplicated sequences ]"
        rm -rf .sorted-huge-file
    else
        >&2 echo "[ W ::: ERR510 - Memdev sorting stage has failed ]"
        cd ../; rm -rf $tp
        exit 1
    fi

    echo "[ M ::: Outputting genes lists ]"

    awk '{print ">smORF_"NR"\t"$1"\t"$2}' .tmp2 | pigz --best > $outfolder/$outtag.ids.tsv.gz

    echo "[ M ::: Outputting fastas ]"

    awk '{print ">smORF_"NR"\n"$1}' .tmp2 > pep.faa
    rm -rf .tmp2

    echo "[ M ::: Collecting statistics ]"

    number_of_peptides=$(tail -2 pep.faa \
        | head -1 \
        | sed 's/|/\t/g' \
        | cut -f2 \
        | sed 's/>smORF_//g')
else
    echo "[ M ::: Skipping clustering ]"
    mv pep.faa.tmp pep.faa
    number_of_peptides=$(grep -c ">" pep.faa)
fi


}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
######################################## 5. computing descriptors ########################################################

descriptors ()
{
echo "[ M ::: Split input into chunks ]"
mkdir splits/
parallel --pipe  -j $j --block "$block" --recstart ">" "cat > splits/small-chunk{#}.faa" < pep.faa >/dev/null
if [[ $? != 0 ]]; then
    >&2 echo "[ Input splitting failed ]"
    clean_temp
    exit 1
fi
rm -rf pep.faa

for i in splits/small-chunk*.faa; do
    echo "Calling $Lib/AMP-features.py $i" >> "$log"
    date >> "$log"
    echo >> "$log"

    python "$Lib/AMP-features.py" "$i" "$i.tabdesc.tsv"
    if [[ $? != 0 ]]; then
        >&2 echo "[ Features computation failed ]"
        clean_temp
        exit 1
    fi
    rm -rf "$i"
done
}

predict ()
{
if [ "$(ls -A splits/)" ]
then

    for i in splits/small-chunk*
    do
        echo "[ M ::: Predicting AMPs -- $i ]" | tee -a "$log"
        finfile="${i/.tabdesc.tsv/.fin}"
        Rscript --vanilla $Lib/Predict_130819.R "$i" "$Lib"/r22_largeTraining.rds "$Lib"/rf_dataset1.rds "$finfile" >"$log"
        if [[ ! -s "$finfile" ]]
        then
            echo "[ W ::: Sample $i did not return any AMP sequences ]"
        fi
        rm -rf "$i" "$Lib"/__pycache__/
    done
else
    cd ../
    rm -rf $tp
    >&2 echo "[ There are no valid sequences to be processed over pipeline, we are sorry -- ERR675 ]"
    exit 1
fi
}

EXTRADATA()
{
awk '{print ">"$1"\n"$2}' .out2.file > .tmp.fasta
cut -f1 .out2.file | sort | uniq  > listall

# To install protein solubility software
# wget --header 'Host: protein-sol.manchester.ac.uk' --referer 'https://protein-sol.manchester.ac.uk/software' --header 'Upgrade-Insecure-Requests: 1' 'https://protein-sol.manchester.ac.uk/cgi-bin/utilities/download_sequence_code.php' --output-document 'protein-sol-sequence-prediction-software.zip'
# unzip protein-sol-sequence-prediction-software.zip
# mv protein-sol-sequence-prediction-software/ $Lib/envs/FACS_env/bin/

# To install emboss:
# wget -m 'ftp://emboss.open-bio.org/pub/EMBOSS/'
# cp emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz ./
# tar -zxf EMBOSS-6.6.0.tar.gz
# cd EMBOSS-6.6.0/
# ./configure
# make
# cd emboss/
# rm -rf *.c *.h
# mv * $Lib/envs/FACS_env/bin/

echo "[ W ::: Just extra moments to assess the predicted AMPs ]"
echo "[ W ::: Predicting solubility ]"
cp $Lib/envs/FACS_env/bin/protein-sol-sequence-prediction-software/* ./
sh multiple_prediction_wrapper_export.sh .tmp.fasta
grep "SEQUENCE PREDICTIONS,>" seq_prediction.txt | sed 's/SEQUENCE PREDICTIONS,>//g' | sed 's/,/\t/g' > protein-sol.res.tsv
rm -rf .tmp* seq_prediction.txt
awk '{print ">"$1"\n"$2}' .out2.file > .tmp.fasta
cat protein-sol.res.tsv | cut -f1,2,3 | sort -k1,1 > sol
rm -rf protein-sol.res.tsv

echo "[ W ::: Predicting antigenicity ]"
antigenic -sequence .tmp.fasta -sprotein1 Y -sformat1 FASTA -minlen 9 -outfile antigenic -rformat excel -rmaxseq2 1 -raccshow2 1 >/dev/null
sed -i '/SeqName/d' antigenic
echo -e "SeqName\tStart\tEnd\tScore\tStrand\tMax_score_pos" > header
cat antigenic | cut -f1 | sort | uniq | awk '{print $1"\t""+"}' > antigen
cat antigenic | cut -f1 | sort | uniq > list0
cat header antigenic | pigz --best > "$outfolder"/"$outtag".antigenic.tsv.gz
rm -rf header antigenic

echo "[ W ::: Creating lists ]"
grep -v -w -f list0 listall > non_antigenic_AMP; rm -rf list0
awk '{print $1"\t""-"}' non_antigenic_AMP > list1
cat antigen list1 | sort -k1,1 > tmp; mv tmp antigen
rm -rf list1 list2 non_antigenic_AMP

echo "[ W ::: Predicting protease sensitivity ]"
echo "[ W ::: This can take a while, since analysis proceed in single sequences ]"
awk '/^>/ {OUT=substr($0,2) ".seq";print " ">OUT}; OUT{print >OUT}' .tmp.fasta

for i in *.seq;
do
    epestfind -sequence $i  -window 9 -order score -outfile ${i/.seq/.epest} -graph none -nopoor -nomap >/dev/null
    cat ${i/.seq/.epest} >> final
    rm -rf ${i/.seq/.epest} $i
done

cat final \
    | sed '/PEST-find/d' \
    | grep "No PEST motif was identified in " \
    | sed 's/No PEST motif was identified in //g' \
    | sed 's/^.   //g' \
    | sed 's/ .*//g' \
    | awk '{print $1"\t""-"}' \
    > nonprotea.list

cat final \
    | sed '/PEST-find/d' \
    | grep -v "No PEST motif was identified in " \
    | sed '/^[[:space:]]*$/d' \
    | pigz --best \
    > "$outfolder"/"$outtag".protealytic_assessment.txt.gz

cat final \
    | grep -v "No PEST motif was identified in " \
    | grep "PEST motif was identified in " \
    | sed 's/^.* PEST motif was identified in //g' \
    | sed 's/ .*//g' \
    | awk '{print $1"\t""+"}' \
    > protea.list

cat protea.list nonprotea.list | sort -k1,1 > protea
rm -rf final *.list

pl=$(zcat "$outfolder"/"$outtag".protealytic_assessment.txt.gz | awk '{print NR}' | tail -1)
if [[ $pl -gt 1 ]]
then
    echo ""
else
    rm -rf zcat "$outfolder"/"$outtag".protealytic_assessment.tsv.gz
    echo "[ W ::: No AMPs returned EPEST motifs ]"
fi

echo "[ W ::: Checking structures ]"
echo -e "antigen\nprotea\nsol" > quicktest.list
fin=$(awk '{print NR}' .out2.file | tail -1)
for i in $(cat quicktest.list);
do
    rown=$(awk '{print NR}' $i | tail -1)
    if [[ "$fin" == "$rown" ]]
    then
        touch $i
    else
        echo "[ W ::: Ending loop ]"
        echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" > header
        cat header .out2.file > tmp; mv tmp "$outfolder"/"$outtag".tsv
        rm -rf header .out2.file protea antigen sol quicktest.list
        echo "[ W ::: Error // Wrong row formatting -- ERR 989 ]"
    fi
    unset rown
done
unset fin

if [ -s .out2.file ]
then
    echo "[ W ::: Formatting lists ]"
    sort -k1,1 .out2.file > tmp1
    sort -k1,1 protea | awk '{$1=""}1' | sed 's/^.* //g'  > tmp4
    sort -k1,1 antigen | awk '{$1=""}1' | sed 's/^.* //g'  > tmp3
    sort -k1,1 sol | awk '{$1=""}1' | sed 's/^ //g' | sed 's/ /\t/g' > tmp2
    rm -rf .out2.file protea antigen sol quicktest.list

    paste -d'\t' tmp1 tmp2 tmp3 tmp4 > .out.file
    rm -rf .out2.file protea antigen sol

    echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability\tPercent_soluble\tScaled_Solubility\tAntigenicity\tSusceptible_to_proteases" > header
    cat header .out.file > tmp; mv tmp "$outfolder"/"$outtag".tsv
    rm -rf header .out.file

    echo "[ W ::: Deciding promising peptides ]"
    awk '$5 == "NonHEMO" && $8 >= 0.45 && $9 == "-" && $10 == "-"' "$outfolder"/"$outtag".tsv > "$outfolder"/"$outtag".promising_subset.tsv
    if [ -s "$outfolder"/"$outtag".promising_subset.tsv ]
    then
        echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability\tPercent_soluble\tScaled_Solubility\tAntigenicity\tSusceptible_to_proteases" > header
        cat header "$outfolder"/"$outtag".promising_subset.tsv > tmp; mv tmp "$outfolder"/"$outtag".promising_subset.tsv
        pigz --best "$outfolder"/"$outtag".promising_subset.tsv
        rm -rf tmp header
        nr=`zcat "$outfolder"/"$outtag".promising_subset.tsv | sed '1,1d' | awk '{print NR}' | tail -1`
        echo "[ W ::: Exists $nr promising AMP candidates saved as "$outfolder"/"$outtag".promising_subset.tsv.gz ]"
    else
        rm -rf "$outfolder"/"$outtag".promising_subset.tsv
        echo "[ W ::: No promising AMP candidates found ]"
    fi
else
    echo "[ W ::: Ending subloop ]"
fi
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#################################### 6. Generating abundance profiles ######################################################

format_results ()
{
echo "[ M ::: Formatting results ]"

if [[ "$(ls -A splits/)" ]]; then
    cat splits/*  \
        | sed 's/\"//g' \
        > .out2.file
    rm -rf splits/
    AMP=$(cat .out2.file | wc -l | awk '{print $1}')
    echo "[ M ::: Calculating statistics ]"
    pepval=$(awk '{print $3"_"$5}' .out2.file | sort | uniq -c | awk '{ print $2"\t"$1 }')
    echo "[ M ::: Exporting results ]"
    if [[ $ep == "0" ]]
    then
        (echo -e "Access\tSequence\tAMP_family\tAMP_probability\tHemolytic\tHemolytic_probability" ; cat .out2.file) > "$outfolder"/"$outtag".tsv
        rm -rf .out2.file
    elif [[ $ep == "1" ]]
    then
        EXTRADATA
    fi
    pigz --best "$outfolder"/"$outtag".tsv
else
    echo "[ We are sorry to inform. Your procedure did not return any AMP sequence. We are closing now ]" | tee -a "$log"
    echo "[ We really tried. However, none of $number_of_peptides tested peptides have a probability higher than 0.5 of being an AMP ]" | tee -a "$log"
    cd ../
    rm -rf $tp
    exit 0
fi
}

summary ()
{
echo "[ M ::: Generating log ]"

# TODO Remove code duplication below:
if [[ $mode == "r" ]]
then
    ## Preparing report
    echo -e "######################################################## FACS report

[ M ::: Variables ]

Threads     $j
read_R1     $read_1
read_R2     $read_2
Folder      $outfolder
Tag         $outtag
Bucket      $block

========================================= Files were treated

*** Deduplicated ORFs were called

A total of $number_of_peptides of smORFs were called.

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
########################################################" >> "$log"
else
    ## Preparing report
    echo -e "######################################################## FACS report

[ M ::: Variables ]

Threads     $j
Contigs     $fasta
Folder      $outfolder
Tag         $outtag
Bucket      $block

========================================= Files were treated

*** Deduplicated ORFs were called

A total of $number_of_peptides of smORFs were called.

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

cat "$outfolder"/.log | tee -a "$log"

}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###########################################      7. Abundances      ########################################################

mapping ()
{
echo "[ M ::: Indexing references ]"
if [[ "$refmode" == "tsv" ]]; then
    pigz -dc "$reference" | sed '1,1d' | awk '{ print ">"$1"\n"$2 }' > .ref.fa
elif [[ "$refmode" == "fasta" ]]; then
    if [[ "$fasta" =~ \.gz$ ]];
    then
        echo "[ M ::: Decompressing files ]"
        pigz -dc "$fasta" > .ref.fa
    else
        ln -s "$fasta" .ref.fa
    fi
else
    >&2 echo "[ Internal error in FACS (refmode is $refmode). Please report a bug at $BUG_URL ]"
    exit 1
fi

paladin index -r 3 .ref.fa

if [[ $? != 0 ]] || [[ ! -s .ref.fa.amb ]]; then
    rm -rf .ref.fa*
    >&2 echo "[ ERR303 - Error in indexing ]"
    clean_temp
    exit 1
fi

if [[ $mode != "mse" ]] && [[ $mode != "mpe" ]]; then
    >&2 echo "[ Internal sanity check failure (expected internal mode to be mse or mpe, is $mode). Please report a bug at $BUG_URL ]"
    clean_temp
    exit 1
fi

echo "[ M ::: Mapping reads against references, be aware it can take a while ]"

echo "[ M ::: Starting paladin ]"
paladin align -t "$j" -T 20 -f 10 -z 11 -a -V -M .ref.fa preproc.pair.1.fq.gz | samtools view -Sb | samtools sort > .m.bam

if [[ -s .m.bam ]]
then
    touch .m.bam
else
    >&2 echo "[ ERR052 - Mapping failed ]"
    cd ../; rm -rf $tp
    exit 1
fi
}


ab_profiling ()
{
echo "[ M ::: Expressing results of abundance ]"
express --no-bias-correct -o "$outfolder"/ .ref.fa .m.bam
if [[ -s "$outfolder"/results.xprs ]] || [[ -s "$outfolder"/params.xprs ]]
then
    mv "$outfolder"/results.xprs "$outfolder"/"$outtag".xprs
    mv "$outfolder"/params.xprs "$outfolder"/"$outtag".params.xprs
else
    echo "[ W ::: ERR054 - Abundance calling failed ]"
fi
rm -rf .ref.* preproc.*.fq.gz .m*
}


sanity_check

if [[ $mode == "pe" || $mode == "se" ]]
then
    assemble
    predict_smorfs
    descriptors
    predict
    format_results
    mode="r"
    summary
elif [[ $mode == "pep" ]]
then
    if [[ "$fasta" =~ \.gz$ ]];
    then
        echo "[ M ::: Decompressing contigs ]"
        pigz -dc "$fasta" > contigs.fna
    else
        ln -s $fasta contigs.fna
    fi
    predict_smorfs
    descriptors
    predict
    format_results
    summary
elif [[ $mode == "c" ]]
then
    if [[ "$fasta" =~ \.gz$ ]];
    then
        echo "[ M ::: Decompressing contigs ]"
        pigz -dc "$fasta" > contigs.fna
    else
        ln -s $fasta contigs.fna
    fi
    predict_smorfs
    descriptors
    predict
    format_results
    summary
elif [[ $mode == "mse" || $mode == "mpe" ]]
then
    if [[ $mode == "mse" ]]
    then
        read_trimming "single"
    elif [[ $mode == "mpe" ]]
    then
        echo "[ W ::: NOTE _ IMPORTANT ]"
        echo "[ W ::: Abundance data is just inferred from R1 file ]"
        read_trimming "paired"
    fi
    mapping
    ab_profiling
else
    >&2 echo "[ W ::: ERR010 - The user needs to specify a valid FACS mode, please review the command line]"
    cd ../; rm -rf $tp/
    exit 1
fi
cd ../
rm -rf $tp/
date
