#!/usr/bin/env bash

input="$1"
Lib="$2"
log="$3"

trap "rm -f $input" EXIT

echo "[ Counting distribution using SA scale -- $input ]"
python3 "$Lib"/CTDDClass.py "$input" CTDDC-SA-full.tsv 'ALFCGIVW' 'RKQEND' 'MSPTHY' #solventaccess
if [[ $? != 0 ]]; then
    >&2 echo "[ CTDDClass failed ]"
    exit 1
fi
awk 'BEGIN { print "SA.G1.residue0\tSA.G2.residue0\tSA.G3.residue0" } (NR > 1) {print $2"\t"$7"\t"$12}' CTDDC-SA-full.tsv > CTDDC-SA.tsv

rm -rf CTDDC-SA-full.tsv

echo "[ Counting distribution using HB scale -- $input ]"
python3 "$Lib"/CTDDClass.py "$input" CTDDC-HB-full.tsv 'ILVWAMGT' 'FYSQCN' 'PHKEDR' # HEIJNE&BLOMBERG1979
if [[ $? != 0 ]]; then
    >&2 echo "[ CTDDClass failed ]"
    exit 1
fi
awk 'BEGIN { print "hb.Group.1.residue0\thb.Group.2.residue0\thb.Group.3.residue0" } (NR > 1) {print $2"\t"$7"\t"$12}' CTDDC-HB-full.tsv > CTDDC-HB.tsv;

rm -rf CTDDC-HB-full.tsv

if [[ ! -s CTDDC-SA.tsv ]] || [[ ! -s CTDDC-HB.tsv ]]
then
    >&2 echo "[ W ::: Error in predictors calculation during CTD steps -- ERR 230 ]"
    exit 1
fi

echo "[ Computing cheminformatics descriptors -- $input ]"
echo -e "header\tseq\tgroup" > .tmp
sed '/>/d' "$input" > .seqs
grep '>' "$input" | sed 's/ .*//g' | sed 's/>//g' > .heade
paste -d'\t' .heade .seqs | awk '{print $1"\t"$2"\t""Unk"}' >> .tmp
rm -rf .heade .seqs

echo "Calling $Lib/features_130819.R" >> "$log"
date >> "$log"
echo >> "$log"
Rscript --vanilla $Lib/features_130819.R .tmp .out.file >>"$log"
if [[ $? != 0 ]]; then
    >&2 echo "[ Calling features_130819.R failed ]"
    exit 1
fi


rm -rf .tmp

if [[ -s .out.file ]]; then
    paste -d'\t' .out.file CTDDC-SA.tsv CTDDC-HB.tsv | sed 's/\"//g' > $input.tabdesc.tsv
    rm -rf .out.file CTDDC-SA.tsv CTDDC-HB.tsv
else
    >&2 echo "[ Error in predictors calculation during cheminformatics steps -- ERR 229 ]"
    rm -rf CTDDC-SA.tsv CTDDC-HB.tsv
    exit 1
fi

