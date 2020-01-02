#!/usr/bin/env bash

input="$1"
Lib="$2"
log="$3"

trap "rm -f $input" EXIT

echo "[ Counting distribution using SA scale -- $input ]"
python3 "$Lib"/features.py "$input" CTDDC-SA.tsv 'SA.G' 'ALFCGIVW' 'RKQEND' 'MSPTHY' #solventaccess
if [[ $? != 0 ]]; then
    >&2 echo "[ CTDDClass failed ]"
    exit 1
fi

echo "[ Counting distribution using HB scale -- $input ]"
python3 "$Lib"/features.py "$input" CTDDC-HB.tsv 'hb.Group.' 'ILVWAMGT' 'FYSQCN' 'PHKEDR' # HEIJNE&BLOMBERG1979
if [[ $? != 0 ]]; then
    >&2 echo "[ CTDDClass failed ]"
    exit 1
fi

if [[ ! -s CTDDC-SA.tsv ]] || [[ ! -s CTDDC-HB.tsv ]]
then
    >&2 echo "[ W ::: Error in predictors calculation during CTD steps -- ERR 230 ]"
    exit 1
fi

echo "[ Computing cheminformatics descriptors -- $input ]"
echo "Calling $Lib/features_R.py" >> "$log"
date >> "$log"
echo >> "$log"
python3 "$Lib/features_R.py" "$input" Rfeatures.out.file >> "$log"
if [[ $? != 0 ]]; then
    >&2 echo "[ Calling features_R.py failed ]"
    exit 1
fi


if [[ -s Rfeatures.out.file ]]; then
    paste -d'\t' Rfeatures.out.file CTDDC-SA.tsv CTDDC-HB.tsv | sed 's/\"//g' > $input.tabdesc.tsv
    rm -rf Rfeatures.out.file CTDDC-SA.tsv CTDDC-HB.tsv
else
    >&2 echo "[ Error in predictors calculation during cheminformatics steps -- ERR 229 ]"
    rm -rf CTDDC-SA.tsv CTDDC-HB.tsv
    exit 1
fi

