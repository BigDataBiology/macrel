#!/usr/bin/env bash

outfolder_mode_a_fasta="out"
outtag_mode_a_fasta="example_abundance"

macrel abundance \
    -1 ../../example_seqs/R1.fq.gz \
    --fasta ../../example_seqs/ref.faa.gz \
    --output out
if [[ $? != 0 ]]; then
    >&2 echo "macrel failed"
    exit 1
fi

mode_a_fasta_result1=$outfolder_mode_a_fasta/$outtag_mode_a_fasta.xprs

if [[ -f "out/macrel.out.abundance.tsv" ]] ; then
    >&2 echo "bad"
    exit 1
fi
