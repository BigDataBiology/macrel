#!/usr/bin/env bash

set -e

macrel peptides \
    --fasta short_seqs.faa.gz \
    --output out \
    --keep-negatives
gunzip out/macrel.out.prediction.gz
