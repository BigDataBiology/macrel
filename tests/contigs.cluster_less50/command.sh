#!/usr/bin/env bash

set -e
macrel contigs \
    --fasta excontigs_with_repeats.fna.gz \
    --cluster \
    -m less50 \
    -o out

gunzip out/macrel.out.prediction.gz

