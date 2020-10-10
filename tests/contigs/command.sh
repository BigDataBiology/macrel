#!/usr/bin/env bash

set -e
macrel contigs \
    --fasta excontigs.fna.gz \
    -o out \
    --log-file log.txt

gunzip out/macrel.out.prediction.gz

