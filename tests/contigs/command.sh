#!/usr/bin/env bash

set -e
macrel contigs \
    --fasta excontigs.fna.gz \
    -o out

gunzip out/macrel.out.prediction.gz

