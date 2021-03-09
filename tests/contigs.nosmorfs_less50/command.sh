#!/usr/bin/env bash

set -e
macrel contigs \
    --fasta no_smorfs.fna.gz \
    -o out \
    -m less50

gunzip out/macrel.out.prediction.gz
