#!/usr/bin/env bash

set -e
macrel peptides \
    --fasta ../../example_seqs/expep.faa.gz \
    -m less50 \
    --output out

gunzip out/macrel.out.prediction.gz
