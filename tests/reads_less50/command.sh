#!/usr/bin/env bash

set -e
macrel reads \
    -1 ../../example_seqs/R1.fq.gz \
    -2 ../../example_seqs/R2.fq.gz \
    -m less50 \
    -o out

gunzip out/macrel.out.prediction.gz
