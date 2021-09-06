#!/usr/bin/env bash

set -e
macrel reads \
    -1 ../../example_seqs/R1.fq.gz \
    -o out

gunzip out/macrel.out.prediction.gz
