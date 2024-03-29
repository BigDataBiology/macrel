#!/usr/bin/env bash

set -e
macrel reads \
    -1 R1.fq.gz \
    -2 R2.fq.gz \
    -o out

gunzip out/macrel.out.prediction.gz
gunzip out/macrel.out.percontigs.gz
