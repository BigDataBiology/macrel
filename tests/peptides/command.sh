#!/usr/bin/env bash

set -e

macrel peptides \
    --fasta expep.faa.gz \
    --output out
gunzip out/macrel.out.prediction.gz
