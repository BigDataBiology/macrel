#!/usr/bin/env bash

set -e

macrel query-ampsphere \
    --fasta pep8.faa \
    --output out
gunzip out/macrel.out.ampsphere_exact.tsv.gz
