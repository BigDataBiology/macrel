#!/usr/bin/env bash

set -e

mkdir -p out

macrel get-smorfs -f excontigs.fna.gz --file-output - --log-file log1.txt > out/macrel.out.smorfs.faa
macrel get-smorfs -f excontigs.fna.gz --file-output - --keep-fasta-headers --log-file > (gzip > log2.txt.gz) > out/macrel.out.smorfs.w_headers.faa

test -f log1.txt
test -f log2.txt.gz

