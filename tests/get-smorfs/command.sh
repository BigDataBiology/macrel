#!/usr/bin/env bash

set -e

mkdir -p out
macrel get-smorfs -f excontigs.fna.gz --file-output - > out/macrel.out.smorfs.faa

