#!/usr/bin/bash
## Generate ini files to produce metagenomes
	
for i in $(ls ../ini_files/); do ./gen-paired-end-reads ini_files/$i; pigz --best *.fastq; done

