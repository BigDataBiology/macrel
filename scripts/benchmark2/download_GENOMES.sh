#!/bin/bash
echo "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/593/425/GCF_001593425.2_ASM159342v2/GCF_001593425.2_ASM159342v2_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/154/225/GCF_000154225.1_ASM15422v1/GCF_000154225.1_ASM15422v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/005/GCF_000008005.1_ASM800v1/GCF_000008005.1_ASM800v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/825/GCF_000012825.1_ASM1282v1/GCF_000012825.1_ASM1282v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/276/455/GCA_002276455.1_ASM227645v1/GCA_002276455.1_ASM227645v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/006/125/GCF_002006125.1_ASM200612v1/GCF_002006125.1_ASM200612v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/685/985/GCF_000685985.1_ASM68598v1/GCF_000685985.1_ASM68598v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/785/GCF_000007785.1_ASM778v1/GCF_000007785.1_ASM778v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/295/525/GCF_004295525.1_HPATCC43504_01/GCF_004295525.1_HPATCC43504_01_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/287/905/GCF_002287905.1_ASM228790v1/GCF_002287905.1_ASM228790v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/035/GCF_000196035.1_ASM19603v1/GCF_000196035.1_ASM19603v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1/GCF_000016525.1_ASM1652v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/805/GCF_000008805.1_ASM880v1/GCF_000008805.1_ASM880v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/472/405/GCF_000472405.1_ASM47240v1/GCF_000472405.1_ASM47240v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/405/GCF_000016405.1_ASM1640v1/GCF_000016405.1_ASM1640v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/645/GCF_000007645.1_ASM764v1/GCF_000007645.1_ASM764v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/265/GCF_000007265.1_ASM726v1/GCF_000007265.1_ASM726v1_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/465/GCF_000007465.2_ASM746v2/GCF_000007465.2_ASM746v2_genomic.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/465/GCF_000007465.2_ASM746v2/GCF_000007465.2_ASM746v2_genomic.fna.gz" > tmp.list

while read a b; do echo $a; curl --remote-name --remote-time $a; done < tmp.list
rm -rf tmp.list

zcat *.gz | pigz --best > repcont.fna.gz

