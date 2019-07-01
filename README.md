# FACS pipeline - Fast AMP Clustering System

FACS is a pipeline to:

1. merge paired-end reads,
2. predict peptides,
3. cluster them at 100% of similarity and 100% coverage,
4. calculate their abundance in peptides per million (ppm), and
5. select those with antimicrobial potential discriminating their hemolytic pontential.

With FACS you can treat a metagenome file of 631.2Mbp as fast as 24 min, using 3 cpus and 100Mb sequence buckets in a Ubuntu v.18 64x bits.

## History

In order to creat a simple 

## Installing

The installation can be performed with downloading the scripts as:

`$ git clone https://github.com/celiosantosjr/FACS`

Performing the decompression:

`$ gunzip FACS-master.gz`

And executing the installation script:

`$ sh install.sh`

## Usage

Basically it can be run using the command line in bash:

FACS.sh "[options]" --fwd <R1.file.gz> --rev <R2.file.gz>

There are few options to make the running of the program a bit customized and speed up process according to the systems settings available.


    Here's a guide for avaiable options. Defaults for each option are showed between brackets: [default].

    Basic options:
    -h, --help	        Show help page
    --fwd               Illumina sequencing file in Fastq format (R1), please leave it compressed and full adress
    --rev		Illumina sequencing file in Fastq format (R2), please leave it compressed and full address
    --outfolder		Folder where output will be generated [./]
    --outtag          	Tag used to name outputs [OUT]
    -t, --threads [N]	Number of threads [90% of avaiable threads]
    --block		Bucket size (take in mind it is measured in bits and also it determines the memory usage). [100MB]
    --adapters		Adapters currently available in Trimmomatic program, if not available you can create one accordingly
		        to its manual. [NexteraPE-PE.fa]
    --log		Save results of FACS run to a log file in output folder.
    
## FACS merger

FACS merger was designed to join results of each different metagenome and return a table of abundances of each detected peptide in ppm to each metagenome.

Usage: FACS_merger.sh "[options]"
	
	
        Basic options:
	-h, --help	        Shows help message	
	--output [file]	        File where output is sent, it needs to be gzipped (ending em .gz)
	-t, --threads [N]	Number of threads [90% of avaiable threads]
	--reference [folder]	Folder where your reference files are located, if none current folder will be used [Reference_seqs]

## Third party softwares

1. To quality trimming of reads and paired-end reads selection and sorting it is used [Trimmomatic]<http://www.usadellab.org/cms/?page=trimmomatic>

  1.1. As complementary to working Trimmomatic needs openjdk-11-jdk-headless.

2. To reads merging it is used [pandaseq]<https://github.com/neufeld/pandaseq> software.

3. To ORFs prediction it is used [ORFm]<https://github.com/wwood/OrfM> thought to be faster than other ORF prediction systems.

4. To produce descriptors and use the AI models to select peptides, we have used the following R packages:
 
  *randomForest
  *caret
  *Peptides
  *data.table
  *dplyr
  *parallel
  *doParallel

5. The library FAST from CPANM to speed up perl.

6. The GNUParallel library to speed up the script.
  
  6.1. Additional libraries needed can include: zlib1g, zlib1g-dev and libpthread-stubs0-dev

7. The [pigz]<https://zlib.net/pigz/> software to speed up the compressing and decompressing of files.

8. The following scripts from the project [iLearn]<https://github.com/Superzchen/iLearn> to calculate some encodings of the sequences:
  
  *CTDDClass.py
  *saveCode.py
  *readFasta.py
