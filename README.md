# FACS pipeline - Fast AMP Clustering System

**NOTE**: This is still a _work in progress_ and, while the results of the tool
should be correct, we are still working on making FACS easier to install and
use.

Fast AMP Clustering System pipeline is a system created by Celio Dias Santos
Jr. and Luis Pedro Coelho, from Fudan University (Shanghai / CN). It is
distributed under MIT license and represents a new way to prospect AMPs in natural environments using metagenomic data or genomic data to generate large datasets of antimicrobial peptides.

(A modified version of the Prodigal software is distributed with FACS under the GPL license)

## Applications

FACS can be used in a wide-ranging of scenarios, such as screening for novel AMPs generating candidates to further testing and patenting, as well as, determination of microbiome quorum sensing mechanisms linking AMPs to health conditions or presence of diseases.

## Installing

Before installation, make sure your system settings are as follows:

1. Linux (preferably Ubuntu 64 bits, version 18+);

2. Dependencies

		- apt

		- git

		- Python Version 3.0 or above

		- The following packages should be already installed in Python environment:
		  sys, os, shutil, scipy, argparse, collections, platform,
		  math, re, numpy (1.13.1), sklearn (0.19.1), matplotlib (2.1.0),
		  pandas (0.20.1).

		- R Version 3.5.2 or above
		
		- Conda/MiniConda


### Third-party software

Also, before start installation, make sure you know the needed third party software list:

1. To quality trimming of reads and selection/sorting of paired-end reads, it is used [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

   - As a complement: to work with Trimmomatic, you will need OpenJDK-11-JDK-headless.
   
2. To assembly reads, we use [MEGAHIT v.1](https://www.ncbi.nlm.nih.gov/pubmed/25609793)

3. To predict smORFs, we use a modified version of the [Prodigal](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/) program.

4. To produce descriptors and use the randomforest models, we use the following R packages:

     - randomForest
     - caret
     - Peptides
     - data.table
     - dplyr
     - parallel
     - doParallel

    4.1. To calculate distribution parameters we used scripts from the project [iLearn](https://github.com/Superzchen/iLearn):

     - CTDDClass.py
     - saveCode.py
     - readFasta.py
    
5. The library FAST from CPANM to speed up Perl.

6. The GNUParallel library to speed up the script.
   - Additional libraries needed can include: zlib1g, zlib1g-dev, and libpthread-stubs0-dev

7. To speed up the compressing and decompressing of files, we used [pigz](https://zlib.net/pigz/) software.



### Install process

>                                  WARNING:                                 
> To install FACS you will need CONDA or miniCONDA installed in your system 

Proceed with the installation by first cloning this repository:

`$ git clone https://github.com/FACS-Antimicrobial-Peptides-Prospection/FACS/`

Then decompress the file:

`$ gunzip FACS-master.gz
$ cd FACS-master/`

And execute the installation script:

`$ sh install.sh`


## Usage

There are few options to make the running of the program a bit customized and speed up process according to the settings of the system available.

    Usage: FACS.sh --mode c/r/p/a [options]
    
    Here's a guide for avaiable options. Defaults for each option are showed inside the
    brackets: [default].
    
    Basic options:
    
    -h, --help            Show this help message
    
    -m          Mode of operation, type:
                "c" to work with contigs,
                "p" to predict AMPs directly from a peptides FASTA file,
                "r" to work with reads, 
                "a" to map reads against AMP output database and generate abundances table
		
    --fasta               Compressed (or not gzipped) contigs or peptides fasta file
    
    --fwd                 Illumina sequencing file in Fastq format (R1), please leave it
                          compressed and full path
    
    --rev                 Illumina sequencing file in Fastq format (R2), please leave it
                          compressed and pass the full path
    
    --ref                 Output of module "c" in its raw format [file type tsv and compressed]
    
    --outfolder           Folder where output will be generated [Default: .]
    
    --outtag              Tag used to name outputs [Default: FACS_OUT]
    
    -t, --threads [N]     Number of threads [Default: 90% of available threads]
    
    --block               Bucket size (take in mind it is measured in Bytes and also it 
                          influences memory usage). [Defeault: 100MB]
    
    --log                 Log file name. FACS will save the run results to this log file
                          in output folder.
    
    --mem                 Memory available to FACS ranging from 0 - 1. [Defult: 0.75]
    
    --tmp                 Temporary folder (address)
    
    --cls                 Cluster peptides: yes (1) or no (0). [Default: 1 - yes]


### Quick running commands:

> FACS supports gzipped inputs

To run FACS on peptides, type directly:

```bash FACS.sh -m p --fasta /path/to/file.fa --outfolder /path/to/output_folder --outtag name_to_your_results -t #_of_cpus --tmp tempdir_name```

To run FACS on contigs, type directly:

```bash FACS.sh -m c --fasta /path/to/file.fa --outfolder /path/to/output_folder --outtag name_to_your_results -t #_of_cpus --tmp tempdir_name```

To run FACS on paired-end reads, type directly:

```bash FACS.sh -m r --fwd /path/to/R1_file.fq.gz --rev /path/to/R2_file.fq.gz --outfolder /path/to/output_folder --outtag name_to_your_results -t #_of_cpus --tmp tempdir_name```

To run FACS on single-end reads, type directly:

```bash FACS.sh -m r --fwd /path/to/R1_file.fq.gz --outfolder /path/to/output_folder --outtag name_to_your_results -t #_of_cpus --tmp tempdir_name```

To run FACS to get abundance profiles, you only need the forward reads file and a reference with peptide sequences. You can run abundance mode using two different types of references:

1. Fasta file with peptide sequences:

```bash FACS.sh -m a --fwd /path/to/R1_file.fq.gz --fasta /path/to/file.fa --outfolder /path/to/output_folder --outtag name_to_your_results -t #_of_cpus --tmp tempdir_name```

2. FACS output with predicted AMPs, it will be the "$outtag.tsv.gz" file generated. To use this file containing the peptide sequences and the probabilities associated to the prediction, you can run the following commmand:

```bash FACS.sh -m a --fwd /path/to/R1_file.fq.gz --ref /path/to/FACS_output.tsv.gz --outfolder /path/to/output_folder --outtag name_to_your_results -t #_of_cpus --tmp tempdir_name```

We have provided to you example files that can help you to test your installed pipeline:
     
     - R1.fq.gz, R2.fq.gz to reads mode
     - excontigs.fa.gz to contigs mode
     - expep.fa.gz to peptides mode
     
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     +++++                    NOTE:                        ++++
     +++++ You also can use expep.fa.gz and R1.fq.gz files ++++
     +++++ to test the abundance mode.                     ++++
     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Pipeline overview

In general, FACS uses the pigz software to compress and decompress files quickly, as well as GNUParallel to make the shell processes faster and parallelized. The 22 descriptors adopted by FACS are hybrid comprising local and global contexts to do the sequence encoding. FACS performs firstly a distribution analysis of three classes of residues in two different features (Solvent accessibility and *Free energy to transfer from water to lipophilic phase*) as shown in Table 1. The novelty in this method is using the *Free energy to transfer from water to lipophilic phase* (FT) firstly described by Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x) to capture the spontaneity of the conformational change that AMPs suffer while their transference from water to the membrane. For more info about the other descriptors used in FACS and the algorithms used to train the classifiers, please refer to the FACS reference.

**Table 1.** Classes adopted to the sequence encoding of the distribution at the first residue of each class. The Solvent Accessibility was adopted as in previous studies (Dubchak et al. [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/), [1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667)), however, the new feature FT was adapted from Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x).

| Properties     | Class I     | Class II     | Class III     |
|:--------------------------------------------------------:    |:------------------------:    |:------------------    |------------------:    |
| Solvent accessibility     | A, L, F, C, G, I, V, W     | R, K, Q, E, N, D     | M, S, P, T, H, Y     |
| FT     | I,L,V,W,A,M,G,T     | F,Y,S,Q,C,N     | P,H,K,E,D,R     |

The other descriptors used in FACS classifiers are widely used in the AMPs description, as follows:

 - tinyAA (A + C + G + S + T)
 - smallAA (A + B + C + D + G + N + P + S + T + V)
 - aliphaticAA (A + I + L + V)
 - aromaticAA (F + H + W + Y)
 - nonpolarAA (A + C + F + G + I + L + M + P + V + W + Y)
 - polarAA (D + E + H + K + N + Q + R + S + T + Z)
 - chargedAA (B + D + E + H + K + R + Z)
 - basicAA (H + K + R)
 - acidicAA (B + D + E + Z)
 - charge (pH = 7, pKscale = "EMBOSS")
 - pI (pKscale = "EMBOSS")
 - aindex (relative volume occupied by aliphatic side chains - A, V, I, and L)
 - instaindex -> stability of a protein based on its amino acid composition
 - boman -> overall estimate of the potential of a peptide to bind to membranes or other proteins as receptor
 - hydrophobicity (scale = "KyteDoolittle") -> GRAVY index
 - hmoment (angle = 100, window = 11) -> quantitative measure of the amphiphilicity perpendicular to theaxis of any periodic peptide structure, such as the alpha-helix or beta-sheet

We opted to use random forests after some tests with alternative algorithms. The training of AMPs classifier used the training and validation data sets used by Bhadra et al. ([2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9)), while the classifier of hemolytic peptides was trained and tested with the data sets previously established by Chaudhary et al. ([2016](https://www.nature.com/articles/srep22843)).

The models here mentioned were implemented in an R script to filter off the non-AMP peptides and classify AMPs into hemolytic or not. After that, this script also submits the predicted AMPs to a decisions tree, classifying AMPs into 4 families:

 - Anionic linear peptides (ALP)
 - Anionic disulfide-bond forming peptides (ADP)
 - Cationic linear peptides (CLP)
 - Cationic disulfide-bond forming peptides (CDP)

Benchmark procedures showed that R22 models are efficient in retrieving AMPs with statistics that overcome the state-of-art methods (Table 2).

**Table 2.** Comparison of FACS and other state-of-art AMP prediction systems. All systems were tested with the benchmark data set from the AMPEP study available in Bhadra et al. ([2018](doi:10.1038/s41598-018-19752-w)).

| Model/Method     | Acc.     | Sp.     | Sn.     | Precision     | MCC     | Refererence     |
|:--------------------:    |:-----:    |:-----:    |:-----:    |:---------:    |:-------:    |:-----:   |
| AMPep     | 0.962     | 0.965     | 0.95     | 0.913     | 0.9     | doi:10.1038/s41598-018-19752-w     |
| Peptide Scanner v2     | 0.755     | 0.686     | 0.943     | 0.523     | 0.557     | doi: 10.1093/bioinformatics/bty179     |
| CAMPr3-NN     | 0.728     | 0.704     | 0.794     | 0.494     | 0.445     | doi:  10.1093/nar/gkv1051     |
| CAMPr3-RF     | 0.584     | 0.461     | 0.923     | 0.384     | 0.354     | doi:  10.1093/nar/gkv1051     |
| CAMPr3-SVM     | 0.6     | 0.506     | 0.858     | 0.388     | 0.328     | doi:  10.1093/nar/gkv1051     |
| CAMPr3-DCA     | 0.617     | 0.542     | 0.821     | 0.396     | 0.324     | doi:  10.1093/nar/gkv1051     |
| iAMP     | 0.548     | 0.413     | 0.918     | 0.363     | 0.313     | doi: 10.1038/srep42362     |
| AMPA     | 0.779     | 0.941     | 0.336     | 0.675     | 0.361     | doi: 10.1093/bioinformatics/btr604     |
| R22_Ctrained     | 0.952     | 0.972     | 0.932     | 0.971     | 0.904     | This study     |
| R22_LargeTrainingset     | 0.967     | 1     | 0.934     | 1     | 0.936     | This study     |

Meanwhile, the hemolytic prediction model implemented in FACS has a comparable performance of the state-of-art methods previously tested by Chaudhary et al., [2016](https://www.nature.com/articles/srep22843) as shown in Table 3.

**Table 3.** Comparison of the performance of different hemolytic activity prediction systems. All the systems were trained and benchmarked with the same data sets (HemoPI) used by Chaudhary et al., [2016](https://www.nature.com/articles/srep22843).

| Methods     | Study     | Sn (%)     | Sp (%)     | Acc (%)     | MCC     |
|:---------------------:    |:----------------------:    |:------:    |:------:    |:-------:    |:----:    |
| SVM     | Chaudhary et al., 2016     | 95.7     | 94.8     | 95.3     | 0.91     |
| IBK     | Chaudhary et al., 2016     | 95.5     | 93.7     | 94.6     | 0.89     |
| Multilayer Perceptron     | Chaudhary et al., 2016     | 93.9     | 92.8     | 93.3     | 0.87     |
| Logistic     | Chaudhary et al., 2016     | 93.4     | 93.7     | 93.6     | 0.87     |
| J48     | Chaudhary et al., 2016     | 89.6     | 88.5     | 89.0     | 0.78     |
| Random Forest     | Chaudhary et al., 2016     | 94.1     | 94.6     | 94.3     | 0.89     |
| Random Forest    | This study     | 95.4     | 93.8     | 94.5     | 0.91     |
