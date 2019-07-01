# FACS pipeline - Fast AMP Clustering System

FACS is a pipeline to:

1. merge paired-end reads,
2. predict peptides,
3. cluster them at 100% of similarity and 100% coverage,
4. calculate their abundance in peptides per million (ppm), and
5. select those with antimicrobial potential discriminating their hemolytic pontential.

With FACS you can treat a metagenome file of 631.2Mbp as fast as 24 min, using 3 cpus and 100Mb sequence buckets in a Ubuntu v.18 64x bits.

## History

Since the discovery of penicillin and its use in the 1940's, the antibiotics resistance developed by the microorganisms during the medical treatment is a recurrent problem in modern medicine. The fact is to each new discovered antibiotic an emerging resistance trait raises right after six months after their release in market. Part of this is due to the huge metabolism diversity presented by prokaryotes as can be observed in Figure 1. There are two main mechanisms of resistance, one obtained via vertical transference; and the other is action of genes in mobile elements, transmitted both vertically and horizontally to other bacteria. These mobile genetic elements such as plasmids, can carry one or more resistance genes. The prevalent and extremelly quick mobility of resistance genes in previously sensitive bacterial populations, now established an world crisis.

![](https://github.com/celiosantosjr/FACS/blob/master/fig1.png)

**Figure 1.** The antibiotics mechanisms and their overcoming (Extracted from: Wright, [2010](https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-123)).

The superbugs risen is faster than the time it takes to develop new antibiotics (Figure 2). However, many of these new antibiotics are just chemical modifications of the molecular structure of the old compounds, which makes them prone to be skipped by bacteria using slightly modificated strategies. But the question is "can anything be done to slow down the emergence of resistance?". Antibiotics represent an evolutionary pressure that eventually is the reason to them become obsolete. So, reducing the exposure of microbes to antibiotics can reduce the opportunity for selection and dissemination of resistance. Despite initiatives such as those taken by European Union and in North America, foccus mainly in surveillance and restriction of use. However, these measures are only able to delay the emergence of antibiotcs resistance. Thus, those strategies are welcome, but new drugs will always be needed, since resistance risen is inevitable.

![](https://github.com/celiosantosjr/FACS/blob/master/fig2.png)

**Figure 2.** Superbugs running against the farmaceutical companies in the antibiotics development and overcoming ([Source](http://2014hs.igem.org/Team:Lethbridge_Canada/project)).

The developing world (Figure 3) usually does not regulate the access to antibiotics and their use in widespread since from agriculture to daily life. This makes antibiotic stewardship an important death cause in those countries, besides the statistics are usually understimated, since the report in those cases is neglected by many health attendants or the diagnosis is not completed before the patient's death. Rapid intercontinental travels also are efficient to bring pathogens that are no longer geographically contained and can transpose countries easily, like cases involving the spread of the severe acute respiratory syndrome (SARS) virus from Guangdong province in China to Hong Kong and then Canada in 2003. 

![](https://github.com/celiosantosjr/FACS/blob/master/fig21.png)

**Figure 3.** Deaths attributable to antimicrobial resistance every year by 2050 (Source: O'Neil, [2014](https://amr-review.org/sites/default/files/AMR%20Review%20Paper%20-%20Tackling%20a%20crisis%20for%20the%20health%20and%20wealth%20of%20nations_1.pdf))

In a prevision using the growing trend of some diseases that are known to be highly mortal, the resistance to antimicrobials overcome them in 2050, becoming more mortal than cancer (Figure 4). This shows the global importante of the mater, besides to evidence the main countries (Figure 3) affected by this. China is among of the most affected countries and the number of deaths can be higher than 4 million people per year.

![](https://github.com/celiosantosjr/FACS/blob/master/fig3.png)

**Figure 4.** Deaths attributable to antimicrobial resistance every year compared to other major causes of death. (Source: O'Neil, [2014](https://amr-review.org/sites/default/files/AMR%20Review%20Paper%20-%20Tackling%20a%20crisis%20for%20the%20health%20and%20wealth%20of%20nations_1.pdf))

Other important side of this problem is the economy. Antibiotic resistance can substantially reduce gross domestic product - but unlike a financial crisis, the damage will last longer and will be greater in the poorest countries (Figure 5). The question here is about the costly treaments and the hospital expenses during treatments of resistance, that can be exhausting and take several months until solved or the patients death.

![](https://github.com/celiosantosjr/FACS/blob/master/fig4.png)

**Figure 5.** The economy of treatments of antimicrobial resistance cases. (Source: World Bank Group)

So, why new drugs are not discovered? 



![](https://github.com/celiosantosjr/FACS/blob/master/fig5.png)

**Figure 6.** The economy of treatments of antimicrobial resistance cases. (Source: World Bank Group)



**Table 1.** Classess adopted to the sequence encoding of the distribution at Residue0. The Solvent Accessibility was adopted as other studies previously Dubchak et al. [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/),[1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667), however the new feature "Free energy to transfer to lipophilic phase" was adopted from Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x).

| Properties 	| Class I 	| Class II 	| Class III 	|
|--------------------------------------------------------	|------------------------	|------------------	|------------------	|
| Solvent accessibility 	| A, L, F, C, G, I, V, W 	| R, K, Q, E, N, D 	| M, S, P, T, H, Y 	|
| Free energy to transfer from water to lipophilic phase 	| I,L,V,W,A,M,G,T 	| F,Y,S,Q,C,N 	| P,H,K,E,D,R 	|


**Table 2.** Models to predict antimicrobial peptides were tested to benchmark the results obtained with the new set of descriptors adopted in this classifier. All models and systems were tested with the benchmark validation dataset from the AMPEP study available [elsewhere](doi:10.1038/s41598-018-19752-w).

| Model/Method 	| Acc. 	| Sp. 	| Sn. 	| Precision 	| F-Score 	| MCC 	| Refererence 	| Descriptors 	|
|:--------------------:	|:-----:	|:-----:	|:-----:	|:---------:	|:-------:	|:-----:	|:----------------------------------:	|:-----------:	|
| AMPep 	| 0.962 	| 0.965 	| 0.95 	| 0.913 	| - 	| 0.9 	| doi:10.1038/s41598-018-19752-w 	| 23 	|
| Peptide Scanner v2 	| 0.755 	| 0.686 	| 0.943 	| 0.523 	| 0.673 	| 0.557 	| doi: 10.1093/bioinformatics/bty179 	| - 	|
| CAMPr3-NN 	| 0.728 	| 0.704 	| 0.794 	| 0.494 	| 0.609 	| 0.445 	| doi:  10.1093/nar/gkv1051 	| - 	|
| CAMPr3-RF 	| 0.584 	| 0.461 	| 0.923 	| 0.384 	| 0.543 	| 0.354 	| doi:  10.1093/nar/gkv1051 	| - 	|
| CAMPr3-SVM 	| 0.6 	| 0.506 	| 0.858 	| 0.388 	| 0.534 	| 0.328 	| doi:  10.1093/nar/gkv1051 	| - 	|
| CAMPr3-DCA 	| 0.617 	| 0.542 	| 0.821 	| 0.396 	| 0.534 	| 0.324 	| doi:  10.1093/nar/gkv1051 	| - 	|
| iAMP 	| 0.548 	| 0.413 	| 0.918 	| 0.363 	| 0.521 	| 0.313 	| doi: 10.1038/srep42362 	| - 	|
| AMPA 	| 0.779 	| 0.941 	| 0.336 	| 0.675 	| 0.449 	| 0.361 	| doi: 10.1093/bioinformatics/btr604 	| - 	|
| R22_Ctrained 	| 0.952 	| 0.972 	| 0.932 	| 0.971 	| 0.951 	| 0.904 	| This study 	| 22 	|
| R22_LargeTrainingset 	| 0.967 	| 1 	| 0.934 	| 1 	| 0.966 	| 0.936 	| This study 	| 22 	|

**Table 3.** Models used to predict hemolytic activity were trained with the same dataset used by Chaudhary et al. [2016](https://www.nature.com/articles/srep22843) and were tested with the same test dataset used by them to benchmark the results obtained with another model used to generate this classifier.

| Methods 	| Study 	| Sn (%) 	| Sp (%) 	| Acc (%) 	| MCC 	|
|:---------------------:	|:----------------------:	|:------:	|:------:	|:-------:	|:----:	|
| SVM 	| Chaudhary et al., 2016 	| 95.7 	| 94.8 	| 95.3 	| 0.91 	|
| IBK 	| Chaudhary et al., 2016 	| 95.5 	| 93.7 	| 94.6 	| 0.89 	|
| Multilayer Perceptron 	| Chaudhary et al., 2016 	| 93.9 	| 92.8 	| 93.3 	| 0.87 	|
| Logistic 	| Chaudhary et al., 2016 	| 93.4 	| 93.7 	| 93.6 	| 0.87 	|
| J48 	| Chaudhary et al., 2016 	| 89.6 	| 88.5 	| 89.0 	| 0.78 	|
| Random Forest 	| Chaudhary et al., 2016 	| 94.1 	| 94.6 	| 94.3 	| 0.89 	|
| ORFsvm 	| This study 	| 95.5 	| 95.5 	| 95.5 	| 0.91 	|


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

1. To quality trimming of reads and paired-end reads selection and sorting it is used [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   - As complementary to working Trimmomatic needs openjdk-11-jdk-headless.

2. To reads merging it is used [pandaseq](https://github.com/neufeld/pandaseq) software.

3. To ORFs prediction it is used [ORFm](https://github.com/wwood/OrfM) thought to be faster than other ORF prediction systems.

4. To produce descriptors and use the AI models to select peptides, we have used the following R packages:
 
 - randomForest
 - caret
 - Peptides
 - data.table
 - dplyr
 - parallel
 - doParallel

5. The library FAST from CPANM to speed up perl.

6. The GNUParallel library to speed up the script.
   - Additional libraries needed can include: zlib1g, zlib1g-dev and libpthread-stubs0-dev

7. The [pigz](https://zlib.net/pigz/) software to speed up the compressing and decompressing of files.

8. The following scripts from the project [iLearn](https://github.com/Superzchen/iLearn) to calculate some encodings of the sequences:
  
 - CTDDClass.py
 - saveCode.py
 - readFasta.py
