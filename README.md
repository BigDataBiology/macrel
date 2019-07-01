# FACS pipeline - Fast AMP Clustering System

Fast AMP Clustering System pipeline is a system created by Celio Dias Santos Jr. and Luis Pedro Coelho, from Fudan University (Shanghai / CN). It is distributed under MIT license and represents a new way to prospect AMPs in natural environments using metagenomic data or genomic data to generate large datasets of antimicrobial peptides.

## Background

Since the discovery of penicillin and its use in the 1940's, the antibiotics resistance developed by the microorganisms during the medical treatment is a recurrent problem in modern medicine. The fact is to each new discovered antibiotic an emerging resistance trait raises right after six months after their release in market. Part of this is due to the huge metabolism diversity presented by prokaryotes as can be observed in Figure 1. There are two main mechanisms of resistance, one obtained via vertical transference; and the other is action of genes in mobile elements, transmitted both vertically and horizontally to other bacteria. These mobile genetic elements such as plasmids, can carry one or more resistance genes. The prevalent and extremelly quick mobility of resistance genes in previously sensitive bacterial populations, now established an world crisis.

![](https://github.com/celiosantosjr/FACS/blob/master/fig1.png)

**Figure 1.** The antibiotics mechanisms and their overcoming (Source: Wright, [2010](https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-8-123)).

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

The growing antimicrobial resistance problem feed a continuous need for new antibiotic drugs, and there are a number of reasons for the scarcity of new antibiotics. Some of them include government regulatory approval adding risk for the pharmaceutical industry, what is usually too criterious, since they will be taken by patients over short periods of time only to cure the disease. Other reason for antibiotic discovery and development decline is scientific. Few compounds discovered with antibiotic properties have had the requisite properties to become drugs. Researchers have argumented that most antibiotics are natural products isolated from soil bacteria, which could suggest the exhaustion of this source now. Many of these 'natural' antibiotics have desirable drug-like qualities: good bioavailability, they can cross the cell membrane and have the ability to evade efflux systems, and chemical structures that favor binding to vital cellular targets. However, there is an increasing difficulty of identifying new chemical compounds with equally suitable drug-like characteristics from natural sources which has caused natural-product-based screening programs diseappear in the past few decades. However, the advantages of synthetic compounds are clear to industry, after decades of emphasis on such molecules and millions of dollars spent, no new synthetic antibiotics have emerged. 

The genomic era was constrasted by the reality of hundreds of available bacterial genomes that have so far failed to deliver the hoped-for new molecular targets for antibiotics. However, so far it always have focused in the active molecules produced by the metabolism, instead searching for active peptides or proteins. The best reason to bet in host defense antimicrobial peptides or AMPs is that they remained potent for millions of years, constituting a useful strategy to develop a new generation of antimicrobials to meet the growing antibiotic resistance problem worldwide. The current informations about AMPs is extended in that regarding eukaryotes' peptides (Figure 6), and their presence in several phyla in that domain. Although well known in eukaryotes, prokaryotes remain under represented and the few information available does not reflect the entire diversity present in that domain. Archaea is another few explored domain, that can contribute to future drugs development. 

![](https://github.com/celiosantosjr/FACS/blob/master/fig5.png)

**Figure 6.** Number of antimicrobial peptides found in different domains of life. (Source: Wang, [2014](https://link.springer.com/protocol/10.1007/978-1-4939-2285-7_3))

AMPs definition comprises peptides with a huge variety of biological activities (Figure 7), and their sequences are key to that activity. AMP producing microbes can limit the growth of other microorganisms and should be considered another normal source of them. AMPs from microbes are quite distinct from those of vertebrates, since they can be obtained from a nonribosomal peptide synthase. Thus, nonribosomal peptides can adopt different structures, such as cyclic or branched structures, and carry modifications like N-methyl and N-formyl groups, glycosylations, acylations, halogenation, or hydroxylation. Some examples of commercial microbial AMPs include polymyxin B and vancomycin, both FDA-approved antibiotics (Zhang and Gallo, [2016](https://www.sciencedirect.com/science/article/pii/S0960982215014098)).

![](https://github.com/celiosantosjr/FACS/blob/master/fig6.png)

**Figure 7.** Antimicrobial peptide biological activities. (Source: Wang, [2014](https://link.springer.com/protocol/10.1007/978-1-4939-2285-7_3))

Most AMPs are peptides 10-50 amino acids long, also being longer until 100 amino acids in some cases, with charge ranging between 2 and 11 (some of them being anionic) and consituted of approximately 50% of hydrophobic residues (Zhang and Gallo, [2016](https://www.sciencedirect.com/science/article/pii/S0960982215014098)). There is a pronounced pH-dependent AMPs charge, mostly resulting in membrane lysis and antibacterial activity at acidic conditions, with many of them not presenting activity at pH higher than 6.0. Thus, the charge seems a key feature in the interaction of AMPs and membranes, where its distribution and nature along the sequence changes the antimicrobial activity (Malmsten, [2014](https://www.tandfonline.com/doi/full/10.3109/03009734.2014.899278); Pasupuleti et al., [2012](www.ncbi.nlm.nih.gov/pubmed/22074402?dopt=Abstract); Ringstad et al., [2006](https://www.ncbi.nlm.nih.gov/pubmed/16700592?dopt=Abstract)). Furthermore, the formation of amphiphilic ordered structures is correlated to peptide-induced membrane disruption. These structures induction, mostly alpha-helices, works as a driving force for membrane binding. Also, the helix destabilization oftenly can reduce the cytotoxicity of AMPs, although this can result in reduction of antimicrobial effects (Malmsten, [2014](https://www.tandfonline.com/doi/full/10.3109/03009734.2014.899278); Borgden, [2005](https://www.ncbi.nlm.nih.gov/pubmed/15703760?dopt=Abstract); Pasupuleti et al., [2012](www.ncbi.nlm.nih.gov/pubmed/22074402?dopt=Abstract); Hancok and Sahl, [2006](https://www.ncbi.nlm.nih.gov/pubmed/17160061?dopt=Abstract); Shai, [2002](https://www.ncbi.nlm.nih.gov/pubmed/12491537?dopt=Abstract); Stromstedt et al., [2006](https://www.ncbi.nlm.nih.gov/pubmed/19029324?dopt=Abstract)).

AMPs can be classified into 5 families accordingly to their origin and composition (Perumal et al., [2013](http://xueshu.baidu.com/usercenter/paper/show?paperid=6033547256f3e45f884306a14bbff34c&site=xueshu_se)):

	1. Anionic peptides: rich in aspartic and glutamic acids. Example: Maxinimin H5 (from amphibians);

	2. Linear alpha-helical cationic peptides: Lack in cysteine. Example: Cecropins (from insects), dermaceptin (from amphibians);

	3. Cationic peptides: rich in proline, arginine, phenylalanine, glycine and thryptophan. Example: Indolicidin (from cattle), prophenin (from frogs);

	4. Anionic and cationic peptides that contain dissulphide bonds: contain cysteine. Examples: 1 dissulphide bond (brevinins), 2 dissulphide bonds (protregrin) and 3 dissulphide bonds (drosomycins and defensins);

	5. Anionic and cationic peptide fragments of larger proteins: contains unusual amounts of tryptophan, lysine, valine, arginine, proline, histidine, leucine. Examples: Haemoglobin (from humans), lysozyme, ovoalbumin and lactoferricin from lactoferrin.

Families 1-3 are largely found in all domains of life, while families 4 and 5 are more related to eukaryotes and their contribution from microbes is very few (Perumal et al., [2013](http://xueshu.baidu.com/usercenter/paper/show?paperid=6033547256f3e45f884306a14bbff34c&site=xueshu_se)).

As previously mentioned, these families can be folded into some structural arrangements (Wang, [2014](https://link.springer.com/protocol/10.1007/978-1-4939-2285-7_3)). The most common are shown in Figure 8. The alpha-helical peptides (Figure 8-a) are usually related to a strong pore-forming activity in bacterial membranes, as well as, the alpha-beta structures (Figure 8-c). The beta-sheet peptides usually change their conformation in apolar environments to an alpha-helical structure that can be refolded into beta-sheet (Figure 8-b) after transposition of the lipophilic phase. The random coiled peptides (Figure 8-d) are usually associated to a mixed function, and usually assume helical structures in the membrane, forming pores and compromising cell functions. Mostly the activities of AMPs are associated to the rupture of cell membrane or promoting the leakage of cell contents, ending in the bacterial cell death. Different biological activites have different mechanisms, however in this review the antibacterial activity will be prioritized.

![](https://github.com/celiosantosjr/FACS/blob/master/fig8.png)

**Figure 8.** Antimicrobial peptide folding groups. (Source: Wang, [2014](https://link.springer.com/protocol/10.1007/978-1-4939-2285-7_3))

There is a dynamic interchange in AMPs structure and topologies along the interaction with the microbial cell membranes (Samson, [1998](https://www.sciencedirect.com/science/article/pii/S1359029498800277)). The outer surface of prokaryotic cells is negatively charged (mainly due to lipopolysaccharides and teichoic acid), what promotes an electrostatic interaction of AMPs with the membrane being the primary mechanism for antimicrobial activity. Other cases are based in the AMP translocation across the cell membrane and the inhibition of essential cellular processes (e.g. protein synthesis, nucleic acid synthesis, enzymatic activities) (Brogden, [2005](https://www.ncbi.nlm.nih.gov/pubmed/15703760)). Based on the mechanisms of action, AMPs are categorized into membrane acting and nonmembrane acting peptides. The first ones are capable of forming transient pores on the membrane, whereas the second ones have the ability to translocate across the cell membrane without permeabilizing it (Pushpanathan et al., [2013](http://dx.doi.org/10.1155/2013/675391)).

Several models have been proposed to describe the mechanism of action of antimicrobial peptides (Figure 9), and can be categorized into energy dependent and energy independent uptake. In barrel-stave mechanism, there is an aggregation of peptide monomers on the surface of the membrane. This aggregated peptides are inserted into the membrane and get such an orientation that the hydrophilic surfaces of peptides point inward and form a water filled transmembrane pore that kills the cell by leakage. In carpet model, AMPs initially get associated on the surface of the membrane, forming a carpet. Once a concentration reaches a threshold, there is a peptide induced membrane permeation. This leads to the cell membrane disruption. In toroidal pore model, peptides get aggregated prior or after binding with the membrane surface. It induces a membrane depolarization and form a toroidal shaped transmembrane pore. The energy independent uptake involves macropinocytosis. Once uptaken in the form of macropinosomes, the AMPs get released into the cytoplasm exerting their antimicrobial action (Pushpanathan et al., [2013](http://dx.doi.org/10.1155/2013/675391)).

![](https://github.com/celiosantosjr/FACS/blob/master/fig22.png)

**Figure 9.** Proposed mechanisms of actions of AMPs: Energy indpendent mechanisms - barrel stave model, carpet model, and toroidal pore model (a); and energy dependent mechanisms (b). (Source: Pushpanathan et al., [2013](http://dx.doi.org/10.1155/2013/675391))

Microbes were thought to be unable to develop resistance towards AMPs. However, recently some resistance mechanisms have been reported, such as upregulation of proteolytic enzymes able to degrade
AMPs, membrane modifications resulting in decreased negative potential of bacterial membranes, and release of glucose aminoglycans, polysaccharides, and other polyanionic species able to
scavenge AMPs (Nizet, [2006](https://www.ncbi.nlm.nih.gov/pubmed/16450883?dopt=Abstract)). Despite having been convincingly demonstrated *in vitro*, resistance development to AMPs *in vivo* needs to be further clarified, since conditions experienced by bacteria in a laboratory setting are likely to differ from those *in vivo*. In the latter case, the microbes are exposed to a cocktail of AMPs, which may reduce or alter the selection pressure underlying resistance development (Malmsten, [2014](https://www.tandfonline.com/doi/full/10.3109/03009734.2014.899278)). 

In summary, the AMPs represent a multidimensional group of molecules with several applications, among them:

 - Drug delivery vectors
 - Mitogenic agent
 - Antitumour agent
 - Signaling molecules
 - Contraceptive agent for vaginal prophylaxis
 - Plant Transgenesis



## Pipeline overview

FACS is a pipeline to:

1. merge paired-end reads,
2. predict peptides,
3. cluster them at 100% of similarity and 100% coverage,
4. calculate their abundance in peptides per million (ppm), and
5. select those with antimicrobial potential discriminating their hemolytic pontential.

With FACS you can treat a metagenome file of 631.2Mbp as fast as 24 min, using 3 cpus and 100Mb sequence buckets in a Ubuntu v.18 64x bits.

![](https://github.com/celiosantosjr/FACS/blob/master/fig11.png)

**Figure 10.** FACS workflow.

![](https://github.com/celiosantosjr/FACS/blob/master/fig12.png)

**Figure 11.** FACS structure.

![](https://github.com/celiosantosjr/FACS/blob/master/fig17.png)

**Figure 12.** Decision tree to classification of peptides into different classes accordingly to their composition and capacity in forming dissulphide bonds.

## Descriptors system: Distribution

![](https://github.com/celiosantosjr/FACS/blob/master/fig13.png)

**Figure 13.** Method of sequence encoding using CTD (Composition, Distribution and Transition).

![](https://github.com/celiosantosjr/FACS/blob/master/fig14.png)

**Figure 14.** Example of CTD application to AMP discovery in AMPEP software. (Source: []())

**Table 1.** Classess adopted to the sequence encoding of the distribution at Residue0. The Solvent Accessibility was adopted as other studies previously Dubchak et al. [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/),[1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667), however the new feature "Free energy to transfer to lipophilic phase" was adopted from Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x).

| Properties 	| Class I 	| Class II 	| Class III 	|
|--------------------------------------------------------	|------------------------	|------------------	|------------------	|
| Solvent accessibility 	| A, L, F, C, G, I, V, W 	| R, K, Q, E, N, D 	| M, S, P, T, H, Y 	|
| Free energy to transfer from water to lipophilic phase 	| I,L,V,W,A,M,G,T 	| F,Y,S,Q,C,N 	| P,H,K,E,D,R 	|


## Datasets and training

![](https://github.com/celiosantosjr/FACS/blob/master/fig15.png)

**Figure 15.** Dataset for training and validation of antimicrobial peptides prediction model. (Source: []())

## Descriptors system: Distribution

![](https://github.com/celiosantosjr/FACS/blob/master/fig18.png)

**Figure 16.** Dataset for training and validation of hemolytic peptides prediction model. (Source: []())

![](https://github.com/celiosantosjr/FACS/blob/master/fig16.png)

**Figure 17.** Measures of model accuracy. (Source: []())

## Models and prediction accuracy

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

## Testing


**Table 4.** FACS assessment of runs performed with two different metagenomes.

| **Access** 	| SRR9016022 	| SRR9097106 	|
|:-----------:	|:----------:	|:-------------:	|
| **Size** 	| 631.2 Mbp 	| 2.9 Gbp 	|
| **User time** 	| 33m 49.703s 	| 28h 14m 37.909s 	|
| **System time** 	| 39m 25.515s 	| 30h 57m 35.064s 	|
| **Real time** 	| 1m 17.210s 	| 44m 29.606s 	|
| **AMP called** 	| 129,398 	| 6,376,290 	|

![](https://github.com/celiosantosjr/FACS/blob/master/fig19.png)

**Figure 18.** Results of test involving metagenome 631.2Mbp.

![](https://github.com/celiosantosjr/FACS/blob/master/fig20.png)

**Figure 19.** Results of test involving metagenome 2.9Gbp.

## Applications

FACS can be used in a wide ranging of scenarios, such as: screening for novel AMPs generating candidates to further testing and patenting, as well as, determination of microbiome quorum sensing mechanisms linking AMPs to health conditions or presence of diseases. 

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
	--reference [folder]	Folder where your reference files are located, if none current folder will be used   
	                        [Reference_seqs]

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
