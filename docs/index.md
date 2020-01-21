---
title: "MACREL - Meta(genomic) AMPs Classification and REtrievaL"
authors: "Celio Dias Santos Jr.; Shaojun Pan, Xing-Ming Zhao, Luis Pedro Coelho"
date: "2020-01-21"
output: pdf_document
---

## Brief

Macrel stands for Meta(genomic) AMPs Classification and REtrievaL. It can be used in a wide-ranging of scenarios, such as screening for novel AMPs, generating candidates to further testing and patenting, as well as, determination of microbiome quorum sensing mechanisms linking AMPs to health conditions or presence of diseases. The application is distributed under GPL v3 license and represents a joint effort of Celio Dias Santos Jr.; Shaojun Pan; Xing-Ming Zhao and Luis Pedro Coelho from the Institute of Science and Technology for Brain-Inspired Intelligence (ISTBI) at Fudan University (Shanghai, China).

If you use this software in a publication please cite

>   FACS: antimicrobial peptide screening in genomes and metagenomes
>   Celio Dias Santos-Junior, Shaojun Pan, Xing-Ming Zhao, Luis Pedro Coelho
>   bioRxiv 2019.12.17.880385; doi:
>   [https://doi.org/10.1101/2019.12.17.880385](https://doi.org/10.1101/2019.12.17.880385)

(The preprint still uses the old name of the tool, _FACS_ and will be updated
soon).

**NOTE**: This is still a _work in progress_ and, while the results of the tool
should be correct, we are still working on making Macrel easier to install and
use.

## Background

AMPs definition comprises peptides with a huge variety of biological activities (such as anticancer, antibacterial, antifungal and insecticidal), and their sequences are key to that activity. AMP producing microbes can limit the growth of other microorganisms and should be considered another normal source of them. Microbial sourced AMPs are quite distinct from those of vertebrates, since they can be obtained from a nonribosomal peptide synthase. Thus, nonribosomal peptides can adopt different structures, such as cyclic or branched structures, and carry modifications like N-methyl and N-formyl groups, glycosylations, acylations, halogenation, or hydroxylation. Some examples of commercial microbial AMPs include polymyxin B and vancomycin, both FDA-approved antibiotics (Zhang and Gallo, [2016](https://www.sciencedirect.com/science/article/pii/S0960982215014098)).

Most AMPs are peptides 10-50 residues long (some reaching 100 amino acids), with charge ranging between 2 and 11 and consituted of approximately 50% of hydrophobic residues (Zhang and Gallo, [2016](https://www.sciencedirect.com/science/article/pii/S0960982215014098)). The formation of amphiphilic ordered structures works as a driving force for membrane binding and disruption, a key AMP feature. The helix destabilization oftenly can reduce the cytotoxicity of AMPs, although this can result in reduction of antimicrobial effects (Malmsten, [2014](https://www.tandfonline.com/doi/full/10.3109/03009734.2014.899278); Borgden, [2005](https://www.ncbi.nlm.nih.gov/pubmed/15703760?dopt=Abstract); Pasupuleti et al., [2012](www.ncbi.nlm.nih.gov/pubmed/22074402?dopt=Abstract); Hancok and Sahl, [2006](https://www.ncbi.nlm.nih.gov/pubmed/17160061?dopt=Abstract); Shai, [2002](https://www.ncbi.nlm.nih.gov/pubmed/12491537?dopt=Abstract); Stromstedt et al., [2006](https://www.ncbi.nlm.nih.gov/pubmed/19029324?dopt=Abstract)). There is a dynamic interchange in AMPs structure and topologies along the interaction with the microbial cell membranes (Samson, [1998](https://www.sciencedirect.com/science/article/pii/S1359029498800277)) and electrostatic interactions between the outer membrane surface of prokaryotic cells (negatively charged) amd AMPs is the primary mechanism for antimicrobial activity. Mostly the activities of AMPs are associated to the rupture of cell membrane, promoting the leakage of cell contents. Other cases are based in the AMP translocation across the cell membrane and the inhibition of essential cellular processes (e.g. protein synthesis, nucleic acid synthesis, enzymatic activities) (Brogden, [2005](https://www.ncbi.nlm.nih.gov/pubmed/15703760)). Based on the mechanisms of action, AMPs are categorized into membrane acting and nonmembrane acting peptides. The first ones are capable of forming transient pores on the membrane, whereas the second ones have the ability to translocate themselves across the cell membrane without permeabilizing it (Pushpanathan et al., [2013](http://dx.doi.org/10.1155/2013/675391)).

The genomic era was constrasted by the reality of hundreds of available bacterial genomes that have so far failed to deliver the hoped-for new molecular targets for antibiotics. However, so far it always have focused in the active molecules produced by the metabolism, instead searching for active peptides or proteins. The best reason to bet in host defense antimicrobial peptides or AMPs is that they remained potent for millions of years, constituting a useful strategy to develop a new generation of antimicrobials to meet the growing antibiotic resistance problem worldwide. However, the prediction of small genes from meta-genomic/transcriptomic sequences and the prediction of active AMP are the main problems with AMP mining from meta- and genomic data sets.

Current methods to small genes prediction tipically lead to unacceptably high rates of false positives (Hyatt et al., [2010](https://www.ncbi.nlm.nih.gov/pubmed/20211023)). Recent smORFs surveys demonstrated that these methods followed by a filtering of false-positives can lead to biologically active smORFs (Miravet-Verde et al., [2019](https://www.ncbi.nlm.nih.gov/pubmed/30796087); Sberro et al., [2019](https://www.ncbi.nlm.nih.gov/pubmed/31402174)). Furthermore, the prediction of AMP activity demands techniques other than homology-based methods, due to the degeneration of searches at smaller sequences. Several machine learning-based methods demonstrated high accuracy in predicting antimicrobial activity in peptides (Xiao et al., [2013](https://www.ncbi.nlm.nih.gov/pubmed/23395824); Meher et al., [2017](https://www.ncbi.nlm.nih.gov/pubmed/28205576); Bhadra et al., [2018](https://www.ncbi.nlm.nih.gov/pubmed/29374199)), although, none of them represented a full pipeline to extract AMPs from genomic data and filter off mispredictions. Our main goal with Macrel is a highthroughput screening system of AMPs, through machine learning, able to retrieve AMP sequences with high confidence from meta(genomic) reads.


## Pipeline overview

Macrel pipeline does:

1. quality trimming of single- and paired-end reads,

2. assembly of reads into contigs

3. small genes prediction,

4. clustering of peptides at 100% of similarity and 100% coverage,

5. calculation of peptides features,

6. classification of peptides into AMPs,

7. classification of AMPs accordingly to their hemolytic activity,

8. calculate AMPs abundance in meta(genomic) samples by reads mapping.

Macrel is fast and works by coordinating [NGLess](https://github.com/ngless-toolkit/ngless), [megahit](https://github.com/voutcn/megahit), [prodigal](https://github.com/hyattpd/Prodigal) and [PALADIN](https://github.com/ToniWestbrook/paladin). It is implemented in Python and R, and its models were trained with [Scikit-Learn](https://github.com/scikit-learn/scikit-learn) python module.

## AMPs classification

The 22 descriptors adopted by Macrel are hybrid comprising local and global contexts to do the sequence encoding. Macrel performs firstly a distribution analysis of three classes of residues in two different features (Solvent accessibility and *Free energy to transfer from water to lipophilic phase*) as shown in Table 1. The novelty in this method is using the *Free energy to transfer from water to lipophilic phase* (FT) firstly described by Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x) to capture the spontaneity of the conformational change that AMPs suffer while their transference from water to the membrane. For more info about the other descriptors used in Macrel and the algorithms used to train the classifiers, please refer to the Macrel reference.

**Table 1.** Classes adopted to the sequence encoding of the distribution at the first residue of each class. The Solvent Accessibility was adopted as in previous studies (Dubchak et al. [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/), [1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667)), however, the new feature FT was adapted from Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x).

| Properties     | Class I     | Class II     | Class III     |
|:--------------------------------------------------------:    |:------------------------:    |:------------------    |------------------:    |
| Solvent accessibility     | A, L, F, C, G, I, V, W     | R, K, Q, E, N, D     | M, S, P, T, H, Y     |
| FT     | I,L,V,W,A,M,G,T     | F,Y,S,Q,C,N     | P,H,K,E,D,R     |

The distribution parameter derivates from The CTD method was firstly described by Dubchak et al. ([1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/),[1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667)) - Figure 13, and it is based in the classification of residues into three different classes accordingly to some specific features, such as hydrophobicity, solvent accessibility or secondary structure. The peptide sequence then is encoded into these three classes and the composition, distribution and transition of classes can be calculated. Mostly the composition and transition are important to other applications than AMP prediction, since they have shown small correlation to AMP peptides, besides not being explicative in some tested models performed by us previously.


![](fig13.png)

**Figure 13.** Method of sequence encoding using CTD (Composition, Distribution and Transition). (Source: Dubchak et al., [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/))


As above mentioned, since distribution of the residues classes seemed to be more effective to explain and classify AMPs, AMPep software (Bhadra et al., [2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9)) was mostly based in the distribution of the classes of the five canonical features (hydrophobicity, normalized van der Waals volume, polarity, polarizability, charge,  secondary structure, and solvent accessibility), as shown in Figure 14.


![](fig14.png)

**Figure 14.** Example of CTD application to AMP discovery in AMPEP software. (Source: Bhadra et al., [2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9))


However, Bhadra et al. ([2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9)) observed that the distribution was more important when taken from the first residue, calculated as:


                             Z = [ R x Y / 100 ], where:
        
            - R is the total number of class residues in the sequence,

            - Y denotes the desired percentage.

&nbsp;
&nbsp;

In this sense, despite the high Accuracy and sensitivity, other works still suggest that methods independent of sequence order and mostly based in cheminformatics have with comparable statistics (Boone et al., [2018](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-018-2514-6)). Thus, Fjell et al. ([2009](https://pubs.acs.org/doi/10.1021/jm8015365)) has shown using a combination of 77 QSAR (quantitative structure-activity relationships) descriptors that artificial neural network models could predict the extension of peptides activity, not only classify them.


Thus, these methods could be joined to achieve a better performance and also economy of computational resources, since sequence encoding represents a considerable cost of processing. The alliance between those two methods can fix their pitfalls, since the sequence order independent methods fail in classify, however are good to describe activity; and sequence encoding is essential to a good classification, but fails when predict activity extent.


The descriptors adopted by FACS are hybrid being partially cheminformatics and sequence encodings, what by itself represents a breakthrough. FACS performs firstly a distribution analysis of three classes of residues in two different features (Solvent accessibility and *Free energy to transfer from water to lipophilic phase*) as shown in Table 1.


The other descriptors used in Macrel classifiers are widely used in the AMPs description, as follows:

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

Finally, in the third step FACS performs a classification using a decisions tree (Figure 11) that classifies the detected AMPs into four different families accordingly their nature (Cationic or Anionic) and structure (linear or dissulphide bond forming). These classifications are then make available in a table where the sequence, random identifiers, abundance in ppm and hemolytic nature is also added. Interestingly, the FACS workflow depends on few third party softwares and some R libraries (Figure 12).

![](fig17.png)

**Figure 11.** Decision tree to classification of peptides into different classes accordingly to their composition and capacity in forming dissulphide bonds (Legend: AcidicAA - Acidic amino acids: B + D + E + Z; BasicAA - Alkaline amino acids: H + K + R).


## Benchmark

Benchmark procedures showed that R22 models are efficient in retrieving AMPs with statistics that overcome the state-of-art methods (Table 2). The AMP prediction model was compared at two levels the first level with it trained with the small training dataset and when trained with the AMPEP complete dataset, presenting a ratio of 1:3 (positives:negatives). The final results (Table 2) shows clearly that R22 models are comparable efficient in retrieving AMPs from the validation dataset reaching accuracies very close to the best system so far (AMPEP). However, R22 trained with the complete dataset (R22_LargeTrainingset) outperforms AMPEP with a better accuracy, specificity and precision. These features also reflects a better global adjustment of this model, since its MCC and F-score were higher than those from AMPEP. Thus, R22_LargeTrainingset seems the best model to be used in the AMPs prediction, besides use less descriptors and memory. Besides, It also is all implemented in R, what gives portability to the process that could be entirely implemented in FACS just by external scripts, composing a bigger workflow.

**Table 2.** Comparison of Macrel and other state-of-art AMP prediction systems. All systems were tested with the benchmark data set from the AMPEP study available in Bhadra et al. ([2018](doi:10.1038/s41598-018-19752-w)).

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

Meanwhile, the hemolytic prediction model implemented in Macrel has a comparable performance of the state-of-art methods previously tested by Chaudhary et al., [2016](https://www.nature.com/articles/srep22843) as shown in Table 3. Meanwhile, the hemolytic prediction model was evaluated using its own datasets and trained as previously informed. The results of this model and the comparisons to the standard system currently available are shown in the Table 4. The model obtained with oblique random forests supported by vector machines was a bit less accurate, however the sensitivity and specificity were higher than that obtained previously. Moreover, the MCC measure also shows our model performance similar to the best model currently implemented in the server, based in support vector machines (SVM). The similarities among their performances are important to make sure our model is reliable, regarding the convenience of being implemented in R with a set of descriptors previously calculated to the AMP prediction model. In this way, the previous tables can be reused in this case, saving time and memory.


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


Our classifiers seems to be extremely interesting in the execution of the filtering off non-AMP peptides and classifying them into hemolytic or non-hemolytic peptides. Also, the models were implemented in the same programming language (R) and used the same set of descriptors, what saved time and memory in the process implemented in FACS.

As a future update to FACS classifier systems, some efforts recently have being done, in order to achieve a model to classify the biological activity presented by AMPs: anti-bacterial, anti-fungal, anti-viral, anti-HIV and anti-tumor. A test was carried out using the same set of descriptors and training with the AMPEP training dataset available [here](doi:10.1038/s41598-018-19752-w). The training procedures were same adopted before to both models, testing the random forest (rf) and the oblique random forests with support vector machines (orfsvm) algorithms. The model here presented is a result of the 10-cross fold validated random forest training and showed a limited capacity of classification (Figure 18), with a specificity and sensitivity ranging to very low values in some classes.

## Install

The recommended method of installation is through
[bioconda](https://anaconda.org/bioconda/macrel):

```bash
conda install -c bioconda macrel
```

## Install from source

If you want to use an unreleased version from Github, for example, we provide a
script which _conda_ (in particular, using
[bioconda](https://bioconda.github.io/) and
[conda-forge](https://conda-forge.org/)) to install all dependencies in a
Macrel-private environment:

```bash
git clone https://github.com/BigDataBiology/Macrel
cd macrel
./install.sh
conda activate envs/Macrel_env
```

Henceforth, to use macrel, activate this environment.

### Examples

> Macrel uses a _subcommand interface_. You run `macrel COMMAND ...` with the
> COMMAND specifying which components of the pipeline you want to use.

To run Macrel on peptides, use the `peptides` subcommand:

```bash
macrel peptides \
    --fasta example_seqs/expep.faa.gz \
    --output out_peptides \
    -t 4
```

In this case, we use `example_seqs/expep.faa.gz` as input sequence. This should
be an amino-acid FASTA file. The outputs will be written into a folder called
`out_peptides`, and Macrel will 4 threads.

To run Macrel on contigs, use the `contigs` subcommand:

```bash
macrel contigs \
    --fasta example_seqs/excontigs.fna.gz \
    --output out_contigs
```

In this example, we use the example file `excontigs.fna.gz` which is a FASTA
file with nucleotide sequences, writing the output to `out_contigs`.

To run Macrel on paired-end reads, use the `reads` subcommand:

```bash
macrel reads \
    -1 example_seqs/R1.fq.gz \
    -2 example_seqs/R2.fq.gz \
    --output out_metag \
    --outtag example_metag
```

The paired-end reads are given as paired files (here, `example_seqs/R1.fq.gz`
and `example_seqs/R2.fq.gz`). If you only have single-end reads, you can omit
the `-2` argument.

To run Macrel to get abundance profiles, you only need the short reads file
and a reference with peptide sequences. Use the `abundance` subcommand:


```bash
macrel abundance \
    -1 example_seqs/R1.fq.gz \
    --fasta example_seqs/ref.faa.gz \
    --output out_abundance \
    --outtag example_abundance
```
