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


The second step is about calculation of descriptors. In FACS a two way system of descriptors was adopted using cheminformatics allied to sequence encoding, since the both methods were previously applied to AMP screening, but we believe there is a sinergy of those informations. We will further discuss about the descriptors used in FACS in detail. The descriptors calculation is made entirely in R using subscripts that run along the FACS pipeline using the sequence buckets and returning only AMP sequences that are then classified into hemolytic or non-hemolytic. The classifiers adopted in FACS are based in random forest algorithms (further discussed in details later) that proved to be more efficient than those previously reported (Gabere and Noble, [2017](https://www.ncbi.nlm.nih.gov/pubmed/28203715); Bhadra et al., [2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9), Meher et al., [2017](https://www.nature.com/articles/srep42362)). Finally, in the third step FACS performs a classification using a decisions tree (Figure 11) that classifies the detected AMPs into four different families accordingly their nature (Cationic or Anionic) and structure (linear or dissulphide bond forming). These classifications are then make available in a table where the sequence, random identifiers, abundance in ppm and hemolytic nature is also added. Interestingly, the FACS workflow depends on few third party softwares and some R libraries (Figure 12).


![](fig17.png)

**Figure 11.** Decision tree to classification of peptides into different classes accordingly to their composition and capacity in forming dissulphide bonds (Legend: AcidicAA - Acidic amino acids: B + D + E + Z; BasicAA - Alkaline amino acids: H + K + R).


6. Descriptors calculations is made by using R scripts that relies on R packages: Peptides, data.table, dplyr, parallel, doParallel;

7. AMPs prediction and Hemolytic activity classification as well as the families identification is basically implemented in R language and uses mostly the following R packages: randomForest, caret, data.table, dplyr;

8. The final formatting of files is performed by shell functions.

&nbsp;
&nbsp;

## Descriptors system: Distribution

&nbsp;
&nbsp;

Spänig and Heider ([2019](https://biodatamining.biomedcentral.com/track/pdf/10.1186/s13040-019-0196-x)) recently released a series of sequence encoding methods and a review of the main machine learning models using them. For a better comprehension of the following topics, we strongly recommend the its reading.


A recent prediction method released by Bhadra et al. ([2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9)), has shown that sequence encoding methods are enough to a good classification of AMPs in a large dataset of peptides. However, the proportion of AMPs to non-AMPs in the dataset influenced the results of the classifier. In this sense, they used a total of **23** different descriptors mostly based in the CTD method (Figure 13).


The CTD method was firstly described by Dubchak et al. ([1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/),[1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667)) and it is based in the classification of residues into three different classes accordingly to some specific features, such as hydrophobicity, solvent accessibility or secondary structure. The peptide sequence then is encoded into these three classes and the composition, distribution and transition of classes can be calculated. Mostly the composition and transition are important to other applications than AMP prediction, since they have shown small correlation to AMP peptides, besides not being explicative in some tested models performed by us previously.


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







**Table 1.** Classess adopted to the sequence encoding of the distribution at Residue0. The Solvent Accessibility was adopted as other studies previously Dubchak et al. [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/), [1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667), however the new feature "Free energy to transfer to lipophilic phase" was adopted from Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x).

| Properties 	| Class I 	| Class II 	| Class III 	|
|:----------:|:----------:|:----------:|:----------:|
| Solvent accessibility 	| A, L, F, C, G, I, V, W 	| R, K, Q, E, N, D 	| M, S, P, T, H, Y 	|
| FT 	| I,L,V,W,A,M,G,T 	| F,Y,S,Q,C,N 	| P,H,K,E,D,R 	|

&nbsp;
&nbsp;

The novelty in this method is use the *Free energy to transfer from water to lipophilic phase* (FT) firstly described by Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x). This measure is based in the estimation of free energy difference for the transfer of a residue from a random coil conformation in water to an alpha-helical conformation in a lipophilic phase (membrane). FT is calculate taking in account hydrophobicity, charge and polarity of the residues what reduces 3 features used to calculate CTD to one. Besides that, it also seems much more credible and important to predict AMPs since their functions are mostly based in the membranes interaction, and FT seems a measure of the potential of peptide insertion into them. To build the three classes, FT measures by residue were normalized as Zeta-Scores and  then sorted into three groups (Table 1). From the distribution of those three classes, FACS uses the first residue measure, getting the following descriptors:


 - SA.G1.residue0
 
 - SA.G2.residue0
 
 - SA.G3.residue0
 
 - hb.Group.1.residue0
 
 - hb.Group.2.residue0
 
 - hb.Group.3.residue0

&nbsp;
&nbsp;

The other cheminformatic descriptors are widely used in the AMPs description, and follows:


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

 - boman -> overall estimate of the potential of a peptide to bind tomembranes or other proteins as receptor

 - hydrophobicity (scale = "KyteDoolittle") -> GRAVY index

 - hmoment (angle = 100, window = 11) -> quantitative measure of the amphiphilicity perpendicular to the axis of any periodic peptide structure, such as the alpha-helix or beta-sheet

&nbsp;
&nbsp;

These descriptors are used to prediction and are calculated to each sequence that was identified as a potential peptide.

## Datasets and training

&nbsp;
&nbsp;

In a recent study, Gabere and Noble ([2017](https://www.ncbi.nlm.nih.gov/pubmed/28203715)) have shown that Random Forests models provide a statistically significant improvement in performance of AMPs detection, as measured by the area under the receiver operating characteristic (ROC) curve in comparison to other methods. Following this trend, many classifiers, such as AMPep and others, also make use of random forests models and presented highly accuracy in AMPs detection. Due to this, we opted to use random forests after some tests using alternative algorithms, such as: treebag, rpart, cpart, adaboost and others (results not shown). In order to be able to trace a comparison among the other classifiers and our model, we opted to use the same training and validation datasets used by Bhadra et al. ([2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9)), shown in Figure 15.


![](fig15.png)

**Figure 15.** Dataset for training and validation of antimicrobial peptides prediction model. (Source: Bhadra et al., [2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9))


The 22 descriptors were generated by using the customized script written for this purpose implemented in FACS. These descriptors were organized into tables containing 26 columns, where the first three standed for: peptide access code, sequence and abundance. Then, two models were generated by traning the AIs with the small and large datasets. Models were trained using random forest R package with a 10-cross fold validation and auto-mtry. The final models were tested by prediction internal system, being validated against the validation dataset containing AMP:non-AMP proportion of 1:1, and a total of 1840 peptides (Figure 15).


The hemolytic activity classifiers was obtained by using the datasets previously established by Chaudhary et al. ([2016](https://www.nature.com/articles/srep22843)). The HemoPI-1 datasets were used both to training and validation (Figure 16). Following the same pre-established logic, the descriptors were calculated as the standard method implemented in FACS by using Peptides R package and CTDDclass.py script from ILearn project, already implemented by the installation procedures bellow further discussed.


![](fig18.png)

**Figure 16.** Dataset for training and validation of hemolytic peptides prediction model. (Source: Chaudhary et al., [2016](https://www.nature.com/articles/srep22843))


Models for hemolytic activity were trained using caret R package, using the 22-descriptors dataset with the algorithm for Oblique Random Forests using Support Vector Machines (orfSVM). Training was performed by using a 5-cross fold validation repeated 3 times. This algorithm is a breakthrough in the predictions, since Chaudhary et al. ([2016](https://www.nature.com/articles/srep22843)) used SVM and got similar results. However, differently from their work, this model works faster using much less descriptors and is implemented in R language, not Java, which was highly dependant of the users settings and not portable after all.


The models here mentioned were implemented in a R script to filter off the non-AMP peptides and classify AMPs into hemolytic or not. After that, this very script also makes the decisions tree sorting presented in Figure 11, classifying AMPs also into the 4 families:


 - Anionic linear peptides (ALP)

 - Anionic dissulphide-bond forming peptides (ADP)

 - Cationic linear peptides (CLP)

 - Cationic dissulphide-bond forming peptides (CDP)
 
&nbsp;
&nbsp;

## Models and prediction accuracy

&nbsp;
&nbsp;

The classifiers of AMP and hemolytic peptides were then assessed and compared to the state of art methods of AMP prediction and hemolytic classification. The results were taken using the same databases of the works cited as reference. Calculations of the statistics of accuracy, sensitivity, McNeymar's Correlation Coefficient (MCC), and F-score were shown in Figure 17.


![](fig16.png)

**Figure 17.** Measures of model accuracy. Legend: FN - False Negatives, TP - True Positives, TN - True Negatives, FP - False Positives. (Source: Bhadra et al., [2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9))


The AMP prediction model was compared at two levels the first level with it trained with the small training dataset and when trained with the AMPEP complete dataset, presenting a ratio of 1:3 (positives:negatives). The final results (Table 2) shows clearly that R22 models are comparable efficient in retrieving AMPs from the validation dataset reaching accuracies very close to the best system so far (AMPEP). However, R22 trained with the complete dataset (R22_LargeTrainingset) outperforms AMPEP with a better accuracy, specificity and precision. These features also reflects a better global adjustment of this model, since its MCC and F-score were higher than those from AMPEP. Thus, R22_LargeTrainingset seems the best model to be used in the AMPs prediction, besides use less descriptors and memory. Besides, It also is all implemented in R, what gives portability to the process that could be entirely implemented in FACS just by external scripts, composing a bigger workflow.


**Table 2.** Models to predict antimicrobial peptides were tested to benchmark the results obtained with the new set of descriptors adopted in this classifier. All models and systems were tested with the benchmark validation dataset from the AMPEP study available in Bhadra et al. ([2018](doi:10.1038/s41598-018-19752-w)).

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
| R22_Large 	| 0.967 	| 1 	| 0.934 	| 1 	| 0.966 	| 0.936 	| This study 	| 22 	|

&nbsp;
&nbsp;

The specific results of the confusion matrix are presented now (Table 3) to the R22_LargeTrainingset. The 100% of specificity does not mean an overfitting since there is a clear misclassification of positive peptides, which also ensures that the model is still reliable to be used in other datasets.


**Table 3.** Confusion Matrix of R22_LargeTrainingset model.

|Prediction/Reference|AMP|Non-AMP|
|:--------:|:-----:|:-----:|
|AMP|859|0|
|Non-AMP|61|920|

&nbsp;
&nbsp;

Meanwhile, the hemolytic prediction model was evaluated using its own datasets and trained as previously informed. The results of this model and the comparisons to the standard system currently available are shown in the Table 4. The model obtained with oblique random forests supported by vector machines was a bit less accurate, however the sensitivity and specificity were higher than that obtained previously. Moreover, the MCC measure also shows our model performance similar to the best model currently implemented in the server, based in support vector machines (SVM). The similarities among their performances are important to make sure our model is reliable, regarding the convenience of being implemented in R with a set of descriptors previously calculated to the AMP prediction model. In this way, the previous tables can be reused in this case, saving time and memory.


**Table 4.** Models used to predict hemolytic activity were trained with the same dataset used by Chaudhary et al. [2016](https://www.nature.com/articles/srep22843) and were tested with the same test dataset used by them to benchmark the results obtained with another model used to generate this classifier.

| Methods 	| Study 	| Sn (%) 	| Sp (%) 	| Acc (%) 	| MCC 	|
|:---------------------:	|:----------------------:	|:------:	|:------:	|:-------:	|:----:	|
| SVM 	| Chaudhary et al., 2016 	| 95.7 	| 94.8 	| 95.3 	| 0.91 	|
| IBK 	| Chaudhary et al., 2016 	| 95.5 	| 93.7 	| 94.6 	| 0.89 	|
| Multilayer Perceptron 	| Chaudhary et al., 2016 	| 93.9 	| 92.8 	| 93.3 	| 0.87 	|
| Logistic 	| Chaudhary et al., 2016 	| 93.4 	| 93.7 	| 93.6 	| 0.87 	|
| J48 	| Chaudhary et al., 2016 	| 89.6 	| 88.5 	| 89.0 	| 0.78 	|
| Random Forest 	| Chaudhary et al., 2016 	| 94.1 	| 94.6 	| 94.3 	| 0.89 	|
| ORFsvm 	| This study 	| 95.5 	| 95.5 	| 95.5 	| 0.91 	|

&nbsp;
&nbsp;

Our classifiers seems to be extremely interesting in the execution of the filtering off non-AMP peptides and classifying them into hemolytic or non-hemolytic peptides. Also, the models were implemented in the same programming language (R) and used the same set of descriptors, what saved time and memory in the process implemented in FACS.


As a future update to FACS classifier systems, some efforts recently have being done, in order to achieve a model to classify the biological activity presented by AMPs: anti-bacterial, anti-fungal, anti-viral, anti-HIV and anti-tumor. A test was carried out using the same set of descriptors and training with the AMPEP training dataset available [here](doi:10.1038/s41598-018-19752-w). The training procedures were same adopted before to both models, testing the random forest (rf) and the oblique random forests with support vector machines (orfsvm) algorithms. The model here presented is a result of the 10-cross fold validated random forest training and showed a limited capacity of classification (Figure 18), with a specificity and sensitivity ranging to very low values in some classes.


![](fig23.png)

**Figure 18.** Confusion matrix of model "Mixedmodel_multiclassifier". This model tried to classify the 8 different AMP biological activities. To train it we have used the large training dataset, previously provided by Bhadra et al., [2018](https://www.nature.com/articles/s41598-018-19752-w#Sec9).


While this classifier works well for prediction of anti-bacterial activity with (sensitivity and specificity of ~97%), other activities, such as anti-viral activity is badly classified with a sensitivity of only 23%. This reflects, in a general way, a non-homogeneous performance, also regarding the biochemical diversity of the peptides present in each one of those classes. These results suggest that a bigger training set, well annotated and maybe another conditions of training could improve the prediction performance.

&nbsp;
&nbsp;

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
