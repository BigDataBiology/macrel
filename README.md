# MACREL: Meta(genomic) AMP Classification and Retrieval

Pipeline to mine antimicrobial peptides (AMPs) from (meta)genomes.

[![Build Status](https://travis-ci.com/BigDataBiology/macrel.svg?branch=master)](https://travis-ci.com/BigDataBiology/macrel)
[![GPL licensed](https://img.shields.io/badge/license-GPL-blue.svg)](https://raw.githubusercontent.com/hyperium/hyper/master/LICENSE)
[![Install with Bioconda](https://anaconda.org/bioconda/macrel/badges/installer/conda.svg)](https://anaconda.org/bioconda/macrel)
[![Install with Bioconda](https://anaconda.org/bioconda/macrel/badges/downloads.svg)](https://anaconda.org/bioconda/macrel)

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

## Applications

Macrel can be used in a wide-ranging of scenarios, such as screening for novel
AMPs, generating candidates to further testing and patenting, as well as,
determination of microbiome quorum sensing mechanisms linking AMPs to health
conditions or presence of diseases.

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

## Pipeline overview

The 22 descriptors adopted by Macrel are hybrid comprising local and global contexts to do the sequence encoding. Macrel performs firstly a distribution analysis of three classes of residues in two different features (Solvent accessibility and *Free energy to transfer from water to lipophilic phase*) as shown in Table 1. The novelty in this method is using the *Free energy to transfer from water to lipophilic phase* (FT) firstly described by Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x) to capture the spontaneity of the conformational change that AMPs suffer while their transference from water to the membrane. For more info about the other descriptors used in Macrel and the algorithms used to train the classifiers, please refer to the Macrel reference.

**Table 1.** Classes adopted to the sequence encoding of the distribution at the first residue of each class. The Solvent Accessibility was adopted as in previous studies (Dubchak et al. [1995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC41034/), [1999](https://www.ncbi.nlm.nih.gov/pubmed/10382667)), however, the new feature FT was adapted from Von Heijne and Blomberg, [1979](https://febs.onlinelibrary.wiley.com/doi/pdf/10.1111/j.1432-1033.1979.tb13100.x).

| Properties     | Class I     | Class II     | Class III     |
|:--------------------------------------------------------:    |:------------------------:    |:------------------    |------------------:    |
| Solvent accessibility     | A, L, F, C, G, I, V, W     | R, K, Q, E, N, D     | M, S, P, T, H, Y     |
| FT     | I,L,V,W,A,M,G,T     | F,Y,S,Q,C,N     | P,H,K,E,D,R     |

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

Benchmark procedures showed that R22 models are efficient in retrieving AMPs with statistics that overcome the state-of-art methods (Table 2).

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

Meanwhile, the hemolytic prediction model implemented in Macrel has a comparable performance of the state-of-art methods previously tested by Chaudhary et al., [2016](https://www.nature.com/articles/srep22843) as shown in Table 3.

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
