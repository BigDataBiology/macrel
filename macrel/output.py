readme_output_abundance_mode = """
--------------------------------------------------------------------------------------------------------------------------
------------------------------------------          READ  ME        ------------------------------------------------------
                  Thank you for using Macrel - the (meta)genome AMPs classification and retrieval. 
                  If you find Macrel useful, please cite:
                  Santos-Junior, C.D. et al. MACREL: antimicrobial peptide screening in genomes and
                  metagenomes. The PeerJ 8:e10555 DOI 10.7717/peerj.10555
--------------------------------------------------------------------------------------------------------------------------
== This is an explanatory file to the outputs given by Macrel.
--------------------------------------------------------------------------------------------------------------------------
Usage example in version 1.0 on the abundance subcommand:
>>> macrel abundance \
    -1 example_seqs/R1.fq.gz \
    --fasta example_seqs/ref.faa.gz \
    --output out_abundance \
    --tag example_abundance
**out_abundance/example_abundance.out.abundance.txt
    > This file contains a table in the following format:
        MACREL
    [1]       [2] 
    = Field 1. Access of smORFs arbitrarily given by Macrel or the name of the proteins in the reference file
    = Field 2. Name of the tag given to Macrel or generically "MACREL" and represents the abundance given in reads mapped
    to the peptide, we recommend to a normalization the user make a pseudo-RPKM calculus as:
    
        pseudo-RPKM =  mapped reads * 1000 / total reads * (length in bp)
"""

readme_output_contigs_mode = """
--------------------------------------------------------------------------------------------------------------------------
------------------------------------------          READ  ME        ------------------------------------------------------
                  Thank you for using Macrel - the (meta)genome AMPs classification and retrieval. 
                  If you find Macrel useful, please cite:
                  Santos-Junior, C.D. et al. MACREL: antimicrobial peptide screening in genomes and
                  metagenomes. The PeerJ 8:e10555 DOI 10.7717/peerj.10555
--------------------------------------------------------------------------------------------------------------------------
== This is an explanatory file to the outputs given by Macrel.
--------------------------------------------------------------------------------------------------------------------------
Usage example in version 1.0 on the contigs subcommand:
>>> macrel contigs \
    --fasta example_seqs/excontigs.fna.gz \
    --output out_contigs
**out_contigs/macrel.out.all_orfs.faa
    > This file contains a prediction for all genes from the contigs given to Macrel independently the size outputted
**out_contigs/macrel.out.smorfs.faa
    > This file contains only those genes obtained after a length filtering. Here, the genes with a length ranging from
    10 to 100 amino acids were kept. This file can be used, for example for independent runs in Macrel under the "peptides"
    mode.
**out_contigs/macrel.out.prediction.gz
    > This file contains a table in the following format:
    Access  Sequence    AMP_family  AMP_probability Hemolytic   Hemolytic_probability
    [1]       [2]             [3]                   [4]                          [5]
    = Field 1. Access of smORFs arbitrarily given by Macrel
    = Field 2. Predicted amino acid sequence 
    = Field 3. AMP family classified into anionic/cationic and cysteine-containing/linear coded as
        A/C and D/L, respectively. (e.g. a cationic cysteine-containing peptide would be CDP)
    = Field 4. Probability of a given peptide be antimicrobial in Macrel's model
    = Field 5. Classification of hemolytic activity of a given peptide into hemolytic (Hemo) or
        non-hemolytic (NonHemo)
    = Field 6. Probability of a given peptide be hemolytic in Macrel's model
"""

readme_output_peptides_mode = """
--------------------------------------------------------------------------------------------------------------------------
------------------------------------------          READ  ME        ------------------------------------------------------
                  Thank you for using Macrel - the (meta)genome AMPs classification and retrieval. 
                  If you find Macrel useful, please cite:
                  Santos-Junior, C.D. et al. MACREL: antimicrobial peptide screening in genomes and
                  metagenomes. The PeerJ 8:e10555 DOI 10.7717/peerj.10555
--------------------------------------------------------------------------------------------------------------------------
== This is an explanatory file to the outputs given by Macrel.
--------------------------------------------------------------------------------------------------------------------------
Usage example in version 1.0 on the peptides subcommand:
>>> macrel peptides \
    --fasta example_seqs/expep.faa.gz \
    --output out_peptides \
    -t 4
**out_peptides/macrel.out.prediction.gz
    > This file contains a table in the following format:
    Access  Sequence    AMP_family  AMP_probability Hemolytic   Hemolytic_probability
    [1]       [2]             [3]                   [4]                          [5]
    = Field 1. Access of smORFs arbitrarily given by Macrel
    = Field 2. Predicted amino acid sequence 
    = Field 3. AMP family classified into anionic/cationic and cysteine-containing/linear coded as
        A/C and D/L, respectively. (e.g. a cationic cysteine-containing peptide would be CDP)
    = Field 4. Probability of a given peptide be antimicrobial in Macrel's model
    = Field 5. Classification of hemolytic activity of a given peptide into hemolytic (Hemo) or
        non-hemolytic (NonHemo)
    = Field 6. Probability of a given peptide be hemolytic in Macrel's model
"""

readme_output_reads_mode = """
--------------------------------------------------------------------------------------------------------------------------
------------------------------------------          READ  ME        ------------------------------------------------------
                  Thank you for using Macrel - the (meta)genome AMPs classification and retrieval. 
                  If you find Macrel useful, please cite:
                  Santos-Junior, C.D. et al. MACREL: antimicrobial peptide screening in genomes and
                  metagenomes. The PeerJ 8:e10555 DOI 10.7717/peerj.10555
--------------------------------------------------------------------------------------------------------------------------
== This is an explanatory file to the outputs given by Macrel.
--------------------------------------------------------------------------------------------------------------------------
Usage example in version 1.0 on the reads subcommand:
>>> macrel reads \
    -1 example_seqs/R1.fq.gz \
    -2 example_seqs/R2.fq.gz \
    --output out_metag \
    --tag example_metag
**out_metag/example_metag.megahit_output/
    > This folder contains the assembly results obtained with megahit. It contains several of subfolders and intermediate
    steps of assembly with important information that allow the user trace back crucial step to infer missassemblies or
    even detect events/binning elements of interest.
**out_metag/example_metag.all_orfs.faa
    > This file contains a prediction for all genes from the contigs given to Macrel independently the size outputted
**out_metag/example_metag.smorfs.faa
    > This file contains only those genes obtained after a length filtering. Here, the genes with a length ranging from
    10 to 100 amino acids were kept. This file can be used, for example for independent runs in Macrel under the "peptides"
    mode.
**out_metag/example_metag.prediction.gz
    > This file contains a table in the following format:
    Access  Sequence    AMP_family  AMP_probability Hemolytic   Hemolytic_probability
    [1]       [2]             [3]                   [4]                          [5]
    = Field 1. Access of smORFs arbitrarily given by Macrel
    = Field 2. Predicted amino acid sequence 
    = Field 3. AMP family classified into anionic/cationic and cysteine-containing/linear coded as
        A/C and D/L, respectively. (e.g. a cationic cysteine-containing peptide would be CDP)
    = Field 4. Probability of a given peptide be antimicrobial in Macrel's model
    = Field 5. Classification of hemolytic activity of a given peptide into hemolytic (Hemo) or
        non-hemolytic (NonHemo)
    = Field 6. Probability of a given peptide be hemolytic in Macrel's model
"""