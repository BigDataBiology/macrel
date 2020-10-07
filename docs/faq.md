# Frequent Asked Questions - FAQ

*[1]. I tested Macrel with some propeptide and it returned a false negative. Do Macrel detect RiPPs?*

    Macrel v.0.5.0 just detects mature peptides. We expect to develop and test approaches to mask the

    leader or tail peptides and get predictions from RiPPs too in the next versions of the program.

*[2]. Do Macrel detect peptides of any length?*

    Macrel v.0.5.0 just detects peptides ranging from 10 to 100 residues.
  
*[3]. Can Macrel work well on isolate genomes or isolate genome reads?* 
  
     Macrel is able to handle isolated genomes, as contigs, short reads (in this case it will trim and assemble

     before starting the gene prediction) or directly on predicted peptide sequences (obtained through some

     alternative pipeline).
  
*[4]. How is Macrel licensed?*

    Macrel is under a hybrid license. Most of its system is under MIT license, however, parts of its system were

    developed from third parties and are licensed under GNU license.
