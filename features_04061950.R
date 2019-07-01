#!/usr/bin env

##########################################################################
# Features extractor and table generator
# Analysis the peptide fasta to generate multiple features table
# Author: Célio D. Santos Jr.
##########################################################################

##########################################################################
# Required libraries
##########################################################################
if(!require(Peptides)){
  install.packages("Peptides")
  library(Peptides)
}

if(!require(data.table)){
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

if(!require(dplyr)){
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if(!require(parallel)){
  install.packages("parallel", dependencies = TRUE)
  library(parallel)
}

if(!require(doParallel)){
  install.packages("doParallel", dependencies = TRUE)
  library(doParallel)
}

##########################################################################
# CMDs
##########################################################################
print("[ M :::: Setting speeding up configs ]")
intervalStart <- Sys.time()
cluster <- makeCluster(detectCores() - 1, # number of cores to use
                         type = "FORK") # type of cluster
registerDoParallel(cluster)
set.seed(95014)
##########################################################################
print("[ M :::: Taking argument ]")
args <- commandArgs(TRUE)
print("[ M :::: Reading converted Tab/multifasta ]")
hs <- fread(file = args[1], sep="\t")
headers <- as.list(as.character(hs$header))
seqs <- as.list(as.character(hs$seq))
group <- as.list(as.character(hs$group))
rm(hs)
print("[ M :::: Generating lists ]")
aalist <- list()
charge <- list()
pI <- list()
aindex <- list()
instaindex <- list()
boman <- list()
hydrophobicity <- list()
hmoment <- list()
print("[ M :::: Making each amino acid types list ]")
print("[ M :::: aaCOMP ]")
aalist <- aaComp(seq = seqs)
print("[ M :::: Formatting ]")
aalist <- as.data.frame(aalist, header = TRUE)
d.edit <- aalist[ , grepl( "Mole." , names( aalist ) ) ]
rm(aalist)
colnames(d.edit) <- NULL
print("[ M :::: Making specific lists ]")
tinyAA <- as.list(d.edit[1,])
smallAA <- as.list(d.edit[2,])
aliphaticAA <- as.list(d.edit[3,])
aromaticAA <- as.list(d.edit[4,])
nonpolarAA <- as.list(d.edit[5,])
polarAA <- as.list(d.edit[6,])
chargedAA <- as.list(d.edit[7,])
basicAA <- as.list(d.edit[8,])
acidicAA <- as.list(d.edit[9,])
rm(d.edit)
print("[ M :::: Charge ]")
charge <- charge(seq = seqs, pH = 7, pKscale = "EMBOSS")
print("[ M :::: pI ]")
pI <- pI(seq = seqs, pKscale = "EMBOSS")
print("[ M :::: AIndex ]")
aindex <- aIndex(seq = seqs)
print("[ M :::: Instaindex ]")
instaindex <- instaIndex(seq = seqs)
print("[ M :::: Boman ]")
boman <- boman(seq = seqs)
print("[ M :::: Hydrophobicity ]")
hydrophobicity <- hydrophobicity(seq = seqs, scale = "Eisenberg")
print("[ M :::: Hmoment ]")
hmoment <- hmoment(seq = seqs, angle = 100, window = 11)
print("[ M :::: Formatting ]")
charge <- as.list(as.numeric(charge))
pI <- as.list(as.numeric(pI))
aindex <- as.list(as.numeric(aindex))
instaindex <- as.list(as.numeric(instaindex))
boman <- as.list(as.numeric(boman))
hydrophobicity <- as.list(as.numeric(hydrophobicity))
hmoment <- as.list(as.numeric(hmoment))
print("[ M :::: Preparing tables and binding lists ]")
t <- rbindlist(list(headers, seqs, group, tinyAA, smallAA, aliphaticAA, aromaticAA, nonpolarAA, polarAA, chargedAA, basicAA, acidicAA, charge, pI, aindex, instaindex, boman, hydrophobicity, hmoment), use.names=FALSE, fill=FALSE, idcol=NULL)
rm(headers, seqs, group, tinyAA, smallAA, aliphaticAA, aromaticAA, nonpolarAA, polarAA, chargedAA, basicAA, acidicAA, charge, pI, aindex, instaindex, boman, hydrophobicity, hmoment)
rownames(t) <- NULL
print("[ M :::: Formating rownames and colnames ]")
df <- c("access", "sequence", "group", "tinyAA", "smallAA", "aliphaticAA", "aromaticAA", "nonpolarAA", "polarAA", "chargedAA", "basicAA", "acidicAA", "charge", "pI", "aindex", "instaindex", "boman", "hydrophobicity", "hmoment") 
df <- as.data.frame(t(df))
t<-as.data.frame(t(t))
Final_table <- rbindlist(list(df,t))
rm(df, t)
print("[ M :::: Exporting ]")
write.table(Final_table, args[2], sep="\t", col.names=FALSE, row.names=FALSE)
print("[ M :::: Stop cluster ]")
stopCluster(cluster)
