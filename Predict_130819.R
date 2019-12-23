#!/usr/bin env

#####################################################
#####################################################

# PEP_gen: analysis the peptide fasta to generate
# multiple features table to AMP and NONAMP prediction

### Generation of profiles and traning - 2nd stage
### Prediction - 3rd stage
### Classification - 4th stage

######################### Author: Célio D. Santos Jr.
#####################################################

#####################################################
# Being feisty
.libPaths("PEPPERIDY")
#####################################################

#####################################################
# Required libraries
#####################################################
if(!require(randomForest)){
  install.packages("randomForest")
  library(randomForest)
}

if(!require(caret)){
  install.packages("caret")
  library(caret)
}

if(!require(dplyr)){
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if(!require(data.table)){
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
#####################################################

# CMDs
#####################################################
# Taking argument
args <- commandArgs(TRUE)

# Setting random seed
set.seed(95014)

# Reading and transforming datasets
hs <- fread(file=args[1], sep="\t")
Testingset <- as.data.frame(Testingset)
Testingset$group <- as.factor(Testingset$group)
str(Testingset)

######## Screening with model1

model <- readRDS(args[2])

predictionsFit <- predict(model, Testingset)
predictionsProb <- predict(model, Testingset, type = "prob")[,1]

# Filtering results

Testingset[["group"]] <- as.character(predictionsFit)
Testingset <- cbind(Testingset, predictionsProb)
Testingset <- Testingset %>% filter(group == "AMP")
Testingset$group <- ifelse (Testingset$acidicAA > Testingset$basicAA, ifelse(grepl("C", Testingset$sequence), "ADP", "ALP") , ifelse(grepl("C", Testingset$sequence), "CDP", "CLP"))

######## Test section

if(dim(Testingset)[1] == 0){
	print("No AMPs were found, cannot keep")
	quit("no")
}else{
	print("Keeping")
}

######## Screening with model2

model <- readRDS(args[3])

predictionsFit <- predict(model, Testingset)
predictionsProb2 <- predict(model, Testingset, type = "prob")[,1]

# Filtering results

predictionsFit <- as.character(predictionsFit)
prob_AMP <- as.numeric(Testingset$predictionsProb)
prob_HEMO <- as.numeric(predictionsProb2)
group <- as.character(Testingset$group)
access <- as.character(Testingset$access)
sequence <- as.character(Testingset$sequence)
final <- cbind(access, sequence, group, prob_AMP, predictionsFit, prob_HEMO)

## Exporting
write.table(final, args[4], sep="\t", col.names=FALSE, row.names=FALSE)
