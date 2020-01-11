##########################################################################
# CMDs
##########################################################################
if(!require(randomForest)){
  install.packages("randomForest")
  library(randomForest)
}

if(!require(caret)){
  install.packages("caret", dependencies = TRUE)
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


#  Setting random seed
set.seed(95014)

# Taking argument
args <- commandArgs(TRUE)

# reading descriptors file
training <- fread(file = args[1], sep="\t")
training <- training[,3:25]
training$group <- as.factor(training$group)

# building model
model <- randomForest(x = training[,-1], y = training$group, ntree = 100, cv.fold = 10) # Testing

# testing model
testing <- fread(file = args[2], sep="\t")
testing <- testing[,3:25]
testing$group <- as.factor(testing$group)

p <- predict(model, testing)
confusionMatrix(p, testing$group, positive = levels(training$group)[1])
varImpPlot(model)

# saving model
saveRDS(model, args[3])

