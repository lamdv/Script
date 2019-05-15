setwd("~/Documents/simus")

library(nnet)
library(bigstatsr)
library(bigsnpr)
library(bigreadr)

# Set parameters
width = 100

# Declare dataset here
data_train <- "data_train"
data_test <- "data_test"

# Read training example
snp_readBed(bedfile = paste(data_train,"bed",sep="."), backingfile = data_train)
obj.bigSNP <- snp_attach(paste(data_train,"rds",sep = "."))

obj.bigSNP$genotypes$show()
head(obj.bigSNP$fam)

# Declare NNet

#########
# Validation step
#########



pred <- big_prodMat(h2.out,beta.out)
length(pred)
length(test.bigSNP$fam$affection)
AUC(pred, test.bigSNP$fam$affection)
