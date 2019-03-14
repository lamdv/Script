setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(bigreadr)

# Set parameters
width = 1000

# Declare dataset here
data_train <- "data_train.bed"
data_test <- "data_test.bed"

# Read training example
tmpfile <- tempfile()
snp_readBed(bedfile = data_train, backingfile = tmpfile)
obj.bigSNP <- snp_attach(paste(tmpfile,"rds",sep = "."))

obj.bigSNP$genotypes$show()
head(obj.bigSNP$fam)

# Generate input weight
# Gaussian input weight; sum = 0
tmpfile_weight <- tempfile()
W <- FBM(obj.bigSNP$genotypes$ncol, width, backingfile = paste(tmpfile,"rds",sep = "."))
W$show()
W [] <- rnorm(length(W))

apply_weight <- function(X, row) {
  X[]
}