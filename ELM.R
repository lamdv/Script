setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(bigreadr)

# Set parameters
width = 1000

# Declare dataset here
data_train <- "data_train"
data_test <- "data_test"

# Read training example
snp_readBed(bedfile = paste(data_train,"bed",sep="."), backingfile = data_train)
obj.bigSNP <- snp_attach(paste(data_train,"rds",sep = "."))

obj.bigSNP$genotypes$show()
head(obj.bigSNP$fam)

# Generate input weight
# Gaussian input weight; sum = 0
W <- FBM(obj.bigSNP$genotypes$ncol, width, backingfile = paste("weight","rds",sep = "."))
W$show()
W [] <- rnorm(length(W))
head(W[])

apply_weight <- function(X, ind) {
  print(min(ind))
  r <- X[ind, ]
  product <- big_apply(W, a.FUN = function(W, ind, vect){
    return vect * W[ind, ]
  })
}

prods <- big_apply(obj.bigSNP$genotypes, a.FUN = apply_weight, 
                   a.combine = 'c', ind = rows_along(fbm),
                   block.size = 100e3, ncores = nb_cores())  