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
W <- FBM(obj.bigSNP$genotypes$ncol, width, backingfile = ("weight"))
#W <- big_attach("weight.rds.bk")
W$show()
W [] <- rnorm(length(W),mean=0,sd = 1/sqrt(ncol(obj.bigSNP$genotypes)))
head(W[])

#apply_weight <- function(X, ind) {
#  print(min(ind))
#  r <- X[ind, ]
#  big_prodVec(W, r.T)
#}
#prods <- big_apply(obj.bigSNP$genotypes, a.FUN = apply_weight, 
#                   a.combine =h1 <- big_prodMat(W, 'c', ind = rows_along(obj.bigSNP$genotypes),
#                   block.size = 100e3, ncores = nb_cores())  

# ELM Steps
h1 <- as_FBM(big_prodMat(obj.bigSNP$genotypes,W))
hist(h1[,1])
h1
h1.out <- matrix(unlist(big_apply(h1, a.FUN = function(X, ind){
  1/(1+exp((-X[,ind])))
})), ncol = 1000, byrow = TRUE)

h1.out <- matrix(unlist(h1.out), ncol = 1000, byrow = TRUE)

h1.out <- as_FBM(h1.out)

# Logistic Regression on the output of hidden layer
test <- big_univLogReg(h1.out, y01.train = obj.bigSNP$fam$affection)

#########
# Validation step
#########

snp_readBed(bedfile = paste(data_test,"bed",sep="."), backingfile = data_test)
obj.bigSNP <- snp_attach(paste(data_test,"rds",sep = "."))

# ELM Steps
h1 <- as_FBM(big_prodMat(obj.bigSNP$genotypes,W))
hist(h1[,1])
h1
h1.out <- matrix(unlist(big_apply(h1, a.FUN = function(X, ind){
  1/(1+exp((-X[,ind])))
})), ncol = 1000, byrow = TRUE)

h1.out <- as_FBM(h1.out)
