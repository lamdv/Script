setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(bigreadr)
library(sigmoid)
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
W <- FBM(obj.bigSNP$genotypes$ncol, setwd("~/Documents/simus")width, backingfile = ("weight"))
W$save()
W <- big_attach("weight.rds")
W$show()
W [] <- rnorm(length(W),mean=0,sd = 1/nrow(obj.bigSNP$genotypes))
hist(W[])
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
big_apply(obj.bigSNP$genotypes, a.FUN = function(X, ind){
  X[, ind] <- X[, ind] - 0.5
})
h1 <- as_FBM(big_prodMat(h1,W))
hist(h1[])
h1
h1.out <- matrix(unlist(big_apply(h1, a.FUN = function(X, ind){
  relu(X[,ind])
})), ncol = width, byrow = TRUE)

h1.out <- matrix(unlist(h1.out), ncol = width, byrow = TRUE)

h1.out <- as_FBM(h1.out)
hist(h1.out[])
h1.out$show()
# Logistic Regression on the output of hidden layer
test <- big_univLogReg(h1.out, y01.train = obj.bigSNP$fam$affection)

#########
# Validation step
#########

snp_readBed(bedfile = paste(data_test,"bed",sep="."), backingfile = data_test)
test.bigSNP <- snp_attach(paste(data_test,"rds",sep = "."))

# ELM Steps
h2 <- as_FBM(big_prodMat(test.bigSNP$genotypes,W))

h2.out <- matrix(unlist(big_apply(h2, a.FUN = function(X, ind){
  1/(1+exp((-X[,ind])))
})), ncol = width, byrow = TRUE)

h2.out <- as_FBM(h2.out)
beta.out <- as_FBM(matrix(unlist(test$estim), nrow = width))

pred <- big_prodMat(h2.out,beta.out)
length(pred)
length(test.bigSNP$fam$affection)
AUC(pred, test.bigSNP$fam$affection)

