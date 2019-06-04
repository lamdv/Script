setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(lassosum)
library(data.table)
library(sigmoid)
library(neuralnet)


sumstats <- bigreadr::fread2("sumstats.txt")
snp_readBed("data_train.bed")
train <- snp_attach("data_train.rds")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)
y.train <- train$fam$affection

# first layer compose of C+T features
all_keep <- snp_grid_clumping(G.train, CHR, POS, lpval, ncores = NCORES)
saveRDS(all_keep, "all_keep.rds")
all_keep <- readRDS("all_keep.rds")
head(all_keep)
PRS <-snp_grid_PRS(G.train,all_keep = all_keep, betas = sumstats$beta,lpval)
head(PRS)
PRS$expo

# Neural net learning on output
# single layer, 50 neuron
shape <- c(1,50)
# Init weights
w <- list ()
w[[1]] <- matrix(runif(len, max = 1/PRS$ncol), nrow = PRS$ncol)
for (i in 2:(shape[1]+1)){
  len <- nrow(w[i - 1]) * shape[2]
  w[[i]] <- matrix(runif(len, max = 1/nrow(w[i - 1])), nrow = nrow(w[i-1]))
}
matrix(unlist(w[1]), nrow = PRS$col)
# Feed forward
h <- list()
w[[1]]
PRS[1]
h[[1]] <- sigmoid(big_prodMat(PRS, as_FBM(w[[1]])))
h[[1]]
for (i in 2:(shape[1]+1)){
  sigmoid(big_prodMat(PRS, as_FBM(w[[1]])))
} 

# 
# M <- snp_grid_stacking(multi_PRS = h1.activation, y.train = train$fam$affection, ncores = NCORES)
saveRDS(h1.out, "output_weight.rds")
# 
snp_readBed("data_test.bed")
test <- snp_attach("data_test.rds")
G.test <- test$genotypes
pred.h1 <- snp_grid_PRS(G.test, all_keep = all_keep, betas = sumstats$beta,lpval)
pred <- big_prodMat(pred.h1, matrix(h1.out$estim))
AUC(pred = pred, test$fam$affection)
h1.out
# Reference
M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
beta <- as_FBM(matrix(M$beta.G))
pred.SCT <- big_prodMat(G.test, beta)
AUC(pred = pred.SCT, test$fam$affection)
