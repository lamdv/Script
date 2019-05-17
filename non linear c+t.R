setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(lassosum)
library(data.table)
library(sigmoid)
library(neuralnet)


sumstats <- bigreadr::fread2("sumstats.txt")

train <- snp_attach("data_train.rds")
train <- snp_readBed("data_train.bed")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)

# first layer compose of C+T features
all_keep <- snp_grid_clumping(G.train, CHR, POS, lpval, ncores = NCORES)
saveRDS(all_keep, "all_keep.rds")
all_keep <- readRDS("all_keep.rds")
head(all_keep)
PRS <-snp_grid_PRS(G.train,all_keep = all_keep, betas = sumstats$beta,lpval)
head(PRS)
data <- PRS$bm()

# # Non linear activation fuction
h1.activation <- as_FBM(matrix(unlist(big_apply(PRS, a.FUN = function(X, ind){
  sigmoid(X[,ind])})),
  nrow = PRS$nrow, ncol = PRS$ncol))
attr(h1.activation, "lpS") <- attr(PRS,"lpS")
train$fam$affection

h1.out <- big_univLinReg(h1.activation, y.train = train$fam$affection)
# h1.out <- big_univLogReg(PRS, y01.train = train$fam$affection)

# Neural net learning on output
head(data)
f <- as.formula(paste("affection",paste() ))

nn <- neuralnet()

# 
# M <- snp_grid_stacking(multi_PRS = h1.activation, y.train = train$fam$affection, ncores = NCORES)
saveRDS(h1.out, "output_weight.rds")
# 
test <- snp_attach("data_test.rds")
G.test <- test$genotypes
pred.h1 <- snp_grid_PRS(G.test, all_keep = all_keep, betas = sumstats$beta,lpval)
nrow(matrix(h1.out$estim))
pred <- big_prodMat(pred.h1, matrix(h1.out$estim))
AUC(pred = pred, test$fam$affection)
h1.out
# Reference
M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
beta <- as_FBM(matrix(M$beta.G))
pred.SCT <- big_prodMat(G.test, beta)
AUC(pred = pred.SCT, test$fam$affection)
