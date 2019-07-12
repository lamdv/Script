setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(bigreadr)
library(sigmoid)
# Set parameters

width = 100000


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

G <- PRS
G$add_columns(1)
G[,G$ncol] <- 1
# Generate input weight
# Gaussian input weight; sum = 0
W <- FBM(G$ncol, width)
W[] <- runif(length(W))
W$save()
W <- big_attach("weight.rds")
W$show()
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

# # ELM Steps
# h1 <- as_FBM(big_apply(obj.bigSNP$genotypes, a.FUN = function(X, ind){
#   X[, ind] <- X[, ind] - 0.5
# }))
h1 <- as_FBM(big_prodMat(G,W))
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
beta.out.ELM <- as_FBM(matrix(unlist(big_univLogReg(h1.out, y01.train = train$fam$affection)$estim), nrow = width))
beta.out <- as_FBM(matrix(unlist(big_univLinReg(G, y.train = train$fam$affection)$estim), nrow = width))

#########
# Validation step
#########

snp_readBed("data_test.bed")
test <- snp_attach("data_test.rds")
G.test <- test$genotypes
pred.h1 <- snp_grid_PRS(G.test, all_keep = all_keep, betas = sumstats$beta,lpval)
pred.h1$add_columns(1)
pred.h1[pred.h1$ncol] <- 1
pred.h1 <- as_FBM(big_prodMat(pred.h1, W[]))
pred <- big_prodMat(pred.h1, beta.out.ELM)
1-AUC(pred = pred, test$fam$affection)
h1.out
# Reference
M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
beta <- as_FBM(matrix(M$beta.G))
pred.SCT <- big_prodMat(G.test, beta)
AUC(pred = pred.SCT, test$fam$affection)
