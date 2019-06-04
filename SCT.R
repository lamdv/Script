setwd("~/Documents/simus")

library(bigsnpr)
library(lassosum)
library(data.table)
library(bigstatsr)
# library(bigsnpr)

sumstats <- bigreadr::fread2("sumstats.txt")

train <- snp_attach("data_train.rds")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)
seeds <- c(1:1000)
test <- snp_attach("data_test.rds")
G.test <- test$genotypes

# G.train[1] <- 0.001

all_keep <- snp_grid_clumping(G.train, CHR, POS, lpval, ncores = NCORES)
PRS <-snp_grid_PRS(G.train,all_keep = all_keep, betas = sumstats$beta,lpval)
PRS

aucs <- list()
for(seed in seeds){
  set.seed(seed)
  M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
  # saveRDS(M, "stack.rds")
  beta <- as_FBM(matrix(M$beta.G))
  pred <- big_prodMat(G.test, beta)
  aucs[seed] <- AUC(pred = pred, test$fam$affection)
  aucs[seed]
}
boxplot(unlist(aucs))

