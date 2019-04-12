setwd("~/Documents/simus")

library(bigsnpr)

sumstats <- bigreadr::fread2("sumstats.txt")

train <- snp_attach("data_train.rds")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)

all_keep <- snp_grid_clumping(G.train, CHR, POS, lpval, ncores = NCORES)
saveRDS(all_keep, "all_keep.rds")
PRS <-snp_grid_PRS(G.train,all_keep = all_keep, betas = sumstats$beta,lpval)
PRS
M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
saveRDS(M, "stack.rds")

beta <- as_FBM(matrix(M$beta.G))

test <- snp_attach("data_test.rds")
G.test <- test$genotypes
pred <- big_prodMat(G.test, beta)
AUC(pred = pred, test$fam$affection)
