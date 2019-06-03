library(lassosum)
library(data.table)
library(bigstatsr)
library(bigsnpr)
library(nnet)

#Lassosum step
setwd("~/Documents/simus")

ss <- fread("sumstats.txt")

test.bfile <- "data_test"
train.bfile <- "data_train"
obj.bigSNP <- snp_attach(paste(train.bfile,"rds",sep = "."))
LDblocks <- "EUR.hg19"

train <- snp_attach("data_train.rds")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(ss$pval)


cor <- p2cor(p = ss$pval, n = 8000, sign=ss$beta)

out <- lassosum.pipeline(cor =cor, chr=ss$chromosome, pos = ss$physical.pos, 
                         A1 = ss$allele1, A2 = ss$allele2,
                         ref.bfile = train.bfile, test.bfile = test.bfile,
                         LDblocks = LDblocks)

###
#Extract result

v <- validate(out)
# Original Result
out2 <- subset(out, s=v$best.s, lambda = v$lambda)
v2 <- validate(out2)
v2$best.validation.result
AUC(v$best.pgs, v$pheno)
# 
###
#Stacking step

nn.inp <- as_FBM(out$beta[[v$best.s[1]]])
keep <- which(out$test.extract == TRUE)
inp <- big_copy(G.train, ind.col = keep)
nn.inp
res <- as_FBM(big_prodMat(inp, nn.inp))
beta.out <- big_univLogReg(res, y01.train = train$fam$affection)

# Validation

test <- snp_attach("data_test.rds")
G.test <- big_copy(test$genotypes, ind.col = keep)
pred.h1 <- as_FBM(big_prodMat(G.test, nn.inp))
pred <- big_prodVec(pred.h1, beta.out$estim)
AUC(pred = pred, test$fam$affection)
err <- pred -test$fam$affection
err <- sqrt(err %*% err) / nrow(pred)
  