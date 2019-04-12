library(lassosum)
library(data.table)
library(bigstatsr)
setwd("~/Documents/simus")

ss <- fread("sumstats.txt")

test.bfile <- "data_test"
train.bfile <- "data_train"
LDblocks <- "EUR.hg19"

cor <- p2cor(p = ss$pval, n = 8000, sign=ss$beta)

out <- lassosum.pipeline(cor =cor, chr=ss$chromosome, pos = ss$physical.pos, 
                         A1 = ss$allele1, A2 = ss$allele2,
                         ref.bfile = train.bfile, test.bfile = test.bfile,
                         LDblocks = LDblocks)

v <- validate(out)

out2 <- subset(out, s=v$best.s, lambda = v$best.lambda)
v2 <- validate(out2)
v2$best.validation.result
AUC(v2$best.pgs, v$pheno)