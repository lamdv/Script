library(lassosum)
library(data.table)
library(bigstatsr)

#Lassosum step
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

###
#Extract result

v <- validate(out)

out2 <- subset(out, s=v$best.s, lambda = v$lambda)
v2 <- validate(out2)
v2$best.validation.result
AUC(v$best.pgs, v$pheno)

###
#Stacking step

beta <- out2$beta['0.5']
keep <- subset(range(0, size))

grid <- beta
