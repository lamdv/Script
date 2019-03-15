library(lassosum)
library(data.table)
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
