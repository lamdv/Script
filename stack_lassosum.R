library(lassosum)
library(data.table)
library(bigstatsr)
library(bigsnpr)

#Lassosum step
setwd("~/Documents/simus")

ss <- fread("sumstats.txt")

test.bfile <- "data_test"
train.bfile <- "data_train"
obj.bigSNP <- snp_attach(paste(train.bfile,"rds",sep = "."))
LDblocks <- "EUR.hg19"

cor <- p2cor(p = ss$pval, n = 8000, sign=ss$beta)

out <- lassosum.pipeline(cor =cor, chr=ss$chromosome, pos = ss$physical.pos, 
                         A1 = ss$allele1, A2 = ss$allele2,
                         ref.bfile = train.bfile, test.bfile = test.bfile,
                         LDblocks = LDblocks)

###
#Extract result

# v <- validate(out)
# 
# out2 <- subset(out, s=v$best.s, lambda = v$lambda)
# v2 <- validate(out2)
# v2$best.validation.result
# AUC(v$best.pgs, v$pheno)
# 
###
#Stacking step

beta <- out2$beta['0.5']
keep <- subset(range(0, size))

keep_ind <- out$also.in.refpanel

keep_ind
# keep <- obj.bigSNP$genotypes$as.FBM()[keep_ind]
keep <- obj.bigSNP$genotypes$bm()

keep

keep <- keep[drop = - keep_ind]
