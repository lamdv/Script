setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(bigreadr)

# Set parameters

# Declare dataset here
data_train <- "data_train"
data_test <- "data_test"

# Read training example
snp_readBed(bedfile = paste(data_train,"bed",sep="."), backingfile = data_train)
obj.bigSNP <- snp_attach(paste(data_train,"rds",sep = "."))

ind.excl <- snp_indLRLDR(obj.bigSNP$map$chromosome, obj.bigSNP$map$physical.pos)

ind.keep <- snp_clumping(obj.bigSNP$genotypes, obj.bigSNP$map$chromosome, thr.r2 = 0.05, exclude = ind.excl)

obj.svd <- big_randomSVD(obj.bigSNP$genotypes, fun.scaling = snp_scaleBinom(), 
                         ind.col = ind.keep,
                         ncores = 4)

plot(obj.svd, type = "score")
y01 <- obj.bigSNP$fam$affection

gwas.train <- big_univLogReg(obj.bigSNP$genotypes, y01.train = y01)

thrs <- seq(0, 4, by = 0.5)
nb.pred <- sapply(thrs, function(thr) sum(lpS.keep > thr))

obj.bigSNP <- snp_attach(paste(data_test,"rds",sep = "."))