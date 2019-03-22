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
G <- obj.bigSNP$genotypes

ind.excl <- snp_indLRLDR(obj.bigSNP$map$chromosome, obj.bigSNP$map$physical.pos)

ind.keep <- snp_clumping(G, obj.bigSNP$map$chromosome, thr.r2 = 0.05, exclude = ind.excl)

obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), 
                         ind.col = ind.keep,
                         ncores = 4)

plot(obj.svd, type = "score")
y01 <- obj.bigSNP$fam$affection

obj.gwas <- big_univLogReg(obj.bigSNP$genotypes, y01.train = y01,
                           ncores = 4)
#thrs <- seq(0, 4, by = 0.5)
#nb.pred <- sapply(thrs, function(thr) sum(lpS.keep > thr))

summary(lpS.keep <- -predict(gwas.train)[ind.keep])



obj.bigSNP <- snp_attach(paste(data_test,"rds",sep = "."))