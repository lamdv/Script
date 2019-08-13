setwd("~/Documents/simus")

library(bigsnpr)
library(lassosum)
library(data.table)
library(bigstatsr)
library(sigmoid)
library(neuralnet)
# library(bigsnpr)

sumstats <- bigreadr::fread2("sumstats.txt")

train <- snp_attach("data_train.rds")
G.train <- train$genotypes
y.train <- train$fam$affection
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)
ths <- 10 # Only take into consideration p-values < 10^-10

head(lpval)
keep <- subset(c(1: nrow(sumstats)),
               lpval > ths)
inp <- big_copy(G.train, ind.col = keep)

# TODO: implement Backprop on neural net
# neuralnet library do not support FBM so new implementation is needed

