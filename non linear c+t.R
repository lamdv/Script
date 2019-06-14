setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(lassosum)
library(data.table)
library(sigmoid)
library(neuralnet)
library(keras)


sumstats <- bigreadr::fread2("sumstats.txt")
snp_readBed("data_train.bed")
train <- snp_attach("data_train.rds")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)
y.train <- train$fam$affection

# first layer compose of C+T features
all_keep <- snp_grid_clumping(G.train, CHR, POS, lpval, ncores = NCORES)
saveRDS(all_keep, "all_keep.rds")
all_keep <- readRDS("all_keep.rds")
head(all_keep)
PRS <-snp_grid_PRS(G.train,all_keep = all_keep, betas = sumstats$beta,lpval)
head(PRS)
PRS_train <- matrix(PRS[,], nrow = PRS$nrow)
head(PRS_train)
PRS[1:PRS$nrow]
# Neural net learning on output
# single layer, 50 neuron

model <- keras_model_sequential() 
model %>%
  layer_dense(units = 256, activation = "relu", input_shape = c(2800)) %>% 
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 10, activation = "softmax")

model %>% compile(
  # loss = 'categorical_crossentropy',
  loss = 'sparse_categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)
PRS$show()
y.train <-  to_categorical(train$fam$affection)
y.train <- train$fam$affection

history <- model$fit(
  PRS[], y.train,
  epochs = 10, batch_size = 64
)

shape(PRS[])
G.train
G.test
# 
# M <- snp_grid_stacking(multi_PRS = h1.activation, y.train = train$fam$affection, ncores = NCORES)
saveRDS(h1.out, "output_weight.rds")
# 
snp_readBed("data_test.bed")
test <- snp_attach("data_test.rds")
G.test <- test$genotypes
pred.h1 <- snp_grid_PRS(G.test, all_keep = all_keep, betas = sumstats$beta,lpval)
pred <- big_prodMat(pred.h1, matrix(h1.out$estim))
AUC(pred = pred, test$fam$affection)
h1.out
# Reference
M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
beta <- as_FBM(matrix(M$beta.G))
pred.SCT <- big_prodMat(G.test, beta)
AUC(pred = pred.SCT, test$fam$affection)
