setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(lassosum)
library(data.table)
library(sigmoid)
library(neuralnet)
library(keras)
library(ranger)


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

# Oversampling
y_case <- which(train$fam$affection %in% c(1))

append(ds,PRS)
ds <- PRS[y_case,]

ds<-PRS[]
y.train <- train$fam$affection
for (i in range(1, 10)){
  y.train <- append(y.train, train$fam$affection[y_case])
  ds <- rbind(ds, PRS[y_case,])
}
# 
# model <- keras_model_sequential() 
# model %>%
#   layer_dense(units = 1000, activation = "softmax", input_shape = c(35707), kernel_regularizer=regularizer_l1(0.02)) %>% 
#   
#   # layer_dense(units = 1000, activation = "softmax", input_shape = c(2800)) %>% 
#   # layer_dropout(rate=0.1) %>%
#   layer_dense(units = 2, activation = "softmax", kernel_regularizer=regularizer_l1(0.01))
# 
# model %>% compile(
#   #loss = 'categorical_crossentropy',
#   loss = 'mean_absolute_percentage_error',
#   optimizer = 'SGD',
#   metrics = c('accuracy')
# )
# PRS$show()
#y.train <-  to_categorical(y.train)

# svd <- big_randomSVD(G.train, snp_scaleBinom())
# pca <- big_prodMat(G.train, svd$scale, ind.col = keeps)
# 
# history <- model$fit(
#   ds, y.train,
#   class_weight = list(1, 10),
#   epochs = as.integer(10), batch_size = as.integer(100)
# )

# keeps <- snp_clumping(G.train, CHR)
# G.train[,keeps]
# history <- model$fit(
#   G.train[,keeps], to_categorical(train$fam$affection),
#   class_weight = list(1, 1000),
#   epochs = as.integer(100), batch_size = as.integer(100)
# )

snp_readBed("data_test.bed")
test <- snp_attach("data_test.rds")
G.test <- test$genotypes
PRS_test <- snp_grid_PRS(G.test, all_keep = all_keep, betas = sumstats$beta,lpval)
pred <- model$predict(PRS_test[,])
pred_prob <- model$predict_proba(PRS_test[,])
AUC(pred = model$predict(PRS_test[,])[,1], test$fam$affection)
h1.out
# Reference
M <- snp_grid_stacking(multi_PRS = PRS, y.train = train$fam$affection, ncores = NCORES)
beta <- as_FBM(matrix(M$beta.G))
pred.SCT <- big_prodMat(G.test, beta)
AUC(pred = pred, test$fam$affection)
