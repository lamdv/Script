setwd("~/Documents/simus")

library(bigstatsr)
library(bigsnpr)
library(lassosum)
library(data.table)
library(sigmoid)
library(xgboost)

# Train data
sumstats <- bigreadr::fread2("sumstats.txt")
snp_readBed("data_train.bed")
train <- snp_attach("data_train.rds")
G.train <- train$genotypes
CHR <- train$map$chromosome
POS <- train$map$physical.pos
NCORES <- nb_cores()
lpval <- -log10(sumstats$pval)
y.train <- train$fam$affection

# Test data
snp_readBed("data_test.bed")
test <- snp_attach("data_test.rds")
G.test <- test$genotypes
PRS_test <- snp_grid_PRS(G.test, all_keep = all_keep, betas = sumstats$beta,lpval)

# first layer compose of C+T features
all_keep <- snp_grid_clumping(G.train, CHR, POS, lpval, ncores = NCORES)
saveRDS(all_keep, "all_keep.rds")
all_keep <- readRDS("all_keep.rds")
head(all_keep)
PRS <-snp_grid_PRS(G.train,all_keep = all_keep, betas = sumstats$beta,lpval)
saveRDS(PRS, "stack.rds")
PRS <- readRDS('stack.rds')
head(PRS)
PRS[1:PRS$nrow]

# Oversampling
y_case <- which(train$fam$affection %in% c(1))

ds<-PRS[]
y.train <- train$fam$affection
for (i in range(1, 10)){
  y.train <- append(y.train, train$fam$affection[y_case])
  ds <- rbind(ds, PRS[y_case,])
}

ds.FBM <- as_FBM(ds)


bstSparse <- xgboost(data = ds[], 
                     label = y.train, 
                     booster="gblinear",
                     max_depth = 6, 
                     eta = 1, 
                     nthread = 2, 
                     nrounds = 2, 
                     lambda = 0.1,
                     objective = "count:poisson")
pred <- predict(bstSparse, PRS_test[])
AUC(pred = pred, test$fam$affection)


AUCs <- matrix(,nrow = 10, ncol = 10)
for (i in seq(2 : 10)){
  for (j in seq(2: 10)) {
    temp <- c()
    for(t in seq(1:10)){
      bstSparse <- xgboost(data = ds[], 
                           label = y.train, 
                           booster="gblinear",
                           max_depth = i, 
                           eta = 1, 
                           nthread = 2, 
                           nrounds = j, 
                           lambda = 0.1,
                           objective = "count:poisson",
                           )
      pred <- predict(bstSparse, PRS_test[])
      temp <- append(temp,AUC(pred = pred, test$fam$affection))
    }
    AUCs[i,j] <- mean(temp)
  }
}
AUCs
boxplot(AUCs)
h1.out
# Reference
M <- snp_grid_stacking(multi_PRS = ds.FBM, y.train = y.train, ncores = NCORES)
M <- snp_grid_stacking(multi_PRS = ds.FBM, y.train = y.train)
beta <- as_FBM(matrix(M$beta.G))
pred.SCT <- big_prodMat(G.test, beta)
AUC(pred = pred, test$fam$affection)



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

# PRS.df <- as.data.frame(PRS[])
# cols <- colnames(PRS.df)
# PRS.df$newcolumn <- train$fam$affection
# colnames(PRS.df) <- append(cols, "y")
# PRS.df$y

