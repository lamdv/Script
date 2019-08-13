library(bigstatsr)
library(bigsnpr)

big_neuralnet <- setRefClass("big_neuralnet", 
                             fields = c ("shape", "data.train", "y.train",
                                         "data.test", "y.test", "func", "weights"))

big_neuralnet$methods(initialize = function(shape, data.train, y.train,
                                         data.test, y.test, func) 
{
  shape <<- shape
  data.train <<- data.train
  y.train <<- y.train
  data.test <<- data.test
  y.test <<- y.test
  func <<- func
  weight <- list()
  weight[1] <- as_FBM(matrix(rnorm(data.train$ncol*shape[2],mean=0,sd = 1/nrow(data.train)), 
                             ncol = shape[2], nrow = data.train$ncol
                             ))
})

data_train <- 

nn <- big_neuralnet()