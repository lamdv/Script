---
title: "Survey on xgboost's objective functions"
author: "VuLamDANG"
date: "August 26, 2019"
output: 
  html_document:
    keep_md: true
editor_options: 
  chunk_output_type: console
---


## Prequisites

The following packages are required to knit this document:

- bigstatsr
- bigsnpr
- data.table
- xgboost



## Dataset

*Please change working direction in `directory`*

I use the "Simus" simulated dataset provided by Florian. This dataset contain ~650,000 SNPs in 2 chromosomes. 20% of them are cases and the rest 80% are controls.




## C+T step
Generating a matrix of C+T with varied clumping radius and $p-value$ thresholds.

$\beta s$ and $p-values$ are aquired in sumstats file (see `sumstats.txt`)



XGBoost only accept dataframe, so all FBM output must be converted to data frame



## Data treatment: Oversampling

To oversampling I extract all positive (*affection = 1*), and append them a number of times. For this simulated dataset, I enhanced the number of cases 10 times, thus the ratio become ~70:30 (from 20:80).



##XGBoost

XGBoost provide slightly better result compare to SCT. However this result heavily depend on booster and objective selection. Experiments shown that **gblinear** booster with **count:poisson** objective give the best results.

For the first experimentation, I run 9 different *max_depth* and 9 *nrounds* with *count:poisson* objective function. Each combination is repeated 10 times and the mean value calculated to rule out the randomess of *gblinear* booster

```r
poisson_AUCs <- matrix(,nrow = 10, ncol = 10)
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
                           objective = "count:poisson"
                           )
      pred <- predict(bstSparse, PRS_test[])
      temp <- append(temp,AUC(pred = pred, test$fam$affection))
    }
    poisson_AUCs[i,j] <- mean(temp)
  }
}
```

```
##            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
##  [1,] 0.7817714 0.7823509 0.7823041 0.7822926 0.7825028 0.7824836
##  [2,] 0.7820513 0.7822373 0.7823633 0.7824287 0.7822630 0.7824896
##  [3,] 0.7814487 0.7821501 0.7823655 0.7826146 0.7825234 0.7825795
##  [4,] 0.7820921 0.7824120 0.7824182 0.7826145 0.7825073 0.7824446
##  [5,] 0.7821936 0.7824070 0.7821748 0.7823153 0.7826239 0.7824142
##  [6,] 0.7819266 0.7822376 0.7822175 0.7824563 0.7823541 0.7824829
##  [7,] 0.7818918 0.7820830 0.7823140 0.7823306 0.7824132 0.7824501
##  [8,] 0.7820375 0.7822164 0.7824588 0.7824749 0.7825382 0.7825095
##  [9,] 0.7822184 0.7819440 0.7825596 0.7823647 0.7824570 0.7823492
## [10,]        NA        NA        NA        NA        NA        NA
##            [,7]      [,8]      [,9] [,10]
##  [1,] 0.7823923 0.7821336 0.7818004    NA
##  [2,] 0.7823814 0.7821682 0.7817430    NA
##  [3,] 0.7823602 0.7821881 0.7817512    NA
##  [4,] 0.7823655 0.7820304 0.7818152    NA
##  [5,] 0.7822578 0.7821342 0.7818768    NA
##  [6,] 0.7824436 0.7822035 0.7817998    NA
##  [7,] 0.7824093 0.7822224 0.7817339    NA
##  [8,] 0.7824048 0.7821033 0.7817913    NA
##  [9,] 0.7823197 0.7822398 0.7817214    NA
## [10,]        NA        NA        NA    NA
```

![Poisson objective](xgboost_fig/result_poisson-1.png)

I omitted the code for the next 2 experiments as they are largely similar, with exception of the option *objective* for *xgboost* function.

For the next experimentation, I run 9 different *max_depth* and 9 *nrounds* with *reg:logistic* objective function. Each combination is repeated 10 times and the mean value calculated to rule out the randomess of *gblinear* booster


![Logistic objective](xgboost_fig/result_log-1.png)

For the final experimentation, I run 9 different *max_depth* and 9 *nrounds* with *binary:logicraw* objective function. This function is interesting since it's specific for binary classification.


![Binary Logistic objective](xgboost_fig/result_bilog-1.png)

## Comparison

Here is a composite boxplot comparing different objective functions. Index 1 is *count:poisson*, 2 is *reg:logistic* and 3 is *binary:logistic*

![Comparision between Objective Functions](xgboost_fig/comp-1.png)
