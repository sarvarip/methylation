# Copyright (C) 2019 Peter Sarvari, University of Southern California
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
stop('regression.R <methylation-data> <age-data> <alpha>',
     call.=FALSE)
}
methylation.filename <- args[1]
metadata.filename <- args[2]
alpha.lasso <- args[3]

library(glmnet)
df <- read.table(methylation.filename, header=T)
rownames(df) <- gsub("X", "", rownames(df))
if (length(grep(".rds", metadata.filename))>0){
df_s <- readRDS(metadata.filename)
} else {
df_s <- read.csv(file=metadata.filename, header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Male.BMI) | df_s$Male.BMI==""), ]
}
if ("Laboratory.ID" %in% colnames(df_s)) {
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL
}

data <- merge(df, df_s, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

get.rsq <- function(y, y.hat) {
  sst <- sum((y - mean(y))^2)
  sse <- sum((y.hat - y)^2)
  return(1 - sse/sst)
}

cross.val.fold <- 10
n.individuals <- nrow(data)
print("Number of samples")
print(n.individuals)
parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))

rsq_list <- rep(0,cross.val.fold)
varexp_train <- rep(0,cross.val.fold)
varexp_list <- rep(0,cross.val.fold)
mse_list <- rep(0,cross.val.fold)
if (alpha.lasso > 0) {
zeroed_list <- rep(0,cross.val.fold)
}

for (i in 1:cross.val.fold) {

training.individuals <- rownames(data)[which(parts != i)]
testing.individuals <- rownames(data)[which(parts == i)]
X_ <- data[training.individuals,]
Y <- as.matrix(X_$Male.BMI)
X_$Male.BMI <- NULL
X <- as.matrix(X_)

cv0 = cv.glmnet(X,Y,type.measure="mse",alpha=alpha.lasso,
                standardize=T)
if (alpha.lasso > 0) {
cs <- coef(cv0, s='lambda.min')
zeroed <- sum(cs == 0)
zeroed_list[i] <- zeroed
}

X.test <- data[testing.individuals,]
y_real <- as.matrix(X.test$Male.BMI)
X.test$Male.BMI <- NULL
x_pred <- as.matrix(X.test)

y_pred <- predict(cv0, newx = x_pred, s="lambda.min")
y.train.pred <- predict(cv0, newx=X, s="lambda.min")
result.lm = lm(y_real ~ y_pred)
print(cbind(y_real, y_pred))
rsq_list[i] <- summary(result.lm)$r.squared
mse_list[i] <- mean((y_real - y_pred)^2)
varexp_list[i] <- get.rsq(y_real, y_pred)
varexp_train[i] <- get.rsq(Y, y.train.pred)

}

print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Std dev of R squared: ", sd(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Std dev of MSE: ", sd(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Std dev of variance explained: ", sd(varexp_list)))
print(paste("Training set: Mean of variance explained: ", 
            mean(varexp_train)))
print(paste("Training set: Std dev of variance explained: ", 
            sd(varexp_train)))
print(paste("Number of features: ", ncol(data)-1))
if (alpha.lasso > 0) {
print(paste("Average number of zeroed features: ",
            mean(zeroed_list)))
}



