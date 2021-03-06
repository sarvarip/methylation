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

# Rscript regression_v2.R ../input/AgingFull_Human_Sperm_hmr.txt
# ../input/Human_Sperm_hmr.txt ../input/AgingFullMeta.rds 
# ../input/agedat.tsv 1

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
stop('regression.R <methylation-data> <independent-data> 
     <age-data> <independent-age-data>' ,
     call.=FALSE)
}
methylation.filename <- args[1]
newdata.filename <- args[2]
metadata.filename <- args[3]
newmetadata.filename <- args[4]

read.meta <- function(metadata.filename){
if (length(grep(".rds", metadata.filename))>0){
df.m <- readRDS(metadata.filename)
} else {
df.m <- read.csv(file=metadata.filename, header=TRUE, sep="\t")
}
if ("Laboratory.ID" %in% colnames(df.m)) {
rownames(df.m) <- df.m$Laboratory.ID
df.m$Laboratory.ID <- NULL
df.m$Age <- df.m$Male.Age
df.m$Male.Age <- NULL
}
return (df.m)
}

library(pls)
df <- read.table(methylation.filename, header=T)
id.df <- read.table(newdata.filename, header=T)
rownames(df) <- gsub("X", "", rownames(df))
rownames(id.df) <- gsub("X", "", rownames(id.df))
df.m <- read.meta(metadata.filename)
id.df.m <- read.meta(newmetadata.filename)

data <- merge(df, df.m, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

id.data <- merge(id.df, id.df.m, by=0)
rownames(id.data) <- id.data$Row.names
id.data$Row.names <- NULL

col.names.data <- colnames(data)
col.names.id.data <- colnames(id.data)
col.names.same <- intersect(col.names.data, col.names.id.data)
id.data <- id.data[, colnames(id.data) %in% col.names.same]
data <- data[, colnames(data) %in% col.names.same]

print("Number of features in D1")
print(ncol(data))
print("Number of features in D2")
print(ncol(id.data))

get.rsq <- function(y, y.hat) {
  sst <- sum((y - mean(y))^2)
  sse <- sum((y.hat - y)^2)
  return(1 - sse/sst)
}

cross.val.fold <- 10
n.individuals <- nrow(data)
print("Number of samples")
print(n.individuals)

set.seed(1)
pls.fit <- plsr(Age~., data=data, ncomp=25, scale=T, validation='CV')
ncomp <- which.min(RMSEP(pls.fit)$val[1,1,])-1

X.test <- id.data
y_real <- as.matrix(id.data$Age)

X.test$Age <- NULL
x_pred <- as.matrix(X.test)

print("Predicting on second dataset using coefficients 
      and optimal number of components")

y_pred <- predict(pls.fit, x_pred, ncomp=ncomp)

dim(y_pred) <- c(nrow(y_pred), 1)

y.train.pred <- predict(pls.fit, data, ncomp=ncomp)

result.lm = lm(y_pred ~ y_real)
print(cbind(y_real, y_pred))

figure.filename <- sprintf('%s_pls_nocenter.pdf', newdata.filename)

rsq_list <- summary(result.lm)$r.squared
mse_list <- mean((y_real - y_pred)^2)
varexp_list <- get.rsq(y_real, y_pred)
varexp_train <- get.rsq(data$Age, y.train.pred)

pdf(figure.filename, width=6.5, height=8, paper='US')
template <- "Independent data: \n %s \n R2=%.3f"
plot(y_real, y_pred, xlab="Actual age", ylab="Predicted age",
     main=sprintf(template, newdata.filename, varexp_list))
abline(result.lm)
abline(a=0, b=1, col="red", lty=2)
dev.off()

print(paste("Mean R squared: ", mean(rsq_list)))
print(paste("Mean MSE: ", mean(mse_list)))
print(paste("Mean of variance explained: ", mean(varexp_list)))
print(paste("Training set: Mean of variance explained: ", 
            mean(varexp_train)))
print(paste("Number of features: ", unname(ncomp)))



