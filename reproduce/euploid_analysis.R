#!/usr/bin/env Rscript

### Copyright 2019 Andrew D Smith, Peter Sarvari and Timothy Jenkins

### To run from the command line (assuming the file is executable):
###
### ./euploid_analysis.R features.txt metadata.txt methylation.txt Euploid.Rate output
###
### The files are as follows:
###
### features.txt: this is a file that contains a data frame with row
### names equal to the names of the probes to be used in the
### model. These must be a subset of, possibly equal to, the full set
### of methylation features provided. Other columns in this file are
### not used, only the rownames.
###
### metadata.txt: this file contains a data frame with one row for
### each sample and has two columns named Euploid.Rate and
### Female.Age. The column 'Euploid.Rate' will be used as the
### outcome. The column 'Female.Age' will be used as a predictive
### feature along with the methylation features.
###
### methylation.txt: this is a data frame with rows corresponding to
### samples, and columns corresponding to methylation probe
### features. These columns must contain all those indicated in
### 'features.txt'. The sample names must have non-zero intersection
### with those in 'metadata.txt'
###
### Euploid.Rate: This is the name of the outcome feature. It is
### hard-coded, so must be specified exactly on the command line. It
### must also exist in the file 'metadata.txt'
###
### output: this is a prefix for the output files, and can be changed
### to anything, but must not contains slashes or spaces.

## needed if script is run from command line
require(methods)

### "x" is the existing data frame
### "lbl" is the label for the new column
### "training" and "testing" are the R2 values
### and "training.age" and "testing.age" are the
### R2 values for the model only using female age
update.rss.info <- function(x, lbl, training, testing,
                            training.age, testing.age) {
  return(cbind(x, setNames(data.frame(c(training, testing,
                                        training.age, testing.age),
                                      row.names=row.names(x)), lbl)))
}

### "x" is the existing data frame
### "lbl" is the label for the new column
### "coefficients" is the values of coefficients to add
update.coef.info <- function(x, lbl, coefficients) {
  return(cbind(x, setNames(as.data.frame(as.matrix(coefficients)), lbl)))
}

### utility function to merge two data frames as safely as possible
brief.merge <- function(a, b) {
  stopifnot(rownames(a) == rownames(b))
  ## keep original sorted order of rows
  saved.row.names <- rownames(a)
  c <- merge(a, b, by='row.names')
  row.names(c) <- c$Row.names
  c$Row.names <- NULL
  ## make sure the rows are properly sorted
  c <- c[saved.row.names, ]
  return(c)
}

### this function returns the unadjusted RSS
get.rsq <- function(y, y.hat) {
  sst <- sum((y - mean(y))^2)
  sse <- sum((y.hat - y)^2)
  return(1 - sse/sst)
}

### this function returns the unadjusted RSS, too liberal, not recommended
get.lm.rsq <- function(y, y.hat) {
  return(summary(lm(as.matrix(y.hat) ~ as.matrix(y)))$r.squared)
}

### this function is to get the linear model without regularization
### 1) used to pass in female age as variable.name to compare the model
### using all the 13 features to the model only using female age
### 2) also used to produce an 'optimistic' bar graph regressing
### the outcome on the 13 features using the whole data
### (might overfit and hence the 'optimistic' adjective)
lm.fit <- function(X, Y, variable.name) {
  resp.name <- colnames(Y)[1]
  if (length(variable.name) > 1) {
    regression.formula <- as.formula(sprintf('%s ~ %s', resp.name,
                                             paste(variable.name,
                                                   collapse='+')))
  }
  else {
    regression.formula <- as.formula(sprintf('%s ~ %s',
                                             resp.name,
                                             variable.name))
  }
  return(lm(regression.formula,
            as.data.frame(cbind(Y, X[, variable.name, drop=F]))))
}

### the function that does all the work
euploid.analysis <- function(features.filename,
                             metadata.filename,
                             meth.filename,
                             outcome.label,
                             output.basename,
                             cross.val.fold = 10) {


  suppressPackageStartupMessages(library(glmnet, quiet=T))

### PART I: REGRESSION ON AGE AND METHYLATION WITH CROSS VALIDATION

  ## READ THE FEATURES FILE
  ## must be separated by tab and contain the features as the first
  ## line (header) example: features.txt
  cat(sprintf('[features file: %s]\n', features.filename))
  features <- read.table(features.filename, sep='\t', header=T)
  cat(sprintf('[N features: %d]\n', nrow(features)))

  ## READ THE METADATA INCLUDING OUTCOMES
  ## must be separated by tabs, contain a header and
  ## have the row names as the first column
  ## example: metadata.txt
  cat(sprintf('[meta data file: %s]\n', metadata.filename))
  metadata <- read.table(metadata.filename, sep='\t', header=T, row.names=1)
  rownames(metadata) <- gsub("X", "", rownames(metadata))

  ## THIS IS THE OUTCOME LABEL
  ## if using the example datasets, then this must be
  ## Euploid.Rate (no need for quotes)
  cat(sprintf('[outcome label: %s]\n', outcome.label))

  ## READ THE METHYLATION DATA
  ## if you're only working with the 13 features, use
  ## 13_methylation.txt
  ## In a more general case, you can use methylation.txt
  ## that contains much more features, provided that they are
  ## tab separated and have the sample names in the first row
  ## Sample names cannot contain X (must start with a number)
  cat(sprintf('[methylation data file: %s]\n', meth.filename))
  methylation <- read.table(meth.filename, sep='\t', header=T, row.names=1)

  n.individuals <- nrow(methylation)
  cat(sprintf('[N individuals with data: %d]\n', n.individuals))
  cat(sprintf('[N probes: %d]\n', ncol(methylation)))

  ## make sure to only use methylation from specified cpg sites
  ## provided as features (first argument)
  methylation <- methylation[,rownames(features)]

  ## make sure to only use individuals with data for both methylation
  ## and the specified outcome label
  shared.rows <- intersect(rownames(methylation), rownames(metadata))
  n.individuals <- length(shared.rows)
  cat(sprintf('[N individuals with data and outcomes: %d]\n', n.individuals))
  methylation <- methylation[shared.rows, ]
  metadata <- metadata[shared.rows, ]
  stopifnot(rownames(metadata) == rownames(methylation)) # sanity check

  ## set the outcomes, they must be contained in metadata
  outcomes <- metadata[, outcome.label, drop=F]
  ## Make sure to remove the outcomes from potential predictors
  metadata[, outcome.label] <- NULL
  ## create the set of potential predictors
  ## Note: this puts the metadata alongside the methylation
  merged.data <- brief.merge(methylation, metadata)

  ## partition the individuals into an approximately equal sized
  ## subsets for cross-validation
  parts <- sample(rep(1:cross.val.fold, length.out=n.individuals))

  ## create a data frame to store summary about use of merged.data in the
  ## models from each fold of cross-validation
  coef.info <- data.frame(row.names=c('Intercept', colnames(merged.data)))
  rss.info.names <- c('R2train', 'R2test', 'R2agetrain', 'R2agetest')
  rss.info <- data.frame(row.names=rss.info.names)

  ## create the directory for saving the estimates
  ## dir.create(file.path('RDS'), showWarnings = FALSE)

  ## open the file for all the plots
  figures.filename <- sprintf('%s.pdf', output.basename)
  pdf(figures.filename, width=6.5, height=8, paper='US')
  par(mfrow=c(2,2))

  for (i in 1:cross.val.fold) {
    cat(sprintf('[cross-validation iteration %d]\n', i))

    ## extract the training/testing individuals and their outcomes
    cat('[extracting training and testing sets]\n')
    training.individuals <- rownames(merged.data)[which(parts != i)]
    training.merged.data <- as.matrix(merged.data[training.individuals, ])
    training.feat.meta <- training.merged.data[, 'Female.Age', drop=F]
    training.outcomes <- as.matrix(outcomes[training.individuals, , drop=F])
    cat(sprintf('[n training individuals: %d]\n', length(training.individuals)))

    testing.individuals <- rownames(merged.data)[which(parts == i)]
    testing.merged.data <- as.matrix(merged.data[testing.individuals, ])
    testing.feat.meta <- testing.merged.data[, 'Female.Age', drop=F]
    testing.outcomes <- as.matrix(outcomes[testing.individuals, , drop=F])
    cat(sprintf('[n testing individuals: %d]\n', length(testing.individuals)))

    ## just a sanity check to make sure the data is consistent
    stopifnot(rownames(training.merged.data) == rownames(training.outcomes))
    stopifnot(rownames(testing.merged.data) == rownames(testing.outcomes))

    cat('[fitting the model using glmnet]\n')
    ## fit the model using the training data and also fit the model
    ## on only female age
    current.model <- cv.glmnet(training.merged.data, training.outcomes)
    current.model.meta <- lm.fit(training.merged.data,
                                 training.outcomes, 'Female.Age')

    cat('[applying the model on training and testing sets]\n')
    ## prediction using only female age, training set
    training.est.meta <- predict(current.model.meta,
                                 newdata=as.data.frame(training.feat.meta))
    training.rsq.meta <- get.rsq(training.outcomes, training.est.meta)

    ## apply the model to the *training* data, get estimates
    training.estimates <- predict(current.model,
                                  newx=training.merged.data, s='lambda.min')
    training.rsq <- get.rsq(training.outcomes, training.estimates)

    ## apply the model to the *testing* data, get estimates
    testing.estimates <- predict(current.model,
                                 newx=testing.merged.data, s='lambda.min')
    testing.rsq <- get.rsq(testing.outcomes, testing.estimates)

    ## prediction using only female age, testing set
    testing.est.meta <- predict(current.model.meta,
                                newdata=as.data.frame(testing.feat.meta))
    testing.rsq.meta <- get.rsq(testing.outcomes, testing.est.meta)

    ## ADS: skpping this for now (saving current model in RDS file)
    ## saveRDS(current.model, sprintf('RDS/model%d.rds', i))

    ## get the coefficients out of the current model
    current.coef <- coef(current.model, s='lambda.min')

    ## update the data frames with info from current fold of cross
    ## validation
    fold.label <- sprintf('model%d', i)
    coef.info <- update.coef.info(coef.info, fold.label, current.coef)
    rss.info <- update.rss.info(rss.info, fold.label,
                                training.rsq, testing.rsq,
                                training.rsq.meta, testing.rsq.meta)

    ## plot the current model
    the.y <- 'predicted'
    the.x <- 'observed'
    title.template.train <- 'model %d training (R2=%.3f)'
    title.template.test <- 'model %d testing (R2=%.3f)'
    plot(training.outcomes, training.estimates,
         main=sprintf(title.template.train, i, training.rsq),
         xlab=the.x, ylab=the.y, pch=19, cex=0.2)
    abline(lm(training.estimates ~ training.outcomes))
    plot(testing.outcomes, testing.estimates,
         main=sprintf(title.template.test, i, testing.rsq),
         xlab=the.x, ylab=the.y, pch=19, cex=0.2)
    abline(lm(testing.estimates ~ testing.outcomes))
  }
  ## dev.off()

  ## get the summaries of how many times each feature was used, and
  ## its mean value when used
  coef.not.zero <- rowSums(coef.info != 0)
  coef.info.used <- as.data.frame(coef.not.zero)
  coef.info.mean <- as.data.frame(rowSums(coef.info)/max(1, coef.not.zero))
  coef.info <- cbind(coef.info, setNames(coef.info.mean, 'mean'))
  coef.info <- cbind(coef.info, setNames(coef.info.used, 'used'))

  ## add dummy NA values to rss.info so it can be merged with
  ## coef.info; calculate mean RSS for the testing and training
  padding <- data.frame(matrix(NA, 4, 2), row.names=rss.info.names)
  for (i in rss.info.names) {
    padding[i, 1] <- mean(as.numeric(rss.info[i, ]))
  }
  rss.info <- cbind(rss.info, setNames(padding, c('mean', 'used')))

  ## now merge the summary data for the rss and for the coefficients
  all.info <- rbind(coef.info, rss.info)

  ## finally, write the summaries for each feature to a text file
  summary.filename <- sprintf('%s.txt', output.basename)
  write.table(format(all.info, digits=4), summary.filename, quote=F)

### PART II: PRODUCING BAR GRAPHS FOR PREDICTION CATEGORIES

  ## fitting the whole data without regularization
  the.model <- lm.fit(merged.data, outcomes, colnames(merged.data))

  ## predict on the whole data
  predictions <- predict(the.model, merged.data)
  all.data <- cbind(merged.data, 'pred'=predictions)
  all.data <- cbind(all.data, 'Euploid.Rate'=outcomes)

  ## define the breaks and labels for the age groups
  age.groups <- c('under 30', '30-35', '35-40', '40 and over')
  age.breaks <- c(-Inf, 30, 35, 40, Inf)
  n.ages <- length(age.groups)

  ## define the labels for the prediction group
  pred.groups <- c('poor', 'borderline', 'normal', 'good', 'ideal')
  pred.breaks <- c(-Inf, 0, 0, 0, 0, Inf) # middle 4 must be learned
  n.preds <- length(pred.groups)

  ## determine the cutoffs for the prediction groups
  mu <- mean(all.data[, 'pred'])
  sigma <- sd(all.data[, 'pred'])
  k <- (n.preds - 1)/2
  for (p in 1:k) {
    pred.breaks[p + 1] <- mu - (k - p + 1)*sigma
    pred.breaks[n.preds + 1 - p] <- mu + (k - p + 1)*sigma
  }

  ## allocate two data frames: one for counting the total samples in
  ## each age group and prediction group; the other for those with no
  ## euploid
  counts.total <- data.frame(matrix(ncol=n.ages, nrow=n.preds))
  colnames(counts.total) <- age.groups
  rownames(counts.total) <- pred.groups
  counts.noeup <- data.frame(matrix(ncol=n.ages, nrow=(n.preds + 1)))
  colnames(counts.noeup) <- age.groups
  rownames(counts.noeup) <- c('age only', pred.groups)

  ## iterate over the age groups
  for (a in 1:n.ages) {
    ag <- age.groups[a]

    ## select by age
    by.age <- subset(all.data, Female.Age >= age.breaks[a] &
                               Female.Age < age.breaks[a + 1])

    ## tabulate for the female age only
    counts.noeup['age only', ag] <- nrow(subset(by.age,
                                                Euploid.Rate == 0))/nrow(by.age)

    ## iterate over the prediction categories
    for (p in 1:n.preds) {
      pg <- pred.groups[p]

      ## select by prediction criteria and tabulate
      x <- subset(by.age, pred >= pred.breaks[p] & pred < pred.breaks[p + 1])
      counts.total[pg, ag] <- nrow(x)
      counts.noeup[pg, ag] <- nrow(subset(x, Euploid.Rate == 0))/nrow(x)
    }
  }

  ## select the colors for the bar chart for the no euploid rate
  colours <- c('blue', 'red', 'green', 'purple', 'lightblue', 'orange')

  ## reset the layout of plots
  par(mfcol=c(2,1))
  barplot(as.matrix(counts.noeup), main='What does the current model add?',
          ylab='no-euploid rate', ylim=c(0, 1), beside=TRUE, col=colours)
  legend('topleft', rownames(counts.noeup), bty='n', fill=colours)

  ## Producing bar chart for the total counts
  colours <- colours[2:(n.preds + 1)]
  barplot(as.matrix(counts.total),
          ylab='N samples', beside=TRUE, col=colours)
  legend('topleft', rownames(counts.total), bty='n', fill=colours)
  garbage <- dev.off()

  ## write the tabulated counts to the terminal
  cat(sprintf('\n[counts of no euploid]\n'))
  print(format(counts.noeup, digits=4))
  cat(sprintf('\n[counts total per group]\n'))
  print(format(counts.total, digits=4))
}


main <- function() {
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) < 5) {
    stop('euploid_analysis.R <features> <outcome-and-female-age> ',
         '<methylation-data> <outcome-label> ',
         '<output-prefix> [cv.fold=10]',
         call.=FALSE)
  }

  features.filename <- args[1]
  metadata.filename <- args[2]
  meth.filename <- args[3]
  outcome.label <- args[4]
  output.basename <- args[5]

  cross.val.fold <- 10
  if (length(args) == 6) {
    cross.val.fold <- as.integer(args[6])
  }

  euploid.analysis(features.filename,
                   metadata.filename,
                   meth.filename,
                   outcome.label,
                   output.basename,
                   cross.val.fold)
}

main()
