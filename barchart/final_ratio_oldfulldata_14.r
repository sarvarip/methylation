## Getting and merging data
df <- read.delim('13_methylation', header = TRUE, sep = '\t')
rownames(df) <- df$X
df$X <- NULL
df.meta <- read.csv(file='OldDataFemale_Age_Euploid_Rate.tsv', header=TRUE, sep='\t')
colnames(df.meta) <- c('laboratoryID', 'femaleAge', 'ratio')
df.meta <- df.meta[!(is.na(df.meta$femaleAge) | df.meta$femaleAge==''), ]
df.meta <- df.meta[!(is.na(df.meta$ratio) | df.meta$ratio==''), ]
rownames(df.meta) <- df.meta$laboratoryID
df.meta$laboratoryID <- NULL
rownames(df.meta) <- gsub('X', '', rownames(df.meta))
data <- merge(df, df.meta, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

## Fitting data
fit1 <- lm(ratio ~ ., data=data)
dataX <- data
dataX$ratio <- NULL
pred <- predict(fit1, dataX)
dat <- cbind(data$ratio, pred)
print(dat)
print('No rows')
length(rownames(data))

## Creating data that contains prediction values
finaldf <- merge(dat, df.meta, by=0)
print(finaldf)
print('Number of rows')
length(rownames(finaldf))

## Define the breaks and labels for the age groups
age.groups <- c('under 30', '30-35', '35-40', '40 and over')
age.breaks <- c(-Inf, 30, 35, 40, Inf)
n.ages <- length(age.groups)
## Define the labels for the prediction group
pred.groups <- c('poor', 'borderline', 'normal', 'good', 'ideal')
pred.breaks <- c(-Inf, 0, 0, 0, 0, Inf) # middle 4 must be learned
n.preds <- length(pred.groups)

## Determine the cutoffs for the prediction groups
mu <- mean(finaldf[, 'pred'])
sigma <- sd(finaldf[, 'pred'])
k <- (n.preds - 1)/2
for (p in 1:k) {
  pred.breaks[p + 1] <- mu - (k - p + 1)*sigma
  pred.breaks[n.preds + 1 - p] <- mu + (k - p + 1)*sigma
}

# allocate two data frames: one for counting the total samples in each
# age group and prediction group; the other for those with no euploid
total.counts <- data.frame(matrix(ncol=n.ages, nrow=n.preds))
colnames(total.counts) <- age.groups
rownames(total.counts) <- pred.groups
female.total <- data.frame(matrix(ncol=n.ages, nrow=1))
colnames(female.total) <- age.groups
rownames(female.total) <- 'female'
ratio.counts <- data.frame(matrix(ncol=n.ages, nrow=n.preds))
colnames(ratio.counts) <- age.groups
rownames(ratio.counts) <- pred.groups
female.ratio <- data.frame(matrix(ncol=n.ages, nrow=1))
colnames(female.ratio) <- age.groups
rownames(female.ratio) <- 'female'


## iterate over the age groups
for (a in 1:n.ages) {
    ## iterate over the prediction categories
    x <- subset(finaldf,
                femaleAge >= age.breaks[a] & femaleAge < age.breaks[a + 1])
    female.total[,age.groups[a]] <- nrow(x)
    female.ratio[,age.groups[a]] <- nrow(subset(x, ratio == 0)) / nrow(x)
  for (p in 1:n.preds) {
    ## select by age and prediction criteria
    x <- subset(finaldf,
                femaleAge >= age.breaks[a] & femaleAge < age.breaks[a + 1] &
                pred >= pred.breaks[p] & pred < pred.breaks[p + 1])
    total.counts[pred.groups[p], age.groups[a]] <- nrow(x)
    ratio.counts[pred.groups[p], age.groups[a]] <- nrow(subset(x, ratio == 0)) / nrow(x)
  }
}

all.samples <- rbind(female.total, total.counts)
all.ratio <- rbind(female.ratio, ratio.counts)

##Producing bar chart
colours <- c('blue', 'red', 'green', 'purple', 'yellow', 'orange')
png('TimReproduceFinal.png')
barplot(as.matrix(all.ratio), main='What does current model add?', ylab='No-Euploid Blasts Rate', cex.lab=1.5, cex.main=1.4, ylim=c(0,1), beside=TRUE, col=colours)
legend('topleft', c('Female Only','Poor','Borderline','Normal','Good','Ideal'), cex=1.3, bty='n', fill=colours) 
dev.off()
png('SampleTimReproduceFinal.png')
barplot(as.matrix(all.samples), main='What does current model add?', ylab='Sample size', cex.lab=1.5, cex.main=1.4, beside=TRUE, col=colours)
legend('topleft', c('Female Only','Poor','Borderline','Normal','Good','Ideal'), cex=1.3, bty='n', fill=colours) 
dev.off()


