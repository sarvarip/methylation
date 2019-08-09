df <- readRDS('MethylationByGene.rds')

df_s <- read.csv(file="ratiodat", header=TRUE, sep="\t")
df_s <- df_s[!(is.na(df_s$Female.Age) | df_s$Female.Age==""), ]
df_s <- df_s[!(is.na(df_s$Ratio) | df_s$Ratio==""), ]
rownames(df_s) <- df_s$Laboratory.ID
df_s$Laboratory.ID <- NULL

data <- merge(df, df_s, by=0)
rownames(data) <- data$Row.names
data$Row.names <- NULL

data[,colnames(data)[colSums(is.na(data)) > 0]] <- NULL

data$euploid <- as.factor(as.integer((data$Ratio==0)))
data$Ratio <- NULL

#library(Rtsne)
library(tsne)
library(irlba)

train <- data
train$euploid <- NULL

train.scaled <- scale(train)
df.svd <- irlba::irlba(train.scaled, nv = 100)

pdf('EigenvalsGeneAvg.pdf', width=6.5, height=8, paper='US')
plot(df.svd$d, main="Eigenvals")
dev.off()

# 40 eigenvectors were chosen based on plot above
df.reduced <- as.matrix(train.scaled) %*% df.svd$v[,1:40]
# same as df.reduced <- df.svd$u[,1:40] %*% diag(df.svd$d[1:40])

## for plotting
colors = rainbow(length(unique(data$euploid)))
names(colors) = unique(data$euploid)

## Executing the algorithm on curated data
#red <- Rtsne(df.reduced, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
red <- tsne(df.reduced)

## Plotting
pdf('tsneGeneAvg.pdf', width=6.5, height=8, paper='US')
plot(red, t='n', main="tsne")
text(red, labels=data$euploid, col=colors[data$euploid])
dev.off()
