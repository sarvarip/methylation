df <- readRDS('methylation.rds')
rownames(df) <- gsub("X", "", rownames(df))
numrows = length(rownames(df))

fi <- file("colnames",open="r")
featnames <- readLines(fi)

fil <- file("colidx", open="r")
featidx <- readLines(fil)
featidx <- as.integer(featidx)

f <- file("coord_lists_per_gene",open="r")
lines <- readLines(f)
num <- length(lines)

reduced <- data.frame(matrix(ncol=num, nrow=numrows))
rownames(reduced) <- rownames(df)
chosen.feat <- featnames[featidx]
colnames(reduced) <- chosen.feat
for (i in 1:num) {
    print(i)
    feat <- unlist(strsplit(lines[i], split='\t'))
    reduced[,i] <- tryCatch(apply(df[,feat, drop='False'], 1, mean), 
             error = function(e){
                 featlist <- strsplit(lines[i], split='\t')
                 featlist <- featlist[[1]]
                 newlist <- character(0)
                 for (j in 1:length(featlist)){
                     if (featlist[j] %in% colnames(df)) {
                         newlist <- c(newlist, featlist[j])
                     }
                 }
                 feat <- unlist(newlist)
                 col <- apply(df[,feat, drop='False'], 1, mean)
                 return(col)
             })
}
saveRDS(reduced, file = "MethylationByGene.rds")
