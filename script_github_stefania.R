#Import libraries
library(gdata)
library(Rtsne)
library(dendextend)
library(circlize)

#Set working directory
setwd("Documents/Projects/zebrafish/")

#Load input data
x <- read.xls("data.xlsx", sheet = 1, header = T)
rownames(x) = x$cell
x = x[, 2:17] #full data table
x_subset = x[, 1:13] #partial data table

#Create metafile
metafile <- matrix(nrow = 64, ncol=2)
row.names(metafile) = rownames(x)
metafile[ , 1] = row.names(metafile)
metafile = as.data.frame(metafile)
colnames(metafile) = c("Cell", "Type")
metafile$Type <- ifelse(grepl("Strong adapting", metafile$Cell, ignore.case = T), "Type I",
                        ifelse(grepl("adapting", metafile$Cell, ignore.case = T), "Type II",
                               ifelse(grepl("not", metafile$Cell, ignore.case = T), "Type III", "Type IV")))
metafile[37:56, 2] = c("Type III")
write.table(metafile, file="metafile_200423.txt", sep="\t", quote = F, row.names = T, col.names = NA)

#Create color code
col_x <- metafile$Type
l2cols <- c("blue", "red", "green", "magenta")[as.integer(factor(col_x, levels = c("Type I", "Type II", "Type III", "Type IV")))]

#Dimensionaly reduction on partial data table
set.seed(123)
tsne_norm <- Rtsne(x_subset, pca_scale=T, perplexity = 5, max_iter = 2000)
saveRDS(tsne_norm, file="scale_included_perp5.RDS")
tsne_norm <- readRDS("scale_included_perp5.RDS")
write.table(tsne_norm$Y, file="tSNE_coordinates.txt", quote=F, sep="\t")
pdf(file="tSNE_partial_data_5perplexity_200424.pdf")
plot(tsne_norm$Y, pch = 19, col=l2cols, main="tSNE colored by Type - partial data", xlab = "tSNE1", ylab="tSNE2")
legend("topright", legend=c("Type I", "Type II", "Type III", "Type IV"), bty = 'n', col=c("blue", "red", "green", "magenta"), pch=c(19))
dev.off()

#Herarchical clustering
x_dend <- as.data.frame(scale(normalize_input(as.matrix(x_subset))))
dist_mat <- dist(x_dend, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')

dend <- as.dendrogram(hclust_avg)
dend <- dend %>% 
  color_branches(k=6) 

pdf(file="circle_200424.pdf")
circlize_dendrogram(dend)
dev.off()
