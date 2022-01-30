library(DESeq2)
library(RColorBrewer)

# read in one example file
test <- read.table(file_1_count.txt)
# read in all count files to make a dataframe
file_list <- list.files(count_folder, "*.txt")
all_lib <- data.frame(matrix(nrow=nrow(test)))
row.names(all_lib) <- test[,1]
for (i in c(1:length(file_list))) {
  dataset <- read.table(file_list[i], header=FALSE, col.names = c("geneID", file_list[i]))
  all_lib[,i] <- dataset[,2]
  colnames(all_lib)[i] <- colnames(dataset)[2]
}

# setup a filter to get genes having at least 5 reads in at least 2 samples
filter <- apply(all_lib, 1, function(x) length(x[x>5])>=2)
filtered <- all_lib[filter,]

# setup the experiments factor:
x <- as.factor(c(rep("ctr_p20", length(grep(colnames(filtered), pattern = "^ctr"))),
                 rep("dcKO_p20", length(grep(colnames(filtered), pattern = "^dcKO"))),
                 rep("dis3cKO_p20", length(grep(colnames(filtered), pattern = "^dis3cKO"))),
                 rep("exo10cKO_adult", length(grep(colnames(filtered), pattern = "^Exo10_GV_cKO")))
                 ))
# separate genes and ERCC
filtered_gene <- filtered[grep(row.names(filtered), pattern = "^ENSMU"),]
filtered_ERCC <- filtered[grep(row.names(filtered), pattern = "^ERCC"),]

coldata <- data.frame(row.names=colnames(filtered_gene), x)
coldataERCC <- data.frame(row.names=colnames(filtered_ERCC), x)
ddsERCC <- DESeqDataSetFromMatrix(countData=filtered_ERCC, colData=coldataERCC, design=~x)
dds <- DESeqDataSetFromMatrix(countData=filtered_gene, colData=coldata, design=~x)
sizeFactor <- estimateSizeFactors(ddsERCC)
sizeFactors(dds) <- sizeFactor
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

# Regularized log transformation
rld <- rlogTransformation(dds)
ncol(assay(rld))
hist(assay(rld))

# PCA visualization and in combination of K-means clustering
res.pca <- pca.data
library(factoextra)
library("scales")
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation

# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

#plot PCA, label all groups
dotSize1 <- 3
alp <- 0.8
p<-ggplot(scores,aes(x=PC1,y=PC2,label=row.names(transpose_df), color=x)) +
  geom_point(size=dotSize1, alpha = alp) +
  coord_fixed(ratio=1)+
  labs(title=paste0('PCA_p20_adult'),
       x=paste0("PC1: ", percent(eig.val$variance.percent[1]/100)),
       y = paste0("PC2: ", percent(eig.val$variance.percent[2]/100)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
p

# K-Means Cluster Analysis
set.seed(2)
fit <- kmeans(res.ind$coord, 3) # 3 cluster solution
# get cluster means
aggregate(res.ind$coord,by=list(fit$cluster),FUN=mean)
# append cluster assignment
cluster.res <- data.frame(res.ind$coord, fit$cluster)
cluster.res
colnames(cluster.res)
fit$cluster

p<-ggplot(scores,aes(x=PC1,y=PC2,label=row.names(transpose_df), color=fit$cluster)) +
  geom_point(size=dotSize1, alpha = alp) +
  coord_fixed(ratio=1)+
  labs(title=paste0('PCA_p20_adult_kmeans3'),
       x=paste0("PC1: ", percent(eig.val$variance.percent[1]/100)),
       y = paste0("PC2: ", percent(eig.val$variance.percent[2]/100)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
p

# differential analysis of MA plots and results
con1 <- c("dcKO_p20","dis3cKO_p20", "exo10cKO_p20")
con2 <- c("ctr_p20", "ctr_p20", "ctr_p20")

for (i in 1:length(con1)) {
  con1[i]
  con2[i]
  res <- results(dds, contrast=c("x", con1[i], con2[i]), alpha=0.1)
  resFilename <- paste0("compare_all_", con1[i], "_vs_", con2[i], ".csv")
  write.csv(res, resFilename)

  # below do MA plot
  library(ggplot2)
useForPlot <- data.frame(log10(res$baseMean),res$log2FoldChange,res$padj)
colnames(useForPlot) <- c("log10baseMean", "log2FC", "padj")
useForPlot
up <- useForPlot[which((useForPlot$padj < 0.01) & (useForPlot$log2FC>0)), ]
down <- useForPlot[which((useForPlot$padj < 0.01) & (useForPlot$log2FC<0)), ]
ncol(down)
down
yuplim = 8.0
ydownlim = -8.0

upout <-useForPlot[which(useForPlot$log2FC>= yuplim), ]
if (length(upout$log2FC) >0) {
  upout$log2FC <- yuplim
}
downout <-useForPlot[which(useForPlot$log2FC<= ydownlim), ]
if (length(downout$log2FC)>0) {
  downout$log2FC <- ydownlim
}

alp=0.1
dotSize <- 2
ggplot()+
  geom_point(data = useForPlot, aes(x = useForPlot$log10baseMean, y = useForPlot$log2FC), size = dotSize, alpha=alp) +
  geom_point(data = up, aes(x = up$log10baseMean, y = up$log2FC), size = dotSize, color="red", alpha=alp) +
  geom_point(data = down, aes(x = down$log10baseMean, y = down$log2FC), size = dotSize, color="blue", alpha=alp) +
  labs(title=paste0(con1[i], "_vs_", con2[i]), x="log10baseMean", y = "log2FC") + ylim(-8,8) +
  geom_point(data = upout, aes(x = upout$log10baseMean, y = upout$log2FC), size = dotSize, color="red", shape=5) +
  geom_point(data= downout, aes(x = downout$log10baseMean, y = downout$log2FC), size = dotSize, color="blue", shape=5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=0, linetype="dashed", color = "gray", size=1)

ggsave(paste0("compare_all_", con1[i], "_vs_", con2[i], ".tiff"))

}
