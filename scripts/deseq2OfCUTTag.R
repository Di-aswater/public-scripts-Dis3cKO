###########################################################
#below for differential analysis with DESeq2 of the polIIs2+s5 as example:
# read_in the SF file
dat <- read.table("polII_SF.txt", sep="\t",header=TRUE,fill=TRUE,stringsAsFactors = FALSE, quote="",check.names=F)
dat
SF <- dat$SF
SF

# use the raw count:
countData <- read.table("polII_repeat_rawCounts.txt",header = TRUE,row.names = 1)
countData
# do a filter of countData: at least 5 lib have more than 10
filter <- apply(countData, 1, function(x) length(x[x>10])>=5)
filter
filtered <- countData[filter,]
length(filtered)
head(filtered)
colnames(filtered)
x <- as.factor(c(rep('pols2s5_adultGV_ctr',3), rep('pols2s5_adultGV_dis3cKO',3),rep('pols2s5_p20_ctr',3)))
library(DESeq2)

countData <- filtered
dds <- DESeqDataSetFromMatrix(countData, DataFrame(x), ~ x)
dds <- estimateSizeFactors(dds)
coldata<- colData(dds)
coldata$sizeFactor <- coldata$sizeFactor *SF
colData(dds) <- coldata
dds <- estimateDispersions(dds) # extremely slow. better if done filter
dds <- nbinomWaldTest(dds) # extremely slow

# below are from previous DESeq2 file:
rld <- rlogTransformation(dds) # extremely slow
ncol(assay(rld))
hist(assay(rld))

library(RColorBrewer)
mycols <- brewer.pal(8, "Dark2")[1:length(unique(x))]

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
tiff("../deseq_polII/sampleDist2.tiff", w=1000, h=1000, pointsize=15)
heatmap.2(as.matrix(sampleDists), key=T, density.info="none", trace="none",
          col=colorpanel(100, "blue", "yellow"),
          ColSideColors=mycols[x], RowSideColors=mycols[x],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

#principal component analysis

#this puts samples as rows, genes as columns 
transpose <- t(assay(rld)) 
transpose_df <- as.data.frame(transpose)
#this is the function that does the PCA
pca.data <- prcomp(transpose_df)
scores = as.data.frame(pca.data$x) 
summary(pca.data)

#2d PCA plot:
pdf("../deseq_polII/pcantop2000_2.pdf")
DESeq2::plotPCA(rld, intgroup="x",ntop=2000)
dev.off()
x
#differential expressed gene lists: not finished
con1 <- c("pols2s5_adultGV_ctr", "pols2s5_adultGV_dis3cKO")
con2 <- c("pols2s5_p20_ctr", "pols2s5_adultGV_ctr")
length(con1)
for (i in 1:length(con1)) {
  con1[i]
  con2[i]
  res <- results(dds, contrast=c("x", con1[i], con2[i]), alpha=0.1)
  resFilename <- paste0("../deseq_polII/",con1[i], "_vs_", con2[i], ".csv")
  write.csv(res, resFilename)
  
  #below do MA plot
  library(ggplot2)
  useForPlot <- data.frame(log10(res$baseMean),res$log2FoldChange,res$padj)
  colnames(useForPlot) <- c("log10baseMean", "log2FC", "padj")
  useForPlot
  up <- useForPlot[which((useForPlot$padj < 0.01) & (useForPlot$log2FC>0)), ]
  down <- useForPlot[which((useForPlot$padj < 0.01) & (useForPlot$log2FC<0)), ]
  ncol(down)
  down
  yuplim = 10.0
  ydownlim = -10.0
  
  upout <-useForPlot[which(useForPlot$log2FC>= yuplim), ]
  if (length(upout$log2FC) >0) {
    upout$log2FC <- yuplim
  } 
  downout <-useForPlot[which(useForPlot$log2FC<= ydownlim), ]
  if (length(downout$log2FC)>0) {
    downout$log2FC <- ydownlim
  }
  dotSize <- 2
  ggplot()+
    geom_point(data = useForPlot, aes(x = useForPlot$log10baseMean, y = useForPlot$log2FC), size = dotSize, alpha=0.1) +
    geom_point(data = up, aes(x = up$log10baseMean, y = up$log2FC), size = dotSize, color="red", alpha=0.1) + 
    geom_point(data = down, aes(x = down$log10baseMean, y = down$log2FC), size = dotSize, color="blue", alpha=0.1) +
    labs(title=paste0(con1[i], "_vs_", con2[i]), x="log10baseMean", y = "log2FC") + ylim(-5,5) + geom_point(data = upout, aes(x = upout$log10baseMean, y = upout$log2FC), size = dotSize, color="red", shape=5) + geom_point(data= downout, aes(x = downout$log10baseMean, y = downout$log2FC), size = dotSize, color="blue", shape=5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) + geom_hline(yintercept=0, linetype="dashed", color = "gray", size=1) 
  
  ggsave(paste0("../deseq_polII/",con1[i], "_vs_", con2[i], ".tiff"))
  
}


