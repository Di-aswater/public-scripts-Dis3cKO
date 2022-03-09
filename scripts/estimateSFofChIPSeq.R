#####################################################
# installation all done.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicAlignments")

# Install this package from GitHub
install.packages("devtools")
library(devtools)
install_github("stjude/ChIPseqSpikeInFree")
packageVersion('ChIPseqSpikeInFree')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")
#####################################################
# below generate SF:
library("ChIPseqSpikeInFree")
library(BiocParallel)
# read in the metadata sheet:
# metadata sheet is in the metadata folder
samples <- read.csv('samplesPolIIK27.csv')
names(samples)
samples$SampleID
# only pick PolII data here. 
df_polII <- samples[c(1:9),]
df_polII
# difine metaFile to have ID	GROUP	ANTIBODY
metaFile <- df_polII[,c(1,3,2)]
colnames(metaFile) <- c('ID','ANTIBODY','GROUP')
bamfiles <- list.files(path="diffBind/reads/",pattern="^pol2*.*bam$")
metaFile$ID <- bamfiles

# must be saved as a file
write.table(metaFile,'metaFile_polII.txt',row.names = FALSE, sep = "\t", quote=FALSE)
bams <- paste0('reads/',bamfiles)
bams
metaFile
#Run ChIPseqSpikeInFree pipeline 
?ChIPseqSpikeInFree
# setwd the output folder
setwd("../scale_repeat")
metaFile <- "/Users/Di/Documents/J_lab/20211108_cutTag/scale_DEseq/scale_repeat/metaFile_polII.txt"
metaFile
# generate all output files, including SF values.
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "mm10", metaFile = metaFile, prefix = "polII")

# obtained values of SF: 
#pols2s5.pols2s5_adultGV_ctr ,	ave.SF = 2.19
#pols2s5.pols2s5_adultGV_dis3cKO ,	ave.SF = 1.503
#pols2s5.pols2s5_p20_ctr ,	ave.SF = 1.943

# use the same code, take samples[c(10:21),] can generate the SF of H3K27me3 samples.
# obtained values of SF: 
#h3k27me3.h3k27me3_adultGV_ctr ,	ave.SF = 7.43
#h3k27me3.h3k27me3_adultGV_dis3cKO ,	ave.SF = 1.44
#h3k27me3.h3k27me3_p20_ctr ,	ave.SF = 7.33
#h3k27me3.h3k27me3_p20_dcKO ,	ave.SF = 1.27

