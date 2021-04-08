library(DESeq2)
library(tidyverse)
library(matrixStats)
library(pheatmap)
source("bio321g_rnaseq_utils.R")

rnaCounts 

### 1.  DESeq


## rnaCounts: data.frame containing RNA-seq counts for each gene in each sample
## (genes are in rows of data.frame, samples in columns):
rnaCounts = read.table("CONVERTED_COAD_GCM.csv",
                       sep=",", header=TRUE, row.names=1, check.names=FALSE)

## sampleAnnotation: data.frame with one row per sample; columns say what
## group (=combination of genotype+time), genotype, and time
## describe each sample:
sampleAnnotation = read.table("COAD_sample_annotation.csv",
                              sep=",", header=TRUE, row.names=1, check.names=FALSE)
## first initialize DESeq object:
dds <- DESeqDataSetFromMatrix(
  ## countData argument should be numeric matrix (or data.frame with
  ## only numeric columns which can be coerced to numeric matrix):
  countData = rnaCounts,
  ## colData argument should be data.frame with grouping information
  ## (and any other non-gene-expression covariates that might be
  ## useful in modeling the counts) whose rows correspond to the
  ## columns (hence name colData) of countData:
  colData = sampleAnnotation,
  ## note tilde preceding the name group of the relevant column in
  ## the data.frame sampleAnnotation provided as the colData argument
  ## in the design argument---design argument must be formula object,
  ## similar to the right-hand side of the formulas used for linear
  ## modeling with lm:
  design = ~ tissue)

## now run DESeq analysis pipeline on object dds in order to perform a hypothesis test specifically for the significance of the interaction term 
dds <-  DESeq(dds)
##Finally, use the results function from the DESeq2 package to extract a table of results for each gene from the DESeqDataSet object.
## sort dds results data.frame by p-value:
ddsresults <- results(dds)
ddsresults <- ddsresults[order(ddsresults$pvalue), ]

ddsresults

## genes from the full rnacounts data are considered to have significant time:geneotype interaction terms to keep the False Discovery Rate at 10%. 
False_Discovered_Genes <- sum(ddsresults$padj < 0.1, na.rm=TRUE ) 
print(False_Discovery_Genes)

## of our our list of genes with evidence of significant interactions would be considered False Positives. 
False_Positive_Genes <- False_Discovery_Genes*.1
print(False_Positive_Genes)


### 2. Extract the Normalized Counts from the DESeqDataSet Object 

## use counts function with normalized arg set to TRUE to extract
## "normalized counts" (which are not integers) from dds
## (note: counts *function* from DESeq2 package is different object
## in R than counts *data.frame* we loaded from tsv file;
## R can tell based on context which one we mean):
normed <- counts(dds, normalized=TRUE)

write.table(data.frame(normed),"COAD_normalized_counts.csv", sep = ",", row.names = FALSE, quote = FALSE)

## save the normalized counts matrix to a tsv file:
lgNorm <- log2(normed + 1)

### 5. Normalizing Genes

lgNorm <- log2(normed + 1)
lgNorm <- data.frame(lgNorm)
## make sure lgNorm is a data.frame
## make the rownames a column name
lgNorm <- rownames_to_column(lgNorm, var = "gene")
## make a new data.frame containing only the rows of normalized data that correspond to genes with a gene ontology id of GO:0006090
lgGo <- lgNorm %>% filter(lgNorm[,1] %in% geneset$geneID)
## get rid of the X's in front of the column names... I know there is a better way to do this, revisit 
names(lgGo) <- c("gene","14BENDDAY2", "14BENDDAY4", "14BEXDARK2", "14BEXDARK3", "14BEXDARK4", "4GENDDAY2", "4GENDDAY3", "4GENDDAY4", "4GEXDARK2", "4GEXDARK3" , "COLENDDAY3" , "COLENDDAY5" , "COLEXDARK2" , "COLEXDARK3" , "COLEXDARK4")
## remove gene column- make it row names
lgGo <- column_to_rownames(lgGo, var = "gene")

### 7. Gene Expression Heatmap of Normalized Gene Expression Data

heatData = lgGo - rowMeans(lgGo)

heatData[heatData > 2] = 2; heatData[heatData < -2] = -2

pheatmap(
  heatData,
  color = heatPalette,
  clustering_method = "average",
  labels_row=geneNamesAndDescriptions[rownames(heatData), "symbol"]
)













