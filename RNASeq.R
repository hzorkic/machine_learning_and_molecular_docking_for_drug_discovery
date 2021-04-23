# In this report, we are primarily interested in understanding how the 
# differences in gene expression between healthy and canerous tissue samples
# from The Cancer Genome Altas (TCGA)


##################
# load libraries #
##################
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(matrixStats)
library(pheatmap)


################################################
# Perform DESeq analysis on our raw RNA counts #
################################################
# A DESeqDataSet Object was constructed on the gene expression data using the 
# DESeq package to perform a hypothesis test for the significance of tissue type
# in the model

# rnaCounts: data.frame containing RNA-seq counts for each gene in each sample
# (genes are in rows of data.frame, samples in columns):
rnaCounts = read.table("CONVERTED_COAD_GCM.csv",
                       sep=",", 
                       header=TRUE, 
                       row.names=1, 
                       check.names=FALSE)

# sampleAnnotation: data.frame with one row per sample; columns say what
# group (=tissue cancerous or not) describe each sample:
sampleAnnotation = read.table("COAD_sample_annotation.csv",
                              sep=",", 
                              header=TRUE, 
                              row.names=1, 
                              check.names=FALSE)

# first initialize DESeq object:
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

# now run DESeq analysis pipeline on object dds in order to perform a hypothesis 
# test specifically for the significance of the interaction term 
dds <-  DESeq(dds)

#Finally, use the results function from the DESeq2 package to extract a table of 
# results for each gene from the DESeqDataSet object.
# sort dds results data.frame by p-value:
ddsresults <- results(dds)
ddsresults <- ddsresults[order(ddsresults$pvalue), ]
ddsresults


##############################################################
# Extract the Normalized Counts from the DESeqDataSet Object #
##############################################################

# use counts function with normalized arg set to TRUE to extract
# "normalized counts" (which are not integers) from dds
# (note: counts *function* from DESeq2 package is different object
# in R than counts *data.frame* we loaded from csv file;
# R can tell based on context which one we mean):
normed <- counts(dds, normalized=TRUE)
# log transform in order to  makes it easy to see proportional changes in 
# expression levels in differential expression. For example we would observe 
# a tumor gene to have positive expression in tumor tissue, but the same gene 
# would have negative proportional expression in a healthy tissue
lgNorm <- log2(normed + 1)

# save the normalized counts matrix to a tsv file:
write.table(data.frame(normed),
            "COAD_normalized_counts.csv", 
            sep = ",", row.names = FALSE, 
            quote = FALSE)


##########################################################
# Generate a PCA Plot of Normalized Gene Expression Data #
##########################################################
# make sure lgNorm is a data.frame
lgNorm <- data.frame(lgNorm)
# make the rownames a column name
lgNorm <- rownames_to_column(lgNorm, var = "gene")
lgGo <- column_to_rownames(lgNorm, var = "gene")

sampleAnnotation = data.frame(
  ## set row.names of sampleAnnotation to match col.names of normed:
  row.names = colnames(normed))

# save
write.csv(lgGo, file="log2transformed_and_normalized_gene_expression_data.csv")


pca = prcomp(t(lgGo))
## perform PCA on Normalized data
pcaFit = rowMeans(lgGo) + t( pca$x %*% t(pca$rotation) )
## have to transpose pca$x %*% t(pca$rotation) above b/c lgNorm is t(z)
## set up data.frame pcaData for ggplot...
pcaData = data.frame(pca$x[ , 1:2])
## add in sample annotation info
pcaData$group = sampleAnnotation[rownames(pcaData), "group"]
## and sample names
pcaData$sample = rownames(pcaData)

pcaData <- pcaData %>% data.frame()
pcaData$group <- ifelse(grepl("Healthy", pcaData$sample, ignore.case = T), "Healthy", 
                  ifelse(grepl("Tumor", pcaData$sample, ignore.case = T), "Tumor", "Other"))

pcaData


ggplot(pcaData, aes(x=PC1, y=PC2, label=sample, color = group)) + 
  geom_point(size=4, shape = 18, alpha = 0.75) + 
  scale_color_manual(values = c("Healthy" = "orange", "Tumor" = "maroon")) +
  ggtitle("Principal Components Analysis (PCA) of TCGA COAD Tissue RnaSeq Data") + 
  theme_light()


##################################################
# Genes with Significant Differential Expression #
##################################################

# order for top gene
top <- ddsresults[order(ddsresults$padj),]

# select only those genes with a padj greater than X (0.05)
top <- data.frame(top) %>% filter(padj < 0.05)
nrow(top)

# select the top X genes (3000)
top <- top[1:3000, ]
# remove any rows that have NAs 
top <- na.omit(top)

# save
write.csv(top, file="top_genes.csv")


##############################################
# Heatmap of Normalized Gene Expression Data #
##############################################

# IMPORTANT::: I would actually just input the 
# log2transformed_and_normalized_gene_expression_data.csv 
# Gene Expression Data into Morpheus: 
# https://software.broadinstitute.org/morpheus/
















