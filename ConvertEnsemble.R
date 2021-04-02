library(readr)
library(dplyr)


gcm <- read.csv('/stor/work/Brock/ag69378/2020TagSeq/results/FusionHybrids/GeneCountMatrix.csv')

#genesymbols <- read_tsv('/stor/work/Brock/ag69378/2020TagSeq/genename_convert/ensgene_to_genesymbol.table', col_names = FALSE)
#colnames(genesymbols) <- c('ensgene', 'gene', 'HGNC')

##################################################################################
# prepare data frame for converting from ensemble to gene symbol
##################################################################################

# pull in data which relates ensemble transcript id to refseq (RNA) and refpro (protein) annotation
              ref <- read_tsv('/stor/work/Brock/ag69378/2020TagSeq/genename_convert/refseq.table', col_names = FALSE)
    colnames(ref) <- c('enstrans','refseq','refpro')

# pull in data which has relationships between gene symbol, ensemble gene, ncbi_id, and refeq (RNA) annotation    
           symbol <- read_tsv('/stor/work/Brock/ag69378/2020TagSeq/genename_convert/genesymbols_refseq.txt', , col_names = TRUE)
 colnames(symbol) <- c('gene_symbol','ncbi_id','refseq','ensgene')

# extract name before the decimal point for refseq annotation and ensemble transcript id
ref$refseq    <- sub('\\.[0-9]*$', '', ref$refseq)
ref$enstrans  <- sub('\\.[0-9]*$', '', ref$enstrans )

# create data frame which relates ensemble transcript id to ensemble gene id and gene_symbol
ens_to_gene <- left_join(symbol, ref)


##################################################################################
# Convert ensemble names in gene count matrix (gcm) to gene symbols
##################################################################################
# DESeq2 requires that all rows (all genes) are uniquely identified. Because our mapping was done to the transcriptome, there 
# are multiple transcripts (isoforms) which map to the same gene which complicates uniqueness of gene symbols in analysis
# to overcome this, we will identify gene symbols which have multiple occurences and for those cases, we will tack on a unique identifier
# (the row index number) such that each isoform of a gene will be considered as different in the analysis

# create new column name of ensemble gene ids before the decimal
gcm$ensgene <- sub('\\.[0-9]*$', '', gcm$Gene)

# join gene count matrix to ensemble conversion matrix
gcm <- left_join(gcm, ens_to_gene)

# remove all rows for which gene names do not exist
## this is fine for early analysis, but may want to rethink this step later or rerun DESeq2 with ensemble ids to see if any
## transcripts with unassigned gene symbols come up
gcm <- gcm[complete.cases(gcm[,'gene_symbol']),]

# count number of occurrences of each gene symbol in the joined matrix
occurrences <- table(unlist(gcm$gene_symbol))

# if number of occurrences of the gene is 1 or less, then assign the new gene name to just be the gene symbol
# if the number of occurrences is multiple, then make the gene symbol unique by tacking on the row index number (column 'X')
# else assign NA
gcm  <- gcm %>% mutate(num_trans = occurrences[gene_symbol]) %>% 
  mutate(new_gene = case_when(num_trans <= 1 ~ gene_symbol, num_trans > 1 ~ paste(gene_symbol, "_", X, sep=""), TRUE ~ NA_character_ ))


# pull out the info you want for DESeq2 -- in this case, we pull the new gene names 'new_gene' and the counts for all of our samples
remove <- c("ensgene", "gene_symbol", "ncbi_id", "refseq", "enstrans", "refpro", "num_trans", "Gene", "X")
gcm_des <- gcm %>% select(-remove) %>% distinct()



#write_csv(gcm_des, 'test_gcm_deseq_table.csv')









