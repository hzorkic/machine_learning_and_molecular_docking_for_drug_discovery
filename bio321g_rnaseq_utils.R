## -----------------------------------------------------------------------------
## load data set (make sure files are in R working directory!
## -----------------------------------------------------------------------------
## rnaCounts: data.frame containing RNA-seq counts for each gene in each sample
## (genes are in rows of data.frame, samples in columns):
rnaCounts = read.table("rna_counts.tsv.gz",
                       sep="\t", header=TRUE, row.names=1, check.names=FALSE)

## sampleAnnotation: data.frame with one row per sample; columns say what
## group (=combination of genotype+time), genotype, and time
## describe each sample:
sampleAnnotation = read.table("rna_sample_annotation.tsv",
                              sep="\t", header=TRUE, row.names=1, check.names=FALSE)

## geneNamesAndDescriptions: data.frame with rownames corresponding to gene
## ids and three columns:
## (1) gene :: gene id (same as rownames,
## (2) symbol :: gene name/symbol
## (3) description :: gene description
geneNamesAndDescriptions = read.table("arabidopsis_thaliana_gene_names.tsv.gz",
                       sep="\t", row.names=1, header=TRUE,
                       quote="", comment.char="")
geneNamesAndDescriptions$gene = rownames(geneNamesAndDescriptions)
geneNamesAndDescriptions =
        geneNamesAndDescriptions[ , c("gene", "symbol", "description")]

## goAssociations: data.frame indicating what genes are associated with what
## gene sets, with four columns:
## (1) gene_ontology_primary_id :: gene set identifier
## (2) gene_ontology_name :: gene set name
## (3) gene_ontology_all_ids :: indicates what other gene ontology groups
##                              have been merged into this gene set; you
##                              don't need to worry about this column!
## (4) gene :: the gene ids associated with the gene set identified by
##             gene_ontology_primary_id column
goAssociations = read.table("gene_sets.tsv.gz",
                            sep="\t", row.names=NULL, header=TRUE,
                            quote="", comment.char="")

## -----------------------------------------------------------------------------
## define colors for plotting:
## -----------------------------------------------------------------------------
heatPalette = colorRampPalette(c("dodgerblue", "lightskyblue", "white",
                                 "lightgoldenrod", "orangered"))(100)

groupColors = c(
    "14BENDDAY" = "orangered",
    "14BEXDARK" = "darkred",
    "4GENDDAY" = "lightseagreen",
    "4GEXDARK" = "dodgerblue3",
    "COLENDDAY" = "lightslategray",
    "COLEXDARK" = "darkslategray"
)

## -----------------------------------------------------------------------------
stripchart321g = function(data,
                          sampleAnnotation,
                          geneNames = geneNamesAndDescriptions,
                          colorValues = groupColors) {
    ## requires dplyr, tidyr, and ggplot2 libraries to be available
    ## (all can be installed via install.packages)
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    ggdata = data %>%
        as.data.frame %>%
        mutate(gene=rownames(data)) %>%
        pivot_longer(-gene, names_to="sample", values_to="expression")
    ggdata[ , c("group", "time", "genotype")] = sampleAnnotation[
        ggdata$sample,
        c("group", "time", "genotype")
    ]
    if (!missing(geneNames)) {
        geneNames = unique(geneNames[ , c("gene", "symbol")])
        geneNames = structure(geneNames$symbol, names=geneNames$gene)
        ggdata$name = paste0(geneNames[ggdata$gene], " (", ggdata$gene, ")")
        ggdata$name = factor(ggdata$name, levels=unique(ggdata$name))        
    } else {
        ggdata$name = factor(ggdata$gene, levels=unique(ggdata$gene))
    }
    gg = ggdata %>% ggplot(aes(x = genotype,
                               y = expression,
                               color = group,
                               shape = time))
    gg = gg + facet_wrap(~ name, scales="free_y")
    gg = gg + geom_point()
    if (length(colorValues) > 0) {
        gg = gg + scale_color_manual(values=colorValues)
    }
    gg = gg + scale_shape_manual(values=c(2, 6))
    invisible(gg)
}
