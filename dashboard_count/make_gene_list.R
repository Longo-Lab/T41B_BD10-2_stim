#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: make_gene_list.R
##
## Version: 0.0.1
##
## Purpose of script: make list of genes for counts dashboard
##
## Author: Robert B Butler III
##
## Date Created: 2023-08-29
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   makes just a table of Ensembl mouse gene IDs and their GeneSymbol for
##   Patricia's app
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(biomaRt)

# fetch from biomaRt
mouse = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version="110")
genes = getBM(mart=mouse, attributes=c("ensembl_gene_id", 
                                       "external_gene_name",
                                       "chromosome_name"))

# use only canonical chromosomes
genes = subset(genes, chromosome_name %in% c(1:19, "X", "Y", "MT"))

# rename
colnames(genes) = c("gene_id", "GeneSymbol", "chr")
rownames(genes) = genes$gene_id

write.csv(genes[, "GeneSymbol", drop=FALSE], file="gene_list.csv")
