#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: adatlas_prep.R
##
## Version: 0.0.1
##
## Purpose of script: Make human gene lists for ADatlas.org
## DE genesets
##
## Author: Robert B Butler III
##
## Date Created: 2021-11-14
##
## Copyright (c) 2021
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Currently hardcodes the input files (good target for options funtion)
##   looks at curated biodomains
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(stringr)
library(biomaRt)

# directory locations
rootdir = '/labs/flongo/t41b_BD10-2_stim'
wk_dir = paste0(rootdir, '/drug_adatlas')
setwd(paste0(rootdir, '/drug_adatlas'))
today = format(Sys.Date(), '%Y%m%d')
# treatdir = "/labs/flongo/TREAT-AD/Biodomain_annotations"
nameset = "human_genes"
# DE cutoff by "pvalue" or "padj"
p = "padj"
pval = 0.1

#####################################################
# fetch lists from directories
f_names = c("2", "9", "12")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
replacement = c("APP", "Drug", "APP + Drug")

####################################################################FUNCTIONS
convertMouseHuman = function(dt, colname, keep.mouse=F){
  # for a given column of mouse ensIDs, replace them with human ensembl gene ids 
  # includes a filter for autosomal genes that are not pseudo or TEC
  # dt: data table to replace names
  # colname: name of the column containing gene symbols
  # keep.symbols: boolean, retain the symbols column?
  # returns dt with a GENE column
  require("biomaRt")
  print(paste(length(unique(dt[[colname]])), "genes to look up"))
  mouse = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version="104")
  human = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version="104")
  # match with ensembl_gene_id across both
  x = data.table(getLDS(mart=mouse, attributes=c('ensembl_gene_id'), martL=human,
                        attributesL=c('ensembl_gene_id', 'chromosome_name', 
                                      'gene_biotype'), 
                        filters='ensembl_gene_id', values=unique(dt[[colname]]),
                        bmHeader=F))
  names(x) = make.unique(names(x))
  print(paste(nrow(x), "total genes returned"))
  # filter out pseudogenes and TEC
  keep_categ = grep("pseudo", unique(x$gene_biotype), value=T, invert=T)
  keep_categ = grep("TEC", keep_categ, value=T, invert=T)
  print(paste("keeping categories:", keep_categ))
  x = x[ gene_biotype %in% keep_categ ]
  print(paste(nrow(x), "genes correct biotype"))
  # use only canonical chromosomes
  use_chrs = c(1:22, "X", "Y")
  x = x[chromosome_name %in% use_chrs]
  print(paste(nrow(x), "genes in chr 1:22 X Y"))
  setnames(x, c("ensembl_gene_id", "ensembl_gene_id.1"), c("GENE", "hGENE"))
  # replace gene names and write to file
  dt = merge(x[, .(GENE, hGENE)], dt, all.y=T, by.x="GENE", by.y=colname)
  dt = dt[!is.na(hGENE)]
  print(paste(length(unique(dt[["hGENE"]])), "unique human genes with ensembl ids"))
  if (keep.mouse == FALSE) dt[, GENE := NULL]
  return(unique(dt))
}

####################################################################MAIN
# read in DE files
all(file.exists(target_files))
# reading in tables, combine and melt by GWAS
file_list = lapply(target_files, fread)
lapply(file_list, function(x) setorder(x, Gene_id))
setattr(file_list, 'names', f_names)

# truncate lists
short_list = lapply(file_list, function(x) {
  x[!(is.na(Gene_id)), c("Gene_id", "GeneSymbol", "Gene_type", "log2FoldChange",
                         "baseMean", "pvalue", "padj")]
})
colnames = c("log2FoldChange", "baseMean", "pvalue", "padj")

# join all sets
result = short_list[[1]]
for (i in head(seq_along(short_list), -1)) {
  result = merge(x=result, y=short_list[[i+1]], by=c("Gene_id", "GeneSymbol", 
                                                     "Gene_type"), all=T, 
                 suffixes=c("", paste0(".", f_names[i+1])))
}
setnames(result, colnames, paste0(colnames, ".", f_names[1]))

#######################
## Fisher exact testing for enrichment of up and downregulated genes
# filter out pseudogenes and TEC
keep_categ = grep("pseudo", unique(result$Gene_type), value=T, invert=T)
keep_categ = grep("TEC", keep_categ, value=T, invert=T)
result = result[ Gene_type %in% keep_categ ]

# DE genes
# generate list of gene sets for AD genes
de_list = lapply(f_names, function(i) {
  x = result[get(paste0(p, ".", i)) < pval, .(Gene_id)]
  if (nrow(x) > 0) {
    fwrite(convertMouseHuman(x, "Gene_id"), col.names=F, 
         file=paste(nameset, i, "geneset", sep="."))
  }
})
# names(de_list) = replacement

# DE genes (split up/down)
# make list
f_names_updown = c(paste(replacement, "Up"), paste(replacement, "Down"))
# up
up_list = lapply(f_names, function(i) {
  x = result[get(paste0(p, ".", i)) < pval 
             & get(paste0("log2FoldChange.", i)) > 0, .(Gene_id)]
  if (nrow(x) > 0) {
    fwrite(convertMouseHuman(x, "Gene_id"), col.names=F, 
           file=paste(nameset, "up", i, "geneset", sep="."))
  }
})
# down
down_list = lapply(f_names, function(i) {
  x = result[get(paste0(p, ".", i)) < pval 
             & get(paste0("log2FoldChange.", i)) < 0, .(Gene_id)]
  if (nrow(x) > 0) {
    fwrite(convertMouseHuman(x, "Gene_id"), col.names=F, 
           file=paste(nameset, "down", i, "geneset", sep="."))
  }
})
# # join and name
# updown_list = c(up_list, down_list)
# names(updown_list) = f_names_updown


