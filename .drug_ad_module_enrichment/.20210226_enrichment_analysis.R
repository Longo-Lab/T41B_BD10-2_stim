#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(data.table)
library(stringi)
library(ComplexHeatmap)
library(circlize)
library(GeneOverlap)

# directory locations
rootdir = '/labs/flongo/2020_T41B_BD10-2_Stimulation'
wk_dir = paste0(rootdir, '/ad_module_enrichment')
setwd(paste0(rootdir, '/ad_module_enrichment'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "_ad_module_"

#####################################################
# fetch lists from directories
f_names = c("2", "9", "12")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
# pattern = c("1", "2", "5", "6", "8", "9", "10")
replacement = c("Genotype Effect", "Drug Effect", "Genotype + Drug Effect")

###############FUNCTIONS
# function runs enrichment and saves to file for p and padj
# takes two set lists, background gene count and output name and pdf width
enrich = function(list1, list2, background, outname, width=7, fdr=TRUE) {
  # enrichment basic
  gom.obj = newGOM(list1, list2, background)
  # geneoverlap heatmap
  pdf(file=paste0(today, nameset, outname, "_pval.pdf"), width=width*2)
  drawHeatmap(gom.obj, grid.col="Blues", cutoff=0.1)
  dev.off()
  # write pvals to file
  a = getMatrix(gom.obj, name="pval")
  b = melt(data.table(a, keep.rownames="sample"), variable.name="group", 
           value.name="P", id.vars="sample")
  b[, FDR := p.adjust(P, method = "fdr") ]
  b[, LOGP := -log10(P) ]
  b[, LOGFDR := -log10(FDR) ]
  # b[!is.finite(b)] = -log10(2.225074e-308)
  fwrite(b, file=paste0(today, nameset, outname, ".txt"), sep="\t", quote=F)
  # log10p heatmap
  if (isTRUE(fdr)) {
    rws = rownames(a)
    cls = colnames(a)
    c = -log10(matrix(p.adjust(as.vector(a), method='fdr'), ncol=length(cls),
                      nrow=length(rws), dimnames=list(rws, cls)))
    pdf(file=paste0(today, nameset, outname, "_logfdr.pdf"), width=width)
    draw(heat(c, "#762a83", "-log FDR"))
    dev.off()
  } else { 
    c = -log10(a)
    pdf(file=paste0(today, nameset, outname, "_logp.pdf"), width=width)
    draw(heat(c, "#762a83", "-log P"))
    dev.off()
  }
}

# plot heatmap given matrix, color, name
heat = function(x, tone, name) {
  # col_fun = colorRamp2(c(0, max(c(ceiling(max(x)),3))), c("#FFFFFF", tone))
  # col_fun = colorRamp2(c(0, 6), c("#FFFFFF", tone))
  col_fun = rev(RColorBrewer::brewer.pal(name="RdYlBu", n=11))
  x[!is.finite(x)] = -log10(2.225074e-308)
  a = Heatmap(x, cluster_rows=F, cluster_columns=F, col=col_fun, 
              name = name, column_names_rot = 60)
          #column_order=tiss$full, row_order=gwas)
  return(a)
}

###############MAIN
# read in files
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

## Fisher exact testing for enrichment of Wan et al 2020 AD modules

# filter out pseudogenes and TEC
keep_categ = grep("pseudo", unique(result$Gene_type), value=T, invert=T)
keep_categ = grep("TEC", keep_categ, value=T, invert=T)
result = result[ Gene_type %in% keep_categ ]

# overall gene list
gs = nrow(result)

# AD
# generate list of gene sets for AD genes
de_list = lapply(f_names, function(i) {
  result[get(paste0("pvalue.", i)) < 0.05, Gene_id]
})
names(de_list) = replacement

# AD (split up/down)
# make list
f_names_updown = c(paste(replacement, "Up"), paste(replacement, "Down"))
# up
up_list = lapply(f_names, function(i) {
  result[get(paste0("pvalue.", i)) < 0.05 
         & get(paste0("log2FoldChange.", i)) > 0, Gene_id]
})
# down
down_list = lapply(f_names, function(i) {
  result[get(paste0("pvalue.", i)) < 0.05 
         & get(paste0("log2FoldChange.", i)) < 0, Gene_id]
})
# join and name
updown_list = c(up_list, down_list)
names(updown_list) = f_names_updown

# Wan et al
# make AD list
ad = fread(paste0(wk_dir, "/AD_modules_full.tsv"))
# drop genes w/o orthologs
ad = ad[Mmus_GeneID != ""]
# order by Consensus Cluster (for figure coloring order)
setorder(ad, Consensus_Cluster, Module)
ad_names = unique(ad$Module)
# build list for AD modules
ad_list = sapply(ad_names, function(k) {
  k = ad[Module == k, Mmus_GeneID]
}, USE.NAMES=T)
# ad_list = list()
# for (k in ad_names) {
  # ad_list[[k]] = ad[Module == k, Mmus_GeneID]
# }

# enrich for all DE genes
enrich(de_list, ad_list, gs, "AD_Modules", width=20)

# enrich for up/down genes
enrich(updown_list, ad_list, gs, "AD_Modules_Up-Down", width=20)





