#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(data.table)
library(stringi)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
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
replacement = c("Genotype - Tg_V vs Wt_V", "Drug - Tg_D vs Tg_V", 
                "Genotype + Drug - Tg_D vs Wt_V")

####################################################################FUNCTIONS
# function runs enrichment and saves to file for p and padj
# takes two set lists, background gene count and output name and pdf width
enrich = function(list1, list2, background, outname, width=7, fdr=TRUE, ...) {
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
  # write ORs to file
  c = getMatrix(gom.obj, name="odds.ratio")
  d = melt(data.table(c, keep.rownames="sample"), variable.name="group", 
           value.name="OR", id.vars="sample")
  b = merge(b, d, by=c("sample", "group"))
  # write intersection to file
  f = getMatrix(gom.obj, name="intersection")
  g = melt(data.table(f, keep.rownames="sample"), variable.name="group", 
           value.name="intersection", id.vars="sample")
  b = merge(b, g, by=c("sample", "group"))
  # write union to file
  h = getMatrix(gom.obj, name="union")
  i = melt(data.table(h, keep.rownames="sample"), variable.name="group", 
           value.name="union", id.vars="sample")
  b = merge(b, i, by=c("sample", "group"))
  # # # # write overalp to file
  # # # j = f/h
  # # # k = melt(data.table(j, keep.rownames="sample"), variable.name="group", 
           # # # value.name="overalp", id.vars="sample")
  # # # b = merge(b, k, by=c("sample", "group"))
  # pct_overlap
  b[, pct_overlap := intersection/union]
  # add fdr
  b[, FDR := p.adjust(P, method = "fdr") ]
  b[, LOGP := -log10(P) ]
  b[, LOGFDR := -log10(FDR) ]
  # b[!is.finite(b)] = -log10(2.225074e-308)
  fwrite(b, file=paste0(today, nameset, outname, ".txt"), sep="\t", quote=F)
  # print(c("before if", (f/h)[1,1]))
  # log10p heatmap
  if (isTRUE(fdr)) {
    rws = rownames(a)
    cls = colnames(a)
    e = -log10(matrix(p.adjust(as.vector(a), method='fdr'), ncol=length(cls),
                      nrow=length(rws), dimnames=list(rws, cls)))
    # filename and heatmap
    # print(c("during if", (f/h)[1,1]))
    pdf(file=paste0(today, nameset, outname, "_logfdr.pdf"), width=width)
    draw(heat(e, e, f, "#762a83", "-log FDR"))
    dev.off()
  } else { 
    e = -log10(a)
    # filename and heatmap
    pdf(file=paste0(today, nameset, outname, "_logp.pdf"), width=width)
    draw(heat(e, e, f, "#762a83", "-log P"))
    dev.off()
  }
  # OR heatmap
  rws = rownames(c)
  cls = colnames(c)
  # filename and heatmap
  pdf(file=paste0(today, nameset, outname, "_logOR.pdf"), width=width)
  draw(heat(c, e, f, "#762a83", "Odds Ratio"))
  dev.off()
}

# plot heatmap given matrix, color, name
heat = function(mat, pval, int, tone, name) {
  # round Inf to real number
  mat[is.infinite(mat)] = -log10(2.225074e-308)
  pval[is.infinite(pval)] = -log10(2.225074e-308)
  # color palette
  col_fun = colorRamp2(seq(min(mat, na.rm=T), max(mat, na.rm=T), length=9), 
                       brewer.pal(name="BuPu", n=9))
  # column annotations
  ad_anno = unique(ad, by=c("Module", "Consensus_Cluster"))[, Consensus_Cluster]
  ha = HeatmapAnnotation("Consensus Cluster"=ad_anno, 
                         col=list("Consensus Cluster"=c("A"="#FAFF34", 
                                                        "B"="#62B0B0", 
                                                        "C"="#7E2221", 
                                                        "D"="#6CC835", 
                                                        "E"="#280808")))
  # find max int buffer (for circles)
  # buffer = (1-max(int, na.rm=T))/2
  # heatmap plot 
  # print(c("in funtion", int[1,1]))
  a = Heatmap(mat, cluster_rows=F, cluster_columns=F, na_col="#FFFFFF", 
              col=col_fun, name=name, column_names_rot=60, top_annotation=ha, 
              column_split=ad_anno, column_gap = unit(5, "mm"), 
              rect_gp=gpar(type="none"),
              cell_fun=function(j, i, x, y, width, height, fill) { 
                grid.circle(x=x, y=y, r=(((int[i, j]/max(int, na.rm=T))/4 + .25)
                                          * min(unit.c(width, height))), 
                            gp=gpar(fill=col_fun(mat[i, j]), 
                                    col=ifelse(pval[i, j] < -log10(0.05), NA, 
                                               "black")))
                if ((pval[i, j] > -log10(0.05)) & (int[i, j] > 200)) {
                  grid.text(sprintf("%i", int[i, j]), x, y, 
                            gp=gpar(fontsize=10, 
                                    col=ifelse(mat[i, j] > mean(mat), "white", 
                                               "black")))
                }
              })
  return(a)
}

# add gom objects to a data.table for use

####################################################################MAIN
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

# DE genes
# generate list of gene sets for AD genes
de_list = lapply(f_names, function(i) {
  result[get(paste0("pvalue.", i)) < 0.05, Gene_id]
})
names(de_list) = replacement

# DE genes (split up/down)
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

# enrich for all DE genes
enrich(de_list, ad_list, gs, "AD_Modules", width=20)

# enrich for up/down genes
enrich(updown_list, ad_list, gs, "AD_Modules_Up-Down", width=20)





