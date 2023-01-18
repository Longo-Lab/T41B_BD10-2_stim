#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 20211014_BD-self_enrichment.R
##
## Version: 0.0.1
##
## Purpose of script: Show Fisher's Exact Enrichment of TREAT-AD Biodomains in 
## DE genesets
##
## Author: Robert B Butler III
##
## Date Created: 2021-10-14
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
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(GeneOverlap)
library(biomaRt)

# directory locations
rootdir = '/labs/flongo/t41b_BD10-2_stim'
wk_dir = paste0(rootdir, '/stim_biodomain_enrichment')
setwd(paste0(rootdir, '/stim_biodomain_enrichment'))
today = format(Sys.Date(), '%Y%m%d')
# genesets to be analyzed (path and prefix)
treatdir = "/labs/flongo/TREAT-AD/Biodomain_annotations"
nameset = "BD-SELF-trace"
# DE cutoff by "pvalue" or "padj"
p = "padj"

#####################################################
# fetch lists from directories
f_names = c("1", "2")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
replacement = c("APP", "APP + Stim")

####################################################################FUNCTIONS
gommat_2_dt = function(go, dt, dtype=c("pval", "odds.ratio", 
                       "intersection", "union", "Jaccard")) {
  # add gom object to a data.table
  # gom.obj: GeneOverlap object
  # dt: data.table to add to
  # dtype: which GO measure to extract
  d = melt(data.table(getMatrix(go, name=dtype), keep.rownames="sample"), 
           variable.name="group", value.name=dtype, id.vars="sample")
  dt = merge(dt, d, by=c("sample", "group"))
  return(dt)
}

get_mcols = function(str, result, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # result: data table with ens_ids in Gene_id
  # col: colname to grab
  require("stringr")
  y = result[[col]][match(unlist(str_split(str, ",")), result$Gene_id)]
  return(paste0(y, collapse=","))
}

gomtable = function(list1, list2, background, outname, width=7, 
                    height=7, ...) {
  # Takes two gene overlap lists and writes a table of stats, makes list object
  # list1: named list object
  # list2: named list object
  # outname: name prefix for file
  # width: width of pdf files to be plotted in inches*2
  # height: height of pdf files
  # returns a list object of table, r_names, c_names
  # enrichment basic
  gom.obj = newGOM(list1, list2, background)
  r_names = names(list1)
  c_names = names(list2)
  # geneoverlap heatmap
  pdf(file=paste(today, nameset, outname, "pval.pdf", sep="."), 
      width=width*2, height=height)
  drawHeatmap(gom.obj, grid.col="Blues", cutoff=0.1)
  dev.off()
  # write DE vals to file
  b = melt(data.table(getMatrix(gom.obj, name="pval"), keep.rownames="sample"), 
           variable.name="group", value.name="P", id.vars="sample")
  # write dtypes to table
  b = gommat_2_dt(gom.obj, b, dtype="odds.ratio")
  b = gommat_2_dt(gom.obj, b, dtype="Jaccard")
  b = gommat_2_dt(gom.obj, b, dtype="intersection")
  b = gommat_2_dt(gom.obj, b, dtype="union")
  # add fdr
  b[, FDR := p.adjust(P, method = "fdr") ]
  b[, LOGP := -log10(P) ]
  b[, LOGFDR := -log10(FDR) ]
  # get intersection genes
  nest = getNestedList(gom.obj, name="intersection")
  nest2 = lapply(seq(nest), function(i) {
    lapply(nest[[i]], function(j) paste0(j, collapse=","))
  })
  nest3 = lapply(seq(nest2), function(i) stack(nest2[[i]]))
  setattr(nest3, 'names', names(nest))
  dt = rbindlist(nest3, idcol=T)
  setnames(dt, c("group", "ens_id", "sample"))
  dt$symbol = lapply(dt$ens_id, get_mcols, result, "GeneSymbol")
  b = merge(b, dt, by=c("sample", "group"), all.x=T, sort=F)
  # write and return
  fwrite(b, file=paste(today, nameset, outname, "txt", sep="."), sep="\t", 
         quote=F)
  return(list("b"=b, "r_names"=r_names, "c_names"=c_names))
}

enrich = function(list1, list2, gs, outname, width=7, height=7, fdr=TRUE, ...) {
  # function runs enrichment and saves to file for odds ratio and padj/p
  # takes two set lists, background gene count and output name and pdf width
  # DE enrich
  de_enrich = gomtable(list1, list2, gs, outname, width=width, height=height)

  # fetch Jaccard index
  f = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="Jaccard"), 
                rownames="sample")
  f = f[de_enrich$r_names, de_enrich$c_names, drop=F]
  # log10p heatmap
  if (isTRUE(fdr)) {
    e = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="LOGFDR"), 
                  rownames="sample")
    e = e[de_enrich$r_names, de_enrich$c_names, drop=F]
    # filename and heatmap
    pdf(file=paste(today, nameset, outname, "logfdr.pdf", sep="."), width=width, 
        height=height)
    heat(e, e, f, "-log FDR")
    dev.off()
  } else { 
    e = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="LOGP"), 
                  rownames="sample")
    e = e[de_enrich$r_names, de_enrich$c_names, drop=F]
    # filename and heatmap
    pdf(file=paste(today, nameset, outname, "logp.pdf", sep="."), width=width, 
        height=height)
    heat(e, e, f, "-log P")
    dev.off()
  }
  # OR heatmap
  c = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="odds.ratio"), 
                rownames="sample")
  c = c[de_enrich$r_names, de_enrich$c_names, drop=F]
  # filename and heatmap
  pdf(file=paste(today, nameset, outname, "logOR.pdf", sep="."), width=width, 
      height=height)
  heat(c, e, f, "Odds Ratio")
  dev.off()
}

heat = function(mat, pval, int, name) {
  # plot heatmap of overlaps colored by value, and sized by overlap (given matrices, name)
  # mat: matrix of heatmap values (fdr, p-val, or odds.ratio)
  # pval: matched matrix of fdr values, to bold significant circles
  # int: matched matrix of Jaccard indexes for intersection sizes
  # name: heatscale name for plot (fdr, p-val, or odds.ratio)
  # returns a draw graph call

  # resolve Inf to real number (round) and check for NAs
  mat[is.infinite(mat)] = -log10(2.225074e-308)
  pval[is.infinite(pval)] = -log10(2.225074e-308)
  if (any(is.na(int))) {
    print("Jaccard index contains NAs!")
    quit("no", 1)
  }

  # color palette
  col_fun = colorRamp2(seq(min(mat, na.rm=T), max(mat, na.rm=T), length=9), 
                       brewer.pal(name="BuPu", n=9))
  col_fun2 = brewer.pal(name="Dark2", n=7)
  # # column annotations
  # ad_anno = unique(ad, by=c("Module", "Consensus_Cluster"))[, Consensus_Cluster]
  # ad_tiss = unique(ad, by=c("Module", "brainRegion"))[, brainRegion]
  # names(col_fun2) = unique(ad_tiss)
  # ha = HeatmapAnnotation("Consensus Cluster"=ad_anno, 
                         # col=list("Consensus Cluster"=c("A"="#FAFF34", 
                                                        # "B"="#62B0B0", 
                                                        # "C"="#7E2221", 
                                                        # "D"="#6CC835", 
                                                        # "E"="#280808")))
  # hb = HeatmapAnnotation("Brain Region"=ad_tiss, col=list("Brain Region"=col_fun2))
  # scale max int (for circle size) to matrix "buffer"
  buffer = ((int/max(int))/4 + .25)
  dot_scalar = unit(1.5, "cm")
  # heatmap plot 
  a = Heatmap(mat, cluster_rows=F, cluster_columns=T, na_col="#FFFFFF", 
              col=col_fun, name=name, column_names_rot=60, 
              # top_annotation=ha, bottom_annotation=hb, column_split=ad_anno, 
              # column_gap = unit(5, "mm"), 
              # column_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp=gpar(type="none"), 
              cell_fun=function(j, i, x, y, width, height, fill) { 
                grid.circle(x=x, y=y, r=(buffer[i, j] * dot_scalar), 
                            gp=gpar(fill=col_fun(mat[i, j]), 
                                    col=ifelse(pval[i, j] < -log10(0.05), NA, 
                                               "black")))
                if (pval[i, j] > -log10(0.05)) {
                  grid.text(sprintf("%.2f", int[i, j]), x, y, 
                            gp=gpar(fontsize=10, fontface="bold", 
                                    col=ifelse(mat[i, j] > quantile(mat, 0.95), 
                                               "white", "black")))
                }
              })
  # heatmap legend
  lgd_list = list(
    Legend(labels=c(rep("", 3)), 
           title = "Jaccard Index",
           grid_height = unit(12, "mm"), 
           grid_width = unit(12, "mm"), 
           graphics = list(
              function(x, y, w, h) {
                grid.circle(x, y, r=min(buffer) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", 0), x, y, 
                          gp=gpar(fontsize=10, fontface="bold", col="white"))
              },
              function(x, y, w, h) {
                grid.circle(x, y, r=mean(buffer) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", mean(int)), x, y, 
                          gp=gpar(fontsize=10, fontface="bold", col="white"))
              },
              function(x, y, w, h) {
                grid.circle(x, y, r=max(buffer) * dot_scalar,
                            gp=gpar(fill="black", col="white"))
                grid.text(sprintf("%.2f", max(int)), x, y, 
                          gp=gpar(fontsize=10, fontface="bold", col="white"))
              }
           )
    )
  )
  draw(a, merge_legend=T, annotation_legend_list=lgd_list)
}

convertHumanMouse = function(dt, colname, keep.human=F){
  # for a given column of human ensIDs, replace them with mouse ensembl gene ids 
  # includes a filter for autosomal genes that are not pseudo or small RNAs
  # dt: data table to replace names
  # colname: name of the column containing gene symbols
  # keep.symbols: boolean, retain the symbols column?
  # returns dt with a GENE column
  require("biomaRt")
  print(paste(length(unique(dt[[colname]])), "genes to look up"))
  mouse = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version="104")
  human = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version="104")
  # match with ensembl_gene_id across both
  x = data.table(getLDS(mart=human, attributes=c('ensembl_gene_id'), martL=mouse,
                        attributesL=c('ensembl_gene_id', 'chromosome_name', 
                                      'gene_biotype'), 
                        filters='ensembl_gene_id', values=unique(dt[[colname]]),
                        bmHeader=F))
  names(x) = make.unique(names(x))
  print(paste(nrow(x), "total genes returned"))
  # restricting to only protein coding 
  remove_categ = c("miRNA", "polymorphic_pseudogene")
  x = x[!(gene_biotype %in% remove_categ)]
  print(paste(nrow(x), "genes correct biotype"))
  # use only canonical chromosomes
  use_chrs = c(1:19, "X", "Y")
  x = x[chromosome_name %in% use_chrs]
  print(paste(nrow(x), "genes in chr 1:19 X Y"))
  setnames(x, c("ensembl_gene_id", "ensembl_gene_id.1"), c("hGENE", "GENE"))
  # replace gene names and write to file
  dt = merge(x[, .(GENE, hGENE)], dt, all.y=T, by.x="hGENE", by.y=colname)
  dt = dt[!is.na(GENE)]
  print(paste(length(unique(dt[["GENE"]])), "unique genes with ensembl ids"))
  if (keep.human == FALSE) dt[, hGENE := NULL]
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
  result[get(paste0(p, ".", i)) < 0.05, Gene_id]
})
names(de_list) = replacement

# DE genes (split up/down)
# make list
f_names_updown = c(paste(replacement, "Up"), paste(replacement, "Down"))
# up
up_list = lapply(f_names, function(i) {
  result[get(paste0(p, ".", i)) < 0.05 
         & get(paste0("log2FoldChange.", i)) > 0, Gene_id]
})
# down
down_list = lapply(f_names, function(i) {
  result[get(paste0(p, ".", i)) < 0.05 
         & get(paste0("log2FoldChange.", i)) < 0, Gene_id]
})
# join and name
updown_list = c(up_list, down_list)
names(updown_list) = f_names_updown


#######################
# TREAT-AD Biodomains
treat_files = list.files(path=treatdir, pattern=paste0(nameset, "_"), full.names=T)
names(treat_files) = gsub(paste0(treatdir, "/", nameset, "_"), "", 
                          gsub("\\.geneset", "", treat_files))
treat_list = lapply(treat_files, fread, header=F)
treat_mlist = lapply(treat_list, convertHumanMouse, colname="V1")
treat_mlist2 = lapply(treat_mlist, function(i) i[, GENE])


##############
# background DE genes list
gs = length(unique(result[,Gene_id], unlist(treat_mlist2)))
fwrite(list(unique(result[,Gene_id], unlist(treat_mlist2))), 
       file=paste(nameset, "background.geneset", sep="."))

# enrich for all DE genes
enrich(de_list, treat_mlist2, gs, "Modules_allDE", width=20, height=6)

# enrich for up/down genes
enrich(updown_list, treat_mlist2, gs, "Modules_Up-Down", width=20, height=6)





