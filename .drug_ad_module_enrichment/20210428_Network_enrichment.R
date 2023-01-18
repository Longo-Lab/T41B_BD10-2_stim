#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(GeneOverlap)

# directory locations
rootdir = '/labs/flongo/2020_T41B_BD10-2_Stimulation'
wk_dir = paste0(rootdir, '/ad_module_enrichment')
setwd(paste0(rootdir, '/ad_module_enrichment'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "_WGCNA_"

#####################################################
# fetch lists from directories
f_names = c("2", "9", "12")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
replacement = c("Genotype - Tg_V vs Wt_V", "Drug - Tg_D vs Tg_V", 
                "Genotype + Drug - Tg_D vs Wt_V")

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

get_mcols = function(str, dt, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # dt: data table with ens_ids in id
  # col: colname to grab
  require("stringr")
  y = dt[[col]][match(unlist(str_split(str, ",")), dt$id)]
  return(paste0(y, collapse=","))
}

gomtable = function(list1, list2, background, outname, matchname=NULL, width=7, 
                    height=7, ...) {
  # Takes two gene overlap lists and writes a table of stats, makes list object
  # list1: named list object
  # list2: named list object
  # outname: name prefix for file
  # matchname: what is the overlap being checked
  # width: width of pdf files to be plotted in inches*2
  # height: height of pdf files
  # returns a list object of table, r_names, c_names
  # enrichment basic
  gom.obj = newGOM(list1, list2, background)
  r_names = names(list1)
  c_names = names(list2)
  # geneoverlap heatmap
  pdf(file=paste0(today, nameset, outname, matchname, "_pval.pdf"), 
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
  dt$symbol = lapply(dt$ens_id, get_mcols, mod, "GeneSymbol")
  b = merge(b, dt, by=c("sample", "group"), all.x=T, sort=F)
  # write and return
  fwrite(b, file=paste0(today, nameset, outname, matchname, ".txt"), sep="\t", 
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
    pdf(file=paste0(today, nameset, outname, "_logfdr.pdf"), width=width, 
        height=height)
    heat(e, e, f, "-log FDR")
    dev.off()
  } else { 
    e = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="LOGP"), 
                  rownames="sample")
    e = e[de_enrich$r_names, de_enrich$c_names, drop=F]
    # filename and heatmap
    pdf(file=paste0(today, nameset, outname, "_logp.pdf"), width=width, 
        height=height)
    heat(e, e, f, "-log P")
    dev.off()
  }
  # OR heatmap
  c = as.matrix(dcast(de_enrich$b, sample ~ group, value.var="odds.ratio"), 
                rownames="sample")
  c = c[de_enrich$r_names, de_enrich$c_names, drop=F]
  # filename and heatmap
  pdf(file=paste0(today, nameset, outname, "_logOR.pdf"), width=width, 
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
  # column annotations
  ad_anno = unique(ad, by=c("Module", "Consensus_Cluster"))[, Consensus_Cluster]
  ad_tiss = unique(ad, by=c("Module", "brainRegion"))[, brainRegion]
  names(col_fun2) = unique(ad_tiss)
  ha = HeatmapAnnotation("Consensus Cluster"=ad_anno, 
                         col=list("Consensus Cluster"=c("A"="#FAFF34", 
                                                        "B"="#62B0B0", 
                                                        "C"="#7E2221", 
                                                        "D"="#6CC835", 
                                                        "E"="#280808")))
  hb = HeatmapAnnotation("Brain Region"=ad_tiss, col=list("Brain Region"=col_fun2))
  # scale max int (for circle size) to matrix "buffer"
  buffer = ((int/max(int))/4 + .25)
  dot_scalar = unit(1.5, "cm")
  # heatmap plot 
  a = Heatmap(mat, cluster_rows=F, cluster_columns=F, na_col="#FFFFFF", 
              col=col_fun, name=name, column_names_rot=60, top_annotation=ha, 
              column_split=ad_anno, column_gap = unit(5, "mm"), 
              column_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp=gpar(type="none"), bottom_annotation=hb, 
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


####################################################################MAIN
## WGCNA genes
# generate list of gene sets for AD WGCNA gene modules
mod = fread(paste0(rootdir, "/wt_veh_50_gene_per_module.csv"), header=T)
mod_names = unique(mod$module)
# drop unclustered
mod_names = mod_names[mod_names != "grey"]
# build list for AD WGCNA modules
mod_list = sapply(mod_names, function(k) {
  k = mod[module == k, id]
}, USE.NAMES=T)

## Wan et al
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

# background WGCNA genes list
gs2 = length(unique(mod[,id], ad[,Mmus_GeneID]))

# enrich for all DE genes
enrich(mod_list, ad_list, gs2, "AD_Modules", width=20, height=20)




