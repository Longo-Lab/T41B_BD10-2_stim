#!/usr/bin/env Rscript

library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(patchwork)


# directory locations
rootdir='/labs/flongo/t41b_BD10-2_stim'
setwd(paste0(rootdir, '/enrichment_plots'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "top-bottom30"

###################################################################FUNCTIONS
heat = function(mat, pval, name, ...) {
  # plot heatmap of overlaps colored by value, and sized by overlap (given matrices, name)
  # mat: matrix of heatmap values (fdr, p-val, or odds.ratio)
  # pval: matched matrix of fdr values, to bold significant circles
  # name: heatscale name for plot (fdr, p-val, or odds.ratio)
  # returns a draw graph call

  # color palette
  col_fun = colorRamp2(c(-6, 0, 7), c("#5570b6", "#FFFFFF", "#f16a6c"))
  # col_fun = colorRamp2(seq(-7, 7, length=9), 
                       # rev(brewer.pal(name="RdBu", n=9)))
  # heatmap plot 
  a = Heatmap(mat, cluster_rows=F, cluster_columns=F, na_col="#FFFFFF", 
              col=col_fun, name=name, show_column_names=F, 
              heatmap_legend_param=list(at=c(-6, 0, 7)), 
              cell_fun=function(j, i, x, y, w=unit(1.5, "cm"), h, col) {
                width=unit(1.5, "cm")
                # add text to each grid
                if (pval[i, j] < 0.05) {
                  grid.text(sprintf("%s", "*"), x, y, 
                            gp=gpar(fontsize=10, fontface="bold", col="black"))
                }
              }, ...)
  draw(a)
}

get_mats = function(dt, rownam, colname1, colname2){
  # gets matricies for complex heatmap to use
  # dt: heatscale name for plot (fdr, p-val, or odds.ratio)
  # colname1: colname for first matrix (often log2fc)
  # colname2: colname for second matrix (often padj), values sorted by this column
  # returns a list of two matrices
  dt1 = dt[get(colname1) > 0][order(get(colname2))][1:30]
  dt2 = dt[get(colname1) < 0][order(get(colname2))][1:30]
  dt = rbindlist(list(dt1, dt2))
  a = as.matrix(dt[, c(mget(rownam), mget(colname1))], rownames=rownam)
  b = as.matrix(dt[, c(mget(rownam), mget(colname2))], rownames=rownam)
  return(list(a, b))
}

########################################################################MAIN
# fetch lists from directories
f_names = c("2", "9", "12")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
replacement = c("APP", "BD10-2", "APP BD10-2")

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

# heatmap
m.list = lapply(seq(f_names), function(i) {
  get_mats(result, "GeneSymbol", 
           paste("log2FoldChange", f_names[i], sep="."),
           paste("padj", f_names[i], sep="."))
})
p.list = lapply(seq(f_names), function(i) {
  grid.grabExpr(heat(m.list[[i]][[1]], m.list[[i]][[2]], "log2FoldChange"))
})

pdf(file=paste(today, nameset, "all.pdf", sep="."), width=7.5, height=14)
plot_grid(p.list[[1]], p.list[[2]], p.list[[3]], nrow = 1, vjust=0.25, 
          labels=replacement, rel_widths=c(1.03, 1.23, 1)) + 
  theme(plot.margin = unit(c(0.25, 0, 0, 0), "in"))
dev.off()
# ggsave(filename=paste(today, nameset, ".pdf", sep="_"), plot=(p1 + p2 + p3), device=cairo_pdf)

