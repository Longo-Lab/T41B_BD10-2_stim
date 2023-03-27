#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 20230326_overlap_correlations_plot.R
##
## Version: 0.0.1
##
## Purpose of script: Correlation and Fold Change comparison for stimulation effect
##
## Author: Robert B Butler III
##
## Date Created: 2023-03-26
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   No longer generalized to x vs y comparisons of fold change, special one off.
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(data.table)
library(stringi)
# library(Cairo)


# directory locations
rootdir='/labs/flongo/t41b_BD10-2_stim'
setwd(paste0(rootdir, '/l2fc_comparison'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "APP_Activity_Dep"

# fetch lists from directories
f_names = data.table(short=c("5", "6", "7"),
                     long=c("TBS Effect (in Wt)", "TBS Effect (in APP)", 
                            "TBS Effect (in APP BD10-2)"))

# Functions -------------------------------------------------------
get_cols = function(str, dds, lookup, return) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # dds: datatable with with GENE col containing ENSGs
  # lookup: colname to lookup list of terms
  # return: colname to return matching values from
  require("stringr")
  y = dds[match(unlist(str_split(str, ",")), dds[, get(lookup)]), get(return)]
  return(paste0(y, collapse=","))
}

# get q2 dataframe and order by the sum of abs L2fc
get_quartile_enrich = function(dt, qtl=c(1, 2, 3, 4), colx="avg_log2FC.stim", 
                               coly="avg_log2FC.geno", gencol="gene", out=NULL) {
  # for a merged datatable of two snRNAseq DEs, looks up a given qudrant in gProfiler
  # dt: merged data table
  # qtl: which quadrant to look at
  # colx: x-axis column in merged table
  # coly: y-axis column in merged table
  # gencol: column containing gene names
  # out: appends an identifier to the filename
  # returns: nothing, writes results to file
  require(gprofiler2)
  # sort which quadrant
  if (qtl == 1) {x = dt[get(colx) > 0 & get(coly) > 0]}
  else if (qtl == 2) {x = dt[get(colx) < 0 & get(coly) > 0]}
  else if (qtl == 3) {x = dt[get(colx) < 0 & get(coly) < 0]}
  else if (qtl == 4) {x = dt[get(colx) > 0 & get(coly) < 0]}
  # sort by greatest sum of abs(L2FC)
  x = x[order(-(abs(get(colx)) + abs(get(coly))))]
  # lookup using ordered subset
  g = gost(query=x[[gencol]], organism="mmusculus", ordered_query=T, evcodes=T,
           significant=T)
  # retrieve gene names
  # gl = data.table(stack(g$meta$genes_metadata$query$query_1$mapping))
  # write to file
  if(is.list(g$result)){
    # get symbols for terms
    g$result$ids = lapply(g$result$intersection, get_cols, dt, "Gene_id", "GeneSymbol")
    # write to file
    fwrite(g$result, paste(today, nameset, out, sprintf("q%s.csv", qtl), sep="."))
    saveRDS(g, paste(today, nameset, out, sprintf("q%s.rds", qtl), sep="."))
  }
}

# correlation plot function (w/ pearson)
plot_corr = function (dt, x, y, labcol="GeneSymbol", title=NULL, ...) {
  # scatterplot with a trendline and repel labels
  # dt: combined data.table for both axes to plot
  # x: column name to plot on x
  # y: column name to plot on y
  # labcol: column name for point labels
  # title: main title, sub will be the correlation equation
  # ...: additional objects passed to ggscatter
  # returns: ggobject
  sub = sprintf("R = %.3f; p < %.3g", 
                cor.test(dt[[x]], dt[[y]])$estimate[[1]], 
                max(cor.test(dt[[x]], dt[[y]])$p.value, .Machine$double.xmin))
  ggscatter(dt, x=x, y=y, add="reg.line", conf.int=T, col="#8073ac",
            add.params=list(color="#053061", fill="#fee0b6"), ...) + 
    ggtitle(title, subtitle=sub) + ### add a graying for close to y=x
    geom_hline(yintercept=0, linetype="dashed", col='gray') +
    geom_vline(xintercept=0, linetype="dashed", col='gray') +
    geom_text_repel(
      aes(label=ifelse(get(x) > 6 | get(y) > 6, 
      get(labcol),"")), 
                    nudge_x=1,
                    nudge_y=1,
                    size=5,
                    fontface=2,
                    force=10, 
                    force_pull=0.1, 
                    max.overlaps=40, 
                    point.size=3, 
                    min.segment.length=0,
                    segment.color="gray") +
    scale_x_continuous(expand=expansion(mult=0.1)) +
    geom_text_repel(
      aes(label=ifelse(get(x) < -2.75 | get(y) < -2.75, 
      get(labcol),"")), 
                    nudge_x=-1,
                    nudge_y=-1,
                    size=5,
                    fontface=2,
                    force=10, 
                    force_pull=0.1, 
                    max.overlaps=40, 
                    point.size=3, 
                    min.segment.length=0,
                    segment.color="gray") +
    scale_y_continuous(expand=expansion(mult=0.1))
}

# Main ---------------------------------------------
# L2FC comparison
fl = lapply(paste0(rootdir, paste0("/res.", f_names$short, ".csv")), fread)
# which are the important test comparison cols
my_cols = c("log2FoldChange", "pvalue", "padj")
# which columns are the keys to join on
key_cols = c("Gene_id", "GeneSymbol", "Gene_type")

# truncate lists, filter out entries with no Gene_id (7 matching errors)
fl = lapply(fl, function(i) {
  cols = c(key_cols, my_cols)
  i[!is.na(get(key_cols[1])), ..cols]
})

# join all sets
fn = c("x", "y", "z")
both = fl[[1]]
for (i in head(seq(fl), -1)) {
  both = merge(x=both, y=fl[[i+1]], by=key_cols, all=T, sort=F,
                 suffixes=c("", paste0(".", fn[i+1])))
}
setnames(both, my_cols, paste0(my_cols, ".", fn[1]))

# replacing NAs with 1
setnafill(both, cols=names(both)[names(both) %like% "padj."], fill=1)

# filter to see what shows up padj in LTP or no LTP (or both)
both = both[(padj.x < 0.05 & padj.z < 0.05) | padj.y < 0.05]
# difference in FC between APP and Wt, activity dependent changes
both[, deltaFC.yx := log2FoldChange.y - log2FoldChange.x]
# difference in FC between APP and APP BD10-2, restored activity dependent changes
both[, deltaFC.yz := log2FoldChange.y - log2FoldChange.z]
# difference in FC between APP BD10-2 and Wt, residual diff between LTP groups
both[, deltaFC.zx := log2FoldChange.z - log2FoldChange.x]
# values for coloring points
both[, xy.padj := ifelse(padj.x < 0.05, ifelse(padj.y < 0.05, "both", "x only"), "y only")]

# scatterplot of L2FCs
p1 = plot_corr(both, "log2FoldChange.x", "log2FoldChange.y", 
               xlab=paste("Log2 Fold Change", f_names$long[1]),
               ylab=paste("Log2 Fold Change", f_names$long[2]))
ggsave(filename=paste(today, nameset, "L2FC.pdf", sep="."), 
       plot=p1, device=cairo_pdf)
ggsave(filename=paste(today, nameset, "L2FC.png", sep="."), 
       plot=p1, type=cairo)
# look up term enrichments
get_quartile_enrich(both, 1, colx="log2FoldChange.x", 
                    coly="log2FoldChange.y", gencol="Gene_id")
get_quartile_enrich(both, 2, colx="log2FoldChange.x", 
                    coly="log2FoldChange.y", gencol="Gene_id")
get_quartile_enrich(both, 3, colx="log2FoldChange.x", 
                    coly="log2FoldChange.y", gencol="Gene_id")
get_quartile_enrich(both, 4, colx="log2FoldChange.x", 
                    coly="log2FoldChange.y", gencol="Gene_id")

# Grab gprofiler results
qs = c("q1", "q2", "q3", "q4")
fh = lapply(qs, function(i) lapply(c(""), function(j){
  paste(today, nameset, j, i, "csv", sep=".")
}))
# read in and fail silently for missing
fl = lapply(unlist(fh), function(i) try(fread(i), silent=T))
fn = lapply(qs, function(i) print(i, sep="."))
setattr(fl, 'names', unlist(fn))
# drop unneeded & combine
fl = fl[unlist(lapply(fl, is.data.table))]
fl = lapply(fl, function(i) i[, "query" := NULL])
tbl = rbindlist(fl, idcol="quadrant")
fwrite(tbl, paste(today, nameset, "all_quadrants.nominalp.csv", sep="."))


# Plot of differential L2FC between TBS in APP and Wt
exclude = c("TEC", "snoRNA", "misc_RNA", "rRNA", "ribozyme", "snRNA", 
            "processed_pseudogene", "transcribed_processed_pseudogene", 
            "transcribed_unprocessed_pseudogene", "unprocessed_pseudogene")
both = both[!(Gene_type %in% exclude)]

labs = data.frame(
  xpos = c(-Inf, -Inf, Inf, Inf),
  ypos = c(-Inf, Inf, -Inf, Inf),
  text = c("bottom-left", "bottom-right", "top-left", "top-right"),
  hj = c(0,0,1,1) ,
  vj = c(0,1,0,1)) #<- adjust

ggplot(both, aes(x=log2FoldChange.x, y=deltaFC.yx, color=xy.padj)) + 
  geom_point() + 
  xlab(f_names$long[1]) + 
  ylab("DeltaFC TBS Effect (APP - Wt groups)") +
  geom_hline(yintercept=0, col='black') +
  geom_text_repel(aes(label=ifelse(abs(deltaFC.yx) > 2, GeneSymbol,"")), 
    force=15,
    size=5,
    fontface=2,
    max.overlaps=40, 
    segment.color="gray",
    show.legend=F) + # set x limits for labels > 2.5 or < -2.5
  # scale_color_brewer(palette="GnBu", direction=-1) + 
  geom_text(data=labs, aes(x=xpos, y=ypos, hjust=hj, vjust=vj, label=text), inherit.aes=F) +
  theme_classic() +
  theme(axis.line.x=element_line(color="white"))




# # dummy gost data object to insert results into
# gost = gost(c("Lrfn5", "Dab1", "Slc17a7"), organism="mmusculus")
# gost$result = rbindlist(fl, idcol="query")
# gost$result = gost$result[term_size < 2000]

# # highlight terms
# terms = unique(gost$result[source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG")
                           # ][order(p_value)
                             # ][, head(.SD, 10), by="query"
                               # ][, term_id])

# CairoPDF(file=file.path(type, paste(today, nameset, type, "gost_plot.pdf", sep=".")), width=14, height=20)
# p = gostplot(gost, capped = FALSE, interactive = FALSE)
# pp = publish_gostplot(p, highlight_terms=terms)
# dev.off()







