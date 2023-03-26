#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(data.table)
library(stringi)
library(Cairo)


# directory locations
rootdir='/labs/flongo/t41b_BD10-2_stim'
setwd(paste0(rootdir, '/l2fc_comparison'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "APP_Activity_Dep"

# fetch lists from directories
f_names = data.table(short=c("5", "6"),
                     long=c("TBS Effect (in Wt)", "TBS Effect (in APP)"))
target_files = unlist(lapply(f_names$short, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))

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
      as.character(get(labcol)),"")), 
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
      as.character(get(labcol)),"")), 
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
both = merge(fl[[1]], fl[[2]], by=c("Gene_id", "GeneSymbol"), 
             suffix=c(".x", ".y"))
both = both[!is.na(Gene_id)]
exclude = c("TEC", "snoRNA", "misc_RNA", "rRNA", "ribozyme", "snRNA", 
            "processed_pseudogene", "transcribed_processed_pseudogene", 
            "transcribed_unprocessed_pseudogene", "unprocessed_pseudogene")
both = both[!(Gene_type.y %in% exclude)]
both = both[padj.x < 0.05 | padj.y < 0.05]

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








