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

#####################################################
# fetch lists from directories
f_names = c("5", "6")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
pattern = c("log2FoldChange.", "baseMean.", "5", "6")
replacement = c("Log2 Fold Change ", "Base Mean Expression ", "Stimulation Effect (in Wt)",
                "Stimulation Effect (in APP)")

###############FUNCTIONS
# plots for AAIC -----------------------------------
# Functions -------------------------------------------------------
# GProfiler2 - loop through each pheno in list and run gprofiler2
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
                               coly="avg_log2FC.geno", gencol="gene", out="") {
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

# Main ---------------------------------------------
# L2FC comparison
lapply(c(""), function(i) {
  fl = lapply(target_files, fread)
  both = merge(fl[[1]], fl[[2]], by=c("Gene_id", "GeneSymbol"), 
               suffix=c(".geno", ".stim"))
  both = both[!is.na(Gene_id)]
  both = both[(pvalue.stim < 0.05 & pvalue.geno < 0.05)]

  # correlation plot function (w/ pearson)
  plot_corr = function (table, x, y, title=NULL) {
    sub = sprintf("R = %.3f; p < %.3g", 
                  cor.test(table[[x]], table[[y]])$estimate[[1]], 
                  max(cor.test(table[[x]], table[[y]])$p.value, .Machine$double.xmin))
    ggscatter(table, x=x, y=y, add="reg.line", conf.int=T, col="#8073ac",
              add.params=list(color="#053061", fill="#fee0b6")) + 
      ggtitle(title, subtitle=sub) +
      geom_hline(yintercept=0, linetype="dashed", col='gray') +
      geom_vline(xintercept=0, linetype="dashed", col='gray') +
      geom_text_repel(
        # aes(label=ifelse(abs(get(x)) > 0.5 | abs(get(y)) > 0.55, 
        aes(label=ifelse(get(x) > 5 | get(y) > 5 | get(x) < -2 | get(y) < -2, 
        as.character(GeneSymbol),"")), max.overlaps=40)
  }
  p1 = plot_corr(both, "log2FoldChange.stim", "log2FoldChange.geno")
  ggsave(filename=paste(today, nameset, "L2FC.pdf", sep="."), 
         plot=p1, device=cairo_pdf)
  ggsave(filename=paste(today, nameset, "L2FC.png", sep="."), 
         plot=p1, type=cairo)
  # look up term enrichments
  get_quartile_enrich(both, 1, colx="log2FoldChange.stim", 
                      coly="log2FoldChange.geno", gencol="Gene_id", out=i)
  get_quartile_enrich(both, 2, colx="log2FoldChange.stim", 
                      coly="log2FoldChange.geno", gencol="Gene_id", out=i)
  get_quartile_enrich(both, 3, colx="log2FoldChange.stim", 
                      coly="log2FoldChange.geno", gencol="Gene_id", out=i)
  get_quartile_enrich(both, 4, colx="log2FoldChange.stim", 
                      coly="log2FoldChange.geno", gencol="Gene_id", out=i)
})

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








