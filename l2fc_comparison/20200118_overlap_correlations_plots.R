#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(data.table)
library(stringi)

# directory locations
rootdir='/labs/flongo/2020_T41B_BD10-2_Stimulation'
setwd(paste0(rootdir, '/l2fc_comparison'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "_overlap_correlations_"

#####################################################
# fetch lists from directories
f_names = c("1", "2", "6", "8", "9", "10")
target_files = unlist(lapply(f_names, function(i){
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
pattern = c("log2FoldChange.", "baseMean.", "1", "\\b2", "6", "8", "9", "10")
replacement = c("Log2 Fold Change ", "Base Mean Expression ", "Genotype Effect",
                "Genotype Effect w/ Stimulation", "Stimulation Effect", 
                "Drug Effect", "Drug Effect w/ Stimulation", 
                "Drug Effect Interaction")

###############FUNCTIONS
# correlation plot function (w/ pearson)
plot_corr = function (table, x, y, title=NULL) {
  sub = sprintf("R = %.3f; p < %.3g; n = %i", 
                cor.test(table[[x]], table[[y]])$estimate[[1]], 
                cor.test(table[[x]], table[[y]])$p.value,
                length(table[[x]]))
  ggscatter(table, x=x, y=y, add="reg.line", conf.int=T,
            add.params=list(color="blue", fill="lightgray")) + 
    ggtitle(title, subtitle=sub) +
    xlab(stri_replace_all_fixed(x, vectorize_all=F, pattern=c(pattern), 
                                replacement=c(replacement))) +
    ylab(stri_replace_all_fixed(y, vectorize_all=F, pattern=pattern, 
                                replacement=replacement)) +
    geom_hline(yintercept=0, linetype="dashed", col='gray') +
    geom_vline(xintercept=0, linetype="dashed", col='gray') +
    theme(text=element_text(size=20),
          plot.title=element_text(size=22, face="bold", hjust=0.5),
          plot.subtitle=element_text(size=18, hjust=0.5)
         )
}

# square a plot with the same x and y lims
squarePlot = function(plt){
    return(plt+coord_equal()+
            expand_limits(x=layer_scales(plt)$y$range$range,
                          y=layer_scales(plt)$x$range$range))
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
  x[!(is.na(Gene_id)), c("Gene_id", "GeneSymbol", "log2FoldChange", "baseMean", "pvalue", "padj")]
})
colnames = c("log2FoldChange", "baseMean", "pvalue", "padj")

# join all sets
result = short_list[[1]]
for (i in head(seq_along(short_list), -1)) {
  result = merge(x=result, y=short_list[[i+1]], by=c("Gene_id", "GeneSymbol"), all=T, 
                 suffixes=c("", paste0(".", f_names[i+1])))
}
setnames(result, colnames, paste0(colnames, ".", f_names[1]))

# plotting graphs
pdf(file=paste0(today, nameset, "plots.pdf"), onefile=T, paper="USr",
    width=11, height=8.5, useDingbats=F )
squarePlot(plot_corr(result[pvalue.8 < 0.05 & pvalue.10 > 0.05], 
          "log2FoldChange.8", "log2FoldChange.1", 
          paste("Relative Differential Expression of Drug Effect Genes\n",
          "Not Modified By Stimulation")))
squarePlot(plot_corr(result[pvalue.10 < 0.05], 
          "log2FoldChange.9", "log2FoldChange.1", 
          paste("Relative Differential Expression of Drug Effect Genes\n",
          "Altered By Stimulation")))
squarePlot(plot_corr(result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"], 
          "log2FoldChange.9", "log2FoldChange.1", 
          paste("Relative Differential Expression of Drug Effect Genes\n",
          "Altered By Stimulation")))
squarePlot(plot_corr(result[pvalue.8 < 0.05 & pvalue.10 > 0.05], 
          "log2FoldChange.6", "log2FoldChange.1", 
          paste("Relative Differential Expression of Stimulation Genes\n",
          "with a Consistent Drug Effect")))
squarePlot(plot_corr(result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"], 
          "log2FoldChange.6", "log2FoldChange.1", 
          paste("Relative Differential Expression of Stimulation Genes\n",
          "with a Drug Effect Altered By Stimulation")))
dev.off()

# extract genesets
# #1 vs #8 plot
one_v_8_q1 = result[pvalue.8 < 0.05 & pvalue.10 > 0.05
                    ][log2FoldChange.8 > 0 & log2FoldChange.1 > 0, Gene_id]
one_v_8_q2 = result[pvalue.8 < 0.05 & pvalue.10 > 0.05
                    ][log2FoldChange.8 < 0 & log2FoldChange.1 > 0, Gene_id]
one_v_8_q3 = result[pvalue.8 < 0.05 & pvalue.10 > 0.05
                    ][log2FoldChange.8 < 0 & log2FoldChange.1 < 0, Gene_id]
one_v_8_q4 = result[pvalue.8 < 0.05 & pvalue.10 > 0.05
                    ][log2FoldChange.8 > 0 & log2FoldChange.1 < 0, Gene_id]
# #1 vs #9 plot
one_v_9_q1 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.9 > 0 & log2FoldChange.1 > 0, Gene_id]
one_v_9_q2 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.9 < 0 & log2FoldChange.1 > 0, Gene_id]
one_v_9_q3 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.9 < 0 & log2FoldChange.1 < 0, Gene_id]
one_v_9_q4 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.9 > 0 & log2FoldChange.1 < 0, Gene_id]

# #1 vs #6 plot
one_v_6_q1 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.6 > 0 & log2FoldChange.1 > 0, Gene_id]
one_v_6_q2 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.6 < 0 & log2FoldChange.1 > 0, Gene_id]
one_v_6_q3 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.6 < 0 & log2FoldChange.1 < 0, Gene_id]
one_v_6_q4 = result[pvalue.10 < 0.05 & GeneSymbol != "Gm4294"
                    ][log2FoldChange.6 > 0 & log2FoldChange.1 < 0, Gene_id]

q_list = c("one_v_8_q1", "one_v_8_q2", "one_v_8_q3", "one_v_8_q4", "one_v_9_q1",
           "one_v_9_q2", "one_v_9_q3", "one_v_9_q4", "one_v_6_q1", "one_v_6_q2",
           "one_v_6_q3", "one_v_6_q4")

# write lists to file
lapply(q_list, function(i){
  f1 = paste0(today, nameset, i, ".txt")
  fwrite(list(get(i)), quote=F, col.names=F, file=f1)
})

# Pathway analyses
# GProfiler2
require(gprofiler2)
# get background of all genes
bckgrd = result$Gene_id
# run enrichments in a loop
lapply(q_list, function(i){
  # lookup using pval and expressed genes only background
  g = gost(query=get(i), organism="mmusculus")
  if(is.list(g$result)){
    fwrite(g$result, paste0(today, nameset, i, ".Gprofiler.csv"))
  }
})





