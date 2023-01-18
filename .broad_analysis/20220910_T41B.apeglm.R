#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)
library(ggplot2)
library(EnhancedVolcano)
library(gprofiler2)
library(stringr)
library(magrittr)
# library(ReactomePA)
# library(clusterProfiler)

# filepaths
wk_dir = '/labs/flongo/t41b_BD10-2_stim/broad_analysis'
data_dir = '/labs/flongo/t41b_BD10-2_stim/broad_analysis'
rdata = '20220910_LTP_broad_dds.Rdata'
# file naming info
today = format(Sys.Date(), '%Y%m%d')
nameset = 'T41B_LTP_sep'
# lfcThreshold choice
lfc = log2(1.10)
# pvalue, padj or svalue
p = 'svalue'
# cutoff threshold
thresh = 0.05
# order ranked query by: log2FoldChange, P ORDERING NOT IMPLEMENTED CURRENTLY
# ordering by l2fc means a cutoff of lfc above will automatically be applied
ord_by = "log2FoldChange" 
# for enrichment, add an additional cutoff of p above?
cut_p = TRUE 

#####################################################################FUNCTIONS
my_volcano = function(x, titl=NULL) {
  # Make a volcano with pval cutoff
  # x = results DESeq Data set
  # titl = Label for graph
  require("ggplot2")
  require("EnhancedVolcano")
  # limit to pval subset
  x <- x[ !is.na(x[[p]]), ]
  # color by up/downregulated
  keyvals = ifelse(x$log2FoldChange < -lfc & x[[p]] < thresh, 'royalblue2',
    ifelse(x$log2FoldChange > lfc & x[[p]] < thresh, 'darkorange2', 'grey45'))
  names(keyvals)[keyvals == 'darkorange2'] <- 'upregulated'
  names(keyvals)[keyvals == 'grey45'] <- 'not significant'
  names(keyvals)[keyvals == 'royalblue2'] <- 'downregulated'
  # plot with EnhancedVolcano
  p1 = EnhancedVolcano(x,
                       lab = x$symbol,
                       x = 'log2FoldChange',
                       y = p,
                       # xlim = c(-9,9),
                       xlab = bquote(~Log[2]~ 'fold change'),
                       ylab = bquote('-' ~Log[10]~ .(p)),
                       title = titl,
                       pCutoff = thresh,
                       FCcutoff = 0,
                       pointSize = 1.5,
                       labSize = 3.5,
                       colCustom = keyvals,
                       colAlpha = 1,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4,
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE,
                       border = 'full',
                       borderWidth = 1,
                       borderColour = 'black')

  p1 + scale_y_continuous(trans = "sqrt")
}

get_mcols = function(str, dds, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # dds: DESeq Dataset with annotations in mcol
  # col: colname to grab
  require("stringr")
  # x = str_split(str, ",")
  y = mcols(dds)[[col]][match(unlist(str_split(str, ",")), rownames(dds))]
  return(paste0(y, collapse=","))
}


##########################################################################MAIN
# Load and drop outliers
require("DESeq2")

# Deseq dataset to run from
load(file.path(wk_dir, rdata))

# # # drop outlier samples from pca/clustering
# # exclude_list = c("308", "332")
# # dds = dds[, !(colnames(dds) %in% exclude_list) ]

##############################################################################
# DESeq2 analysis
# set up naming for easier metadata
dds$group = str_replace_all(dds$Group, c("TRANS"="Tg", "VEH"="V", "DRUG"="D", "_STIM"=""))
dds$group = factor(dds$group, levels=c("WT_V", "WT_D", "Tg_V", "Tg_D"))
dds$Genotype = str_replace_all(dds$Genotype, c("TRANS"="Tg"))
dds$Genotype = factor(dds$Genotype, levels=c("WT", "Tg"))
dds$Treatment = str_replace_all(dds$Treatment, c("VEH"="V", "DRUG"="D"))
dds$Treatment = factor(dds$Treatment, levels=c("V", "D"))
dds$BatchID = factor(dds$BatchID)

# # # Using pairwise group comparisons ------------------------------
# # # design options
# # # # mm = model.matrix(~ MouseID + group, colData(dds))
# # # # mm[duplicated(mm)|duplicated(mm, fromLast=TRUE), ]
# # design(dds) = ~ group # run with dds1-dds4 separate comparisons below

# # # run deseq
# # dds = DESeq(dds)
# # # extract results
# # resultsNames(dds)

# # # result df names for lapply
# # res_list = c("dds1", "dds2", "dds3", "dds4")
# # # matched title list for plots (same order as above)
# # titl_list = c("Genotype Effect", "Treatment in Wt", 
              # # "Treatment in Tg", "Treatment + Genotype (TgD vs WtV)")

# # # Wald test for lfc thresh significance
# # res.dds1 = results(dds,contrast=c("group", "Tg_V","WT_V"), lfcThreshold=lfc)
# # res.dds2 = results(dds,contrast=c("group", "WT_D","WT_V"), lfcThreshold=lfc)
# # res.dds3 = results(dds,contrast=c("group", "Tg_D","Tg_V"), lfcThreshold=lfc)
# # res.dds4 = results(dds,contrast=c("group", "Tg_D","WT_V"), lfcThreshold=lfc)


# # # lfc shrink for accurate rank ordering
# # res.dds1.lfc = lfcShrink(dds,coef="group_Tg_V_vs_WT_V", lfcThreshold=lfc)
# # res.dds2.lfc = lfcShrink(dds,coef="group_WT_D_vs_WT_V", lfcThreshold=lfc)
# # res.dds4.lfc = lfcShrink(dds,coef="group_Tg_D_vs_WT_V", lfcThreshold=lfc)

# # # relevel for coef comparison to use apeglm
# # dds$group %<>% relevel("Tg_V","WT_V","WT_D","Tg_D")
# # dds = nbinomWaldTest(dds)
# # res.dds3.lfc = lfcShrink(dds,coef="group_Tg_D_vs_Tg_V", lfcThreshold=lfc)

# # # change to interaction model -----------------------------------
# # design(dds) = ~ Genotype + Treatment + Genotype:Treatment # interaction model
# # # cannot add BatchID or MouseID b/c of rank overlap with groups

# # # run deseq
# # dds = DESeq(dds)
# # # extract results
# # resultsNames(dds)

# # # result df names for lapply
# res_list = c("wt_drug", "app_drug", "int_drug", "genotype")
# # matched title list for plots (same order as above)
# titl_list = c("WT Drug Effect", "App Drug Effect", 
              # "Drug Effect Interaction Term", "Genotype Effect")

# # # extract DE genes
# # # the drug effect for genotype WT(the main effect)
# # res.wt_drug = results(dds, name="Treatment_D_vs_V", lfcThreshold=lfc)
# # res.wt_drug.lfc = lfcShrink(dds, coef="Treatment_D_vs_V", lfcThreshold=lfc)

# # # the interaction term, answering: is the drug effect *different* across genotypes?
# # res.int_drug = results(dds, name="GenotypeTg.TreatmentD")
# # res.int_drug.lfc = lfcShrink(dds, coef="GenotypeTg.TreatmentD", lfcThreshold=lfc)

# # # the genotype effect (in vehicle)
# # res.genotype = results(dds, name="Genotype_Tg_vs_WT")
# # res.genotype.lfc = lfcShrink(dds, coef="Genotype_Tg_vs_WT", lfcThreshold=lfc)

# # # the drug effect for genotype Tg
# # # relevel genotype for correct coef comparison to use apeglm
# # dds$Genotype %<>% relevel("Tg","WT")
# # dds = nbinomWaldTest(dds)
# # res.app_drug = results(dds, name="Treatment_D_vs_V", lfcThreshold=lfc)
# # res.app_drug.lfc = lfcShrink(dds, coef="Treatment_D_vs_V", lfcThreshold=lfc)

# Separating pairwise group comparisons ------------------------------
# this attempts to correct for significant dispersion difference between 
# Genotype comparisons and the Treatment comparisons
# design options
# # mm = model.matrix(~ MouseID + group, colData(dds))
# # mm[duplicated(mm)|duplicated(mm, fromLast=TRUE), ]
design(dds) = ~ group # run with dds1-dds4 separate comparisons below

# result df names for lapply
res_list = c("dds1", "dds2", "dds3", "dds4")
# matched title list for plots (same order as above)
titl_list = c("Genotype Effect", "Treatment in Wt", 
              "Treatment in Tg", "Treatment + Genotype (TgD vs WtV)")

# extract results
# resultsNames(dds)

# run deseq separately for each
dds1 = dds[, dds$group %in% c("Tg_V","WT_V")]
dds1$group = droplevels(dds1$group)
dds1 = DESeq(dds1)
# Wald test for lfc thresh significance
res.dds1 = results(dds1, name="group_Tg_V_vs_WT_V", lfcThreshold=lfc)
# lfc shrink for accurate rank ordering
res.dds1.lfc = lfcShrink(dds1, coef="group_Tg_V_vs_WT_V", lfcThreshold=lfc)

# run deseq separately for each
dds2 = dds[, dds$group %in% c("WT_D","WT_V")]
dds2$group = droplevels(dds2$group)
dds2 = DESeq(dds2)
# Wald test for lfc thresh significance
res.dds2 = results(dds2, name="group_WT_D_vs_WT_V", lfcThreshold=lfc)
# lfc shrink for accurate rank ordering
res.dds2.lfc = lfcShrink(dds2, coef="group_WT_D_vs_WT_V", lfcThreshold=lfc)

# run deseq separately for each
dds3 = dds[, dds$group %in% c("Tg_D","Tg_V")]
dds3$group = droplevels(dds3$group)
dds3 = DESeq(dds3)
# Wald test for lfc thresh significance
res.dds3 = results(dds3, name="group_Tg_D_vs_Tg_V", lfcThreshold=lfc)
# lfc shrink for accurate rank ordering
res.dds3.lfc = lfcShrink(dds3, coef="group_Tg_D_vs_Tg_V", lfcThreshold=lfc)

# run deseq separately for each
dds4 = dds[, dds$group %in% c("Tg_D","WT_V")]
dds4$group = droplevels(dds4$group)
dds4 = DESeq(dds4)
# Wald test for lfc thresh significance
res.dds4 = results(dds4, name="group_Tg_D_vs_WT_V", lfcThreshold=lfc)
# lfc shrink for accurate rank ordering
res.dds4.lfc = lfcShrink(dds4, coef="group_Tg_D_vs_WT_V", lfcThreshold=lfc)


##############################################################################
# MA plots for initial analysis
pdf(paste0(wk_dir, "/", today, ".MA_plots.", nameset, ".pdf"))
# visualize results
drawLines = function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
lapply(seq(res_list), function(i) {
    plotMA(get(paste0("res.", res_list[i])), ylim=c(-5,5),
           main=paste(titl_list[i], "-", res_list[i]))
    drawLines()
})
# visualize lfcshrink threshold
lapply(seq(res_list), function(i) {
    plotMA(get(paste0("res.", res_list[i], ".lfc")), ylim=c(-5,5), alpha=0.1, 
           main=paste(titl_list[i], "- apeglm -", res_list[i]))
    drawLines()
})
dev.off()

##############################################################################
# annotate pval and padj and from dds for results
lapply(res_list, function(i){
  x = get(paste0("res.", i, ".lfc"))
  y = get(paste0("res.", i))
  x$pvalue = y$pvalue[match(rownames(y), rownames(x))]
  x$padj = y$padj[match(rownames(y), rownames(x))]
  x$raw_L2FC = y$log2FoldChange[match(rownames(y), rownames(x))]
  assign(paste0("res.", i, ".lfc"), x, inherits=T)
})

# annotate symbol and  entrez id from dds for l2fc results
lapply(res_list, function(i){
  x = get(paste0("res.", i, ".lfc"))
  x$symbol = mcols(dds)$symbol[match(rownames(x), rownames(dds))]
  x$entrez_id = mcols(dds)$entrez_id[match(rownames(x), rownames(dds))]
  assign(paste0("res.", i, ".lfc"), x, inherits=T)
})


##############################################################################
# write results to file (either l2fc or regular)
lapply(res_list, function(i) {
  write.csv(get(paste0("res.", i, ".lfc")), 
            paste0(wk_dir, "/res.", i, ".lfc.txt"))
  # write.csv(get(paste0("res.", i)), paste0(wk_dir, "/res.", i, ".txt"))
})


##############################################################################
# plot each volcano with function above
pdf(paste0(wk_dir, "/", today, ".EnhancedVolcano.", nameset, ".pdf"))
lapply(seq(length(res_list)), function(i){
  my_volcano(get(paste0("res.", res_list[i], ".lfc")), titl=titl_list[i])
})
graphics.off()


##############################################################################
# Pathway analyses
# GProfiler2
require(gprofiler2)
lapply(res_list, function(i){
  # get dataframe based on keyword list
  x = get(paste0("res.", i, ".lfc"))
  # drop NA pvalues
  x = x[ !is.na(x[[ord_by]]), ]
  # sort by pvalue
  x = x[with(x, order(-abs(get(ord_by)))), ]
  # define subset
  if(cut_p){
    q = rownames(x)[abs(x[[ord_by]]) > lfc & x[[p]] < thresh]
  } else {
    q = rownames(x)[abs(x[[ord_by]]) > lfc]
  }
  # lookup using ordered subset
  g = gost(query=q, organism="mmusculus", ordered_query=T, evcodes=T)
  # write to file
  if(is.list(g$result)){
    # get symbols for terms
    g$result$symbols = lapply(g$result$intersection, get_mcols, dds, "symbol")
    # write to file
    fwrite(g$result, paste0(wk_dir, "/res.", i, ".lfc.Gprofiler.csv"))
  }
})


