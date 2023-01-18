#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


# directory locations
rootdir = '/oak/stanford/scg/lab_flongo/2020_T41B_BD10-2_Stimulation'
wk_dir = paste0(rootdir, '/updown_ontology_plots')
setwd(paste0(rootdir, '/updown_ontology_plots'))

# file naming info
today = format(Sys.Date(), '%Y%m%d')
nameset = '_ontology_'

# get file info
tiss = data.table("short"=c("2", "9", "12"), 
                  "full"=c("Genotype (Tg_V vs Wt_V)", 
                           "Drug (Tg_D vs Tg_V)", 
                           "Genotype + Drug (Tg_D vs Wt_V)")) 
# "10", "Drug Effect Interaction", 
cats = c("CORUM","GO:BP","GO:CC","GO:MF","HP","KEGG","MIRNA","REAC","TF","WP")
target_files = unlist(lapply(tiss$short, function(i){
  c(paste0(rootdir, paste0("/gost.", i, ".csv")))
}))

# get list of terms to plot from file
sht_lst = fread('20210405_terms.txt', sep="\t", header=F)

# biomaRt database to draw from
ensembl = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version="103")


#####################################################################FUNCTIONS
# string wrap function for long text
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

# make genlist table for plot (for a given term)
gen_tbl = function(term, tbl) {
  x = tbl[Title == term, intersection]
  y = data.table("gene"=unique(unlist(strsplit(x, ",", fixed=T))))
  y$symbol = gene_info$symbol[match(y$gene, gene_info$gene)]
  return(y)
}

# make heatmap plot (vst(norm) or tpm)
mk_heat = function(term, tbl, heat, scl="none") {
  genlist = gen_tbl(term, tbl)
  # brewer colors
  colors = colorRampPalette( brewer.pal(9, "RdYlBu") )(255)
  # setting up gene matrix
  mat = heat[as.character(genlist$gene) , ]
  rownames(mat) = genlist$symbol
  # regular heatmap
  gg = pheatmap(mat, color=rev(colors), show_rownames=T, cluster_cols=F, 
                main=wrapper(term, width=29), scale=scl)
  return(gg[[4]])
}

# test_tbl = gen_tbl("GO:0005930 axoneme", tbl=fd)
# test = mk_heat("GO:0005930 axoneme", tbl=fd, heat=mean_tpm)
# plot(test)

# make up/down table for plot (for a given term negative downreg LogBH)
plt_tbl = function(term, tbl) {
  x = tbl[Title == term]
  x[DGE == "down", LogBH := -(LogBH)]
  x[, Tiss := tiss$full[match(Tiss, tiss$short)] ]
  return(x)
}

# make plot given a list of terms
mk_plt = function(term, tbl) {
  a = plt_tbl(term, tbl)
  lim = max(abs(a$LogBH))
  gg = ggplot(a, aes(x=Tiss, y=LogBH, fill=DGE)) + ylim(-lim, lim) +
    geom_bar(stat = "identity", position="identity") + theme_classic() +
    geom_text(aes(label=intersection_size, vjust=1*(sign(LogBH))), color="white", size=5) +
    ggtitle(wrapper(unique(a$Title), width=29)) + scale_x_discrete(limits=tiss$full) + 
    geom_hline(yintercept=c(-log10(0.05),log10(0.05)), linetype="dashed", col='gray') +
    geom_hline(yintercept=0, col='black') +
    scale_fill_manual(values=c("#762a83","#1c793d")) +
    theme(legend.position="none", 
          text=element_text(size=20),
          plot.title=element_text(size=16, face="bold", hjust=0.5),
          axis.title.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.line.x=element_line(color="white"),
          axis.text.x=element_text(angle=45, hjust=1)) + 
    ylab(expression("down" %<-% "-Log10 FDR" %->% "up"))
  return(gg)
}

# test = mk_plt("hsa04730:Long-term depression", fd)  
# plot(test)

grid_pdf = function(plist, nrow=3, ncol=4, suffix="plots.pdf") {
  # Function to generate gridded multi-page pdf with multiple ggplots in a list
  # plist: list of ggobjects
  # nrow : number of rows per page
  # ncol: number of columns per page
  # suffix: filename suffix
  # returns none
  ml = marrangeGrob(plist, nrow=nrow, ncol=ncol, respect=F, top=NULL, 
                    padding=unit(15, "mm"), 
                    layout_matrix=t(matrix(1:(nrow*ncol), nrow=ncol, ncol=nrow)))
  ggsave(paste0(today, nameset, suffix), ml, width=420, height=594, unit="mm",
         useDingbats=F)
}


##########################################################################MAIN
## Load data
# read in files
all(file.exists(target_files))
# reading in tables, combine and melt by GWAS
fl = lapply(target_files, fread)
setattr(fl, 'names', tiss$short)
fd = rbindlist(fl, use.names=T, idcol="Tiss")

# clean up for parsing
fd[, c("V1", "evidence_codes") := NULL]
fd[, DGE := gsub("regulated", "", DGE)]
fd[, LogBH := -log10(p_value)]
fd[, Title := paste(term_id, term_name)]

# gene info for gene symbols
if(!file.exists("gene_info.txt")){
  all = data.table("gene"=unique(unlist(strsplit(as.character(fd$intersection),
                                                 ",", fixed=T))))
  gene_info = getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                    filters="ensembl_gene_id", mart=ensembl, values=all$gene)
  names(gene_info) = c("gene", "symbol")
  gene_info = merge(all, gene_info, all.x=T, by="gene")
  gene_info[is.na(symbol), symbol := gene]
  write.table(gene_info, file="gene_info.txt", sep="\t", quote=F, row.names=F)
}
gene_info = fread("gene_info.txt")


## Expression Data
# reading in counts
metadata = fread(paste0(rootdir, "/metadata_samples_relabel.csv"), header=T)
load(paste0(rootdir, "/T41B_BD10_2_RSEM_Quant.genes.counts.RData"))
counts = as.data.frame(counts)
load(paste0(rootdir, "/T41B_BD10_2_RSEM_Quant.genes.tpm.RData"))
tpm = as.data.frame(tpm)


# drop outlier samples from pca/clustering
exclude_list = c("BD10_2_30")
metadata = metadata[!(metadata$V1 %in% exclude_list) ]
rownames(metadata) = metadata$V1
counts = round(counts[, !(colnames(counts) %in% exclude_list) ])
tpm = tpm[, !(colnames(tpm) %in% exclude_list) ]

# reorder metadata rows to match counts, tpm for analysis
counts = counts[, metadata$V1]
tpm = tpm[, metadata$V1]

## TPM analysis
tpm = log2(tpm + 1)
mean_tpm = data.frame("Wt_V"=rowMeans(tpm[,metadata$Subgroup == "A_STIM"]),
                      "Wt_D"=rowMeans(tpm[,metadata$Subgroup == "B_STIM"]),
                      "Tg_V"=rowMeans(tpm[,metadata$Subgroup == "C_STIM"]),
                      "Tg_D"=rowMeans(tpm[,metadata$Subgroup == "D_STIM"]))

## DESeq2 analysis
# checking if counts and metadata are set up right
if (!all(names(counts) == rownames(metadata))) {
  print("counts file and sampleTable do not match!")
  quit("no", 1)
}

#import to DESeq
dds = DESeqDataSetFromMatrix(counts, metadata, ~ Subgroup)

# just estimateSizeFactors, don't need to run DE analysis
dds = estimateSizeFactors(dds)

# vst transform
vst = vst(dds)

# generate mean values per group
mean_vst = data.frame("Wt_V"=rowMeans(assay(vst)[,colData(vst)$Subgroup == "A_STIM"]),
                      "Wt_D"=rowMeans(assay(vst)[,colData(vst)$Subgroup == "B_STIM"]),
                      "Tg_V"=rowMeans(assay(vst)[,colData(vst)$Subgroup == "C_STIM"]),
                      "Tg_D"=rowMeans(assay(vst)[,colData(vst)$Subgroup == "D_STIM"]))


## Generating plots
# lookup terms w/ term_ids
match_terms = unique(fd$Title[fd$term_name %in% sht_lst$V1])

# plot barplots for relevant tissues
plist = lapply(match_terms, mk_plt, tbl=fd)
grid_pdf(plist, suffix="enrich_bars.pdf")

# # plot tpm heatmap
# tlist = lapply(match_terms, mk_heat, tbl=fd, heat=mean_tpm, scl="row")
# grid_pdf(tlist, suffix="enrich_tpm_rs.pdf")

# plot vst heatmap
tlist = lapply(match_terms, mk_heat, tbl=fd, heat=mean_vst)
grid_pdf(tlist, suffix="enrich_vst.pdf")

# plot rowscaled vst heatmap
vlist = lapply(match_terms, mk_heat, tbl=fd, heat=mean_vst, scl="row")
grid_pdf(vlist, suffix="enrich_vst_rs.pdf")

# one plot 3 x 3
nlist = list()
for (i in seq(plist)) {
  nlist = c(nlist, plist[i])
  nlist = c(nlist, vlist[i])
  nlist = c(nlist, tlist[i])
}
grid_pdf(nlist, nrow=3, ncol=3, suffix="enrich_all.pdf")
