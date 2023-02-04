#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(biomaRt)
library(magrittr)
library(BiocParallel)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(PoiClaClu)
library(jsonlite)


# filepaths
rootdir = '/labs/flongo/t41b_BD10-2_stim'
wk_dir = paste0(rootdir, '/.broad_analysis')
setwd(paste0(rootdir, '/.broad_analysis'))

# file naming info
today = format(Sys.Date(), '%Y%m%d')
nameset = '_APP_Activity_Dep_'

# biomaRt database to draw from
ensembl = useEnsembl("ensembl", dataset="mmusculus_gene_ensembl", version="108")

# Load sample data
sampleTable = fread(paste0(rootdir, "/metadata_samples_relabel.csv"), header=T)

# what annotation columns are the sample names and file names coming from?
samp_names = "Sample_ID"
# secondary name
secon_name = "Group"
# don't forget to define pca functions at the bottom


#################################
# Building dds from salmon counts
# reading in counts
load(paste0(rootdir, "/T41B_BD10_2_RSEM_Quant.genes.counts.RData"))
counts = as.data.frame(counts)

# drop outlier samples from initial pca/clustering
exclude_list = c("BD10_2_30")
sampleTable = sampleTable[!(sampleTable$V1 %in% exclude_list) ]
rownames(sampleTable) = sampleTable$V1
counts = round(counts[, !(colnames(counts) %in% exclude_list) ])

# reorder metadata rows to match counts, tpm for analysis
counts = counts[, sampleTable$V1]

# checking if counts and metadata are set up right
if (!all(names(counts) == rownames(sampleTable))) {
  print("counts file and sampleTable do not match!")
  quit("no", 1)
}

dds = DESeqDataSetFromMatrix(counts, sampleTable, ~ Subgroup)


##########################
# Prep and Normalization

# relevel factors to base value
dds = dds[, dds$Genotype == "TRANS" ]
dds$MouseID = factor(dds$MouseID)
dds$Subgroup = droplevels(dds$Subgroup)

# set model
# design(dds) = ~ MouseID + Subgroup


# gene info for filtering by gene biotype
if(!file.exists("gene_info.txt")){
  gene_info = getBM(attributes=c("ensembl_gene_id", "external_gene_name",
                                 "gene_biotype", "chromosome_name", 
                                 "entrezgene_id"),
                    filters="ensembl_gene_id", values=rownames(dds),
                    mart=ensembl)
  write.table(gene_info, file="gene_info.txt", sep="\t", quote=F, row.names=F)
}
gene_info = fread("gene_info.txt")

# ensembl transcript types to keep, excluding pseudogenes and TECs
keep_categ = grep("pseudo", unique(gene_info$gene_biotype), value=T, invert=T)
keep_categ = grep("TEC", keep_categ, value=T, invert=T)

# filter low count genes (n-1 of smallest group)
print(paste(nrow(dds), "unfiltered genes"))
# n-1 of smallest group (minimum of 3)
smallest_n = max(min(table(colData(dds)$Subgroup)) - 1, 3)
# filter low count genes to at least 10 in n-1 of smallest group
dds = dds[ rowSums( counts(dds) >= 1 ) >= smallest_n, ]
print(paste(nrow(dds), "filtered genes, n-1 of smallest group", smallest_n))

# just estimateSizeFactors, don't need to run DE analysis
dds = estimateSizeFactors(dds)

# filter by gene type
use_genes = gene_info$ensembl_gene_id[ gene_info$gene_biotype %in% keep_categ ]
dds = dds[ rownames(dds) %in% use_genes ]

# use only canonicanl chromosomes
use_chrs = gene_info$ensembl_gene_id[ 
  grep("^\\d{1,2}|X|Y|MT" , gene_info$chromosome_name, perl=T) ]
dds = dds[ rownames(dds) %in% use_chrs ]

# add gene names annotation
mcols(dds) = DataFrame(mcols(dds), 
        symbol=gene_info$external_gene_name[ match(rownames(dds), 
                                                   gene_info$ensembl_gene_id) ],
        entrez_id=gene_info$entrezgene_id[ match(rownames(dds), 
                                                 gene_info$ensembl_gene_id) ])

# vst transform
vst = vst(dds)

# checking normalization
# this gives log2(n + 1) transformation
ntd = normTransform(dds)
# plot variances for both transformations
plot_1 = meanSdPlot(assay(ntd))
plot_2 = meanSdPlot(assay(vst))
ml = marrangeGrob(list(plot_1$gg, plot_2$gg), nrow=1, ncol=2, respect=T)
ggsave(paste0(today, nameset, "normalization.pdf"), ml, width=11, height=8.5, 
       useDingbats=F)

remove(ntd, ml, plot_1, plot_2, ensembl, keep_categ, use_chrs, use_genes)


#############################
# Poisson Distance Pairwise
poisd = PoissonDistance(t(counts(dds)))
samplePoisDistMatrix = as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) = paste(dds[[samp_names]], dds[[secon_name]], 
                                       sep=" - ")
colnames(samplePoisDistMatrix) = paste(dds[[samp_names]], dds[[secon_name]], 
                                       sep=" - ")

# set colors and draw
colors = colorRampPalette(brewer.pal(11, "RdYlBu"))(255)
pheatmap(samplePoisDistMatrix, col=colors, clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd, width=11, height=8.5,
         filename=paste0(today, nameset, "poisson.pdf"))

remove(poisd, samplePoisDistMatrix)


###############################
# PCA top X genes by variance
# c_by = column to color by
# c_lab = Label for color legend
# c_type = "disc" discrete or "cont" continuous coloring
# num = number of genes, default is all
my_pca_plot = function(vst, c_by, c_lab, c_type="disc", num=1000000) {
  rv = rowVars(assay(vst))
  top_genes = order(rv, decreasing=T)[seq_len(min(num, length(rv)))]
  pca = prcomp(t(assay(vst)[top_genes, ]))
  percent_pc = round((pca$sdev^2 / sum(pca$sdev^2))*100, 2)
  if (c_type == "disc") {
    col_by = as.factor(vst[[c_by]])
  } else if (c_type == "cont") {
    col_by = vst[[c_by]]
  } else {
    print("Incorrect value for c_type!")
    quit("no", 1)
  }

  # PC1 v PC2
  ggPCA12 = ggplot(data.frame(pca$x), aes(x=PC1, y=PC2, col=col_by)) +
    labs(x=sprintf("PC1 - %.2f%%", percent_pc[1]),
         y=sprintf("PC2 - %.2f%%", percent_pc[2]), col=c_lab) +
    geom_point(size=5, alpha=0.5) + stat_ellipse() +
    # geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
    ggtitle(sprintf("top %s genes", length(top_genes))) + theme_classic() + {
      if(c_type == "disc") scale_color_brewer(type="qual", palette="Paired") 
      else scale_color_distiller(type="div", palette="RdYlBu")
    }

  # PC1 v PC3
  ggPCA13 = ggplot(data.frame(pca$x), aes(x=PC1, y=PC3, col=col_by)) +
    labs(x=sprintf("PC1 - %.2f%%", percent_pc[1]),
         y=sprintf("PC3 - %.2f%%", percent_pc[3]), col=c_lab) +
    geom_point(size=5, alpha=0.5) + stat_ellipse() +
    # geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
    ggtitle(sprintf("top %s genes", length(top_genes))) + theme_classic() + {
      if(c_type == "disc") scale_color_brewer(type="qual", palette="Paired") 
      else scale_color_distiller(type="div", palette="RdYlBu")
    }
  plot(ggPCA12)
  plot(ggPCA13)
}


#################################
# PCA plotting 
pdf(file=paste0(today, nameset, "pca.pdf"), onefile=T, paper="USr", width=11, 
    height=8.5)
my_pca_plot(vst, c_by="Group", c_lab="Group", c_type="disc")
my_pca_plot(vst, c_by="Stimulation", c_lab="Stimulation", c_type="disc")
my_pca_plot(vst, c_by="Treatment", c_lab="Treatment", c_type="disc")
my_pca_plot(vst, c_by="Stimulation", c_lab="Top 500 Stimulation", c_type="disc", num=500)
my_pca_plot(vst, c_by="Treatment", c_lab="Top 500 Treatment", c_type="disc", num=500)
my_pca_plot(vst, c_by="Group", c_lab="Top 500 Group", c_type="disc", 
            num=500)
my_pca_plot(vst, c_by="Stimulation", c_lab="Top 1000 Stimulation", c_type="disc", num=1000)
my_pca_plot(vst, c_by="Treatment", c_lab="Top 1000 Treatment", c_type="disc", num=1000)
my_pca_plot(vst, c_by="Group", c_lab="Top 1000 Group", c_type="disc", 
            num=1000)
dev.off()


# save DESeqDataSet for future use
save(dds, file=paste0(today, nameset, 'dds.Rdata'), compress=T)

