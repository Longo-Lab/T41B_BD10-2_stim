#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(data.table)

# directory locations
rootdir='/labs/flongo/2020_T41B_BD10-2_Stimulation'
setwd(paste0(rootdir, '/enrichment_plots'))
today = format(Sys.Date(), '%Y%m%d')
nameset = "_enrichment_"
# sensi_dir = paste0(rootdir, '')

####FUNCTIONS
# correlation plot function log2 enrichment ratio vs -log10
plot_corr = function (table, x, y, col, lbl, reg, title=NULL, ...) {
  gradcol = colorRampPalette(brewer.pal(9, "BuPu"), bias=1.5)(255)
  ggscatter(table, x=x, y=y, color=col, size=col, shape=reg, label=lbl, repel=T, ...) + 
    gradient_color(gradcol[30:255]) + scale_shape_manual(values=c(25, 24))+
    ggtitle(title)
}

####MAIN
# read in files
# files = list.files(path=rootdir, pattern="gost.*.csv", full.names=T)
# f_names = gsub(".csv", "", list.files(path=rootdir, pattern="gost.*.csv"))
# fl = lapply(files, fread)
# setattr(fl, 'names', f_names)
# fdata = rbindlist(fl, use.names=T, idcol="group")
fdata = fread("selected_terms.csv")
fdata = fdata[, -c("V1", "intersection", "Gene_name", "query", "evidence_codes")]

# calculate plot values
fdata[, expected := (query_size * (term_size / effective_domain_size))]
fdata[, L2ER := log2(intersection_size / expected)]
fdata[, L10FDR := -log10(p_value)]

# filter genes 
fdata = fdata[term_size <= 1000 & significant == "TRUE" ]

# plot
options(ggrepel.max.overlaps = Inf)
pdf(file=paste0(today, nameset, "plots.pdf"), onefile=T, paper="USr",
    width=11, height=8.5, useDingbats=F )
gname = "T41B + BD10-2 vs T41B"
labs = fdata[ group == gname ][ L10FDR >= 7 | L2ER >= 4 ][, term_name ]
plot_corr(fdata[ group == gname ], x="L2ER", y="L10FDR", 
          col="intersection_size", lbl="term_name", reg="DGE", title=gname, 
          label.select=labs)
dev.off()

