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

#######################################MODIFY to run on theses subsets
a = fread("20220922.T41B_BD10-2_stim.Wan_etal.Modules_Up-Down.full.txt")
# neuronal APPL/S BD10-2 unique genes
down = c("APPL/S Down", "APPL/S_BD10-2 Down")
neu = c("CBEyellow", "DLPFCyellow", "FPyellow", "IFGbrown", "PHGbrown", "STGbrown", "TCXgreen")
b = a[sample %in% down & group %in% neu, .(sample, group, symbol)]
# split string and explode 
b[, c("symbol") := lapply(.SD, strsplit, split=","), .SDcols="symbol"]
c = b[, lapply(.SD, unlist), by=1:nrow(b)]
c[, nrow := NULL]
# count up occurences of each gene in APPL/S_BD10-2 but not APPL/S
f <- function(sam,sym) setdiff(sym[sam!="APPL/S DOWN"], sym[sam=="APPL/S DOWN"])
d = c[,.(symbol = f(sample,symbol)),group][, .(N_diff = uniqueN(group)),symbol]
g = gost(query=d[N_diff >= 6, symbol], organism="mmusculus", ordered_query=T, 
         evcodes=T, significant=T)

# now do APPL/S BD10-2 only oligo genes
up = c("APPL/S Up", "APPL/S_BD10-2 Up")
olig = c("CBEbrown", "DLPFCbrown", "FPblue", "IFGblue", "PHGgreen", "STGyellow", "TCXyellow")
b = a[sample %in% up & group %in% olig, .(sample, group, symbol)]
# split string and explode 
b[, c("symbol") := lapply(.SD, strsplit, split=","), .SDcols="symbol"]
c = b[, lapply(.SD, unlist), by=1:nrow(b)]
c[, nrow := NULL]
# count up occurences of each gene in APPL/S_BD10-2 but not APPL/S
f <- function(sam,sym) setdiff(sym[sam!="APPL/S UP"], sym[sam=="APPL/S UP"])
d = c[,.(symbol = f(sample,symbol)),group][, .(N_diff = uniqueN(group)),symbol]
g = gost(query=d[N_diff >= 6, symbol], organism="mmusculus", ordered_query=T, 
         evcodes=T, significant=T)




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

