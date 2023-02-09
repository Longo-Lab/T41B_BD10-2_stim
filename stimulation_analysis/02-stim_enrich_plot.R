#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-stim_enrich_plot.R
##
## Version: 0.0.1
##
## Purpose of script: plot top ontology enrichments across three stim conditions
##
## Author: Robert R Butler III
##
## Date Created: 2023-02-08
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Using padj < 0.05 enrichment from gprofiler, limiting to term size < 1000
##
##   Usage:
##     sbatch -J gos --mem=32G -c 4 -t 01:00:00 -p interactive -A default \
##       -o %A_stim_enrich_plot.log \
##       --wrap "ml R/4.0; Rscript 02-stim_enrich_plot.R"
##
##   or interactive session:
##     sdev -m 32 -c 4 -t 01:00:00 -p interactive -a default
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(data.table)
# library(biomaRt)
# library(magrittr)
library(ggplot2)
# library(gridExtra)
# library(DESeq2)
# library(pheatmap)
# library(RColorBrewer)


# filepaths
rootdir = '/labs/flongo/t41b_BD10-2_stim'
wk_dir = paste0(rootdir, '/stimulation_analysis')
setwd(wk_dir)

# file naming info
today = format(Sys.Date(), '%Y%m%d')
nameset = 'APP_Activity_Dep'

# get file info
tiss = data.table("short"=c("5", "6", "7"), 
                  "full"=c("Wt", "APPL/S", "APPL/S + BD10-2")) 
target_files = unlist(lapply(tiss$short, function(i){
  c(paste0(rootdir, paste0("/gost.padj.", i, ".csv.gz")))
}))
if (!all(file.exists(target_files))) {
  stop("One or more input files are missing .n", call.=FALSE)
}


# Load data ------
fl = lapply(target_files, fread)
setattr(fl, 'names', tiss$full)
fd = rbindlist(fl, use.names=T, idcol="Tiss")

# clean up for parsing --------------------------------
fd[, c("intersection", "evidence_codes", "symbols", "source_order", "parents", "precision", "recall", "significant", "effective_domain_size") := NULL]
fd[, LogBH := -log10(p_value)]
fd = fd[term_size < 1000]
# cats = c("GO:BP", "GO:CC", "KEGG", "REAC", "WP")
cats = c("GO:BP")
fd = fd[source %in% cats]
plt = fd[, head(.SD, 5), by=c("query", "Tiss")][, term_name]
# plt = fd[Tiss == "Wt", head(.SD, 5), by=c("source", "query")][,term_name]
fd[query == "downregulated", LogBH := -(LogBH)]
a = fd[term_name %in% plt]
# shuffle terms so microtubule based processes isn't at the bottom
term_order = rev(unique(a$term_name))
term_order = append(term_order[-1], term_order[1], after=11)
a$term_name = factor(a$term_name, levels=term_order)
a$Tiss = factor(a$Tiss, levels=rev(tiss$full))

# fill in missing values so plot works ------------------
a = rbindlist(list(a, data.table(
  "Tiss"=tiss$full[1], 
  "LogBH"=0, 
  "term_name"=term_order[12]
                      )), fill=T)
a = rbindlist(list(a, data.table(
  "Tiss"=rep(tiss$full[2], length(term_order[13:18])),
  "LogBH"=rep(0, length(term_order[13:18])),
  "term_name"=term_order[13:18]
                      )), fill=T)
a = rbindlist(list(a, data.table(
  "Tiss"=tiss$full[3], 
  "LogBH"=0, 
  "term_name"=term_order[12]
                      )), fill=T)


# Plot -----------------------------------------------
p = ggplot(a, aes(x=LogBH, y=term_name, fill=Tiss)) + 
  geom_bar(stat = "identity", position="dodge") + 
  guides(fill = guide_legend(reverse=TRUE)) +
  geom_vline(xintercept=0) +
  scale_fill_brewer(type="qual", palette="Set2") +
  xlab(expression("Downregulated" %<-% "-Log10 FDR" %->% "Upregulated")) +
  geom_text(
    data = subset(a, LogBH < 0 & Tiss == tiss$full[1]),
    aes(0, y=term_name, label=term_name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "black",
    family = "Arial",
    size = 4
  ) +
  geom_text(
    data = subset(a, LogBH > 0 & Tiss == tiss$full[1]),
    aes(0, y=term_name, label=term_name),
    hjust = 1,
    nudge_x = -0.3,
    colour = "black",
    family = "Arial",
    size = 4
  ) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
        axis.ticks.length = unit(0, "mm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(family = "Arial", size = 16, face="bold"),
        axis.line.y.left = element_line(color = "#A8BAC4"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = "Arial", size = 16),
        plot.title=element_text(size=20, hjust=0.5, face="bold"),
        plot.subtitle=element_text(size=18, hjust=0.5)) +
  ggtitle("Stimulation effect enrichment of GO:BP terms", 
          sub="top 5 up/downregulated in each condition")
ggsave(filename=paste(today, nameset, "stimulation_enrich_bar.pdf", sep="."),
       plot=p, device=cairo_pdf, width=15, height=10)






