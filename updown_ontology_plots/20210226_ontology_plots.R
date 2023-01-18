#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(gridExtra)


# directory locations
rootdir = '/oak/stanford/scg/lab_flongo/2020_T41B_BD10-2_Stimulation'
wk_dir = paste0(rootdir, '/updown_ontology_plots')
setwd(paste0(rootdir, '/updown_ontology_plots'))

# file naming info
today = format(Sys.Date(), '%Y%m%d')
nameset = '_ontology_'

# get file info
tiss = data.table("short"=c("2", "9", "12"), 
                  "full"=c("Genotype Effect", "Drug Effect", 
                           "Genotype Effect + Drug Effect")) 
# "10", "Drug Effect Interaction", 
cats = c("CORUM","GO:BP","GO:CC","GO:MF","HP","KEGG","MIRNA","REAC","TF","WP")
target_files = unlist(lapply(tiss$short, function(i){
  c(paste0(rootdir, paste0("/gost.", i, ".csv")))
}))

# get list of terms to plot from file
sht_lst = fread('20210127_terms.txt', sep="\t", header=F)


########FUNCTIONS
# string wrap function for long text
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

# make table for plot (for a given term negative downreg LogBH)
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

########MAIN
# read in files
all(file.exists(target_files))
# reading in tables, combine and melt by GWAS
fl = lapply(target_files, fread)
setattr(fl, 'names', tiss$short)
fd = rbindlist(fl, use.names=T, idcol="Tiss")

# clean up for parsing
fd[, c("V1", "evidence_codes", "Gene_name") := NULL]
fd[, DGE := gsub("regulated", "", DGE)]
fd[, LogBH := -log10(p_value)]
fd[, Title := paste(term_id, term_name)]

# list of unique terms with significant result
terms = unique(fd[LogBH > -log10(0.05) & source %in% cats
                  ][order(-LogBH), Title])

# lookup terms w/ term_ids
match_terms = unique(fd$Title[fd$term_name %in% sht_lst$V1])

# plot for relevant tissues
plist = lapply(match_terms, mk_plt, tbl=fd)

# generate pdf
ml <- marrangeGrob(plist, nrow=3, ncol=4, respect=F, padding=unit(15, "mm"),
                   layout_matrix=t(matrix(1:12, nrow=4, ncol=3)))
ggsave(paste0(today, nameset, "plots.pdf"), ml, width=420, height=594, unit="mm",
       useDingbats=F)


