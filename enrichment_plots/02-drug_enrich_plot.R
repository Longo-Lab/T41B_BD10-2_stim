#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-drug_enrich_plot.R
#'
#' Version: 0.0.2
#'
#' Purpose of script: plot top ontology enrichments across three drug conditions
#'
#' Author: Robert R Butler III
#'
#' Date Created: 2024-02-12
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   Using padj < 0.05 enrichment from gprofiler, limiting to term size < 1000
#'
#'   Usage:
#'     sbatch -J gos --mem=32G -c 4 -t 01:00:00 -p interactive -A default \
#'       -o %A_stim_enrich_plot.log \
#'       --wrap "ml R/4.0; Rscript 02-stim_enrich_plot.R"
#'
#'   or interactive session:
#'     sdev -m 32 -c 4 -t 01:00:00 -p interactive -a default
#'
#' ---------------------------

#' load up the packages we will need:  (uncomment as required)

library(data.table)
library(ggplot2)


# filepaths
rootdir <- "/labs/flongo/t41b_BD10-2_stim"
wk_dir <- paste0(rootdir, "/enrichment_plots")
setwd(wk_dir)

# file naming info
today <- format(Sys.Date(), "%Y%m%d")
nameset <- "BD10-2"

# get file info
tiss <- data.table(
  "short" = c("2", "9", "12"),
  "full" = c("APP", "BD10-2", "APP-BD10-2")
)
target_files <- unlist(lapply(tiss$short, function(i) {
  c(paste0(rootdir, paste0("/gost.padj.", i, ".csv.gz")))
}))
if (!all(file.exists(target_files))) {
  stop("One or more input files are missing .n", call. = FALSE)
}


# Load data ------
fl <- lapply(target_files, fread)
setattr(fl, "names", tiss$full)
fd <- rbindlist(fl, use.names = T, idcol = "Tiss")

# clean up for parsing --------------------------------
fd[, c("intersection", "evidence_codes", "symbols", "source_order", "parents", "precision", "recall", "significant", "effective_domain_size", "Gene_name") := NULL]
fd[, LogBH := -log10(p_value)]
fd[query == "downregulated", LogBH := -(LogBH)]
fd <- fd[term_size < 1000]
# cats = c("GO:BP", "GO:CC", "KEGG", "REAC", "WP")
# cats = c("GO:BP")
# fd = fd[source %in% cats]
# plt = fd[, head(.SD, 10), by=c("DGE", "Tiss")][, term_name]
# plt = fd[Tiss == "Wt", head(.SD, 5), by=c("source", "query")][,term_name]
plt <- c(
  "cytokine production",
  "inflammatory response",
  "cytokine receptor activity",
  "regulation of neuron death",
  "intrinsic apoptotic signaling pathway",
  "neuroinflammatory response",
  "complement activation",
  "cell aging",
  "Alzheimers Disease",
  "MAP kinase activation",
  "positive regulation of amyloid-beta formation",
  "axoneme",
  "regulation of long-term synaptic potentiation",
  "long-term synaptic potentiation",
  "glutamate receptor activity",
  "EEG with spike-wave complexes",
  "EEG abnormality",
  "EEG with generalized epileptiform discharges",
  "vesicle-mediated transport in synapse",
  "signal release from synapse",
  "neuronal action potential",
  "Glutamatergic synapse",
  "GABA-ergic synapse",
  "postsynaptic membrane",
  "synaptic signaling",
  "postsynapse",
  "presynapse",
  "Complement cascade",
  "Complement activation, classical pathway",
  "PI3K-Akt signaling pathway",
  "Focal adhesion: PI3K-Akt-mTOR signaling pathway",
  "JNK cascade",
  "ERK1 and ERK2 cascade",
  "Wnt signaling pathway and pluripotency",
  "Wnt signaling pathway"
)
a <- fd[term_name %in% plt]
# # shuffle terms so microtubule based processes isn't at the bottom
# term_order = rev(unique(a$term_name))
# term_order = append(term_order[-1], term_order[1], after=11)
a$term_name <- factor(a$term_name, levels = plt)
a$Tiss <- factor(a$Tiss, levels = rev(tiss$full))

# fill in missing values so plot works ------------------
# a = rbindlist(list(a, data.table(
# "Tiss"=tiss$full[2],
# "LogBH"=0,
# "term_name"=term_order[23]
# )), fill=T)
# a = rbindlist(list(a, data.table(
# "Tiss"=rep(tiss$full[2], length(plt[c(15, 21, 23)])),
# "LogBH"=rep(0, length(plt[c(15, 21, 23)])),
# "term_name"=plt[c(15, 21, 23)]
# )), fill=T)
# a = rbindlist(list(a, data.table(
# "Tiss"=tiss$full[3],
# "LogBH"=0,
# "term_name"=term_order[12]
# )), fill=T)


# Plot -----------------------------------------------
p <- ggplot(a, aes(x = LogBH, y = term_name, fill = Tiss)) +
  geom_bar(stat = "identity", position = "dodge") +
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  xlab(expression("Downregulated" %<-% "-Log10 FDR" %->% "Upregulated")) +
  geom_text(
    data = subset(a, LogBH <= 0 & Tiss == tiss$full[3]),
    aes(0, y = term_name, label = term_name),
    hjust = 0,
    nudge_x = 10,
    colour = "black",
    family = "Arial",
    size = 4
  ) +
  geom_text(
    data = subset(a, LogBH > 0 & Tiss == tiss$full[3]),
    aes(0, y = term_name, label = term_name),
    hjust = 0,
    nudge_x = 10,
    colour = "black",
    family = "Arial",
    size = 4
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    axis.ticks.length = unit(0, "mm"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(family = "Arial", size = 16, face = "bold"),
    axis.line.y.left = element_line(color = "#A8BAC4"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Arial", size = 16),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 18, hjust = 0.5)
  ) +
  ggtitle("Selected LTP-relevant Term Enrichment")
ggsave(
  filename = paste(today, nameset, "drug_enrich_bar.pdf", sep = "."),
  plot = p, device = cairo_pdf, width = 13, height = 12
)
