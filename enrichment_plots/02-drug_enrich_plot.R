#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 02-drug_enrich_plot.R
#'
#' Version: 0.0.3
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
library(patchwork)
library(stringr)

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
  c(paste0(rootdir, paste0("/gost.", i, ".csv")))
}))
if (!all(file.exists(target_files))) {
  stop("One or more input files are missing .n", call. = FALSE)
}


# Load data --------------------------------------------
fl <- lapply(target_files, fread)
setattr(fl, "names", tiss$full)
fd <- rbindlist(fl, use.names = T, idcol = "Tiss")

# clean up for parsing --------------------------------
drops = c(
  "intersection", 
  "evidence_codes", 
  "source_order", 
  "precision", 
  "recall", 
  "significant", 
  "effective_domain_size", 
  "Gene_name"
) 
fd[, (drops) := NULL]
fd[, LogBH := -log10(p_value)]
fd[DGE == "downregulated", LogBH := -(LogBH)]
fd <- fd[term_size < 1000]

# prepping terms of interest ---------------------------
plt <- c(
  "cytokine production", 
  "inflammatory response",
  "regulation of neuron death",
  "intrinsic apoptotic signaling pathway", 
  "cytokine receptor activity", 
  "neuroinflammatory response",
  "Complement cascade",
  "cell aging", 
  "Alzheimers Disease", 
  "MAP kinase activation",
  "regulation of amyloid-beta formation",
  "regulation of neurogenesis",
  "cell projection membrane",
  "axoneme",
  "regulation of ERK1 and ERK2 cascade",
  "positive regulation of ERK1 and ERK2 cascade", 
  "negative regulation of ERK1 and ERK2 cascade",
  "PI3K-Akt signaling pathway",
  "JNK cascade",
  "regulation of Wnt signaling pathway",
  "positive regulation of Wnt signaling pathway",
  "negative regulation of Wnt signaling pathway",
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
  "presynapse"
)

# Build short list
a <- fd[term_name %in% plt]

# Make terms sentence case with abbreviations
abbr <- c("ERK", "Wnt", "MAP", "JNK", "PI3K-Akt", "EEG", "GABA")
repl <- rep(abbr, 2)
names(repl) <- c(str_to_sentence(abbr), str_to_lower(abbr))
plt <- str_replace_all(str_to_sentence(plt), repl)
a$term_name <- str_replace_all(str_to_sentence(a$term_name), repl)

# Assign factors for ordering of plots
a$term_name <- factor(a$term_name, levels = plt)
a$Tiss <- factor(a$Tiss, levels = rev(tiss$full))

# fill in missing values so plot works 
miss <- setdiff(plt, a[LogBH <=0 & Tiss == tiss$full[2], term_name])
a <- rbindlist(list(a, 
  data.table(
    "Tiss" = rep(tiss$full[2], length(miss)),
    "LogBH" = rep(0, length(miss)),
    "term_name" = miss
  )
), fill = TRUE)


# Plot -----------------------------------------------
a %>% 
  ggplot(aes(x = LogBH, y = term_name, fill = Tiss)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  xlab(expression("Downregulated" %<-% "-Log10 FDR" %->% "Upregulated")) +
  geom_text(
    data = a[LogBH <= 0 & Tiss == tiss$full[2]],
    aes(0, y=term_name, label=term_name),
    hjust = 0,
    nudge_x = 10,
    colour = "black",
    family = "Arial",
    size = 4
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
    axis.ticks.length = unit(0, "mm"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(family = "Arial", size = 16, face = "bold"),
    axis.line.y.left = element_line(color = "#A8BAC4"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Arial", size = 16),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 18, hjust = 0.5)
  ) +
  ggtitle("Selected LTP-relevant Term Enrichment") -> p1

ggsave(
  filename = paste(today, nameset, "drug_enrich_bar.pdf", sep = "."),
  plot = p1, 
  device = cairo_pdf, 
  width = 13, 
  height = 18
)


