#!/usr/bin/env Rscript --vanilla
# WAJ 2021-11-24

# Drink the tidyverse Kool-Aid and load gProfiler2:
library(tidyverse)
library(data.table)
setDTthreads(0)
library(dtplyr)
library(lubridate)
library(gprofiler2)

# Setup, Preferences, and Data Import --------------------------------------------------
# Let's make sure the script knows where to find everything it needs:
root_dir <- "/oak/stanford/scg/lab_flongo/t41b_BD10-2_stim/Correlation_Plots_and_Pathway_Enrichments/"
data_dir <- paste0(root_dir, "data/")
plots_dir <- paste0(root_dir, "plots/")
output_dir <- paste0(root_dir, "gprofiler/")

# Change your term search and threshold preferences!
term_name_regex <- "neur|dendr|glia|p75|NFKB|glutamate|GABA|gaba|synap"
padj_threshold <- 0.05
gost_threshold <- 0.06
l2fc_threshold <- log2(1.2)

# I have the script set to read from the working directory, so lets set the
# working directory to data_dir here:
setwd(data_dir)

# Read in files:
fileNames <- list.files(pattern = ".csv$", recursive = T)
files <- lapply(fileNames, FUN = function(i) {
  read_csv(i) %>%
    select(2:11) %>%
    filter(log2FoldChange < 20)
})

# Make the file names human readable; they become figure axis labels later on
fileNames <- fileNames %>%
  str_remove("results_") %>%
  str_remove("_annotated.csv") %>%
  str_replace_all("_", " ") %>%
  str_replace_all("vs", " vs. ") %>%
  str_replace_all("  ", " ") %>%
  str_replace_all("TG", "Transgenic") %>%
  str_replace_all("WT", "Wild Type") %>%
  str_replace_all("VEH", "Vehicle") %>%
  str_replace_all("NOSTIM", "noSTIM") %>%
  str_replace_all("noSTIM", "No Stimulation") %>%
  str_to_title() 

# This list will help form  output filenames; replace " " with "_" to reduce
# future headaches:
pdfNames <- str_replace_all(fileNames, " ", "_")


# Functions --------------------------------------------------------------


joinLM <- function(i, j) {
  # Joins two datasets and fits a linear regression to their L2FC's
  results.total <- inner_join(files[[i]], files[[j]], by = "Gene_id",
                              suffix = c("_X", "_Y")) %>%
    filter(padj_X < padj_threshold | padj_Y < padj_threshold) %>%
    filter(!is.na(padj_X) & !is.na(padj_Y))
  
  return(lm(log2FoldChange_Y ~ log2FoldChange_X, data = results.total))
}

corPlot <- function(i, j, pdf=FALSE) {
  # Does what joinLM does, plus creates a pretty l2fc correlation plot
  results.total <- inner_join(files[[i]], files[[j]], by = "Gene_id",
                              suffix = c("_X", "_Y")) %>%
    filter(padj_X < padj_threshold | padj_Y < padj_threshold,) %>%
    filter(abs(log2FoldChange_X) > l2fc_threshold & abs(log2FoldChange_Y) > l2fc_threshold) %>%
    filter(!is.na(padj_X) & !is.na(padj_Y))
  
  lreg <- lm(log2FoldChange_Y ~ log2FoldChange_X, data = results.total)

  ggplot(results.total, aes(x = log2FoldChange_X, y = log2FoldChange_Y)) +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    geom_point(size = 0.5) + 
    geom_smooth(method = "lm", formula = "y ~ x", se = F, size = 0.5, color = "orange") +
    geom_hline(yintercept = c(l2fc_threshold,-(l2fc_threshold)), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(l2fc_threshold,-(l2fc_threshold)), linetype = "dashed", color = "blue") +
    theme_bw() +
    labs(x = fileNames[i],
         y = fileNames[j],
         title = "L2FC Correlation Plots",
         subtitle = paste0("y = ", round(lreg$coefficients[2], 3), "x + ",
                           round(lreg$coefficients[1], 3), " || R2 = ",
                           round(summary(lreg)$adj.r.squared, 3)))
  if (pdf == TRUE) {
    ggsave(paste0(plots_dir, round(summary(lreg)$adj.r.squared, 3),
                  pdfNames[i], "__", pdfNames[j],
                  "_correlation_Plots.pdf"),
           plot = last_plot(), width = 9, height = 6.5)
  } else suppressWarnings(print((last_plot())))
}


# Make a tibble of all relevant data + Create Plots ----------------------
tb <- tibble(
  x = character(),
  y = character(),
  comparison = character(),
  r2Value = numeric()
)
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    lreg <- joinLM(i, j)
    tb <- add_row(tb,
                  x = paste0(i),
                  y = paste0(j),
                  comparison = paste0(fileNames[[i]], " || ", fileNames[[j]]),
                  r2Value = round(summary(lreg)$adj.r.squared, 3))
    corPlot(i, j, pdf = T)
  }
}

tb <- tb %>%
  arrange(desc(r2Value)) %>%
  mutate(rValue = round(sqrt(r2Value), 3))

tb
fileNames


# Find Interesting Pathway Enrichments -------------------------------------------------------------


## TO DO: I repeat myself to obscene degrees in creating "results.total.a-d". Need to fix this. 
for (i in 1:(length(files) - 1)) {
  for (j in (i + 1):length(files)) {
    results.total.a <- inner_join(files[[i]], files[[j]], by = "Gene_id",
                                suffix = c("_X", "_Y")) %>%
      filter(padj_X < padj_threshold | padj_Y < padj_threshold) %>%
      filter(!is.na(padj_X) & !is.na(padj_Y)) %>%
      filter((log2FoldChange_Y > l2fc_threshold & log2FoldChange_X > l2fc_threshold)) %>%
      select("GeneSymbol_X", "padj_X", "padj_Y", everything()) %>%
      arrange(padj_X)

    results.total.b <- inner_join(files[[i]], files[[j]], by = "Gene_id",
                                  suffix = c("_X", "_Y")) %>%
      filter(padj_X < padj_threshold | padj_Y < padj_threshold) %>%
      filter(!is.na(padj_X) & !is.na(padj_Y)) %>%
      filter((log2FoldChange_X < -(l2fc_threshold) & log2FoldChange_Y > l2fc_threshold)) %>%
      select("GeneSymbol_X", "padj_X", "padj_Y", everything()) %>%
      arrange(padj_X)

    results.total.c <- inner_join(files[[i]], files[[j]], by = "Gene_id",
                                  suffix = c("_X", "_Y")) %>%
      filter(padj_X < padj_threshold | padj_Y < padj_threshold) %>%
      filter(!is.na(padj_X) & !is.na(padj_Y)) %>%
      filter((log2FoldChange_X < -(l2fc_threshold) & log2FoldChange_Y < -(l2fc_threshold))) %>%
      select("GeneSymbol_X", "padj_X", "padj_Y", everything()) %>%
      arrange(padj_X)

    results.total.d <- inner_join(files[[i]], files[[j]], by = "Gene_id",
                                  suffix = c("_X", "_Y")) %>%
      filter(padj_X < padj_threshold | padj_Y < padj_threshold) %>%
      filter(!is.na(padj_X) & !is.na(padj_Y)) %>%
      filter((log2FoldChange_X > l2fc_threshold & log2FoldChange_Y < -(l2fc_threshold))) %>%
      select("GeneSymbol_X", "padj_X", "padj_Y", everything()) %>%
      arrange(padj_X)

    if (nrow(results.total.a) > 0 | nrow(results.total.b) > 0 |
        nrow(results.total.c) > 0 | nrow(results.total.d) > 0) {
      gostres <- gost(list("BothUp" = results.total.a$Gene_id,
                           "XDown_YUp" = results.total.b$Gene_id,
                           "BothDown" = results.total.c$Gene_id,
                           "XUp_YDown" = results.total.d$Gene_id),
                      organism = "mmusculus", significant = F,
                      correction_method = "fdr", evcodes = TRUE)

        if (is.list(gostres$result)) {
          results <- as_tibble(gostres$result) %>%
            filter(p_value < gost_threshold,
                   query %in% c("XDown_YUp", "XUp_YDown"),
                   intersection_size > 3,
                   term_size < 1000) %>%
            arrange(p_value)
          
          # results <- results[grep(term_name_regex, results$term_name), ]

          if (nrow(results) > 0) {
            write_csv(results, file = paste0(output_dir, fileNames[i], " || ",
                                             fileNames[j], ".csv"))
            print(paste0("Written ", fileNames[i], " || ", fileNames[j], ".csv to file."))
          }

        }
    }
  }
}
