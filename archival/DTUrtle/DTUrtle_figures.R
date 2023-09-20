#!/usr/bin/env Rscript --vanilla
# WAJ 2022-01-25

library(tidyverse)
library(DTUrtle)
biocpar <- BiocParallel::MulticoreParam(12)

# Options
run <- "local"
gencode <- "no"
# gencode <- "M28"

condition <- "genotype"

root_prefix <- switch(run,
                      "scg" = "/oak/stanford/scg/",
                      "local" = "/Volumes/",
                      stop("Where is the script running from?"))
root_dir <- paste0(root_prefix, "lab_flongo/t41b_BD10-2_stim/DTUrtle/", condition, "/")
salmon_dir <- paste0(root_prefix, "lab_flongo/t41b_BD10-2_stim/DTUrtle/salmon/")
species <- "Mus_musculus"
build <- "GRCm39"
version <- "105"
setwd(root_dir)
getwd()

dturtle <- readRDS(paste0("dturtle_object_", condition, ".rds"))

dturtle <- create_dtu_table(dturtle = dturtle,
                            add_gene_metadata = list("chromosome" = "seqnames"),
                            add_tx_metadata = list("tx_expr_in_max" = c("exp_in", max)))
tb <- tibble::as_tibble(dturtle$dtu_table)

dturtle <- plot_proportion_barplot(dturtle = dturtle,
                                   meta_gene_id = "gene_id.1",
                                   savepath = "./images/barplots",
                                   add_to_table = "barplot",
                                   BPPARAM = biocpar)

dturtle <- plot_proportion_pheatmap(dturtle = dturtle,
                                    include_expression = TRUE,
                                    treeheight_col = 20,
                                    savepath = "./images/heatmaps",
                                    add_to_table = "pheatmap",
                                    BPPARAM = biocpar)

gtf <- switch(gencode, "no" = import_gtf(paste0(salmon_dir,
                                                species, ".",
                                                build, ".",
                                                version, ".", "gtf"),
                                         feature_type = NULL,
                                         out_df = FALSE),
              "M28" = import_gtf(paste0(salmon_dir,
                                        "gencode.v", gencode,
                                        ".annotation.gtf"),
                                 feature_type = NULL,
                                 out_df = FALSE),
              stop("Unknown gencode parameter!"))

dturtle <- plot_transcripts_view(dturtle = dturtle,
                      gtf = gtf,
                      genome = 'mm39',
                      one_to_one = TRUE,
                      savepath = "./images/transcripts",
                      add_to_table = "transcript_view",
                      BPPARAM = biocpar)

column_formatter_list <- list(
  "gene_qvalue" = table_pval_tile("white", "orange", digits = 3),
  "minimal_tx_qvalue" = table_pval_tile("white", "orange", digits = 3),
  "number_tx" = formattable::color_tile('white', "lightblue"),
  "number_significant_tx" = formattable::color_tile('white', "lightblue"),
  "max(Dex2hr-EtOH)" = table_percentage_bar('lightgreen', "#FF9999",
                                            digits = 2),
  "tx_expr_in_max" = table_percentage_bar('white', "lightblue",
                                          color_break = 0, digits = 2)
)

plot_dtu_table(dturtle = dturtle,
               savepath = paste0("DTUrtle_results_", condition, ".html"),
               column_formatters = column_formatter_list)

saveRDS(dturtle, file = paste0("dturtle_results_", condition, ".rds"))
