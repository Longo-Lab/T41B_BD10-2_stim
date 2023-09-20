library(tidyverse)
library(lubridate)

run <- "scg"

if (run == "scg") options("readr.num_threads" = 16L)


change_local <- function(condition) {
  root_dir <- switch(run,
                     "samba" = paste0(
                       "/Volumes/lab_flongo/t41b_BD10-2_stim/DTUrtle/",
                       condition
                     ),
                     "scg" = paste0(
                       "/labs/flongo/t41b_BD10-2_stim/DTUrtle/",
                       condition
                     ),
                     stop("Unknown run parameter!"))
  setwd(root_dir)
}


wrangle <- function(condition, dge=FALSE){
  change_local(condition)
  tb <- read_csv(paste0("dtu_table_", condition, ".csv"))

  if (condition == "genotype") {
    tb <- tb %>%
      filter(
        number_significant_tx > 0,
        minimal_tx_qvalue < 0.05,
        abs(`max(WT-TRANS)`) > 0.1
        ) %>%
      arrange(desc(abs(`max(WT-TRANS)`)))
  } else if (condition == "stim") {
      tb <- tb %>%
        filter(
          number_significant_tx > 0,
          minimal_tx_qvalue < 0.05,
          abs(`max(noSTIM-STIM)`) > 0.1
        ) %>%
        arrange(desc(abs(`max(noSTIM-STIM)`)))
  } else stop("Unknown condition parameter!")

  if (dge == TRUE) {
    deseq_res <- read_csv(paste0("DGE_results_", condition, ".csv"))

    tb <- tb %>%
      left_join(deseq_res, by = c("gene_ID" = "gene"))
  } else warning("DGE results NOT included in this output.")

  write_csv(tb, file = paste0(today(),
                              "_filtered_dtu_table_",
                              condition,
                              ".csv"
                              ))

  return(tb)
}


profiling_those_Gs <- function(query, condition) {
  change_local(condition)
  gprof <- gprofiler2::gost(query = query,
                            organism = "mmusculus",
                            ordered_query = TRUE,
                            significant = FALSE,
                            evcodes = TRUE,
                            correction_method = "fdr",
                            user_threshold = 0.05)

  saveRDS(gprof, file = paste0(today(),
                               "_gProfiler_object_",
                               condition,
                               ".rds"))

  tb <- gprof$result %>%
    as_tibble() %>%
    filter(
      p_value < 0.06,
      intersection_size > 1,
      term_size < 1000
    ) %>%
    select(-query) %>%
    select(term_name, term_id, everything()) %>%
    arrange(p_value)

  write_csv(tb, file = paste0(today(),
                              "_gProfiler_results_",
                              condition,
                              ".csv"))

  return(tb)
}


stim <- wrangle("stim", dge = TRUE)
genotype <- wrangle("genotype", dge = TRUE)

stim_go <- profiling_those_Gs(stim$gene_ID, "stim")
genotype_go <- profiling_those_Gs(genotype$gene_ID, "genotype")
