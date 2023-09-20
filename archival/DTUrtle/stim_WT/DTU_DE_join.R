#!/usr/bin/env Rscript --vanilla
# WAJ 20220318

library(tidyverse)

setwd("/labs/flongo/t41b_BD10-2_stim/DTUrtle/stim_WT")

(dtu <- read_csv("dtu_table_stim_WT.csv") %>%
  filter(number_significant_tx >= 1,
         minimal_tx_qvalue < 0.05) %>%
  arrange(desc("max(noSTIM-STIM)")))

dge <- read_csv("DGE_results_stim_WT.csv") %>%
  mutate(fdr = p.adjust(svalue, method = "fdr")) %>%
  mutate(de_significant = fdr < 0.05)

dge.significant <- filter(dge, de_significant == TRUE)
dge.boring <- filter(dge, de_significant == FALSE)

dtu.sig_dge <- inner_join(dtu, dge.significant, by = c("gene_ID" = "gene")) %>%
  arrange(fdr)
dtu.boring_dge <- inner_join(dtu, dge.boring, by = c("gene_ID" = "gene")) %>%
  arrange(desc("max(noSTIM-STIM"))
dtu.dge <- inner_join(dtu, dge, by = c("gene_ID" = "gene")) %>%
  arrange(fdr)

write_csv(dtu.sig_dge, file = paste0(lubridate::today(), "_DTU_table_with_signifant_DE.csv"))
write_csv(dtu.boring_dge, file = paste0(lubridate::today(), "_DTU_table_with_insignificant_DE.csv"))
write_csv(dtu.dge, file = paste0(lubridate::today(), "_DTU_table_with_all_DE.csv"))


# Ranked gProfiler Query -------------------------------------------------

run_all_gProfilers <- function() {
  require(gprofiler2)

  general_res <<- gost(dtu.dge$gene_ID, organism = "mmusculus",
                      ordered_query = TRUE, significant = TRUE,
                      evcodes = TRUE, correction_method = "fdr")

  general_result <<- as_tibble(general_res$result) %>%
    select(-1) %>%
    select(term_name, everything())

  write_csv(general_result, file = paste0(lubridate::today(), "_gProfiler_ranked_query_all_DE.csv"))

  boring_res <<- gost(dtu.boring_dge$gene_ID, organism = "mmusculus",
              ordered_query = TRUE, significant = TRUE,
              evcodes = TRUE, correction_method = "fdr")

  boring_result <<- as_tibble(boring_res$result) %>%
    select(-1) %>%
    select(term_name, everything())

  write_csv(boring_result, file = paste0(lubridate::today(), "_gProfiler_ranked_query_boring_DE.csv"))

  res_with_significance <<- gost(dtu.dge$gene_ID, organism = "mmusculus",
                                ordered_query = TRUE, significant = TRUE,
                                evcodes = TRUE, correction_method = "fdr")

  result_with_significance <<- as_tibble(res_with_significance$result) %>%
    select(-1) %>%
    select(term_name, everything())

  write_csv(result_with_significance, file = paste0(lubridate::today(), "_gProfiler_ranked_query_significant_DE.csv"))
}

run_all_gProfilers()

# Do it all over again with genotype_noSTIM ------------------------------



setwd("/labs/flongo/t41b_BD10-2_stim/DTUrtle/genotype_noSTIM")

(dtu <- read_csv("dtu_table_genotype_noSTIM.csv") %>%
    filter(number_significant_tx >= 1,
           minimal_tx_qvalue < 0.05) %>%
    arrange(desc("max(WT-TRANS)")))

dge <- read_csv("DGE_results_genotype_noSTIM.csv") %>%
  mutate(fdr = p.adjust(svalue, method = "fdr")) %>%
  mutate(de_significant = fdr < 0.05)

dge.significant <- filter(dge, de_significant == TRUE)
dge.boring <- filter(dge, de_significant == FALSE)

dtu.sig_dge <- inner_join(dtu, dge.significant, by = c("gene_ID" = "gene")) %>%
  arrange(fdr)
dtu.boring_dge <- inner_join(dtu, dge.boring, by = c("gene_ID" = "gene")) %>%
  arrange(desc("max(WT-TRANS"))
dtu.dge <- inner_join(dtu, dge, by = c("gene_ID" = "gene")) %>%
  arrange(fdr)

write_csv(dtu.sig_dge, file = paste0(lubridate::today(), "_DTU_table_with_signifant_DE.csv"))
write_csv(dtu.boring_dge, file = paste0(lubridate::today(), "_DTU_table_with_insignificant_DE.csv"))
write_csv(dtu.dge, file = paste0(lubridate::today(), "_DTU_table_with_all_DE.csv"))

run_all_gProfilers()
