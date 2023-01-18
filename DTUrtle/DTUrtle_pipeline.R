#!/usr/bin/env Rscript --vanilla
# WAJ 2022-01-25


# Parse Arguments ---------------------------------------------------------

library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-d", "--run_dge"), action="store_true",
                     default=FALSE, help="Shall I run a DGE Analysis?")
parser <- add_option(parser, c("-b", "--check_bias"), action="store_true",
                     default=FALSE, help="Shall I check for priming bias?")
parser <- add_option(parser, c("-c", "--condition"), action="store", type="character",
                     default="stim", help="Which condition are you comparing?",
                     metavar='string')
parser <- add_option(parser, c("-g", "--genotype_group"), action="store", type="character",
                     default="WT",
                     help='ONLY WHEN condition="stim": which genotype should be included in the comparison?',
                     metavar='"WT" or "TG"')
parser <- add_option(parser, c("-s", "--stimulation_group"), action="store", type="character",
                     default="noSTIM",
                     help='ONLY WHEN condition="genotype": which stim status should be included in the comparison?',
                     metavar='"STIM" or "noSTIM"')
parser <- add_option(parser, c("-j", "--jobs"), action="store", type="integer",
                     default=32, help='How many cores to use?',
                     metavar='INTEGER')
arguments <<- parse_args(parser)



# Prep and Options ----------------------------------------------------------------

library(DTUrtle)
library(tidyverse)
options("readr.num_threads" = arguments$jobs)
biocpar <- BiocParallel::MulticoreParam(arguments$jobs)

run <- "scg"
gencode <- "no"
# gencode <- "M28"
# gProfiler = TRUE
gProfiler = FALSE
gProfiler_threshold = 0.06

# Functions ---------------------------------------------------------------


set_directory <- function(condition="stim", geno="WT", stim="STIM") {
  root_prefix <- switch(run,
                        "scg" = "/oak/stanford/scg/",
                        "samba" = "/Volumes/",
                        stop("Where is the script running from?"))
  root_dir <<- paste0(root_prefix, "lab_flongo/t41b_BD10-2_stim/DTUrtle/")
  salmon_dir <<- paste0(root_dir, "salmon/")
  output_dir <<- paste0(root_dir, condition, "_", switch(condition,
                                                        "stim" = geno,
                                                        "genotype" = stim,
                                                        stop("Check lines 29-32!")),
                       "/")
  species <<- "Mus_musculus"
  build <<- "GRCm39"
  ensembl_version <<- "105"

  condition <<- condition
  geno <<- geno
  stim <<- stim

  setwd(root_dir)
  }


fetch_metadata <- function() {

  pd <- read_csv(paste0(root_dir, "RNAseq_metadata.csv"))[-(1:2)]
  if (condition == "stim") {
    pd <- pd %>%
      filter(Genotype == switch(geno,
                                "WT" = "WT",
                                "TG" = "TRANS",
                                "TRANS" = "TRANS",
                                "t41b" = "TRANS",
                                "tg" = "TRANS",
                                stop("Wrong Genotype Paramater!")))
  } else if (condition == "genotype") {
    pd <- pd %>%
      filter(Stimulation == switch(stim,
                                   "noSTIM" = "noSTIM",
                                   "stim" = "STIM",
                                   "STIM" = "STIM",
                                   stop("Unknown stimulation parameter!")))
  }
  head(pd, 5)
  return(pd)
  }


fetch_tx2gene <- function() {

  tx2gene <- switch(gencode,
                    "no" = import_gtf(gtf_file = paste0(salmon_dir,species,
                                                        ".",build,
                                                        ".",ensembl_version,
                                                        ".gtf")),
                    "M28" = import_gtf(gtf_file = paste0(salmon_dir,
                                                         "gencode.v", gencode,
                                                         ".annotation.gtf")),
                    stop("Unknown gencode parameter!"))

  tx2gene$gene_name[is.na(tx2gene$gene_name)] <-
    tx2gene$gene_id[is.na(tx2gene$gene_name)]

  tx2gene$transcript_name[is.na(tx2gene$transcript_name)] <-
    tx2gene$transcript_id[is.na(tx2gene$transcript_name)]

  tx2gene <- move_columns_to_front(df = tx2gene,
                                   columns = c("transcript_name", "gene_name"))

  tx2gene$gene_name[is.na(tx2gene$gene_name)] <-
    tx2gene$gene_id[is.na(tx2gene$gene_name)]

  tx2gene$transcript_name[is.na(tx2gene$transcript_name)] <-
    tx2gene$transcript_id[is.na(tx2gene$transcript_name)]

  tx2gene$gene_name <-
    one_to_one_mapping(name = tx2gene$gene_name,
                       id = tx2gene$gene_id)

  tx2gene$transcript_name <-
    one_to_one_mapping(name = tx2gene$transcript_name,
                       id = tx2gene$transcript_id)
  return(tx2gene)
  }


files_and_counts <- function() {

  samples <- pd$Sample_ID
  files <- Sys.glob(paste0(salmon_dir, samples, "_trimmed*/quant.sf"))
  names(files) <- files %>%
    str_remove_all(salmon_dir) %>%
    str_remove_all("_trimmed.fq.gz/quant.sf")
  names(files)
  files <<- files

  cts <- import_counts(files = files,
                       type = 'salmon',
                       tx2gene = tx2gene[, c('transcript_id', 'gene_name')],
                       ignoreTxVersion = TRUE)
  rownames(cts) <- str_remove_all(rownames(cts), "\\...") %>%
    str_remove_all("\\..$")
  rownames(cts) <- tx2gene$transcript_name[match(rownames(cts),
                                                 tx2gene$transcript_id)]
  cts <<- cts
}


dturtle_init <- function() {

  dturtle <- run_drimseq(counts = cts,
                         tx2gene = tx2gene,
                         pd = pd,
                         id_col = "Sample_ID",
                         cond_col =
                           switch(condition,
                                  "stim" = "Stimulation",
                                  "genotype" = "Genotype",
                                  "drug" =  "Treatment",
                                  stop(paste0("'Condition' argument must be in ",
                                              "c('stim', 'genotype', 'drug')."))),
                         filtering_strategy = "bulk",
                         BPPARAM = biocpar)

  dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05, posthoc = 0.1)

  return(dturtle)
}


dge_analysis <- function(dturtle=dturtle,
                         lfc=log2(1.2),
                         p_thresh=1) {

  cts_dge <<- import_dge_counts(files, type = "salmon",
                                tx2gene = tx2gene[, c('transcript_id', 'gene_name')],
                                ignoreTxVersion = TRUE)

  dturtle$dge_analysis <- run_deseq2(counts = cts_dge, pd = pd,
                                     id_col = "Sample_ID",
                                     cond_col =
                                       switch(condition,
                                              "stim" = "Stimulation",
                                              "genotype" = "Genotype",
                                              "drug" =  "Treatment",
                                              stop(paste0("'Condition' argument must be in ",
                                                          "c('stim', 'genotype', 'drug')."))),
                                     lfc_threshold = lfc, # Rob's threshold, may change for future projects
                                     sig_threshold = p_thresh, # filtered during post-hoc wrangling
                                     dge_calling_strategy = "bulk",
                                     BPPARAM = biocpar)

  dge <- as_tibble(dturtle$dge_analysis$results_sig) %>%
    arrange(svalue, desc(abs(log2FoldChange)))

  write_csv(dge, paste0(output_dir,
                        "DGE_results_",
                        condition, "_",
                        switch(condition,
                               "stim" = geno,
                               "genotype" = stim,
                               stop("Check lines 29-32!")),
                        ".csv"))
  return(dturtle)
}


check_priming_bias <- function(dturtle=dturtle) {

  priming_bias_df <- priming_bias_detection_probability(
    counts = cts,
    gtf = paste0(salmon_dir,species,
                 ".",build,
                 ".",ensembl_version,
                 ".gtf"),
    tx2gene = tx2gene,
    one_to_one = TRUE,
    priming_enrichment = "3",
    BPPARAM = biocpar
  )

  priming_bias_df %>%
    as_tibble() %>%
    write_csv(paste0(output_dir,
                     "priming_bias_",
                     condition, "_",
                     switch(condition,
                            "stim" = geno,
                            "genotype" = stim,
                            stop("Check lines 29-32!")),
                     ".csv"))
}


prep_dtu_table <- function(dturtle=dturtle) {

  dturtle <- create_dtu_table(dturtle = dturtle,
                              add_gene_metadata = list("chromosome" = "seqnames"),
                              add_tx_metadata =
                                list("tx_expr_in_max" = c("exp_in", max)))

  tb <- tibble::as_tibble(dturtle$dtu_table)

  write_csv(tb, file = paste0(output_dir,
                              "dtu_table_",
                              condition, "_",
                              switch(condition,
                                     "stim" = geno,
                                     "genotype" = stim,
                                     stop("Check lines 29-32!")),
                              ".csv"))
}


save_work <- function(dturtle=dturtle) {

  save.image(paste0(output_dir,
                    "DTUrtle_Analysis_",
                    condition, "_",
                    switch(condition,
                           "stim" = geno,
                           "genotype" = stim,
                           stop("Check lines 29-32!")),
                    ".RData"))

  saveRDS(dturtle, paste0(output_dir,
                          "dturtle_object_",
                          condition, "_",
                          switch(condition,
                                 "stim" = geno,
                                 "genotype" = stim,
                                 stop("Check lines 29-32!")),
                          ".rds"))
  write_tsv(pd, paste0(output_dir,
                       "salmon_metadata_",
                       condition, "_",
                       switch(condition,
                              "stim" = geno,
                              "genotype" = stim,
                              stop("Check lines 29-32!")),
                       ".txt"))
}


main <- function(condition="stim",
                 geno="WT",
                 stim="STIM",
                 run_dge=FALSE,
                 check_bias=FALSE) {
  set_directory(condition = condition,
                geno = geno,
                stim = stim)

  pd <<- fetch_metadata()
  head(pd, 5)

  tx2gene <<- fetch_tx2gene()

  files_and_counts()

  dturtle <<- dturtle_init()

  if (run_dge == TRUE) {
    dturtle <<- dge_analysis(dturtle = dturtle,
                             lfc = log2(1.2),
                             p_thresh = 1)
  } else if (run_dge == FALSE) {
      warning("No DGE analysis was performed.")
    } else stop("Should I run DGE Analysis or not?")

  if (check_bias == TRUE) {
    check_priming_bias(dturtle = dturtle)
  } else if (check_bias == FALSE) {
      warning("Did not check for priming bias.")
    } else stop("Should I check for priming bias or not?")

  dturtle <<- prep_dtu_table(dturtle = dturtle)

  save_work(dturtle = dturtle)
}



# Actual Run ------------------------------------------------------------

main(condition = arguments$condition,
     geno = arguments$genotype_group,
     stim = arguments$stimulation_group,
     run_dge = arguments$run_dge,
     check_bias = arguments$check_bias)
