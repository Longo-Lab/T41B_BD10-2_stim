#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 01-stim_get_GO.R
##
## Version: 0.0.1
##
## Purpose of script: get the missing stim data for the 5, 6, 7 stim DEs
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
##   originals not in the OneDrive folder
##
##   Usage:
##     sbatch -J gos --mem=32G -c 4 -t 01:00:00 -p interactive -A default \
##       -o %A_stim_get_GO.log \
##       --wrap "ml R/4.0; Rscript 01-stim_get_GO.R"
##
##   or interactive session:
##     sdev -m 32 -c 4 -t 01:00:00 -p interactive -a default
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)

library(data.table)
library(ggplot2)
library(gprofiler2)
library(stringr)
library(magrittr)

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
  c(paste0(rootdir, paste0("/res.", i, ".csv")))
}))
if (!all(file.exists(target_files))) {
  stop("One or more input files are missing .n", call.=FALSE)
}

# lfc col (for up/down split)
lfc = "log2FoldChange"
# pvalue, padj or svalue
p = 'padj'
# cutoff threshold
thresh = 0.05
# ensg id/symbol cols
ensg = "Gene_id"
sym = "GeneSymbol"

# Functions ----------------------------------------------

get_mcols = function(str, dds, col) {
  # for a string of ens ids: explode, look up symbols, return symbols as string
  # string: comma separated ens_ids
  # dds: DESeq Dataset with annotations in mcol
  # col: colname to grab
  require("stringr")
  # x = str_split(str, ",")
  y = dds[[col]][match(unlist(str_split(str, ",")), dds[[ensg]])]
  return(paste0(y, collapse=","))
}


# Main ----------------------------------------------------
# Pathway analyses ----------------------------------------
require(gprofiler2)
lapply(tiss$short, function(i){
  # get dataframe based on keyword list (also bg list)
  x = fread(paste0(rootdir, paste0("/res.", i, ".csv")))
  x = x[ !is.na(get(ensg)) ]
  bg = x[[ensg]]
  # drop NA pvalues
  x = x[ !is.na(get(p)) ]
  # sort by pvalue
  x = x[with(x, order(get(p))), ]
  # define up/down subsets
  up = x[get(lfc) > 0 & get(p) < thresh, get(ensg)]
  down = x[get(lfc) < 0 & get(p) < thresh, get(ensg)]
  # lookup using ordered subset
  g = gost(query=list("upregulated"=up, "downregulated"=down), 
           organism="mmusculus", ordered_query=T, evcodes=T, custom_bg=bg)
  # write to file
  if(is.list(g$result)){
    # get symbols for terms
    g$result$symbols = lapply(g$result$intersection, get_mcols, x, sym)
    # write to file
    fwrite(g$result, paste0(rootdir, "/gost.padj.", i, ".csv.gz"))
  }
})


