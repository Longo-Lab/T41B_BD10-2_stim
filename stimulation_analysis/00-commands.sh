#!/usr/bin/env bash

Rscript 01-stim_get_GO.R
Rscript 02-stim_enrich_plot.R

round_num="R1"
nameset="APP_Activity_Dep"


# T41B_BD10-2_stim (data files obtained from parent directory)
mkdir -p $nameset
gzip -c ../res.5.csv > "${nameset}/${round_num}.${nameset}.RNA.Wt-StvsWt-Un.csv.gz"
gzip -c ../res.6.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-StvsTg-Un.csv.gz"
gzip -c ../res.7.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-D-StvsTg-D-Un.csv.gz"

# analysis names will be screwed up for the third, which will be "Wt_APPL/S" instead of "APPL/S + BD10-2"
biodomain_correlation.R \
  -n $nameset \
  -g Wt \
  -d APPL/S \
  -c Gene_id,log2FoldChange \
  -a Wt-StvsWt-Un,Tg-StvsTg-Un,Tg-D-StvsTg-D-Un

# analysis names will be screwed up for the third, which will be "Wt_APPL/S" instead of "APPL/S + BD10-2"
biodomain_enrichment.R \
  -n $nameset \
  -g Wt \
  -d APPL/S \
  -c Gene_id,log2FoldChange,padj \
  -a Wt-StvsWt-Un,Tg-StvsTg-Un,Tg-D-StvsTg-D-Un

