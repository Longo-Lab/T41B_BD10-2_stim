#!/usr/bin/env bash

round_num="R1"
nameset="T41B_stim"

# T41B_stim (data files obtained from parent directory)
mkdir -p $nameset
gzip -c ../res.5.csv > "${nameset}/${round_num}.${nameset}.RNA.Wt-StvsWt-Un.csv.gz"
# gzip -c ../res.X.csv > "${nameset}/${round_num}.${nameset}.RNA.Wt-D-StvsWt-D-Un.csv.gz"
gzip -c ../res.6.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-StvsTg-Un.csv.gz"
gzip -c ../res.7.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-D-StvsTg-D-Un.csv.gz"

# biodomain_correlation.R
sbatch -J $nameset --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_correlation_%x.log \
  --wrap "ml R/4.2.2; biodomain_correlation.R -n $nameset -g APPL/S -d BD10-2 -i Gene_id -t stim"

# biodomain_enrichment.R
sbatch -J $nameset --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_enrichment_%x.log \
  --wrap "ml R/4.2.2; biodomain_enrichment.R -n $nameset -g APPL/S -d BD10-2 -i Gene_id -t stim"

# save_stim_dashboard_files.R
sbatch -J $nameset --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_save_stim_dashboard_files_%x.log \
  --wrap "ml R/4.2.2; save_stim_dashboard_files.R -n $nameset -g APPL/S -d BD10-2 -i Gene_id"
