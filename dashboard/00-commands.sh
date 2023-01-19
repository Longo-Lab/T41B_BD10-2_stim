#!/usr/bin/env bash

# round_num="R1"
# nameset="T41B_BD10-2_stim"

# T41B_BD10-2_stim (data files obtained from parent directory)
# gzip -c ../res.2.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-VvsWt-V.csv.gz"
# gzip -c ../res.12.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-DvsWt-V.csv.gz"
# gzip -c ../res.9.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-DvsTg-V.csv.gz"

# biodomain_correlation.R
sbatch -J T41B_BD10-2_stim --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_correlation_%x.log \
  --wrap "ml R/4.0; biodomain_correlation.R -n T41B_BD10-2_stim -g APPL/S -d BD10-2 -c Gene_id,log2FoldChange"

# biodomain_enrichment.R
sbatch -J T41B_BD10-2_stim --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_enrichment_%x.log \
  --wrap "ml R/4.0; biodomain_enrichment.R -n T41B_BD10-2_stim -g APPL/S -d BD10-2 -c Gene_id,log2FoldChange,padj"


# round_num="R1"
# nameset="T41B_BD10-2_unstim"

# T41B_BD10-2_unstim (data files obtained from parent directory)
# gzip -c ../res.1.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-VvsWt-V.csv.gz"
# gzip -c ../res.11.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-DvsWt-V.csv.gz"
# gzip -c ../res.8.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-DvsTg-V.csv.gz"

# biodomain_correlation.R
sbatch -J T41B_BD10-2_unstim --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_correlation_%x.log \
  --wrap "ml R/4.0; biodomain_correlation.R -n T41B_BD10-2_unstim -g APPL/S -d BD10-2 -c Gene_id,log2FoldChange"

# biodomain_enrichment.R
sbatch -J T41B_BD10-2_unstim --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_enrichment_%x.log \
  --wrap "ml R/4.0; biodomain_enrichment.R -n T41B_BD10-2_unstim -g APPL/S -d BD10-2 -c Gene_id,log2FoldChange,padj"


# round_num="R1"
# nameset="T41B_stim"

# T41B_stim (data files obtained from parent directory)
# gzip -c ../res.1.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-VvsWt-V.csv.gz"
# gzip -c ../res.14.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-DvsWt-V.csv.gz"
# gzip -c ../res.6.csv > "${nameset}/${round_num}.${nameset}.RNA.Tg-DvsTg-V.csv.gz"
# gzip -c ../res.5.csv > "${nameset}/${round_num}.${nameset}.RNA.Wt-DvsWt-V.csv.gz"

# biodomain_correlation.R
sbatch -J T41B_stim --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_correlation_%x.log \
  --wrap "ml R/4.0; biodomain_correlation.R -n T41B_stim -g APPL/S -d TBS -c Gene_id,log2FoldChange"

# biodomain_enrichment.R
sbatch -J T41B_stim --mem=5G -c 2 -t 01:00:00 -p interactive \
  -o %x/%A_biodomain_enrichment_%x.log \
  --wrap "ml R/4.0; biodomain_enrichment.R -n T41B_stim -g APPL/S -d TBS -c Gene_id,log2FoldChange,padj"
