#!/usr/bin/env Rscript

library(data.table)
library(biomaRt)

# directory locations
rootdir = '/labs/flongo/2020_T41B_BD10-2_Stimulation'
wk_dir = paste0(rootdir, '/ad_module_enrichment')
today = format(Sys.Date(), '%Y%m%d')
nameset = "_ad_module_"


###############FUNCTIONS
# Basic function to convert human to mouse gene names ####convert to human
convertGeneList = function(x){
  require("biomaRt")
  mrt = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version="102")
  # match with ensembl_gene_id across both
  genes = getBM(mart=mrt, filters="ensembl_gene_id", values=x, 
                attributes=c("ensembl_gene_id", 
                             "mmusculus_homolog_ensembl_gene",
                             "mmusculus_homolog_orthology_type",
                             "mmusculus_homolog_perc_id",
                             "mmusculus_homolog_perc_id_r1",
                             "mmusculus_homolog_goc_score",
                             "mmusculus_homolog_wga_coverage",
                             "mmusculus_homolog_orthology_confidence"))
  # Print the first 6 genes found to the screen
  print(head(unique(genes)))
  return(unique(genes))
}


###############MAIN
dt = fread(paste0(wk_dir, "/AD_consensus_clusers.tsv"))

# look up consensus clusters for coloring
# from figure S1 legend in Wan et al
string = list("A"=c("TCXblue", "PHGyellow", "IFGyellow"),
              "B"=c("FPturquoise", "TCXturquoise", "STGblue", "PHGturquoise",
                    "IFGturquoise", "DLPFCblue", "CBEturquoise"),
              "C"=c("CBEyellow", "PHGbrown", "TCXgreen", "DLPFCyellow",
                    "STGbrown", "FPyellow", "IFGbrown"),
              "D"=c("CBEbrown", "DLPFCbrown", "PHGgreen", "STGyellow", "FPblue",
                    "TCXyellow", "IFGblue"),
              "E"=c("FPbrown", "TCXbrown", "STGturquoise", "DLPFCturquoise",
                    "PHGblue", "CBEblue"))

# map consensus cluster to dt
ccs = reshape2::melt(string) 
setnames(ccs, c("Module", "Consensus_Cluster"))
dt = merge(dt, ccs, sort=F)

# map mouse ENSM to dt
converts = data.table(convertGeneList(unique(dt$GeneID)))
# sort one to many duplicates by Mouse Gene-order conservation score
converts = converts[order(-mmusculus_homolog_goc_score)]
converts = converts[!duplicated(ensembl_gene_id)]
# rename and merge
setnames(converts, new=c("GeneID", "Mmus_GeneID"), 
         old=c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"))
dt = merge(dt, converts[, .(GeneID, Mmus_GeneID)], sort=F, all.x=T)

# print to file
fwrite(dt, file=paste0(wk_dir, "/AD_modules_full.tsv"))

