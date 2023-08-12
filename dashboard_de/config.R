geno <- 'APPL/S'
drug <- 'BD10-2'
de_cols <- c('Gene_id', 'log2FoldChange', 'pvalue', 'padj')
de_names <- c('Log2FC', 'Pval', 'Pval (adj)')
pval_col <- 'pvalue'
lfc_unshrunken_col <- 'log2FoldChange'
lfc_cutoff_geno <- 2.5
lfc_cutoff_drug <- 1
is_sc <- F
