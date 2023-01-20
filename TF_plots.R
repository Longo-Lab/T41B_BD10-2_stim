library(data.table)
library(ggplot2)
a = fread("gost.2.csv")
b = fread("gost.12.csv")
c = merge(a[source == "TF", .(p_value, term_id, term_name, DGE)], b[source == "TF", .(p_value, term_id, term_name, DGE)], by=c("term_id", "term_name", "DGE"), suffix=c(".T41B", ".T41B.BD10-2"))
c
c[term_name %like% "Sp1"]
c[term_name %like% "S[pP]1"]
c[term_name %like% "S[pP]1[^0]"]
c[, sp := ifelse(term_name %like% "S[pP]1[^0]", "SP1", "other")]
c[-log10(p_value.T41B) > 10 | -log10(`p_value.T41B.BD10-2`) > 30]
c[(-log10(`p_value.T41B.BD10-2`)) - (-log10(`p_value.T41B`)) > 5]
c[(-log10(`p_value.T41B.BD10-2`)) > 40]
c[term_name %like% "ZF5", sp := "ZF5"]
c[term_name %like% "GKLF", sp := "GKLF"]
c[term_name %like% "E2F", sp := "E2F"]
c[term_name %like% "Pax"]
c[term_name %like% "Spi-B" | term_name %like% "SPIB" | term_name %like% "PU\\.1", sp := "PU.1"]
c[(-log10(`p_value.T41B.BD10-2`)) > 30]
c[term_name %like% "ZIC" | term_name %like% "Zic", sp := "Zic"]
c[sp == "other", sp := NA]
pdf(file="TF_plots.pdf", onefile=T, paper="USr")
ggplot(c[DGE == "upregulated"][order(sp)], aes(x=-log10(p_value.T41B), y=-log10(`p_value.T41B.BD10-2`), col=sp)) + geom_point() + ggtitle("Upregulated TF motifs")
ggplot(c[DGE == "downregulated"][order(sp)], aes(x=-log10(p_value.T41B), y=-log10(`p_value.T41B.BD10-2`), col=sp)) + geom_point() + ggtitle("Downregulated TF motifs")
dev.off()
