source("data_load.R")


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "www.ensembl.org")

# annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name","ensembl_gene_id_version", "ensembl_transcript_id_version", "description", "hgnc_symbol", "entrezgene_id", "start_position", "end_position"))
annotations2 <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"))

annotations <- dplyr::mutate(annotations2, gene_length = end_position - start_position)

ensembl_list = gsub("[.][0-9]*", "", rownames(counts))
# Filter and re-order gene.annotations to match the order in your input genes list
final.genes <- annotations %>% dplyr::filter(ensembl_gene_id %in% ensembl_list)
final.genes <- final.genes[order(match(final.genes$ensembl_gene_id, ensembl_list)),]; rownames(final.genes) <-NULL

head(final.genes)

keep = rowSums(counts) > 0
table(keep)
# emat = Ens2Symbols(emat = counts[keep,], mapping = mapping)
# saveRDS(emat, "../data_out/raw_counts_geneSymb.rds")
emat = readRDS("../data_out/raw_counts_geneSymb.rds")

final.genes <- annotations %>% dplyr::filter(external_gene_name %in% rownames(emat))
final.genes <- final.genes[order(match(final.genes$external_gene_name, rownames(emat))),]; rownames(final.genes) <-NULL
final.genes = final.genes[!duplicated(final.genes$external_gene_name),]
rownames(final.genes) = final.genes$external_gene_name

emat = emat / final.genes$gene_length
tpm_counts <- t( t(emat) * 1e6 / colSums(emat) )


# Permorm cell type enrichment analysis ----
xcell = xCellAnalysis(expr = tpm_counts,
                      rnaseq = TRUE,
                      parallel.sz = 3,
                      parallel.type = "FORK")
pval = xCellSignifcanceBetaDist(scores = xcell, rnaseq = TRUE)
xcell = xcell[!rownames(xcell) %in% c("ImmuneScore", "StromaScore", "MicroenvironmentScore"),]
colnames(pval) = colnames(xcell)

write.csv(xcell, file = "../data_out/xcell.csv")
write.csv(pval, file = "../data_out/xcell_pval.csv")
