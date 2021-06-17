library("tidyverse")

message("Loading the data...")
meta <- read.csv(file = "../data_out/metadata.csv", row.names = 1, stringsAsFactors = F)
counts <- read.csv(file = "../data_out/raw_counts.csv", row.names = 1) %>% as.matrix()

message("Removing outlier RAZZ_247_P3...\n")
meta <- meta[!rownames(meta) %in% "RAZZ_247_P3",]
counts <- counts[,rownames(meta)]


# mapping = t2g()
# write.csv(mapping, file = "../data_in/Ens_gene_mapping.csv")
mapping <- read.csv(file = "../data_in/Ens_gene_mapping.csv", stringsAsFactors = F, row.names = 1)