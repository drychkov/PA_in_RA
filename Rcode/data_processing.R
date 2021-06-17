library("tidyverse")


# Metadata preparation ----
## Short metadata -----
meta <- xlsx::read.xlsx(file = "../data_in/RAZZ-metadata_12.5.20.xlsx", sheetIndex = 1)
rownames(meta) <- paste0(meta$ID, "_", meta$timepoint)

meta2 = xlsx::read.xlsx(file = "../data_in/RAZZ_metadata_12-22-2020_super_simplified.xlsx", sheetIndex = 1)
meta2$sample = paste0("RAZZ_", meta2$subjectid)
meta2 = meta2[meta2$sample %in% meta$ID,]
rownames(meta2) = paste0(meta2$sample, "_P1")
all(meta2$activity_moderate_tertile_v1 == meta[rownames(meta2), "activity_moderate_tertile_v1"])
all(meta2$activity_sedentary_tertile_v1 == meta[rownames(meta2), "activity_sedentary_tertile_v1"])
all(meta2$age == meta[rownames(meta2), "age"])
all(meta2$female == meta[rownames(meta2), "female"])
all(meta2$race == meta[rownames(meta2), "race"])

tmp <- setNames(object = meta2$hispanic, meta2$subjectid)
# all(recode(.x = meta2$subjectid, !!!hisp) == meta2$hispanic)
meta$hispanic <- recode(.x = gsub("RAZZ_", "", meta$ID), !!!tmp)
tmp <- setNames(object = meta2$biologicIS, meta2$subjectid)
meta$biologicIS <- recode(.x = gsub("RAZZ_", "", meta$ID), !!!tmp)


keys <- c(`0` = "zero", `1` = "one", `2` = "two", `3` = "three")
meta[, c(3:7, 9, 12, 13, 17, 18)] <- meta %>%
    dplyr::select(c(3:7, 9, 12, 13, 17, 18)) %>%
    map(~ recode(., !!!keys)) %>%
    as.data.frame()

write.csv(meta, file = "../data_out/metadata.csv")
# meta = read.csv("~/Box/3. RAZZ Transcriptomic Study/Dmitry/data_out/metadata.csv", row.names = 1, stringsAsFactors = F)

rm(tmp, meta2)

## Long metadata ----
metaFull <- openxlsx::read.xlsx(xlsxFile = "../data_in/RAZZ_metadata_12-2-2020_complete.xlsx", na.strings = ".r")


##### 
metaFull$DMARDv1 = 1*(rowSums(metaFull[,c("methonowv1", "plaquenilnowv1", "leflunonowv1", "sulfasalazinenowv1", "imurannowv1")], na.rm=T) >= 1)
metaFull$DMARDv2 = 1*(rowSums(metaFull[,c("methonowv2", "plaquenilnowv2", "leflunonowv2", "sulfasalazinenowv2", "imurannowv2")], na.rm=T) >= 1)
metaFull$DMARDv3 = 1*(rowSums(metaFull[,c("methonowv3", "plaquenilnowv3", "leflunonowv3", "sulfasalazinenowv3", "imurannowv3")], na.rm=T) >= 1)

metaFull$biologicISv1 = 1*(rowSums(metaFull[,c("enbrelnowv1", "remicadenowv1", "humiranowv1", "kineretnowv1", "orencianowv1", "simponinowv1", "cimzianowv1", "xeljanznowv1", "actemranowv1", "rituxanv1")], na.rm=T) >= 1)
metaFull$biologicISv2 = 1*(rowSums(metaFull[,c("enbrelnowv2", "remicadenowv2", "humiranowv2", "kineretnowv2", "orencianowv2", "simponinowv2", "cimzianowv2", "xeljanznowv2", "actemranowv2", "rituxanv2")], na.rm=T) >= 1)
metaFull$biologicISv3 = 1*(rowSums(metaFull[,c("enbrelnowv3", "remicadenowv3", "humiranowv3", "kineretnowv3", "orencianowv3", "simponinowv3", "cimzianowv3", "xeljanznowv3", "actemranowv3", "rituxanv3")], na.rm=T) >= 1)

metaFull$ID = paste0("RAZZ_", metaFull$subjectid)
rownames(metaFull) = metaFull$ID

cbind(meta[meta$timepoint == "P1", "biologicIS"],
metaFull[match(meta[meta$timepoint == "P1", "ID"], metaFull$ID), "biologicISv1"]) %>% View()

smp = meta[is.na(meta[, "biologicIS"]), "ID"] %>% unique() %>% sort()
meta[as.vector(outer(smp, c("P1", "P2", "P3"), paste, sep="_")),]

meta[paste0(smp, "_", "P1"), "biologicIS"] = ifelse(metaFull[smp, "biologicISv1"] == 1, "one", "zero")
meta[paste0(smp, "_", "P2"), "biologicIS"] = ifelse(metaFull[smp, "biologicISv2"] == 1, "one", "zero")
meta[paste0(smp, "_", "P3"), "biologicIS"] = ifelse(metaFull[smp, "biologicISv3"] == 1, "one", "zero")
meta = meta[complete.cases(meta),]

meta[meta$timepoint == "P1", "DMARD"] = ifelse(metaFull[match(meta[meta$timepoint == "P1", "ID"], metaFull$ID), "DMARDv1"] == 1, "one", "zero")
meta[meta$timepoint == "P2", "DMARD"] = ifelse(metaFull[match(meta[meta$timepoint == "P2", "ID"], metaFull$ID), "DMARDv2"] == 1, "one", "zero")
meta[meta$timepoint == "P3", "DMARD"] = ifelse(metaFull[match(meta[meta$timepoint == "P3", "ID"], metaFull$ID), "DMARDv3"] == 1, "one", "zero")

write.csv(meta, file = "../data_out/metadata.csv")

###
tmp = colnames(metaFull)[grepl("v1", colnames(metaFull))]
tmp = gsub("[_]?v1", "", tmp)
metaFullx = data.frame(matrix(nrow = nrow(metaFull)*3, ncol = 71))
colnames(metaFullx) = c("ID", "timepoint", gsub("[_]?v1", "", tmp))
metaFullx$ID = paste0("RAZZ_", metaFull$subjectid)
metaFullx$timepoint = rep(c("P1", "P2", "P3"), each = nrow(metaFull))

metaFullx[metaFullx$timepoint == 'P1', gsub("[_]?v1", "", tmp)] = metaFull[,tmp]
tmp = colnames(metaFull)[grepl("v2", colnames(metaFull))]
metaFullx[metaFullx$timepoint == 'P2', gsub("[_]?v2", "", tmp)] = metaFull[,tmp]
tmp = colnames(metaFull)[grepl("v3", colnames(metaFull))]
metaFullx[metaFullx$timepoint == 'P3', gsub("[_]?v3", "", tmp)] = metaFull[,tmp]


write.csv(metaFullx, file = "../data_out/metadataFull_tidy.csv")
write.csv(metaFull, file = "../data_out/metadataFull.csv")

# Reading and combining raw data from STAR outputs -----
files <- list.files(path = "../data_in/counts/", pattern = "*reads_per_gene.star.tab$", full.names = TRUE)
counts.files <- lapply(files, read.table, skip = 4)
counts <- as.data.frame(sapply(counts.files, function(x) x[, 2]))

files <- gsub("../data_in/counts/", "", files)
files %in% meta$idseq_filename
counts <- counts[, match(meta$idseq_filename, files)]
colnames(counts) <- rownames(meta)
row.names(counts) <- counts.files[[1]]$V1

rm(counts.files, files)

write.csv(counts, file = "../data_out/raw_counts.csv")