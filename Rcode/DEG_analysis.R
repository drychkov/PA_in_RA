# load environments
library("DESeq2")
library("ggplot2")
library("IHW")
library("BiocParallel")

source("plot_functions.R")

# Data load ----
source("data_load.R")

# ----------------
# Focus on P1
meta <- meta[meta$timepoint =='P1',]

counts <- counts[,rownames(meta)]
metaChar = meta

flevels <- c("zero", "one", "two", "three")
meta[, c(3:7, 9, 12, 13, 17, 18)] <- meta %>%
    dplyr::select(c(3:7, 9, 12, 13, 17, 18)) %>%
    map(~ factor(., levels = flevels)) %>%
    map(~droplevels(.)) %>% 
    as.data.frame()


# meta$timepoint = as.factor(meta$timepoint)
meta$race = as.factor(meta$race)
meta$age = (meta$age - min(meta$age))/(max(meta$age) - min(meta$age)) # normalize the age
metaTmp = meta
countsTmp = counts[, rownames(metaTmp)]

# DE analysis ----

## Design the analysis for activity_moderate_tertile_v1 -----
dds_AS <- DESeqDataSetFromMatrix(
    countData = countsTmp,
    colData = metaTmp,
    design = ~ 1#female + age + race + hispanic + biologicIS + activity_sedentary_tertile_v1
)

mcols(dds_AS) = DataFrame(
    ensembl_gene_id_version = rownames(dds_AS), 
    mapping[match(gsub("[.][0-9]*", "", rownames(dds_AS)), mapping$ensembl_gene_id), 
            c("hgnc_symbol", "entrezgene_id", "description")], 
    row.names = rownames(dds_AS)
)

# Build the model
mod = model.matrix(~female + age + race + hispanic + activity_moderate_tertile_v1, metaTmp)

# Checking and removing columns that are unnecessary (hispanic = 2 and biologicIS = 2) and lead to redundant model matrix 
caret::findLinearCombos(mod)  # Harry does not know what this is mean


# Each group should have at least 30% of samples expressing a gene
keep = list()
for (group in c("one", "two", "three")) {
    smp = metaTmp[metaTmp$activity_moderate_tertile_v1 == group,] %>% rownames()
    k <- rowSums(counts(dds_AS[,smp]) > 0) > length(smp) * 0.3
    keep[[group]] <- rownames(counts(dds_AS))[k]
    
}

keep = Reduce(intersect, keep)
length(keep)
dds_AS <- dds_AS[keep, ]

## Batch Correction with SVA ----
dds_AS <- estimateSizeFactors(dds_AS)
dat  <- counts(dds_AS, normalized = TRUE)

mod0 <- model.matrix(~female + age + race + hispanic, colData(dds_AS))


num_sv = sva::num.sv(dat = dat, mod = mod, method = "be")
svseq <- sva::svaseq(dat = dat, mod = mod, mod0 = mod0, n.sv = num_sv)

surVars = svseq$sv
colnames(surVars) = paste0("SV", 1:ncol(surVars))
colData(dds_AS) = cbind(colData(dds_AS), surVars)


design(dds_AS) <- formula(paste0("~", paste(c(colnames(surVars), "female", "age", "race", "hispanic",   "activity_moderate_tertile_v1"), collapse = "+")))

## Run the analysis -----
dds_AS <- DESeq(dds_AS, 
                full = cbind(mod,surVars),
                betaPrior = FALSE, 
                parallel = T, 
                BPPARAM = MulticoreParam(workers = 3)
)

## Extract results -----

### Three vs One -----
library(IHW)
res <- results(
    dds_AS, 
    contrast = list("activity_moderate_tertile_v1three",
                    "activity_moderate_tertile_v1one"),
    filterFun = ihw
)
res <- res[order(-res$log2FoldChange), ]
sum(res$padj < 0.1, na.rm = TRUE)
summary(res)

resAnno_1vs3 = data.frame(
    as.data.frame(subset(res)), 
    mapping[match(gsub("[.][0-9]*", "", rownames(res)), mapping$ensembl_gene_id),
            c("hgnc_symbol", "entrezgene_id", "description")]
)
resSig_1vs3 <- subset(res, padj < 0.1) %>% rownames()
resFilter_1vs3 = data.frame(
    as.data.frame(subset(res, padj < 0.1)), 
    mapping[match(gsub("[.][0-9]*", "", resSig_1vs3), mapping$ensembl_gene_id), 
            c("hgnc_symbol", "entrezgene_id", "description")]
)

resSig_1vs2 <- subset(res_1vs2, padj < 0.1) %>% rownames()
resFilter_1vs2 = data.frame(
    as.data.frame(subset(res_1vs2, padj < 0.1)), 
    mapping[match(gsub("[.][0-9]*", "", resSig_1vs2), mapping$ensembl_gene_id), 
            c("hgnc_symbol", "entrezgene_id", "description")]
)

saveRDS(dds_AS, file = "../data_out/ddsAS_P1_activity_moderate.rds")
saveRDS(res, file = "../data_out/res_P1_activity_moderate_3vs1.rds")

dds_AS <- readRDS("../data_out/ddsAS_P1_activity_moderate.rds")
res <- readRDS("../data_out/res_P1_activity_moderate_3vs1.rds")

res <- results(
    dds_AS, 
    contrast = list("activity_moderate_tertile_v1three",
                    "activity_moderate_tertile_v1one"),
    filterFun = ihw
)
# compare level one and level two
res_1vs2 <- results(
    dds_AS, 
    contrast = list("activity_moderate_tertile_v1two",
                    "activity_moderate_tertile_v1one"),
    filterFun = ihw
)

saveRDS(res_1vs2, file = "../data_out/res_P1_activity_moderate_2vs1.rds")
#### Plots -----

vsd <- vst(dds_AS, blind = T, nsub = 2000, fitType = "parametric")
#plotPCA final version



plotPCA(vsd[resSig_1vs3, rownames(metaChar[metaChar$activity_moderate_tertile_v1 %in% c('one','three'),])], intgroup = "activity_moderate_tertile_v1")+ 
    theme_bw() +  scale_colour_discrete("Activity Level")+
    scale_color_manual(labels = c("one", "three"), values = c("blue", "red")) +
    geom_point(size = 5) +
    coord_fixed(ratio = 1) +theme(axis.title.x = element_text(size = 10),axis.title.y = element_text(size = 10),
                                  axis.text.x = element_text(face="bold" ,size=7),
                                  axis.text.y = element_text(face="bold",size=7),
                                  legend.title = element_text( size=7,  face="bold"),
                                  legend.text = element_text( size=7, 
                                                              face="bold"))

# Pathways analysis
Upreg <- resFilter_1vs3[resFilter_1vs3$log2FoldChange > 0, "hgnc_symbol"]
write.table(Upreg,'../data_out/Upreg.txt', col.names = F,row.names = F, quote = F)
Downreg <- resFilter_1vs3[resFilter_1vs3$log2FoldChange < 0, "hgnc_symbol"]
write.table(Downreg,'../data_out/Downreg.txt', col.names = F,row.names = F, quote = F)

# Make figures
Up <- read.table('../data_out/BioPlanet_2019_table_Up.txt', header = T, sep ='\t')
Up$'-log10(p)'<- -log10(Up$P.value)
Up$P.value <- Up$'-log10(p)'
Up <- Up[order(-Up$P.value),]
ggplot(head(Up,5), aes(x=Term, y=P.value,label = Term, fill = Term) ) + 
    geom_bar(stat = "identity", width=0.7)   + geom_text(size=4.5,hjust=1) +
    labs(x = "Pathways", y = "-log10(p-value)") +
    scale_y_continuous(expand=c(0,0)) + coord_flip() +#labs(x = "") +
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.y = element_blank(),axis.text.x =  element_text(size = 10) ,
          axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10))

Down <- read.table('../data_out/BioPlanet_2019_table_Down.txt', header = T, sep ='\t')
Down$'-log10(p)'<- -log10(Down$P.value)
Down$P.value <- Down$'-log10(p)'
Down <- Down[order(-Down$P.value),]
Down$P.value <- -(Down$P.value)
Merge <- rbind(head(Up,5), head(Down,5))
Merge <- Merge[order(Merge$P.value),]
positions <- Merge$Term
ggplot(Merge, aes(x=Term, y=P.value,label = Term, fill = Term) ) + 
    geom_bar(stat = "identity", width=0.7)   + geom_text(size=4.5,hjust=1) +
    labs(x = "Pathways", y = "-log10(p-value)") +
    scale_x_discrete(limits = positions) +
    coord_flip()  +
    theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),axis.text.x =  element_text(size = 10) ,
          axis.title.y = element_blank(),axis.title.x = element_text(size = 10))

immune_genes_1 <- union(strsplit(Down$Genes[1],';')[[1]],strsplit(Down$Genes[2],';')[[1]])#,
immune_genes_2 <- union(strsplit(Down$Genes[3],';')[[1]], strsplit(Down$Genes[4],';')[[1]])
immune_genes <- union(immune_genes_1, immune_genes_2)

translation_genes <- strsplit(Up$Genes[3],';')[[1]]
#scale_y_continuous(expand=c(0,0)) ++ expand_limits(x = 0， y =0)
# Volcano plot
library("EnhancedVolcano")

volcanoLabs = resAnno_1vs3$hgnc_symbol
#volcanoLabs[abs(resAnno_1vs3$log2FoldChange) < 4 & resAnno_1vs3$padj > 0.01] = ""
volcanoLabs[!(resAnno_1vs3$hgnc_symbol%in% immune_genes | resAnno_1vs3$hgnc_symbol%in% translation_genes)] = ""

EnhancedVolcano(res,
                lab = volcanoLabs,
                pCutoff = 0.1,
                x = 'log2FoldChange',
                y = 'padj',
                title =  '',
                labSize = 6,
                ylim = c(0, 2.3),
                colAlpha = 1,
                legendPosition = "bottom",
                drawConnectors = T,
                legendLabSize = 6,
                legendIconSize = 2) 

EnhancedVolcano(res,
                lab = volcanoLabs,
                pCutoff = 0.1,
                x = 'log2FoldChange',
                y = 'padj',
                title =  '',
                labSize = 4,
                colAlpha = 1,
                legendPosition = "bottom",
                drawConnectors = T,
                legendLabSize = 6,
                col = c('seashell3','aquamarine1', 'chocolate1','tomato1'),
                ylim = c(0, 3.5),
                legendIconSize = 2) +
    theme(axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
          axis.text.x = element_text(face="bold",size=15),axis.text.y = element_text(face="bold",size=15),
          legend.text = element_text( size=15))


# Vennplot
Upreg_1vs3 <- rownames(resFilter_1vs3[resFilter_1vs3$log2FoldChange > 0, ])
Downreg_1vs3 <- rownames(resFilter_1vs3[resFilter_1vs3$log2FoldChange < 0, ])
Upreg_1vs2 <- rownames(resFilter_1vs2[resFilter_1vs2$log2FoldChange > 0, ])
Downreg_1vs2 <- rownames(resFilter_1vs2[resFilter_1vs2$log2FoldChange <0, ])

genelist <- list(Upreg_1vs3,Upreg_1vs2)
names(genelist) <- c('1vs3','1vs2')

plotVenn(genelist)

genelist <- list(Downreg_1vs3,Downreg_1vs2)
names(genelist) <- c('1vs3','1vs2')
plotVenn(genelist)


# barplot
vsd <- vst(dds_AS, blind = T, nsub = 2000, fitType = "parametric")
vsd_mtx <- assay(vsd)
description <- as.data.frame(resAnno_1vs3[,'hgnc_symbol'])

colnames(description) <- "hgnc_symbol"


vsd_cleanY = cleanY(vsd_mtx, mod, surVars)

immune_response <- c('JAK2','RAF1','IL6R','PTPN6','IL1RN', 'HGF', 'NOD2','MX2', 'XAF1','IFIT1','MX1', 'STAT2','EIF4A3')
immune_id = rownames(resAnno_1vs3[resAnno_1vs3$hgnc_symbol %in%immune_response, ])



vsd_immnue <- vsd_cleanY[immune_id,]

   
immune_genes = as.data.frame(resAnno_1vs3[resAnno_1vs3$hgnc_symbol %in%immune_response, ]$hgnc_symbol)
rownames(immune_genes) <- rownames(resAnno_1vs3[resAnno_1vs3$hgnc_symbol %in%immune_response, ])
colnames(immune_genes) <- 'Gene' 
    
vsd_immnue <- cbind(vsd_immnue, immune_genes) 
   
rownames(vsd_immnue) <- vsd_immnue$Gene
vsd_immnue <- vsd_immnue[,-dim(vsd_immnue)[2]]   

vsd_immnue <- t(vsd_immnue)
vsd_immnue <- as.data.frame(vsd_immnue)
vsd_immnue$Level <- metaChar$activity_moderate_tertile_v1 

p = 0
for (gene in colnames(vsd_immnue)[-14]){
    print(gene)
    if (p == 0){
        df <- as.data.frame(vsd_immnue[,c(gene,'Level')])
        df$Gene <- gene
        colnames(df) <- c('Counts','Level','Gene')
        p = p +1
    }else{
        df_ <- as.data.frame(vsd_immnue[,c(gene,'Level')])
        df_$Gene <- gene
        colnames(df_) <- c('Counts','Level','Gene')
        df <- rbind(df,df_)
    }
}
df = df[df$Level != 'zero',]


# plot1
plot_box <- list()
for (gene in c('IL1RN','IFIT1','IL6R','MX1','NOD2','MX2')){
    df <- df[order(df$Level,decreasing = TRUE), ]

    plot_box[[gene]] <-ggboxplot(df[df$Gene == gene,], x = "Level", y = "Counts",
                                 fill = 'Level', palette = c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
        stat_compare_means(label.y = 10)  +
        #scale_color_manual(values=c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
        labs(x=gene, y = "Expression")
    
}
plot_sum <- ggarrange(plot_box[[1]],plot_box[[2]],plot_box[[3]],
                      plot_box[[4]],plot_box[[5]],plot_box[[6]],
                      nrow = 2,ncol=3,  common.legend = TRUE, legend="bottom")
plot_sum


#plot2
plot_box <- list()
for (gene in c('IL1RN','IFIT1','IL6R','MX1','NOD2','MX2')){
    df <- df[order(df$Level,decreasing = TRUE), ]
    
    plot_box[[gene]] <-ggboxplot(df[df$Gene == gene,], x = "Level", y = "Counts",
                                  palette = c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
        geom_jitter(color=c( "#E7B800", "#BC3C29E5","#0072B5E5"), size=0.4, alpha=0.9)+
        stat_compare_means(method = "anova")  +
        #scale_color_manual(values=c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
        labs(x=gene, y = "Expression")
    
}
plot_sum <- ggarrange(plot_box[[1]],plot_box[[2]],plot_box[[3]],
                      plot_box[[4]],plot_box[[5]],plot_box[[6]],
                      nrow = 2,ncol=3,  common.legend = TRUE, legend="bottom")
ggboxplot(df[df$Gene == gene,], x = "Level", y = "Counts",
          fill = 'Level', palette = c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
  #scale_color_manual(values=c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
    labs(x=gene, y = "Expression")  +stat_compare_means(label.y = 10)
ggboxplot(df[df$Gene == gene,], x = "Level", y = "Counts",
          fill = 'Level', palette = c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
    stat_compare_means(method = "anova")  +
      #scale_color_manual(values=c( "#E7B800", "#BC3C29E5","#0072B5E5")) +
    labs(x=gene, y = "Expression") #+  scale_x_discrete(limits = c("one", "two", "three"))





    
