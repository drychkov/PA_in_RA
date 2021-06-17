library("limma")
library("multcomp")
library('dplyr')
library("DESeq2")
library("ggplot2")
library("IHW")
library("BiocParallel")
library('purrr')
source("plot_functions.R")

xcell <- read.csv(file = "../data_out/xcell.csv", row.names = 1)
pval <- read.csv(file = "../data_out/xcell_pval.csv", row.names = 1)

meta <- read.csv(file = "../data_out/metadata.csv", row.names = 1, stringsAsFactors = F)
meta <- meta[meta$timepoint =='P1',]



flevels <- c("zero", "one", "two", "three")
meta[, c(3:7, 9, 12, 13, 17, 18,19)] <- meta %>%
    dplyr::select(c(3:7, 9, 12, 13, 17, 18,19)) %>%
    map(~ factor(., levels = flevels)) %>%
    map(~droplevels(.)) %>% 
    as.data.frame()


# meta$timepoint = as.factor(meta$timepoint)
meta$race = as.factor(meta$race)
meta$age = (meta$age - min(meta$age))/(max(meta$age) - min(meta$age)) 
# Analysis of xCell results for moderate activity at P1 -----
xcell_P1 = xcell[,rownames(meta[meta$timepoint == 'P1', ])]
pval_P1 = pval[,rownames(meta[meta$timepoint == 'P1', ])]
xcell_P1 <- as.matrix(xcell_P1)
## Identify surrogate variables using SVA ----
mod = model.matrix(~female + age + race + hispanic + activity_moderate_tertile_v1, meta[meta$timepoint == 'P1', ])
caret::findLinearCombos(mod)

#mod = mod[,-11]
colnames(mod)
# add SV1 abd do remove batch effect
mod0 <- model.matrix(~female + age + race + hispanic , meta[meta$timepoint == 'P1', ])
#num_sv <-  sva::num.sv(dat = xcell_P1, mod = mod, method = "be")
#svs <- sva::sva(dat = xcell_P1, mod = mod, mod0 = mod0, n.sv = num_sv)
#surVars = svs$sv
num_sv <-  sva::num.sv(dat = as.matrix(xcell_P1), mod = mod, method = "be")
svs <- sva::sva(dat = as.matrix(xcell_P1), mod = mod, mod0 = mod0, n.sv = num_sv)
surVars = svs$sv
colnames(surVars) = paste0("SV", 1:ncol(surVars))
meta2 = cbind(meta, surVars)
# 
# 
# mod = model.matrix(formula(paste0("~", paste(c(colnames(surVars), "female", "age", "race", "hispanic", "activity_sedentary_tertile_v1"), collapse = "+"))), meta2[meta2$timepoint == 'P1', ])
caret::findLinearCombos(mod)
# # mod = mod[,-11]

# Analyze cell by cell removing samples that have insignificant cell detection p-value
keep = list()
p.thres = 0.2
for (group in c("one",'two', "three")) {
    smp = meta[meta$timepoint == 'P1' & meta$activity_moderate_tertile_v1 == group,] %>% rownames()
    k <- rowSums(pval_P1[,smp] < p.thres) > length(smp) * 0.3
    keep[[group]] <- rownames(xcell_P1)[k]
}



keepCells = Reduce(intersect, keep)
length(keepCells)

# cell="B-cells"
# cell = "CD8+ naive T-cells"
p.thres = 0.2 # Threshold for cell detection
xcell.stat_byCells = list() #rownames(xcell_P1)
for (cell in keepCells){
    print(cell)
    # if (sum(pval_P1[cell,] < p.thres) < ncol(xcell_P1)*0.7) next
    xcelldata = xcell_P1[cell,]#[cell,pval[cell,] < p.thres]
    #xcelldata = (xcelldata - min(xcelldata)) / (max(xcelldata) - min(xcelldata))
    # might be cooment it later for comparing
    mod1 = cbind(mod, surVars)
    xcelldata = unlist(xcelldata )
    s.weights = 1*(pval_P1[cell,]< p.thres)[cell,]
    q = lm(xcelldata ~ mod1 + 0, weights = s.weights, data = meta2[meta2$timepoint == 'P1', ])
    # summary(q)
    
    K <- matrix(c(0,0,0,0,0,0,0,0,-1,0,1,0)[!summary(q)$aliased],1)
    # K <- matrix(c(0,0,0,0,0,0,0,0,0,-1,0,1)[!summary(q)$aliased],1)
    if (sum(K)) next
    if (any(summary(q)$aliased)) {
        mod1 = mod1[, !summary(q)$aliased]
        q = lm(xcelldata ~ mod1 + 0, weights = s.weights) #1*(pval_P1[cell,]< p.thres))
    }
    
    tt <- glht(q, linfct = K)
    if (!coef(tt) | !sum(residuals(q))) next
    xcell.stat_byCells[[cell]] = data.frame(
        cell = cell,
        mean.three = coef(q)["mod1activity_moderate_tertile_v1three"] + coef(q)[1], 
        mean.one = coef(q)["mod1activity_moderate_tertile_v1one"] + coef(q)[1],
        FC = (coef(q)["mod1activity_moderate_tertile_v1three"] + coef(q)[1])/ 
            (coef(q)["mod1activity_moderate_tertile_v1one"] + coef(q)[1]), 
        pvalue = summary(tt)$test$pvalues[1],
        coef = summary(tt)$test$coefficients %>% unname()
    )
    print(colnames(mod1))
    
}
xcell.stat  = do.call(rbind, xcell.stat_byCells) %>% as.data.frame() 
xcell.stat$padj = p.adjust(p = xcell.stat$pvalue, method = "BH", n = 26)
# version fpr the heatmap
xcell_P1_0.1 = xcell_P1[rownames(xcell.stat[xcell.stat$pvalue < 0.1,]),]
for (cell in rownames(xcell_P1_0.1)){
    xcelldata_ = xcell_P1_0.1[cell,]
    xcell_P1_0.1[cell,] = (xcelldata_ - min(xcelldata_ )) / (max(xcelldata_ ) - min(xcelldata_ )) 
}
                     
plotHeatmap(xcell_P1_0.1, metaChar[metaChar$timepoint == "P1" & metaChar$activity_moderate_tertile_v1 %in% c("one", "three"),], classes = "activity_moderate_tertile_v1", rowNames = T, colNames = F, method = "ward.D", scale = "none", colorPalette = 1)



