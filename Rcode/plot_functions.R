library("ggpubr")
library('RColorBrewer')
library('colorspace')
library('pheatmap')
library('annotate')
library("ggsci")
library("scales")


plotHeatmap = function(emat, sampleMeta, colorPalette = 1, rowNames = FALSE, colNames = FALSE, 
                       colCluster = TRUE, rowCluster = TRUE, classes = "class", title = "", 
                       rowFont = 10, fontsize = 10, bars = NA, colSplit = 2, rowSplit = 2, 
                       method = "ward.D", scale = c("row", "none")[1]) {
  emat = emat[, colnames(emat) %in% rownames(sampleMeta)]
  sampleMeta = sampleMeta[colnames(emat),]
  annotation_col = sampleMeta[, classes, drop = FALSE]
  
  
  cl = sort(unique(unlist(sampleMeta[, classes])))
  cols = rep("", length(cl))
  col_palette = pal_startrek("uniform", alpha = 0.65)(7)
  cols = pal_igv("default", alpha = 0.75)(50)[1:length(cl)]
  names(cols) = cl
  
  ann_colors = lapply(classes, function(x) cols[sort(unique(sampleMeta[, x]))])
  names(ann_colors) = classes
  
  colors = colorRampPalette(rev(brewer.pal(n = 11, name = 'RdBu')))(100)
  
  if (class(bars) == "data.frame") {
    for (i in 1:length(bars)) {
      if (colnames(bars)[i]  %in% c("log2_FC", "log10_FC")) {
        cols.extra = diverge.color(start.color = "#008837", end.color = "#7b3294", 
                                   min.value = min(bars[[i]]*20), max.value = max(bars[[i]]*20), 
                                   mid.value = 0, mid.color = "white")
      } else if (colnames(bars)[i] == "correlation") {
        cols.extra = diverge.color(start.color = "#67a9cf", end.color = "#fc8d59", 
                                   min.value = min(bars[[i]]*10), max.value = max(bars[[i]]*10), 
                                   mid.value = 0, mid.color = "white")
      } else if (colnames(bars)[i] == "fold change") {
        cols.extra = diverge.color(start.color = "#008837", end.color = "#7b3294", 
                                   min.value = min(bars[[i]]*10), max.value = max(bars[[i]]*10), 
                                   mid.value = 1, mid.color = "white")
      } else {
        cols.extra = pal_material(palette = c("pink", "purple", "indigo", "cyan", "green", "lime",
                                              "amber", "deep-orange", "brown", "blue-grey")[i],
                                  n = length(bars[[i]]), alpha = 0.65)(length(bars[[i]]))
      }
      
      ann_colors = c(ann_colors,  list(cols.extra))
      names(ann_colors)[length(ann_colors)] = names(bars)[i]
    }
  }
  
  if (scale == "row") {
    breaks = seq(from = -3, to = 3, length.out = 100)
    breaks[length(breaks)] <- max(max(emat),max(breaks))
    breaks[1] <- min(min(emat),min(breaks))
    
  } else {
    breaks = NA
    if (colorPalette == 2) {
      colors = colorRampPalette(brewer.pal(n = 9, name = 'OrRd'))(100)
    } else if (colorPalette == 3) {
      colors = colorRampPalette(brewer.pal(n = 9, name = 'Blues'))(100)
    }
  }
  
  featureVars = sapply(rownames(emat), function(x) var(emat[x,], na.rm = TRUE))
  if (any(featureVars == 0)) {
    message("Zero variance genes: ", names(featureVars[featureVars == 0]), ", were removed.")
  }
  emat = emat[featureVars != 0, ]
  
  pheatmap(emat, 
           color = colors,
           scale = scale,
           breaks = breaks, 
           cluster_rows = rowCluster, 
           cluster_cols = colCluster, 
           show_colnames = colNames, 
           show_rownames = rowNames, 
           border_color = NA, 
           clustering_method = method,
           cutree_rows = rowSplit,
           cutree_cols =  colSplit,
           fontsize_row = rowFont,
           fontsize = fontsize,
           annotation_col = annotation_col,
           annotation_row = bars,
           annotation_colors = ann_colors,
           treeheight_row = 30,
           treeheight_col = 50,
           angle_col = "90",
           # cex = 0.9,
           silent = FALSE)
}

diverge.color <- function(start.color, end.color, min.value, max.value, 
                          mid.value = 0, mid.color = "ivory") {
  # based on ideas from Maureen Kennedy, Nick Povak, and Alina Cansler
  
  # creates a palette for the current session for a divergent-color
  # graphic with a non-symmetric range
  # "cuts" = the number of slices to be made in the range above and below "mid.value"
  
  ramp1 <- colorRampPalette(c(start.color, mid.color))
  ramp2 <- colorRampPalette(c(mid.color, end.color))
  
  # now specify the number of values on either side of "mid.value"
  
  max.breaks <- round(max.value - mid.value)
  min.breaks <- round(mid.value - min.value)
  
  num.breaks <- max(max.breaks,min.breaks)
  
  low.ramp <- ramp1(num.breaks)
  high.ramp <- ramp2(num.breaks)
  
  # now create a combined ramp from the higher values of "low.ramp" and 
  # the lower values of "high.ramp", with the longer one using all values 
  # high.ramp starts at 2 to avoid duplicating zero
  
  myColors <- c(low.ramp[(num.breaks - min.breaks):num.breaks], high.ramp[2:max.breaks])
  
  return(myColors)
}


plotPCA = function(emat, sampleMeta, title = "", interest = "class", alpha = 1, ellipse = TRUE, pointSize = 2, PCobj = NULL) {
  sampleMeta = sampleMeta[rownames(sampleMeta) %in% colnames(emat), ]
  emat = emat[, rownames(sampleMeta)]
  classes = gsub("_[A-Z]", "", sampleMeta[, interest])
  names = factor(classes, levels = sort(unique(classes)))
  
  cols = rep("", length(classes))
  names(cols) = classes
  
  col_palette = pal_startrek("uniform", alpha = alpha)(7)
  cols = rainbow(length(unique(names)))[as.integer(names)]
  names(cols) = classes
  
  emat = emat[complete.cases(emat),]
  a = sapply(rownames(emat), function(x) var(emat[x, ], use = "everything", na.rm = T))
  if (sum(a == 0) > 0) {
    emat = emat[!rownames(emat) %in% names(a[a == 0]),]
  }
  
  if (is.null(PCobj)) {
    res.svd = svd(scale(t(emat)))
  } else {
    res.svd = PCobj
  }

  ggbiplot(res.svd, ellipse = ellipse, circle = TRUE, 
            obs.scale = 1, var.scale = 1,
            var.axes = FALSE, labels = NULL, groups = classes, pointSize = pointSize, alpha = alpha) + 
    scale_fill_manual(name = interest, values = cols) +
    scale_color_manual(name = interest, values = cols) +
    theme(text = element_text(size = 15)) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5)) + 
    guides(size = "none", color = "none") +
    # coord_fixed(1) +
    labs(title = title)

}


t2g = function() {
  library(biomaRt)
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host = "www.ensembl.org")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name","ensembl_gene_id_version", "ensembl_transcript_id_version", "description",
                                       "hgnc_symbol", "entrezgene_id"), mart = mart )
  return(t2g)
}

Ens2Symbols = function(emat, mapping = NA, fromTo = c("ensembl_gene_id", "hgnc_symbol")) 
  {
  library("parallel")
  library("doParallel")
  
  emat = as.matrix(emat)
  if (!is.data.frame(mapping)) {
    mapping = t2g()
  }

  
  # if(fromTo[1] == "ensembl_gene_id_version"){
  #   features = gsub("[.][0-9]", "", rownames(emat))
  #   fromTo[1] = "ensembl_gene_id"
  # } else {
  #   features = rownames(emat)
  # }
  
  rownames(emat) = gsub("[.][0-9]*", "", rownames(emat))
  mapping = mapping[, fromTo]
  mapping = mapping[mapping[, fromTo[2]] != "",]
  commEns = intersect(rownames(emat), mapping[,fromTo[1]])
  mapping = mapping[mapping[, fromTo[1]] %in% commEns, ]
  genesSym = unique(mapping[, fromTo[2]])
  
  exprsByGene = data.frame(matrix(nrow = length(genesSym), 
                       ncol = ncol(emat), 
                       dimnames = list(genesSym, colnames(emat))))
  
  cl <-  makeCluster(3, type = "FORK")
  registerDoParallel(cl)
  exprsByGene = foreach::foreach(geneId = genesSym) %dopar% {
    ens = mapping[mapping[, fromTo[2]] == geneId, fromTo[1]]
    ens = ens[ens %in% rownames(emat)]
    if (length(ens) == 0) {
      next
    }
    exprsTmp = emat[ens,, drop = FALSE]
    if (nrow(exprsTmp) == 1) {
      #exprsByGene[geneId,]
      tmp = exprsTmp[1,]
    } else {
      #exprsByGene[geneId,]
      tmp = rowMeans(t(exprsTmp), na.rm = TRUE)
    }
    tmp = c(gene = geneId, tmp)
  }
  stopCluster(cl)

  q = do.call("rbind", exprsByGene)
  rownames(q) = q[,1]
  exprsByGene = as.matrix(q[,2:ncol(q)])
  exprsByGene  = apply(exprsByGene, 2, as.numeric)
  rownames(exprsByGene) = rownames(q)

  return(exprsByGene)
  # return(mapping)
}

#  Copyright 2011 Vincent Q. Vu.
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                      obs.scale = 1 - scale, var.scale = scale, 
                      groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                      labels = NULL, labels.size = 3, alpha = 1, 
                      var.axes = TRUE, pointSize = 1,
                      circle = FALSE, circle.prob = 0.69, 
                      varname.size = 3, varname.adjust = 1.5, 
                      varname.abbrev = FALSE, ...)
{
    library(ggplot2)
    library(plyr)
    library(scales)
    library(grid)
    
    stopifnot(length(choices) == 2)
    
    # Recover the SVD
    if(class(pcobj) == 'prcomp'){
        nobs.factor <- sqrt(nrow(pcobj$x) - 1)
        d <- pcobj$sdev
        u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- pcobj$rotation
    } else if(inherits(pcobj, 'princomp')) {
        nobs.factor <- sqrt(pcobj$n.obs)
        d <- pcobj$sdev
        u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- pcobj$loadings
    } else if(inherits(pcobj, 'PCA')) {
        nobs.factor <- sqrt(nrow(pcobj$call$X))
        d <- unlist(sqrt(pcobj$eig)[1])
        u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
    } else if(class(pcobj) == 'list') {
        nobs.factor <- sqrt(nrow(pcobj$u) - 1)
        d <- pcobj$d
        u <- sweep(pcobj$u, 2, 1 / nobs.factor, FUN = '*')
        # u <- sweep(u, 2, 1 / (d * nobs.factor), FUN = '*') #predict(pcobj)$x/nobs.factor
        v <- pcobj$v
        d.total <- sum(d^2)
    } else {
        stop('Expected a object of class prcomp, princomp, PCA, or lda')
    }
    
    # Scores
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
    
    # Directions
    v <- sweep(v, 2, d^var.scale, FUN='*')
    df.v <- as.data.frame(v[, choices])
    
    names(df.u) <- c('xvar', 'yvar')
    names(df.v) <- names(df.u)
    
    if(pc.biplot) {
        df.u <- df.u * nobs.factor
    }
    
    # Scale the radius of the correlation circle so that it corresponds to 
    # a data ellipse for the standardized PC scores
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
    
    # Scale directions
    v.scale <- rowSums(v^2)
    df.v <- r * df.v / sqrt(max(v.scale))
    
    # Change the labels for the axes
    if(obs.scale == 0) {
        u.axis.labs <- paste('standardized PC', choices, sep='')
    } else {
        u.axis.labs <- paste('PC', choices, sep='')
    }
    
    # Append the proportion of explained variance to the axis labels
    u.axis.labs <- paste(u.axis.labs, 
                         sprintf('(%0.1f%% explained var.)', 
                                 100 * d[choices]^2/sum(d^2)))
    
    # Score Labels
    if(!is.null(labels)) {
        df.u$labels <- labels
    }
    
    # Grouping variable
    if(!is.null(groups)) {
        df.u$groups <- groups
    }
    
    # Variable Names
    if(varname.abbrev) {
        df.v$varname <- abbreviate(rownames(v))
    } else {
        df.v$varname <- rownames(v)
    }
    
    # Variables for text label placement
    df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
    
    # Base plot
    g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
        xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) #+ coord_fixed(1)
    
    if(var.axes) {
        # Draw circle
        if(circle) 
        {
            theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
            circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
            g <- g + geom_path(data = circle, color = muted('white'), 
                               size = 1/2, alpha = 1/3)
        }
        
        # Draw directions
        g <- g +
            geom_segment(data = df.v,
                         aes(x = 0, y = 0, xend = xvar, yend = yvar),
                         arrow = arrow(length = unit(1/2, 'picas')), 
                         color = muted('red'))
    }
    
    # Draw either labels or points
    if(!is.null(df.u$labels)) {
        if(!is.null(df.u$groups)) {
            g <- g + geom_text(aes(label = labels, color = groups), 
                               size = labels.size)
        } else {
            g <- g + geom_text(aes(label = labels), size = labels.size)      
        }
    } else {
        if(!is.null(df.u$groups)) {
            g <- g + geom_point(aes(fill = groups), size = pointSize, alpha = alpha, shape = 21, color = "white")
        } else {
            g <- g + geom_point(size = pointSize, alpha = alpha, shape = 21, color = "white", fill = "red")      
        }
    }
    
    # Overlay a concentration ellipse if there are groups
    if(!is.null(df.u$groups) && ellipse) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
        circle <- cbind(cos(theta), sin(theta))
        
        ell <- ddply(df.u, 'groups', function(x) {
            if(nrow(x) <= 2) {
                return(NULL)
            }
            sigma <- var(cbind(x$xvar, x$yvar))
            mu <- c(mean(x$xvar), mean(x$yvar))
            ed <- sqrt(qchisq(ellipse.prob, df = 2))
            data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                       groups = x$groups[1])
        })
        names(ell)[1:2] <- c('xvar', 'yvar')
        g <- g + geom_path(data = ell, aes(color = groups, group = groups))
    }
    
    # Label the variable axes
    if(var.axes) {
        g <- g + 
            geom_text(data = df.v, 
                      aes(label = varname, x = xvar, y = yvar, 
                          angle = angle, hjust = hjust), 
                      color = 'darkred', size = varname.size)
    }
    # Change the name of the legend for groups
    # if(!is.null(groups)) {
    #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
    #                               palette = 'Dark2')
    # }
    
    # TODO: Add a second set of axes
    
    return(g)
}

