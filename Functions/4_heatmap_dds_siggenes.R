# Heatmap for relatedness of samples in most significant genes

#Libraries circlize, complexHeatmap, dplyr

#Resulting plot is heatmap of gene expression levels based on either input genes list, or top n genes ordered by increasing padj of results table
#"ranks" annotation shows relative total count of gene across all samples, such that purple annotations have higher overall expression, and white/gray annotations have lower overall expression (default is to display, but can turn off with rank.annot = FALSE)

# Inputs required: 
# > dds_obj with sample data and gene counts
# > dds_results with padj for subsetting
# > colname of variable being used for colour annotations in sample clustering
# Optional inputs:
# > second colname of sample data for second sample cluster annotation
# > sample_sub, a vector of which samples to include
# > default n=200, (maximum) number of genes to include in table, provided they meet the provided padj threshold
# > default p.adj=0.5, maximum padj value for inclusion in chart
# > Alternate to using top 200 genes that meet padj threshold, a specific list of genes to include can be provided
#If results table is provided, it will override given list of significant genes

draw_ch_TopGenes <- function(dds_obj, dds_results, col.name, col.name2, sig_genes = NULL, sample_sub, n=200, p.adj=0.5, rank.annot = TRUE) {
  
  #Prep relatedness tables
  if (missing(dds_results)) {
    sigGenes <- sig_genes
  } else {
    sigGenes <- dds_results[which(dds_results$padj < p.adj),]
    sigGenes <- rownames(head(sigGenes[order(sigGenes$padj),], n))
    cat(length(sigGenes), " genes meet padj threshold\n")
  }
  if (missing(sample_sub)) {
    cts <- assay(dds_obj)
  } else {
    dds_obj <- dds_obj[,sample_sub]
    cts <- assay(dds_obj)
  }
  sigGenes_keep <- sigGenes[which(sigGenes %in% rownames(cts))]
  if (length(sigGenes_keep) != length(sigGenes)) {
    rm_genes <- length(sigGenes) - length(sigGenes_keep)
    cat(rm_genes, " genes of ", length(sigGenes), " in sig_genes not in gene counts provided.\n")
  }
  plotDat <- vst(cts)[sigGenes_keep,]
  z.mat <- t(scale(t(plotDat), center=TRUE, scale=TRUE))
  
  #Prep gene count rank
  counts.gene.sort <- dds_obj %>%
    assay() %>%
    rowSums() %>%
    sort()
  names.gene.sort <- names(counts.gene.sort)
  ind <- which(names.gene.sort %in% sigGenes)
  row_ha = rowAnnotation(ranks = ind, col = list(ranks=colorRamp2(c(1, median(ind), max(ind)), c("white", "grey", "purple"))))
  
  #Prep treatment group annotations
  group.col <- which(colnames(colData(dds_obj)) == col.name)
  groups <- factor(colData(dds_obj)[,group.col])
  gp_level <- levels(groups)
  colsA <- c("darkblue", "magenta", "forestgreen", "chocolate1", "gray48", "cadetblue2", "bisque", "darkmagenta", "gray15")
  gp_cols <- colsA[1:length(gp_level)]
  names(gp_cols) <- gp_level
  
  if (!missing(col.name2)) {
    group.col2 <- which(colnames(colData(dds_obj)) == col.name2)
    groups2 <- factor(colData(dds_obj)[,group.col2])
    gp_level2 <- levels(groups2)
    colsB <- c("cornflowerblue", "coral1", "blueviolet", "darkgoldenrod1", "brown4", "springgreen", "slategray1", "seagreen", "navy", "gray26")
    gp_cols2 <- colsB[1:length(gp_level2)]
    names(gp_cols2) <- gp_level2
    
    col_list <- list(group=gp_cols, group2=gp_cols2)
    ha_trt = HeatmapAnnotation(group = groups, group2=groups2,
                               col = col_list,
                               annotation_label = c(col.name, col.name2))
  } else {
    col_list <- list(group=gp_cols)
    ha_trt = HeatmapAnnotation(group = groups,
                               col = col_list,
                               annotation_label = col.name)
  }

  
  #Make Heatmaps
  
  col_fun = colorRamp2(c(-2, 0, 2), c("red", "white", "blue"))
  
  h <- Heatmap(z.mat, name = "z-score",
          show_row_name = FALSE,
          show_row_dend = FALSE,
          col = col_fun,
          row_names_gp = gpar(fontsize = 6),
          cluster_columns = TRUE,
          top_annotation = ha_trt,
  )
  
  if (rank.annot) {
    draw(h + row_ha)
  } else {
    draw(h)
  }
  
}

