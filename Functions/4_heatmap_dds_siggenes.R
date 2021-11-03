# Heatmap for relatedness of samples in most significant genes

#Libraries circlize, complexHeatmap, dplyr

#Treatment annotation - colData column must be labelled "treatment" and can handle up to 5 levels

# Inputs required: 
# > dds_obj with sample data and gene counts
# > dds_results with padj for subsetting
# Optional inputs:
# > sample_sub, a vector of which samples to include
# > default n=200, (maximum) number of genes to include in table, provided they meet the provided padj threshold
# > default p.adj=0.5, maximum padj value for inclusion in chart

draw_ch_TopGenes <- function(dds_obj, dds_results, sig_genes = NULL, sample_sub, n=200, p.adj=0.5) {
  
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
  groups <- factor(colData(dds_obj)[,"treatment"])
  gp_level <- levels(groups)
  cols5 <- c("darkblue", "magenta", "forestgreen", "chocolate1", "gray48")
  gp_cols <- cols5[1:length(gp_level)]
  names(gp_cols) <- gp_level
  
  dates <- factor(colData(dds_obj)[,"date"])
  date_level <- levels(dates)
  date_cols <- rainbow(length(date_level))
  names(date_cols) <- date_level
  
  col_list <- list(group=gp_cols, date=date_cols)
  ha_trt = HeatmapAnnotation(group = groups, date=dates,
                          col = col_list)
  
  #Make Heatmaps
  
  col_fun = colorRamp2(c(-2, 0, 2), c("red", "white", "blue"))
  
  Heatmap(z.mat, name = "z-score",
          show_row_name = FALSE,
          col = col_fun,
          row_names_gp = gpar(fontsize = 6),
          cluster_columns = TRUE,
          top_annotation = ha_trt,
          right_annotation = row_ha
  )
  
}

