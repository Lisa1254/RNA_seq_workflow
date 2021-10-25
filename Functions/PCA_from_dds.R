## Function to make PCA plot from DESeq Data Set object
#Option to use only a subset of genes for PCA
#Option to return PCA loadings instead of plotting
#Default to plot PC1 & PC2

PCs_plot_dds <- function(dds, color.var, shape.var, gene_subset, PCX=1, PCY=2, main = "PCA from VST data", returnLoadings = FALSE) {
  if (missing(gene_subset)) {
    vst_pc <- vst(dds)
  } else {
    vst_pc <- vst(dds)[gene_subset,]
  }
  pca <- prcomp(t(assay(vst_pc)))
  pca_mat <- as.data.frame(pca$x)
  pca_mat <- cbind(pca_mat, colData(dds))
  pca <- summary(pca)
  percentVar <- round(100 * pca$importance[2,], 2)
  if (returnLoadings == FALSE) {
    color.col <- which(colnames(pca_mat) == color.var)
    shape.col <- which(colnames(pca_mat) == shape.var)
    ggplot(pca_mat, aes(x = pca_mat[,PCX], y = pca_mat[,PCY], color = pca_mat[,color.col], shape = pca_mat[,shape.col])) +
      geom_point(size = 3) +
      #scale_shape_manual(values=c(15, 1, 2, 16)) +
      xlab(paste0("PC", PCX, ": ", percentVar[PCX], "% variance")) +
      ylab(paste0("PC", PCY, ": ", percentVar[PCY], "% variance")) +
      labs(color = color.var, shape = shape.var) +
      coord_fixed() +
      ggtitle(main)
  } else if (returnLoadings == TRUE) {
    loadings <- as.data.frame(pca$rotation)
    return(loadings)
  }
}
