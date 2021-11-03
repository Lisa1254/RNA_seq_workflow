## Functions to adjust count matrices based on surrogate variable analysis

#Note these functions are intended for data exploration, and not to modify counts for use in analyses.

# SVs_plot ----
#Plotting function for checking clustering of SVs
#This is useful to check if estimated surrogate variables align with an expected batch variable.
#Default is to plot SV1 and SV2, color and shape should be specified as the column name of sample information in the DESeq data set object
SVs_plot <- function(dds_obj, SVx = "SV1", SVy = "SV2", color, shape) {
  svdata <- as.data.frame(colData(dds_obj))
  SVx.n <- which(colnames(svdata) == SVx)
  SVy.n <- which(colnames(svdata) == SVy)
  color.n <- which(colnames(svdata) == color)
  shape.n <- which(colnames(svdata) == shape)
  ggplot(data=svdata, aes(x = svdata[,SVx.n], y = svdata[,SVy.n], color = svdata[,color.n], shape = svdata[,shape.n])) +
    geom_point(size = 3) +
    xlab(SVx) +
    ylab(SVy) +
    labs(color = color, shape = shape) +
    coord_fixed() +
    ggtitle(paste0("Surrogate Variable Plot: ", SVx, ", ", SVy))
}

# PCA plot of SV-cleaned data ----

#Functions for adjusting count data and plotting are coded separately for clarity of the process, but adjusted counts are recommended for data exploration only. Adjusted values can be negative and non-integers. They should not be considered count data.


#To adjust count data with surrogate variables through matrix multiplication.
#Inputs are matrix of normalized count data and full model matrix that were used as inputs to svaseq function, and sva fit object returned from svaseq function
rem_SVs <- function(norm.counts, model.mat, svaobj) {
  X = cbind(model.mat, svaobj$sv)
  Hat = solve(t(X)%*%X)%*%t(X)
  beta = (Hat%*%t(norm.counts))
  P = ncol(model.mat)
  cleany = norm.counts-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

#Since this PCA does not include a method for normalization, such as the variance stabilizing transformation, it is recommended to use pre-normalized data. It is being presented here for use with the SVA-adjusted data, but can be used with other inputs
PCs_plot <- function(cts, sample.data, color.var, shape.var, gene_subset, PCX=1, PCY=2, main = "PCA from Count Data", returnLoadings = FALSE) {
  if (!missing(gene_subset)) {
    cts <- cts[gene_subset,]
  } 
  pca <- prcomp(t(cts))
  pca_mat <- as.data.frame(pca$x)
  pca_mat <- cbind(pca_mat, sample.data)
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
