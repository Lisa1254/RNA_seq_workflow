## Functions to adjust count matrices based on surrogate variable analysis

#Note these functions are intended for data exploration, and not to modify counts for use in analyses.

#To adjust count data with surrogate variables through matrix multiplication. This function can be followed up with prcomp to get a PCA of data with SVs accounted for. 
#This is recommended for data exploration only, and can return values that are negative values and non-integers. It should not be considered count data.
#Inputs are matrix of normalized count data and full model matrix that were used as inputs to svaseq function, and sva fit object returned from svaseq function
rem_SVs <- function(norm.counts, model.mat, svaobj) {
  X = cbind(model.mat, svaobj$sv)
  Hat = solve(t(X)%*%X)%*%t(X)
  beta = (Hat%*%t(norm.counts))
  P = ncol(model.mat)
  cleany = norm.counts-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}

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
    coord_fixed() +
    ggtitle(paste0("Surrogate Variable Plot: ", SVx, ", ", SVy))
}