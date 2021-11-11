## Exploratory Sample heatmap

# Method must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski", as per dist() function. Default is "euclidean"
#Expected colData column within dds object is "group", which describes important experimental grouping for identifying relevant sample clustering

#Uses libraries ComplexHeatmap and circlize

sample_dist_heatmap <- function(dds, method = "euclidean", gp_name) {
  vsd <- vst(dds, blind = FALSE)
  sampleDists <- dist(t(assay(vsd)), method = method)
  sampleDistMatrix <- as.matrix( sampleDists )
  
  if (!missing(gp_name)) {
    gp_ind <- which(colnames(colData(dds)) == gp_name)
    gp_vector <- colData(dds)[,gp_ind]
    rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix), gp_vector, sep = "_")
  }
  
  colnames(sampleDistMatrix) <- NULL
  col_fun = circlize::colorRamp2(c(0, max(sampleDistMatrix)), c("blue", "white"))
  Heatmap(sampleDistMatrix,
          name = paste(method, "distance"),
          col = col_fun,
          column_title = paste0("Sample distance (", method, ")"))
}


##
