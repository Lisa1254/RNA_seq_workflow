#Function to plot heatmap of genes & enriched go terms from topGO results

#required libraries: "dendextend", "ComplexHeatmap", "circlize"

draw_topGO_heatmap <- function(topGO.table, results.table, annot.table) {
  ##SET UP DATA
  
  #Define genes of interest
  gene.names <- rownames(results.table)
  
  #Initialize data frame
  annGSEA <- data.frame(row.names = gene.names)
  #Populate dataframe for which genes are mapped with which GO terms
  for (j in 1:(length(gene.names))) {
    # create a matching pattern to ensure genes match exactly
    #  '^GENE,'  --> Match at beginning of matching string
    #  ', GENE$'  --> Match at end of matching string
    #  'GENE,'  --> Match between first and last gene in matching string
    gene <- gene.names[j]
    pattern <- paste('^', gene, ',|,', gene, '$|', gene, ',|', gene,  sep = '')
    for (k in 1:(nrow(topGO.table))) {
      if (any(grepl(pattern, topGO.table$genes[k]))) {
        annGSEA[j,k] <- 1
      } else {
        annGSEA[j,k] <- 0
      }
    }
  }
  
  #Assign descriptive colnames
  colnames(annGSEA) <- topGO.table[,"Term"]
  #Ensure no empty rows/cols
  annGSEA <- annGSEA[,apply(annGSEA, 2, mean)!=0]
  annGSEA <- annGSEA[apply(annGSEA, 1, mean)!=0,]
  
  #Make sure results table is ordered same as mapping dataframe
  topTableAligned <- as.data.frame(results.table[which(rownames(results.table) %in% rownames(annGSEA)),])
  topTableAligned <- topTableAligned[match(rownames(annGSEA), rownames(topTableAligned)),]
  if (all(rownames(topTableAligned) != rownames(annGSEA))) {
    cat("Error aligning results table with enriched genes.")
  }
  
  #Assign symbols to rownames
  symbol_matched <- annot.table$external_gene_name[match(rownames(annGSEA), annot.table$ensembl_gene_id)]
  rownames(topTableAligned) <- symbol_matched
  rownames(annGSEA) <- symbol_matched
  
  #ANNOTATIONS
  
  # colour bar for -log10(adjusted p-value) for sigGenes
  dfMinusLog10FDRGenes <- data.frame(-log10(topTableAligned[, 'padj']))
  #My set does not have Inf values, but this is the code if necessary:
  dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
  
  # colour bar for fold changes for sigGenes
  dfFoldChangeGenes <- data.frame(topTableAligned[, 'log2FoldChange'])
  
  #Merge
  dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
  colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
  
  #If changing logfc to categorical, uncomment this
  #dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated', ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated, 'Unchanged'))
  #colours <- list('Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow'))
  
  # colour bar for -log10(Fisher-p values) for topGO results
  dfMinusLog10Fisher<- data.frame(-log10(
    topGO.table[which(topGO.table$Term %in% colnames(annGSEA)), 'weightFisher']))
  colnames(dfMinusLog10Fisher) <- "Weighted fisher (-log10)"
  
  ## HEATMAP PARTS
  
  haGenes <- HeatmapAnnotation(
    df = dfGeneAnno,
    #col = colours,
    #width = unit(1,'cm'),
    annotation_name_side = 'left' #or "top"?
  )
  
  haTerms <- rowAnnotation(
    df = dfMinusLog10Fisher,
    Term = anno_text(
      colnames(annGSEA),
      rot = 0,
      just = 'left',
      gp = gpar(fontsize = 8)),
    annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
    annotation_name_side = 'bottom')
  
  annGSEA.mat <- t(as.matrix(annGSEA))
  
  hmapGSEA <- Heatmap(annGSEA.mat,
                      name = 'Enriched Terms',
                      
                      col = c('0' = 'white', '1' = 'forestgreen'),
                      rect_gp = gpar(col = 'grey85'),
                      border_gp = gpar(col = "black"),
                      
                      cluster_rows = TRUE,
                      show_row_dend = TRUE,
                      row_title = 'Top Terms',
                      row_title_side = 'right',
                      row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                      row_title_rot = 90,
                      show_row_names = TRUE,
                      row_names_rot = 0,
                      row_names_gp = gpar(fontsize = 7, just = 'left'),
                      row_names_side = 'right',
                      
                      cluster_columns = TRUE,
                      show_column_dend = FALSE,
                      column_title = 'Enriched terms',
                      column_title_side = 'top',
                      column_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                      column_title_rot = 0,
                      show_column_names = TRUE,
                      
                      show_heatmap_legend = FALSE,
                      
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      #row_split = term.split,
                      
                      bottom_annotation = haGenes)
  
  draw(hmapGSEA + haTerms,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right') 
  
  
}
