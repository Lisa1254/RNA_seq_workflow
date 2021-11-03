#Function to plot heatmap of genes & enriched go terms from topGO results

#required libraries: "dendextend", "ComplexHeatmap", "circlize"
#Input is results table from topGO wrapper function, results table from either DESeq2 or DRIMSeq (specify), and annotations table with gene id as provided through Ensembl and corresponding external symbol, as mapped with biomaRt

draw_topGO_heatmap <- function(topGO.table, results.table, annot.table, res.style = c("DESeq2", "DRIMSeq")) {
  ##SET UP DATA
  
  #Define genes of interest
  gene.names <- rownames(results.table)
  
  #Initialize data frame
  go.gene.table <- data.frame(row.names = gene.names)
  #Populate dataframe for which genes are mapped with which GO terms
  for (j in 1:(length(gene.names))) {
    # pattern specifies to match format of genes in topGO results chart produced with provided wrapper function
    gene <- gene.names[j]
    pattern <- paste('^', gene, ',|,', gene, '$|', gene, ',|', gene,  sep = '')
    for (k in 1:(nrow(topGO.table))) {
      if (any(grepl(pattern, topGO.table$genes[k]))) {
        go.gene.table[j,k] <- 1
      } else {
        go.gene.table[j,k] <- 0
      }
    }
  }
  
  #Assign descriptive colnames
  colnames(go.gene.table) <- topGO.table[,"Term"]
  #Ensure no empty rows/cols
  go.gene.table <- go.gene.table[,apply(go.gene.table, 2, mean)!=0]
  go.gene.table <- go.gene.table[apply(go.gene.table, 1, mean)!=0,]
  
  #Make sure results table is ordered same as mapping dataframe
  topTableAligned <- as.data.frame(results.table[which(rownames(results.table) %in% rownames(go.gene.table)),])
  topTableAligned <- topTableAligned[match(rownames(go.gene.table), rownames(topTableAligned)),]
  if (all(rownames(topTableAligned) != rownames(go.gene.table))) {
    cat("Error aligning results table with enriched genes.")
  }
  
  #Assign symbols to rownames
  symbol_matched <- annot.table$external_gene_name[match(rownames(go.gene.table), annot.table$ensembl_gene_id)]
  rownames(topTableAligned) <- symbol_matched
  rownames(go.gene.table) <- symbol_matched
  
  #ANNOTATIONS
  
  res.style <- match.arg(res.style)
  
  if (res.style == "DESeq2") {
    # colour bar for -log10(adjusted p-value) for sigGenes
    dfMinusLog10FDRGenes <- data.frame(-log10(topTableAligned[, 'padj']))
    #Fix infinity values if present:
    dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
    
    # colour bar for fold changes for sigGenes
    dfFoldChangeGenes <- data.frame(topTableAligned[, 'log2FoldChange'])
    
    #Merge
    dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
    colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
  }
  
  if (res.style == "DRIMSeq") {
    # colour bar for -log10(adjusted p-value) for sigGenes
    dfMinusLog10FDRGenes <- data.frame(-log10(topTableAligned[, 'adj_pvalue']))
    #Fix infinity values if present:
    dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
    
    colnames(dfMinusLog10FDRGenes) <- c('Gene score')
    
    #Name annotation same as with DESeq2 annotations
    dfGeneAnno <- dfMinusLog10FDRGenes
    colnames(dfGeneAnno) <- c('Gene score')
  }
  
  # colour bar for -log10(Fisher-p values) for topGO results
  dfMinusLog10Fisher<- data.frame(-log10(
    topGO.table[which(topGO.table$Term %in% colnames(go.gene.table)), 'weightFisher']))
  colnames(dfMinusLog10Fisher) <- "Weighted fisher (-log10)"
  
  ## HEATMAP PARTS
  
  haGenes <- HeatmapAnnotation(
    df = dfGeneAnno,
    annotation_name_side = 'left'
  )
  
  haTerms <- rowAnnotation(
    df = dfMinusLog10Fisher,
    Term = anno_text(
      colnames(go.gene.table),
      rot = 0,
      just = 'left',
      gp = gpar(fontsize = 8)),
    annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
    annotation_name_side = 'bottom')
  
  go.gene.table.mat <- t(as.matrix(go.gene.table))
  
  hmapGSEA <- Heatmap(go.gene.table.mat,
                      name = 'Enriched Terms',
                      
                      col = c('0' = 'white', '1' = 'forestgreen'),
                      rect_gp = gpar(col = 'grey85'),
                      border_gp = gpar(col = "black"),
                      
                      cluster_rows = TRUE,
                      show_row_dend = TRUE,
                      row_title = 'GO Terms',
                      row_title_side = 'left',
                      row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                      row_title_rot = 90,
                      show_row_names = TRUE,
                      row_names_rot = 0,
                      row_names_gp = gpar(fontsize = 7, just = 'left'),
                      row_names_side = 'right',
                      
                      cluster_columns = TRUE,
                      show_column_dend = FALSE,
                      column_title = 'Genes',
                      column_title_side = 'bottom',
                      column_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                      column_title_rot = 0,
                      show_column_names = TRUE,
                      column_names_gp = gpar(fontsize = 7),
                      
                      show_heatmap_legend = FALSE,
                      
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      
                      bottom_annotation = haGenes)
  
  draw(hmapGSEA + haTerms,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right') 
  
  
}
