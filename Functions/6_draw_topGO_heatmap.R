#Function to plot heatmap of genes & enriched go terms from topGO results

#required libraries: "dendextend", "ComplexHeatmap", "circlize"

#Input is results table from topGO wrapper function, results table from either DESeq2 or DRIMSeq (specify), and annotations table with gene id as provided through Ensembl and corresponding external symbol, as mapped with biomaRt
#Can specify minimum number of genes per term (gene.min) and terms per gene (term.min) that must be represented
#Can specify title for plot with main=

draw_topGO_heatmap <- function(topGO.table, results.table, annot.table, res.style = c("DESeq2", "DRIMSeq"), gene.min=1, term.min=1, main = "default", gene.font.size=7) {
  
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
  #Ensure no empty rows/cols, or if different min set, that each meets threshold
  go.gene.table <- go.gene.table[,apply(go.gene.table, 2, sum)>=term.min]
  go.gene.table <- go.gene.table[apply(go.gene.table, 1, sum)>=gene.min,]
  
  #Make sure results table is ordered same as mapping dataframe
  topTableAligned <- as.data.frame(results.table[which(rownames(results.table) %in% rownames(go.gene.table)),])
  topTableAligned <- topTableAligned[match(rownames(go.gene.table), rownames(topTableAligned)),]
  if (all(rownames(topTableAligned) != rownames(go.gene.table))) {
    cat("Error aligning results table with enriched genes.")
  }
  
  #Assign symbols to rownames
  symbol_matched <- annot.table$external_gene_name[match(rownames(go.gene.table), annot.table$ensembl_gene_id)]
  
  #Check for & resolve duplicates
  dups <- unique(symbol_matched[duplicated(symbol_matched)])
  if (length(dups) != 0) {
    for (dup in dups) {
      n <- 1
      for (sym in 1:length(symbol_matched)) {
        if (symbol_matched[sym] == dup) {
          symbol_matched[sym] <- paste(dup, n, sep = ".")
          n <- n+1
        }
      }
    }
  }
  
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
    
    #Adjust color ramp depending on whether both directions of change are included, or if only up/down 
    if (sum(dfGeneAnno$Log2FC < 0) == 0) {
      col_fun_l2fc <- colorRamp2(c(0, max(dfGeneAnno[,2])), c("white", "darkslategray3"))
    } else if (sum(dfGeneAnno$Log2FC > 0) == 0) {
      col_fun_l2fc <- colorRamp2(c(min(dfGeneAnno[,2]), 0), c("coral1", "white"))
    } else {
      col_fun_l2fc <- colorRamp2(c(min(dfGeneAnno[,2]), 0, max(dfGeneAnno[,2])), c("coral1", "white", "darkslategray3"))
    }
    
    
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
  
  col_fun_gene = colorRamp2(c(min(dfGeneAnno[,1]), max(dfGeneAnno[,1])), c("white", "firebrick"))
  
  # colour bar for -log10(Fisher-p values) for topGO results
  dfMinusLog10Fisher<- data.frame(-log10(
    topGO.table[which(topGO.table$Term %in% colnames(go.gene.table)), 'weightFisher']))
  colnames(dfMinusLog10Fisher) <- "Weighted fisher (-log10)"
  col_fun_term = colorRamp2(c(min(dfMinusLog10Fisher[,1]), max(dfMinusLog10Fisher[,1])), c("white", "darkblue"))
  
  ## HEATMAP PARTS
  
  #Title
  if (main == "default") {
    plot_title <- paste0("Enriched term annotation for ", res.style, " results")
  } else {
    plot_title <- main
  }
  
  #Gene annotations
  if (res.style == "DESeq2") {
    if ((sum(dfGeneAnno$Log2FC < 0) != 0) & (sum(dfGeneAnno$Log2FC > 0) != 0)) {
      haGenes <- HeatmapAnnotation(
        df = dfGeneAnno,
        annotation_name_side = 'left',
        annotation_name_gp = gpar(fontsize = 8, fontface = 'bold'),
        col = list(`Gene score` = col_fun_gene, Log2FC = col_fun_l2fc),
        annotation_legend_param = list(
          `Gene score` = list(title = "Gene Score",
                              col_fun = col_fun_gene),
          Log2FC = list(title = "Log2FC",
                        break_dist = 1,
                        at = c(round(min(dfGeneAnno[,2]),1), 0, round(max(dfGeneAnno[,2]),1)),
                        col_fun = col_fun_l2fc)
        )
      )
      
    } else {
      haGenes <- HeatmapAnnotation(
        df = dfGeneAnno,
        annotation_name_side = 'left',
        annotation_name_gp = gpar(fontsize = 8, fontface = 'bold'),
        col = list(`Gene score` = col_fun_gene, Log2FC = col_fun_l2fc)
      )
      
    }
    
  } else {
    haGenes <- HeatmapAnnotation(
      df = dfGeneAnno,
      annotation_name_side = 'left',
      annotation_name_gp = gpar(fontsize = 8, fontface = 'bold'),
      col = list(`Gene score` = col_fun_gene)
    )
  }
  
  #Term annotations
  haTerms <- rowAnnotation(
    df = dfMinusLog10Fisher,
    Term = anno_text(
      colnames(go.gene.table),
      rot = 0,
      just = 'left',
      gp = gpar(fontsize = 8)),
    annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
    annotation_name_side = 'bottom',
    annotation_name_gp = gpar(fontsize = 8, fontface = 'bold'),
    col = list(`Weighted fisher (-log10)` = col_fun_term))
  
  
  #Main heatmap
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
                      column_names_gp = gpar(fontsize = gene.font.size),
                      
                      show_heatmap_legend = FALSE,
                      
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      
                      bottom_annotation = haGenes)
  
  #DRAW
  draw(hmapGSEA + haTerms,
       column_title = plot_title) 
  
  
}
