#Function to make heatmap style visual of genes overlapping between different datasets or different comparisons within the same dataset

#Required libraries: ComplexHeatmap

#Inputs: vectors of gene sets to plot. All genes must be expressed in same naming convention for function to locate matches (i.e. all symbols or Ensembl IDs belonging to same species)
#Optional to add gene descriptions by providing annotation table. If adding descriptions to the gene names, input gene sets should be provided in symbol form
#Columns are titled by names of gene sets provided, substituting "_" for a space

draw_overlap_table <- function(..., desc.annots) {
  #Define gene lists for comparison
  if (missing(desc.annots)) {
    gene.sets <- list(...)[1:(length(as.character(match.call()))-1)]
    names.in <- as.character(match.call())[-1]
  } else {
    gene.sets <- list(...)[1:(length(as.character(match.call()))-2)]
    names.in <- as.character(match.call())[-1]
    names.in <- names.in[-length(names.in)]
  }
  
  #Prep table of gene overlap
  all_genes_any_set <- unique(unlist(gene.sets))
  overlap.table <- matrix(nrow=length(all_genes_any_set), ncol=length(names.in))
  rownames(overlap.table) <- all_genes_any_set
  colnames(overlap.table) <- gsub("_", " ", names.in)
  
  #Define 1 as presence of gene in a set, and 0 as missing from set
  for (gs in 1:length(names.in)) {
    for (gene in 1:nrow(overlap.table)) {
      if (rownames(overlap.table)[gene] %in% gene.sets[[gs]]) {
        overlap.table[gene,gs] <- 1
      } else {
        overlap.table[gene,gs] <- 0
      }
    }
  }
  
  
  
  #Keep only genes present in more than 1 set
  overlap.table <- overlap.table[which(rowSums(overlap.table)>1),]
  #Remove sets if no genes have overlap with other sets
  overlap.table <- overlap.table[,which(colSums(overlap.table)>0)]

  #End script if there are no overlapping genes in provided set
  if ((nrow(overlap.table) == 0) | ncol(overlap.table) == 0) {
    stop("No genes overlapping in sets provided")
    }
  
  #Prep heatmap
  hmap <- Heatmap(overlap.table,
                  col = c('0' = 'white', '1' = 'forestgreen'),
                  rect_gp = gpar(col = 'grey85'),
                  border_gp = gpar(col = "black"),
                  
                  show_row_dend = FALSE,
                  show_column_dend = FALSE, 
                  
                  heatmap_legend_param = list(
                    at = c(1,0),
                    labels = c("In Set", "Not in Set"),
                    title = "Presence of Gene",
                    legend_gp = gpar(fill = c("forestgreen", "white")),
                    border = 'grey85'
                  ),
                  
                  column_title = 'Genes in Common')
  
  
  #Draw, and add gene description annotation if including
  if(missing(desc.annots)) {
    
    #Plot map without gene annotations
    draw(hmap)
    
  } else {
    #Get annotations from annotation table; citation removed for visual simplicity
    gene_annot <- desc.annots[match(rownames(overlap.table), desc.annots$external_gene_name),"description"]
    gene_annot <- strsplit(gene_annot, " \\[")
    gene_annot <- unlist(gene_annot)
    gene_annot <- gene_annot[seq(1,length(gene_annot),2)]
    
    #Make annotations for gene descriptions
    ha_gene <- rowAnnotation(
      Gene = anno_text(paste(rownames(overlap.table), 
                             gene_annot, sep = ": "),
                       rot = 0,
                       just = 'left',
                       gp = gpar(fontsize = 8))
    )
    
    #Draw Heatmap
    draw(hmap + ha_gene)
  }
  
  
}