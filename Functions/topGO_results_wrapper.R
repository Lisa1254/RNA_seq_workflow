#Functions for facilitating topGO exploration

# Function to subset significance table for only genes that have provided GO annotations. 
#If provided table is from DESeq2 and has column for log2FoldChange, can also specify if subset should only include genes of increased or decreased expression.
#Rownames of provided significance table should be gene names 

keep_annot <- function(sig.table, go_annot.table, up.or.down) {
  if (missing(up.or.down)) {
    sig_gene <- rownames(sig.table)
  } else if (up.or.down == "up") {
    sig_gene <- rownames(sig.table[which(sig.table$log2FoldChange>0),])
  } else if (up.or.down == "down") {
    sig_gene <- rownames(sig.table[which(sig.table$log2FoldChange<0),])
  }
  keep <- sig_gene %in% go_annot.table[,2]
  keep <- which(keep==TRUE)
  sig_gene <- sig_gene[keep]
  return(sig_gene)
}


# Wrapper function for topGO

#Inputs include vector of significant genes, vector of background genes, ontology of interest, gene annotation table, as constructed in script 2 via biomaRt
#gene annotation table provided is expected to have first column with GO Ids, second column with gene names, and a column labeled "name_1006" which is the biomaRt label for GO Term name
#Term name required in annotation input because topGO cuts off GO term after a specified number of characters, so this is used to complete term names returned
#ontology should be specified as one or more of BP, CC, MF
#Default statistic is fisher, with weight01 algorithm, and pvalue threshold of 0.001, but these can be changed.


GO_2_results <- function(sig_genes, bg.genes, ont, go_annots, thresh_p=0.001, algorithm="weight01", statistic="fisher") {
  gene2go <- unstack(go_annots[,c(1,2)])
  geneList <- factor(as.integer(bg.genes %in% sig_genes))
  names(geneList) <- bg.genes
  keep_res <- data.frame()
  for (ONT in ont) {
    cat("Setting up topGO data for ONT=", ONT, "\n")
    GO_data <- new('topGOdata', ontology=ONT, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2go)
    
    cat("Preparing topGO results for ONT=", ONT, "\n")
    wt_fres <- runTest(GO_data, algorithm=algorithm, statistic=statistic)
    
    cat("Preparing results table for p < ", thresh_p, "for ONT=", ONT, "\n")
    allGO <- usedGO(GO_data)
    all_res <- GenTable(GO_data, weightFisher=wt_fres, orderBy='weightFisher', topNodes=length(allGO))
    all_res$weightFisher <- as.numeric(all_res$weightFisher)
    keep_res_temp <- cbind(ONT=ONT, all_res)[which(all_res$weightFisher < thresh_p),]
    
    cat("Adding gene annotation column for ONT=", ONT, "\n")
    gene_retrieve <- genesInTerm(GO_data, keep_res_temp$GO.ID)
    gene_col <- data.frame()
    for (GO in keep_res_temp$GO.ID) {
      gene_go <- vector()
      for (gene in gene_retrieve[[GO]]) {
        if (gene %in% sig_genes) {
          gene_go <- c(gene_go, gene)
        }
      }
      gene_go <- paste(gene_go, collapse=",")
      gene_col <- rbind(gene_col, gene_go)
      colnames(gene_col) <- "genes"
    }
    keep_res_temp <- cbind(keep_res_temp, gene_col)
    keep_res <- rbind(keep_res, keep_res_temp)
  }
  for (row in 1:length(keep_res$Term)) {
    if (grepl("\\.\\.\\.", keep_res$Term[row])) {
      full.term <- go_annots[which(go_annots$go_id == keep_res[row,2]),"name_1006"]
      keep_res$Term[row] <- full.term[1]
    }
  }
  return(keep_res)
}
