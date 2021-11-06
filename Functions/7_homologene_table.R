## Function for making comparison chart of homologs

#In progress still

## Data formatting and homologene function wrapper
# Library: homologene
#Required inputs:
# > Vector of gene names, and their numeric taxonomic identity for 2 different gene sets
# > annots1 = annotation table for gene1 descriptions. If genes are provided as ensembl IDs, the external gene name should also be in this table
#Optional inputs:
# > If gene set 2 is provided as ensembl IDs, annot2 should be specified with the corresponding external gene name mappings. Genes are assumed to already be given as symbols, unless annot2 is provided
# If genes1 is provided as symbols, specify symbols1 = TRUE
#If no target tax provided, default behaviour is to use tax2 as target tax for homologene. This would be the use, for example, if genes1 are mouse genes, and genes2 are human genes. If a different target tax for homologene is desired, for example if genes1 are mouse genes and genes2 are zebrafish genes, can specify human to be the target tax 

homol_table <- function(genes1, tax1, genes2, tax2, annots1, annots2, symbols1 = FALSE, target.tax) {
  
  #Translate Ensembl genes to symbols
  if (symbols1) {
    gene.symbols1 <- genes1
  } else {
    gene.symbols1 <- annots1[which(annots1$ensembl_gene_id %in% genes1),"external_gene_name"]
  }
  
  if (missing(annots2)) {
    gene.symbols2 <- genes2
  } else {
    gene.symbols2 <- annots2[which(annots2$ensembl_gene_id %in% genes2),"external_gene_name"]
  }
  
  #Check for duplicated symbols
  dup.gene1 <- unique(gene.symbols1[duplicated(gene.symbols1)])
  if (length(dup.gene1) != 0) {
    cat("The following symbols are duplicated in gene set 1: ", dup.gene1)
    gene.symbols1 <- unique(gene.symbols1)
  }
  
  dup.gene2 <- unique(gene.symbols2[duplicated(gene.symbols2)])
  if (length(dup.gene2) != 0) {
    cat("The following symbols are duplicated in gene set 2: ", dup.gene2)
    gene.symbols2 <- unique(gene.symbols2)
  }
  
  #Use homologene to get homologs
  if(missing(target.tax)) {
    target.tax <- tax2
  } else {
    homologs2 <- homologene(gene.symbols2, inTax = tax2, outTax = target.tax,  db = homologene::homologeneData2)
  }
  
  homologs1 <- homologene(gene.symbols1, inTax = tax1, outTax = target.tax,  db = homologene::homologeneData2)
  
  #Check for overlap (break if none)
  if (target.tax == tax2) {
    overlap.n <- length(which(homologs1[,2] %in% gene.symbols2))
  } else {
    overlap.n <- length(which(homologs1[,2] %in% homologs2[,2]))
  }
  
  if (overlap.n == 0) {
    cat("No homologous genes in provided sets")
    break
  }
  
  #Build table
  if (target.tax == tax2) {
    overlap.table <- homologs1[which(homologs1[,2] %in% gene.symbols2),c(1,2)]
    colnames(overlap.table) <- c("Gene Set 1", "Gene Set 2")
  } else {
    overlap.table <- homologs1[which(homologs1[,2] %in% homologs2[,2]),c(1,2)]
  }


  #Pasted stuff for reference in building function
  #ortho_overlap <- merge(test_homol_shdr_human, test_homol_andr_human, by = "9606", all.x = FALSE, all.y = FALSE)
  #ortho_overlap <- ortho_overlap[,c(1,2,5)]
  #colnames(ortho_overlap) <- c("Human", "Mouse", "Zebrafish")
  #homol_overlap_drim$Description <- gene_symbols[match(homol_overlap_drim$Mouse, gene_symbols$external_gene_name), 3]
  
}




#