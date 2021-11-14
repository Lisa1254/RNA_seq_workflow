#Function to aid construction of results table with comparison of interest identified
#Input requires prefix used for naming convention of significant gene dataframes produced from DRIMSeq results object - with rownames of gene_id
#Input also requires table of gene symbol annotations to translate from input Ensembl IDs
#All input significance tables must be listed BEFORE sig.prefix and symbol.annots in order to call the tables and their names for use in the function.
#Option to include adjusted p-value by stating inc.padj=TRUE or inc.padj=FALSE

annot_table_full_dtu <- function(..., sig.prefix, symbol.annots, inc.padj = TRUE) {
  comp.table <- data.frame()
  dfs <- list(...)[1:(length(as.list(match.call()))-4)]
  names.in <- as.list(match.call())[-1]
  for (sig in seq(1:length(dfs))) {
    temp_sig <- dfs[[sig]]
    temp_sig <- data.frame(ensembl_gene_id = rownames(temp_sig), comparison = temp_sig[,"adj_pvalue"])
    prefix <- names.in[sig]
    prefix <- gsub(sig.prefix, "", prefix)
    if (inc.padj) {
      temp_sig$comparison <- paste0(prefix, " (", round(temp_sig$comparison, 3), ")")
    } else {
      temp_sig$comparison <- prefix
    }
    comp.table <- rbind(comp.table, temp_sig)
  }
  dup <- duplicated(comp.table$ensembl_gene_id)
  annot_all_pwise <- merge(symbol.annots, comp.table[!dup,], by = "ensembl_gene_id", all.x = FALSE)
  for (gene in seq(1:length(dup))) {
    if (dup[gene]) {
      temp_gene <- comp.table[gene,1]
      table_ind <- which(annot_all_pwise$ensembl_gene_id == temp_gene)
      prev_comp <- annot_all_pwise[table_ind,"comparison"]
      annot_all_pwise[table_ind,"comparison"] <- paste(prev_comp, comp.table[gene,2], sep = ",")
    }
  }
  return(annot_all_pwise)
}



##


