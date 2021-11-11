## Functions to facilitate multiple pairwise comparisons 

#This script contains 2 functions to streamline the process of making multiple compaisons in a dataset with several variables of interest
#For example, if there are two treatment types, and a control group, each pairwise comparison will be of interest
#The first function sets up a vector of contrasts in the format 
#"Var1.Var2"
#The second function wraps DESeq2's results function with the provided p-adjusment method according to the contrast specified in each
#Using assign() function, all pairwise comparisons identified in the first function can be supplied to the DESeq2 results wrapper with a simple command
#See script 4_DGE_DESeq2 for example usage


#Function to define all pairwise contrasts in defined variable of sample information data frame
#Note: order of pairing will follow order of group in the data frame. Be mindful of this when interpreting direction of change in results table, or define contrasts manually above.

pwise_gps <- function(sample_frame, variable) {
  groups <- unique(sample_frame[,which(colnames(sample_frame) == variable)])
  k <- length(groups)
  n <- k*(k-1)/2
  cat("Constructing ", n, " pairwise comparisons\n")
  pairs <- vector()
  for (g in 1:(k-1)) {
    for (g2 in (g+1):k) {
      pairs <- c(pairs, paste0(groups[g], ".", groups[g2]))
    }
  }
  return(pairs)
}

#Function to use defined contrast (format Var1.Var2 as above) to produce results table from DESeq2 object
#Default p-adjust method is Benjamini and Hochberg ("BH")
#p-adj methods include c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
# see ?p.adjust for more details
results_contrast <- function(des, variable, contrast, pAdjM = "BH") {
  cons <- strsplit(contrast, split = "\\.")
  res <- results(des, contrast = c(variable, cons[[1]][1], cons[[1]][2]), pAdjustMethod = pAdjM)
  return(res)
}

#Function to aid construction of results table with comparison of interest identified
#Input requires prefix used for naming convention of significant gene dataframes produced from DESeq2 results object
#Input also requires table of gene symbol annotations to translate from input Ensembl IDs
#All input significance tables must be listed BEFORE sig.prefix and symbol.annots in order to call the tables and their names for use in the function.

annot_table_full <- function(..., sig.prefix, symbol.annots, lfc_style = c("numeric", "categorical", "none")) {
  lfc.table <- data.frame()
  dfs <- list(...)[1:(length(as.list(match.call()))-4)]
  names.in <- as.list(match.call())[-1]
  for (sig in seq(1:length(dfs))) {
    temp_sig <- dfs[[sig]]
    temp_sig <- data.frame(ensembl_gene_id = rownames(temp_sig), comparison = temp_sig[,"log2FoldChange"])
    prefix <- names.in[sig]
    prefix <- gsub(sig.prefix, "", prefix)
    if (lfc_style == "numeric") {
      temp_sig$comparison <- paste0(prefix, " (", round(temp_sig$comparison, 3), ")")
    } else if (lfc_style == "categorical") {
      temp_sig$comparison <- ifelse(temp_sig$comparison < 0, paste0(prefix, " (DOWN)"), paste0(prefix, " (UP)"))
    } else if (lfc_style == "none") {
      temp_sig$comparison <- prefix
    }
    lfc.table <- rbind(lfc.table, temp_sig)
  }
  dup <- duplicated(lfc.table$ensembl_gene_id)
  annot_all_pwise <- merge(symbol.annots, lfc.table[!dup,], by = "ensembl_gene_id", all.x = FALSE)
  for (gene in seq(1:length(dup))) {
    if (dup[gene]) {
      temp_gene <- lfc.table[gene,1]
      table_ind <- which(annot_all_pwise$ensembl_gene_id == temp_gene)
      prev_comp <- annot_all_pwise[table_ind,"comparison"]
      annot_all_pwise[table_ind,"comparison"] <- paste(prev_comp, lfc.table[gene,2], sep = ",")
    }
  }
  return(annot_all_pwise)
}



##
