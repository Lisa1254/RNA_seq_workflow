##
## Differential Gene Expression Analysis with DESeq2 ----
##

#Data required for this script:
##   > DESeq data set object (output of script 3_Data_Exploration)
### If using output after collapsing technical replicates:
load("Output/gene_dds.Rdata")
### If using output with surrogate variable information from batch correction:
load("Output/gene_dds_sva.Rdata")

##
# Libraries & Functions----
##

#DESeq2 is used for differential gene expression analysis
library(DESeq2)

#Source import functions

##
# Define Key Variables ----
##

#Output directory for saving files
out_dir = paste0(getwd(),"/Output/")

#Define p-value adjustment method for multiple hypothesis testing 
#Options include: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
# see ?p.adjust for more details
pAdjM <- "BH"

#Define significance thresholds for maximum p-value, maximum adjusted p-value, and minimum log2FoldChange magnitude
pval_th <- 0.05
padj_th <- 0.1
l2fc_th <- 1
#To have no limits, set pval_th and padj_th to 1, and l2fc_th to 0

##
# DESeq2 Analysis ----
##

#Define design matrix

#Notes about choosing design:
## Using -1 removes the reference level intercept, making comparisons between multiple groups easier
## Simple format can be expressed as a variation of design = ~ -1 + group + batch
## If there is more than one relevant independent factor to the design, an additive model may be used (e.g. design = ~ -1 + treatment + tissue + batch), if the factors are not independent, interaction factors may be needed in the design formula. See DESeq2 documentation for more information.
## In this workflow, I'll be using the combined factor method in place of the additive model, which uses a single variable ("group") to express all combinations of the multiple design elements (e.g. group is the combined factor of tissue, treatment, and donor). 

#If using gene_dds from SVA batch correction, design was defined already
design(gene_dds)
#If using gene_dds that does not have design yet specified, use or modify line below
design(gene_dds) <- ~ -1 + group

#Apply DESeq2 method for differential gene expression
#Consider if adding modifiers for test as "Wald" vs. "LRT", or for fitType
gene_dds <- DESeq(gene_dds)

#Consider removing this line, as it won't be necessary in most workflows.
#If keeping, describe when a gene might not converge in beta
#Good explanation at https://support.bioconductor.org/p/65091/
#Keep only genes which converge in beta
gene_dds <- gene_dds[which(mcols(gene_dds)$betaConv),]

#Check resultsNames to help determine contrast coefficients for results
resultsNames(gene_dds)

#This dataset had a nested design of multiple donors to each treatment type (and tissue type in full dataset), expressed as a combined factor, so we'll need to define contrasts of interests manually to each comparison of interest

#Locate coefficients for groups of interest
#Using experiment specific names of "control" and "asd" for groups of interest. 
coef.control <- ifelse(grepl("control", resultsNames(gene_dds)), 1, 0)
coef.asd <- ifelse(grepl("ASD", resultsNames(gene_dds)), 1, 0)

#Define contrast coefficients with appropriate weighting. 
#NOTE: The negative coefficients indicate the "reference" level, and the positive coefficients indicate the treatment that corresponds with the directional sign of the log fold change.
coef_asd.control <- coef.asd/sum(coef.asd) - coef.control/sum(coef.control)


#Determine results stats for each comparison
res_asd.control <- results(gene_dds, contrast = coef_asd.control, pAdjustMethod = pAdjM)

#Subset results table for genes that meet defined criteria for significance
sig_asd.control <- as.data.frame(res_asd.control[which(
          (res_asd.control$pvalue < pval_th) &
          (res_asd.control$padj < padj_th) & 
          (abs(res_asd.control$log2FoldChange) > l2fc_th)),])


#Save results table for significant genes, or all genes if desired
save(sig_asd.control, file = paste0(out_dir, "sig_deg_asd_control.Rdata"))
save(res_asd.control, file = paste0(out_dir, "res_deg_asd_control.Rdata"))


## 
# Generalized code for contrasts & results ----
##

#In progress; ideas & pasted code from previous work

#Function to define all pairwise contrasts in group of sample data frame
#Note: order of pairing will follow order of group in the data frame. Be mindful of this when interpreting direction of change in results table, or define contrasts manually above.
#This was from a script where "group" had several factor levels, each of which had a relevant comparison
pwise_gps <- function(sample_frame) {
  groups <- unique(sample_frame[,which(colnames(sample_frame) == "group")])
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

#Function to use defined contrast (format GF.CV as above) to produce results table from DESeq2 object
#Default p-adjust method is Benjamini and Hochberg ("BH")
#p-adj methods include c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
# see ?p.adjust for more details
res_contrast <- function(des, contrast, pAdjM = "BH") {
  cons <- strsplit(contrast, split = "\\.")
  res <- results(des, contrast = c("group", cons[[1]][1], cons[[1]][2]), pAdjustMethod = pAdjM)
  return(res)
}

#Produce results tables for each comparison of interest
#This applies the res_contrast function defined above to each contrast indicated manually or by the contrast function.
for (con in contrasts_0) {
  assign(paste0("res_g0_", con),
         res_contrast(des_0, con))
}

