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
#Default in DESeq2 is BH; for more information see help page for p.adjust
#?p.adjust
p.adj <- "BH"

## Consider adding here intended thresholds for p.value, p.adj, l2fc. 


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
#If using gene_dds that does not have design yet specified, define here
design(gene_dds) <- ~ -1 + group

#Apply DESeq2 method for differential gene expression
#Consider if adding modifiers for test as "Wald" vs. "LRT", or for fitType
gene_dds <- DESeq(gene_dds)

#Consider keeping this line, and describing when a gene might not converge in beta. This won't be necessary in most workflows.
#Keep only genes which converge in beta (11 rows to remove)
gene_dds <- gene_dds[which(mcols(gene_dds)$betaConv),]