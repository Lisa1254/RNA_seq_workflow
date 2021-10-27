##
## Gene Ontology Analysis with topGO ----
##

#Data required for this script:
##   > Significant genes for exploration from DGE (Script 4) & DTU (Script 5) analyses. Just a list of gene names is sufficient for analysis, but full table of significance statistics is used for annotations in plots.
load("Output/sig_deg_asd_control.Rdata")
load("Output/sig_drim_asd_control.Rdata")
##   > Gene to GO mappings (Script 2)
load("Output/go_annotations.Rdata")


##
# Libraries & Functions----
##

#Functional enrichment analysis of gene ontology done with topGO
library(topGO)

#Functions for ensuring all genes used in topGO have GO annotations, and for wrapping the topGO run tests into results table
source("Functions/topGO_results_wrapper.R")

##
# Define Key Variables ----
##

#Output directory for saving files
out_dir = paste0(getwd(),"/Output/")


##
# topGO analysis ----
##

#List of background genes is required to compare for enrichment
#Since GO annotation dataframe was constructed from the genes expressed in the original dataset, they can be used for the background genes
bg_genes <- unique(go_annotations$ensembl_gene_id)

# build the gene 2 GO annotation list in the format required for topGO
gene_2_GO <- unstack(go_annotations[,c(1,2)])



##