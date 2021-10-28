##
## Gene Ontology Analysis with topGO ----
##

#Data required for this script:
##   > Significant genes for exploration from DGE (Script 4) & DTU (Script 5) analyses. Just a list of gene names is sufficient for analysis, but full table of significance statistics is used for annotations in plots.
load("Output/sig_deg_asd_control.Rdata")
load("Output/sig_drim_asd_control.Rdata")
##   > Gene to GO mappings (Script 2)
load("Output/go_annotations.Rdata")
##   > Gene symbol mappings used in plots (Script 2)
load("Output/gene_symbols.Rdata")


##
# Libraries & Functions----
##

#Functional enrichment analysis of gene ontology done with topGO
library(topGO)
#Plot of enriched GO terms and genes done with complex Heatmap, and associated packages
library(dendextend)
library(ComplexHeatmap)
library(circlize)

#Functions for ensuring all genes used in topGO have GO annotations, and for wrapping the topGO run tests into results table
source("Functions/topGO_results_wrapper.R")
#Function to plot topGO results as heatmap with genes & GO terms
source("Functions/draw_topGO_heatmap.R")

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

#Use provided subset function to return significant genes that have GO annotations, and separate DEG with increased and descreased expression
sig.names_up_asd.control <- keep_annot(sig_asd.control, go_annotations, "up")
sig.names_down_asd.control <- keep_annot(sig_asd.control, go_annotations, "down")
sig.names_drim_asd.control <- keep_annot(sig_drim_asd.control, go_annotations)

#Use provided topGO wrapper function to return table of GO with enrichment using default fisher statistic and weight01 algorithm for each significant gene set. Using all three ontologies.

TG_dge_up_asd.control <- GO_2_results(sig.names_up_asd.control, bg_genes, c("BP", "CC", "MF"), go_annotations)
TG_dge_down_asd.control <- GO_2_results(sig.names_down_asd.control, bg_genes, c("BP", "CC", "MF"), go_annotations)
TG_drim_asd.control <- GO_2_results(sig.names_up_asd.control, bg_genes, c("BP", "CC", "MF"), go_annotations)


#Save results tables
save(TG_dge_down_asd.control, TG_dge_up_asd.control, TG_drim_asd.control, file = paste0(out_dir, "TG_res_tables.Rdata"))


##
# Plot results ----
## 

#In progress, updating plotting function with complex heatmap

draw_topGO_heatmap(TG_dge_down_asd.control, sig_asd.control[sig.names_down_asd.control,], gene_symbols, res.style = "DESeq2")

draw_topGO_heatmap(TG_drim_asd.control, sig_drim_asd.control, gene_symbols, res.style = "DRIMSeq")
