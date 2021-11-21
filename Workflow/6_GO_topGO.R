##
## Gene Ontology Analysis with topGO ----
##

#This script demonstrates the use of topGO for gene ontology enrichment analysis, along with results visualization


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
source("Functions/6_topGO_results_wrapper.R")
#Function to plot topGO results as heatmap with genes & GO terms
source("Functions/6_draw_topGO_heatmap.R")

#Data required for this script:
##   > Significant genes for exploration from DGE (Script 4) & DTU (Script 5) analyses. Just a list of gene names is sufficient for analysis, but full table of significance statistics is used for annotations in plots.
load("Output/sig_deg_asd_control.Rdata")
load("Output/sig_drim_asd_control.Rdata")
##   > Gene to GO mappings (Script 2)
load("Output/go_annotations.Rdata")
##   > Gene symbol mappings used in plots (Script 2)
load("Output/gene_symbols.Rdata")


##
# Define Key Variables ----
##

#Output directory for saving files
out_dir <- paste0(getwd(),"/Output/")


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
TG_drim_asd.control <- GO_2_results(sig.names_drim_asd.control, bg_genes, c("BP", "CC", "MF"), go_annotations)


#Save results tables
save(TG_dge_down_asd.control, TG_dge_up_asd.control, TG_drim_asd.control, file = paste0(out_dir, "TG_res_tables.Rdata"))


##
# Plot results ----
## 

#Provided function acts to organize results information into complex heatmap

#Input is results table from topGO wrapper function, results table from either DESeq2 or DRIMSeq (specify), and annotations table with gene id as provided through Ensembl and corresponding external symbol, as mapped with biomaRt

#For the genes with decreased expression (subset the results table for computational efficiency, but it isn't necessary for functionality)
draw_topGO_heatmap(topGO.table = TG_dge_down_asd.control, 
                   results.table = sig_asd.control[sig.names_down_asd.control,], 
                   annot.table = gene_symbols, 
                   res.style = "DESeq2")

#Same for genes with increased expression (not many results in this example dataset)
draw_topGO_heatmap(topGO.table = TG_dge_up_asd.control, 
                   results.table = sig_asd.control[sig.names_up_asd.control,], 
                   annot.table = gene_symbols, 
                   res.style = "DESeq2")

#Can include both increased and decreased expression in same plot
#This plot includes a lot of genes, so I've added argument to decrease gene font size in the plot, but for presentation consider decreasing the number of terms/genes in the plot
draw_topGO_heatmap(topGO.table = rbind(TG_dge_down_asd.control,TG_dge_up_asd.control), 
                   results.table = sig_asd.control, 
                   annot.table = gene_symbols, 
                   res.style = "DESeq2",
                   gene.font.size = 5.5)

#Saved as 6_tg_deseq

#Also works for DRIMSeq style results table
#This example had many genes with only "Cytoplasm" as the enriched term, so I've added arguments for gene.min=2 to only include in the plots genes that have at least 2 terms, and term.min=3 just for further  simplification
draw_topGO_heatmap(topGO.table = TG_drim_asd.control, 
                   results.table = sig_drim_asd.control, 
                   annot.table = gene_symbols, 
                   res.style = "DRIMSeq",
                   gene.min = 2, 
                   term.min = 3)

#Saved as 6_tg_drim
