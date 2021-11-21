##
## Comparisons of Multiple Datasets ----
##

#This script provides an example for how to approach comparisons between different contrasts in a single dataset, and results between different datasets

##
# Libraries & Functions----
##

#homologene is used to map genes to known orthologs in other species for comparison
library(homologene)
#ComplexHeatmap and circlize are used for the heatmap table style of visual for showing which genes are represented in multiple comparisons or datasets
library(ComplexHeatmap)
#ggVennDiagram can quickly make Venn Diagrams using ggplot2
library(ggplot2)
library(ggVennDiagram)

#Function used to format and pass data of 2 gene sets to homologene, returning table of homologs present in both sets
source("Functions/7_homologene_table.R")
#Function to produce heatmap style table of genes in multiple comparisons or datasets
source("Functions/7_gene_overlap_table.R")

#Data required for this script:
##   > significant genes determined in scripts 4 & 5
load("Output/sig_deg_asd_control.Rdata")
load("Output/sig_drim_asd_control.Rdata")
##   > other gene lists for comparison (can use SPARK ASD-risk genes and provided example gene lists as specified below)
##   > Gene to symbol & description mappings (Script 2)
load("Output/gene_symbols.Rdata")



##
# Define Key Variables ----
##

#Output directory for saving files
out_dir = paste0(getwd(),"/Output/")

#Taxonomic ID number for species being compared
human_taxID <- 9606
mouse_taxID <- 10090
zebrafish_taxID <- 7955


#For comparison, will be using list of genes with demonstrated strong association with ASD-risk, according to Simons Foundation Powering Autism Research for Knowledge
#https://www.sfari.org/grant/2021-genomics-of-asd-pathways-to-genetic-therapies-request-for-applications/
#Including both the full list of genes with proposed ASD-risk, and high-priority list
#Also for providing list of zebrafish genes to demonstrate obtaining homolog table as positioning experimental results within field of established research of comparable experiments
load("Data/ex_compare_genes.Rdata")

ex_dr_geneset_ens <- ex_annots_dr$ensembl_gene_id
ex_dr_geneset_symbol <- ex_annots_dr$external_gene_name

#Including for comparison Striatum results for DGE and DTU analysis, and PFC results for DTU analysis from passing full data set to pipeline
load("Data/full_data_results.Rdata")

##
# Compare known homologs ----
##

#The source data was exploring genes differing between GF mice colonized with samples derived from ASD and typically-developing (TD) people, with a focus on alternative splicing. Here I'll compare the results from this experiment with the human genes on the SPARK high ASD-risk for homologs.

#Using provided homologene wrapper function for formatting inputs and creating table of homologous genes in different sets, a few examples are below

#If there are no overlapping genes, the function should return an error message for having no homologous genes in provided sets
#None of the genes significant to the DGE analysis are in the SPARK priority ASD list
test_none <- homol_table(genes1 = rownames(sig_asd.control), 
                    tax1 = mouse_taxID, 
                    annots1 = gene_symbols, 
                    genes2 = human_SPARK_priority, 
                    tax2 = human_taxID)

#Comparison of all DGE analysis results with full set of human SPARK ASD-risk genes (3 genes)
deg_spark_all <- homol_table(genes1 = rownames(sig_asd.control), 
                         tax1 = mouse_taxID, 
                         annots1 = gene_symbols, 
                         genes2 = human_SPARK_all, 
                         tax2 = human_taxID)

#This shows the DTU analysis results that are homologous to the human genes in the full list of SPARK ASD-risk genes (1 gene)
dtu_spark_all <- homol_table(genes1 = rownames(sig_drim_asd.control), 
                    tax1 = mouse_taxID, 
                    annots1 = gene_symbols, 
                    genes2 = human_SPARK_all, 
                    tax2 = human_taxID)

#If comparing to other another experiment that uses a different species but comparable experimental design, it may be interesting to identify if homologous genes were identified in both.
dtu_ex_dr <- homol_table(genes1 = rownames(sig_drim_asd.control), 
                    tax1 = mouse_taxID, 
                    annots1 = gene_symbols, 
                    genes2 = ex_dr_geneset_symbol, 
                    tax2 = zebrafish_taxID)

#Can add third taxonomy as target, if like in the above example would like to also know the human homologous gene to the mouse and zebrafish gene, so that other comparisons can be done downstream, such as to the human SPARK ASD-risk genes
#The second gene set can also be provided as Ensembl gene ids if a second annotation table is provided
dtu_ex_dr_human <- homol_table(genes1 = rownames(sig_drim_asd.control), 
                         tax1 = mouse_taxID, 
                         annots1 = gene_symbols, 
                         genes2 = ex_dr_geneset_ens, 
                         tax2 = zebrafish_taxID,
                         annots2 = ex_annots_dr,
                         target.tax = human_taxID)

#Using DTU results from running complete dataset from Sharon et al. through pipe with all human SPARK ASD-risk genes (6 results)

#Can also use symbols for gene set input instead of Ensembl IDs
#Using only Striatum results for comparison
names_full_drim <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_full_drim_Striatum)), "external_gene_name"]

dtu_fulldata_spark_all <- homol_table(genes1 = names_full_drim, 
                    tax1 = mouse_taxID, 
                    annots1 = gene_symbols, 
                    genes2 = human_SPARK_all, 
                    tax2 = human_taxID, 
                    symbols1 = TRUE)

#Output tables saved for reference
save(deg_spark_all, dtu_ex_dr_human, dtu_fulldata_spark_all, dtu_spark_all, file = paste0(out_dir, "homol_tables.Rdata"))

##
# Comparison Plots (Heatmap) ----
##

#Use provided function for making heatmap style table of overlapping results

#Genes in example data that had homolog in the SPARK ASD-risk sets:
#Define gene sets of interest using descriptive names; using gene symbols of Human homologs for each set to compare overlap in sets
#First three sets here are from the homologene comparison with the SPARK genes
DGE_Subset_Data <- deg_spark_all[,2]
DTU_Subset_Data <- dtu_spark_all[,2]
DTU_Full_Data <- dtu_fulldata_spark_all[,2]
#Including additional genes from comparison between this pipe's data and the provided zebrafish genes
DTU_Example_DR <- dtu_ex_dr_human[,3]
#Will also rename the SPARK genes for use as column label in plot
All_SPARK_Genes <- human_SPARK_all


draw_overlap_table(DGE_Subset_Data, DTU_Example_DR, DTU_Subset_Data, DTU_Full_Data, All_SPARK_Genes)
#The DTU_Example_DR set is not on the plot because none of the genes were also in the SPARK Genes set, which was the basis for defining the other gene lists
#Saved as 7_common_genes_no_annot

#Can also add descriptive annotations to gene names, if an annotation chart is provided that uses the same symbol type as the provided gene sets
#Will demonstrate with showing the overlap between resulting differentially expressed genes from the provided data subset, and full dataset of Striatum samples
#Adding the genes that were identified as overlapping between the subsample data and the SPARK genes

DGE_Subset_Striatum <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_asd.control)),"external_gene_name"]
DGE_Full_Striatum <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_full_deg_Striatum)),"external_gene_name"]
All_SPARK_Genes <- deg_spark_all[,1]

draw_overlap_table(DGE_Subset_Striatum, DGE_Full_Striatum, All_SPARK_Genes,
                           desc.annots = gene_symbols)

#Saved as 7_common_genes_with_annot

#Can use Ensembl or other gene names instead of symbols, but not with annotations
DTU_PFC_Full_Dataset <- rownames(sig_full_drim_PFC)
DTU_Striatum_Full_Dataset <- rownames(sig_full_drim_Striatum)
DTU_Striatum_Sample_Data <- rownames(sig_drim_asd.control)

draw_overlap_table(DTU_PFC_Full_Dataset, DTU_Striatum_Full_Dataset, DTU_Striatum_Sample_Data)

#If there are no overlapping genes in the provided sets, a descriptive error will return
#If you were expecting genes to be in common, double check that gene names are of the same type (e.g. here human gene symbols and mouse Ensembl IDs are being mixed)
draw_overlap_table(DTU_Example_DR, DTU_Striatum_Sample_Data)


##
# Comparison Plots (Venn Diagram) ----
##

##ggVennDiagram uses gradient colors to help visualize size of sets
#Takes list of sets as input
#Can customize fill colour and category names easily

#3-Dimensional plot showing relationship between DTU in sample dataset (Striatum), and full dataset for each of Striatum and PFC
DTU_venn_data <- list(
  A = DTU_PFC_Full_Dataset, 
  B = DTU_Striatum_Full_Dataset, 
  C = DTU_Striatum_Sample_Data
)
DTU_venn_names <- c("PFC (Full)","Striatum (Full)","Striatum (Sample)")

#scale_x_continuous argument is added to increase the space allocated along the x axis for the descriptive labels
ggVennDiagram(DTU_venn_data, category.names = DTU_venn_names) +
  ggplot2::scale_fill_gradient(low="white",high = "darkred") +
  ggplot2::scale_x_continuous(expand = expansion(mult = .2))

#Saved as 7_3D_Venn

#Can also make 4-Dimensional Venn Diagrams
#Will demonstrate with showing the overlap of increased and decreased expression of ASD to CON samples in the subset Striatum data vs. the full set of Striatum data

#Define gene sets
#Can paste the rowname subsetting directly into the creation of the venn data list instead of creating first, but keeping separate here for clarity
DGE_Down_Sample <- rownames(sig_asd.control[which(sig_asd.control$log2FoldChange<0),])
DGE_Up_Sample <- rownames(sig_asd.control[which(sig_asd.control$log2FoldChange>0),])
DGE_Down_Full <- rownames(sig_full_deg_Striatum[which(sig_full_deg_Striatum$log2FoldChange<0),])
DGE_Up_Full <- rownames(sig_full_deg_Striatum[which(sig_full_deg_Striatum$log2FoldChange>0),])

#Define list of gene sets
DGE_venn_data <- list(
  A = DGE_Down_Sample,
  B = DGE_Up_Sample,
  C = DGE_Down_Full,
  D = DGE_Up_Full)

#Define category names
DGE_venn_names <- c("Down (Sample)","Up (Sample)","Down (Full)", "Up (Full)")

#Plot Venn Diagram
#label_alpha=0 removes the lighter box around the count labels in the plot
ggVennDiagram(DGE_venn_data, label_alpha = 0,
              category.names = DGE_venn_names) +
  ggplot2::scale_fill_gradient(low="white",high = "deepskyblue3") +
  ggplot2::scale_x_continuous(expand = expansion(mult = .2))

#This shows that the included sample data appears to have far more significantly differentially expressed genes than the full dataset, and only 3 of the 13 genes in the full set are also in the subset data, the samples selected are not entirely representative of the group.

#Saved as 7_4D_Venn

## Bonus Venn Diagram, Experimental Usage ----

#Since the overlap sections of the Venn Diagram are calculated by having exact character matches for gene names in each set, can plot the overlap of different experiment results (mouse gene symbol) with the SPARK ASD-risk genes (human gene symbol) by replacing the relevant mouse symbols with human symbols in each list to generate a character match

#This code is left rough, as this is just an experiment to demonstrate different usage of the Venn Diagram package. There may be genes misclassified if they overlap with SPARK in one set (and get attributed human gene symbol) but not in a different set, so the match would be missed.

#Prep Sample data DEG
deg_sub_spark <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_asd.control)),"external_gene_name"]
for (gene in 1:length(deg_sub_spark)) {
  if (deg_sub_spark[gene] %in% deg_spark_all[,1]) {
    new_gene <- deg_spark_all[which(deg_spark_all[,1] == deg_sub_spark[gene]),2]
    deg_sub_spark[gene] <- new_gene
  }
}


#Prep Sample data DTU
dtu_sub_spark <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_drim_asd.control)),"external_gene_name"]
for (gene in 1:length(dtu_sub_spark)) {
  if (dtu_sub_spark[gene] %in% dtu_spark_all[,1]) {
    new_gene <- dtu_spark_all[which(dtu_spark_all[,1] == dtu_sub_spark[gene]),2]
    dtu_sub_spark[gene] <- new_gene
  }
}

#Prep Full dataset DTU 
dtu_full_sub_spark <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_full_drim_Striatum)),"external_gene_name"]
for (gene in 1:length(dtu_full_sub_spark)) {
  if (dtu_full_sub_spark[gene] %in% dtu_fulldata_spark_all[,1]) {
    new_gene <- dtu_fulldata_spark_all[which(dtu_fulldata_spark_all[,1] == dtu_full_sub_spark[gene]),2]
    dtu_full_sub_spark[gene] <- new_gene
  }
}


#Define gene list & names for plot
test_venn_list <- list(
  A = deg_sub_spark,
  B = dtu_sub_spark,
  C = dtu_full_sub_spark,
  D = human_SPARK_all
)

test_venn_names <- c("DEG (Sample)", "DTU (Sample)", "DTU (Full)", "SPARK")

#Draw Venn Diagram
ggVennDiagram(test_venn_list, label_alpha = 0,
              category.names = test_venn_names) +
  ggplot2::scale_fill_gradient(low="white",high = "deepskyblue3") +
  ggplot2::scale_x_continuous(expand = expansion(mult = .2))

#Saved as 7_venn_spark_ex

##