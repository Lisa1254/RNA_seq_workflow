##
## Comparisons of Multiple Datasets ----
##

#Data required for this script:
##   > significant genes determined in scripts 4 & 5
load("Output/sig_deg_asd_control.Rdata")
load("Output/sig_drim_asd_control.Rdata")
##   > other gene lists for comparison (can use SPARK ASD-risk genes and provided example gene lists as specified below)
##   > Gene to symbol & description mappings (Script 2)
load("Output/gene_symbols.Rdata")

##
# Libraries & Functions----
##

#homologene is used to map genes to known orthologs in other species for comparison
library(homologene)
#ComplexHeatmap and circlize are used for the heatmap table style of visual for showing which genes are represented in multiple comparisons or datasets
library(ComplexHeatmap)
library(circlize)

#Function used to format and pass data of 2 gene sets to homologene, returning table of homologs present in both sets
source("Functions/7_homologene_table.R")
#Function to produce heatmap style table of genes in multiple comparisons or datasets
source("Functions/7_gene_overlap_table.R")


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

#This part considering including
#Loaded genes from separately passing the full Sharon dataset to the pipeline
#Saving here for use when re-loading script
save(sig_full_deg_Striatum, sig_full_drim_PFC, sig_full_drim_Striatum, file = "Data/full_data_results.Rdata")
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

#Comparison of all DGE analysis results with full set of human SPARK ASD-risk genes
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

#Can also use symbols for gene set input insead of Ensembl IDs
names_full_drim <- c(rownames(sig_full_drim_PFC), rownames(sig_full_drim_Striatum))
names_full_drim <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% names_full_drim), "external_gene_name"]

dtu_fulldata_spark_all <- homol_table(genes1 = names_full_drim, 
                    tax1 = mouse_taxID, 
                    annots1 = gene_symbols, 
                    genes2 = human_SPARK_all, 
                    tax2 = human_taxID, 
                    symbols1 = TRUE)

#Saving here the ones I plan to use in the plots below for safety, but consider deleting, and only keeping figures.
save(deg_spark_all, dtu_ex_dr_human, dtu_fulldata_spark_all, dtu_spark_all, file = paste0(out_dir, "homol_tables.Rdata"))

##
# Comparison Plots ----
##

#Use provided function for making heatmap style table of overlapping results

#Genes in example data that had homolog in the SPARK ASD-risk sets:
#Define gene sets of interest using descriptive names
DGE_Subset_Data <- deg_spark_all[,2]

colnames(deg_spark_all)[1:2] <- c("Mouse", "Human")
colnames(dtu_ex_dr_human)[1:3] <- c("Mouse", "Zebrafish", "Human")
colnames(dtu_spark_all)[1:2] <- c("Mouse", "Human")
colnames(dtu_fulldata_spark_all)[1:2] <- c("Mouse", "Human")

#OR, maybe more useful to have standardized name for each gene set as vector
DGE_Subset_Data <- deg_spark_all$Human
DTU_Example_DR <- dtu_ex_dr_human$Human
DTU_Subset_Data <- dtu_spark_all$Human
DTU_Full_Data <- dtu_fulldata_spark_all$Human

DGE_Subset_Striatum <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_asd.control)),"external_gene_name"]
DGE_Full_Striatum <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_full_deg_Striatum)),"external_gene_name"]
SPARK_All <- deg_spark_all$Mouse

draw_overlap_table(DGE_Subset_Striatum, DGE_Full_Striatum, SPARK_All,
                           desc.annots = gene_symbols)
draw_overlap_table(DGE_Subset_Data, DTU_Example_DR, DTU_Subset_Data, DTU_Full_Data, human_SPARK_all)


# Would like to add Venn diagram style comparison








##