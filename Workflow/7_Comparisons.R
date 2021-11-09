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

#Function used to 
source("Functions/7_homologene_table.R")

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

##
# Comparison Plots ----
##

# Would like to add Venn diagram style comparison

# This is a paste from my other script of creating a heatmap style of overlap in significant genes in several sets
# This is more useful with more than 2 comparison sets
# Consider modifying the code to keep as function or process so that it can accommodate desired number of gene sets for comparison and/or adding another set of relevant genes to the import for comparison

any_asdSet_sym <- annot_sig_dge$external_gene_name
any_vanSet_sym <- annot_sig_dge_vanco$external_gene_name

top50_risk_drSym <- annot_ens_dr$external_gene_name

gene_sets <- c("any_asdSet_sym", "any_vanSet_sym", "top50_risk_drSym", "resc_genes_update")

all_genes_any_set <- unique(c(any_asdSet_sym, any_vanSet_sym, top50_risk_drSym, resc_genes_update))

SigGene_table_sym <- matrix(nrow=length(all_genes_any_set), ncol=length(gene_sets))
rownames(SigGene_table_sym) <- all_genes_any_set
colnames(SigGene_table_sym) <- c("Any ASD-NT", "VANCO", "top 50 Risk", "Prev. Rescue Study")

for (gl in 1:length(gene_sets)) {
  temp_list <- get(gene_sets[gl])
  for (gene in 1:nrow(SigGene_table_sym)) {
    if (rownames(SigGene_table_sym)[gene] %in% temp_list) {
      SigGene_table_sym[gene,gl] <- 1
    } else {
      SigGene_table_sym[gene,gl] <- 0
    }
  }
}

SigGene_table_sym2 <- SigGene_table_sym[which(rowSums(SigGene_table_sym)>1),]

gene_annot <- gene_symbol_all[match(rownames(SigGene_table_sym2),gene_symbol_all$external_gene_name),"description"]
gene_annot <- strsplit(gene_annot, " \\[")
gene_annot <- unlist(gene_annot)
gene_annot <- gene_annot[seq(1,length(gene_annot),2)]

ha_gene <- rowAnnotation(
  Gene = anno_text(paste(rownames(SigGene_table_sym2), gene_annot, sep = ": "),
                   rot = 0,
                   just = 'left',
                   gp = gpar(fontsize = 8))
)

hmap <- Heatmap(SigGene_table_sym2,
                name = 'Genes in Common',
                col = c('0' = 'white', '1' = 'forestgreen'),
                rect_gp = gpar(col = 'grey85'),
                border_gp = gpar(col = "black"),
                
                show_row_dend = FALSE)
draw(hmap + ha_gene)





##