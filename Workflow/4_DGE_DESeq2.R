##
## Differential Gene Expression Analysis with DESeq2 ----
##

#Data required for this script:
##   > DESeq data set object (output of script 3_Data_Exploration)
### If using output after collapsing technical replicates:
load("Output/gene_dds.Rdata")
### If using output with surrogate variable information from batch correction:
load("Output/gene_dds_sva.Rdata")
##   > annotation table with gene symbols and descriptions to identify results (output of script 2_Data_Import)
load("Output/gene_symbols.Rdata")

##
# Libraries & Functions----
##

#DESeq2 is used for differential gene expression analysis
library(DESeq2)
#ggplot2 is used for count and volcano plots
library(ggplot2)
#circlize, complexHeatmap, dplyr used in heatmap of expression of significant genes
library(ComplexHeatmap)
library(circlize)
library(dplyr)

#Source import functions
#Used as example in generalized code for streamlining multiple pairwise comparisons, but not specifically applicable to this dataset
source("Functions/4_multiple_comparisons_deseq.R")
#Includes function for drawing gene count plot
source("Functions/4_dge_count_plot.R")
#Function to make heatmap of gene expression of significant genes
source("Functions/4_heatmap_dds_siggenes.R")

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
gene_dds <- DESeq(gene_dds)

#If running the differential expression analysis returns the warning that some rows did not converge in beta (trouble converging the general linear model), the rows that did not converge can be filtered out with the following:
gene_dds <- gene_dds[which(mcols(gene_dds)$betaConv),]

#This may happen when the complex model encounters genes with very low counts. If a lot of rows do not converge, consider re-evaluating if the correct model is being used, or adding additional filtering for low count genes prior to running the differential expression analysis 
#Good explanation at https://support.bioconductor.org/p/65091/


#Check resultsNames to help determine contrast coefficients for results
resultsNames(gene_dds)

#This dataset had a nested design of multiple donors to each treatment type (and tissue type in full dataset), expressed as a combined factor, so we'll need to define contrasts of interests manually to each comparison of interest

#Locate coefficients for groups of interest
#Using experiment specific names of "control" and "asd" for groups of interest. 
coef.control <- ifelse(grepl("CON", resultsNames(gene_dds)), 1, 0)
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

#To save a tsv of significant genes & descriptions
#Included for demonstration, but not saved with Outputs
annot_asd.control <- gene_symbols[which(gene_symbols$ensembl_gene_id %in% rownames(sig_asd.control)),]

#To add column with log2fold change
annot_asd.control$direction <- sig_asd.control[match(annot_asd.control$ensembl_gene_id, rownames(sig_asd.control)),"log2FoldChange"]
#If categorical "DOWN"/"UP" is prefered to numeric value
annot_asd.control$direction <- ifelse(annot_asd.control$direction < 0, "DOWN", "UP")

write.table(annot_asd.control, file=paste0(out_dir, "sig_gene_table_deg.tsv"), quote=FALSE, sep='\t', col.names = NA)

## 
# Generalized code for contrasts & results ----
##

#If dataset has multiple comparisons of interest, can use provided functions to facilitate and streamline process
#For example, if there are two treatment types, and a control group, each pairwise comparison will be of interest

#NOTE: the group variable being supplied to the pairwise comparison functions MUST be the additive variable in the design used for differential expression analysis.
## If multiple comparisons are desired from a complex design matrix, apply the contrast coefficient method above multiple times instead

#Pairwise comparisons are not especially useful for the grouped factors in this dataset, but will be used here as an example
#Specify colname of factor of interest using variable argument
#If only a subset of comparisons are of interest can define manually in same format of Var1.Var2

all_pwise <- pwise_gps(sample_frame = colData(gene_dds), variable = "group")

#Produce results tables for each comparison of interest
#Specify colname of factor of interest using variable argument
#p-adjustment method is specified above, or remove argument to use default method "BH"
#This applies the res_contrast function provided to each contrast indicated manually or by the contrast function (format Var1.Var2).
#Output is results object for 

for (con in all_pwise) {
  assign(paste0("res_", con),
         results_contrast(des = gene_dds, 
                      variable = "group", 
                      contrast = con, 
                      pAdjM = pAdjM))
}

#Subset each results table for genes that meet defined criteria for significance
for (res in ls()[grepl("res_", ls())]) {
  temp_res <- get(res)
  assign(gsub("res", "sig", res),
         as.data.frame(temp_res[which(
           (temp_res$pvalue < pval_th) &
             (temp_res$padj < padj_th) & 
             (abs(temp_res$log2FoldChange) > l2fc_th)),]))
  rm(temp_res)
}

#Save any results and significance tables of interest

#To produce annotations table of significant genes in all comparisons, can use included function, which returns subset of gene symbol annotation table, with additional column specifying which comparison and direction of lfc
#Code makes use of consistent naming conventions attributed above, where tables are named sig_Var1.Var2 and direction of log2FoldChange is specific to Var1
#Input as many significance tables as desired, followed by specifying the prefix used for significance tables ("sig_"), and the gene symbols complete annotations table
#Options for lfc_style are c("numeric", "categorical", "none")
#In this demonstration with provided dataset, the comparisons between individual donor groups is not biologically relevant, and produces many genes considered statistically significant. So I will input 3 of the resulting tables for the example.

test_multi_annotations <- annot_table_full(sig_ASD_9.ASD_3, sig_CON_5.ASD_3, sig_CON_5.ASD_9, sig.prefix = "sig_", symbol.annots = gene_symbols, lfc_style = "categorical")

#To save as tsv:
write.table(test_multi_annotations, file=paste0(out_dir, "sig_gene_compare_deg.tsv"), quote=FALSE, sep='\t', col.names = NA)
#Will save subset as R file for reference
test_multi_annotations <- test_multi_annotations[1:10,]
save(test_multi_annotations, file = paste0(out_dir, "test_multi_annot.Rdata"))



##
# Plots for Significant genes ----
##

#Using provided function, can visualize the differences in counts of significant genes
#Demonstrating with subset of first 10 genes in significance table constructed
#Note: sample data plots well with default binwidth=1/30; can add binwidth argument if points are too big/small.

count.plot(gene_dds, 
           rownames(sig_asd.control)[1:10], 
           gene_symbols, 
           "treatment")

#Saved as 4_gene_counts


## Volcano plot

#Input is complete results table from DESeq2, as dataframe, excluding any rows that have NA as value in padj column
#Significance is coloured based on padj and log2FoldChange thresholds defined above, or modify for different visualization

sub_res_df <- as.data.frame(res_asd.control)[!is.na(res_asd.control$padj),]
ggplot(sub_res_df, mapping=aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = 
        ifelse((padj>padj_th) | (abs(log2FoldChange)<l2fc_th), 'black',
               ifelse(log2FoldChange>l2fc_th, "blue", "red")))) +
  scale_colour_manual(labels = c("No Change", "Sig Up", "Sig Down"),
                      values=c('black', 'blue', 'red')) + 
  labs(color = "DGE", title = "Volcano Plot of DGE")

#Saved as 4_volcano

#To add annotations to the genes with the largest log2fc

#Define threshold here, or if different thresholds for increased and decreased exprssion are desired, input directly into plot code below
l2fc_lab_th <- 5

#Add symbols column to results dataframe
sub_res_df$symbol <- gene_symbols[match(rownames(sub_res_df), gene_symbols$ensembl_gene_id), "external_gene_name"]

#First part of volcano plot is identical to above, but with an additional geom_text line for each up and down expressed genes to allow independent adjustment of position (through hjust and vjust)
ggplot(sub_res_df, mapping=aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = 
                   ifelse((padj>padj_th) | (abs(log2FoldChange)<l2fc_th), 'black',
                          ifelse(log2FoldChange>l2fc_th, "blue", "red")))) +
  scale_colour_manual(labels = c("No Change", "Sig Up", "Sig Down"),
                      values=c('black', 'blue', 'red')) + 
  labs(color = "DGE", title = "Volcano Plot of DGE") +
  geom_text(aes(label=ifelse((log2FoldChange>l2fc_lab_th1) & (padj<padj_th), symbol,'')),hjust=0.5,vjust=-0.5) +
  geom_text(aes(label=ifelse((log2FoldChange<(-l2fc_lab_th1)) & (padj<padj_th), symbol,'')),hjust=0.25,vjust=-0.5)


# Heatmap of how gene expression clusters samples
#Using provided function. Several modifiers are describe in function script, a few examples are shown here

#With only ASD/CON group annotation & top 200 genes by increasing padj in results table:
draw_ch_TopGenes(dds_obj = gene_dds, 
                 dds_results = res_asd.control, 
                 col.name = "treatment", 
                 rank.annot = FALSE)

#Same, but with annotations added for sample donor group, and rank of gene expression (ordered by total number of counts to gene across samples):
draw_ch_TopGenes(dds_obj = gene_dds, 
                 dds_results = res_asd.control, 
                 col.name = "treatment", 
                 col.name2 = "group")

#This one is saved as 4_gene_heatmap

#With CON_5 samples removed:
rem_con5 <- rownames(colData(gene_dds)[-(which(colData(gene_dds)$group == "CON_5")),])
draw_ch_TopGenes(dds_obj = gene_dds, 
                 dds_results = res_asd.control, 
                 col.name = "treatment", 
                 col.name2 = "group", 
                 sample_sub = rem_con5)

#With all genes in the sig_asd.control table produced:
draw_ch_TopGenes(dds_obj = gene_dds, 
                 sig_genes = rownames(sig_asd.control), 
                 col.name = "treatment", 
                 col.name2 = "group")

##