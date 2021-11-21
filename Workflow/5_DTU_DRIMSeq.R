##
## Differential Transcript Usage Analysis with DRIMSeq ----
##

#This script demonstrates the use of DRIMSeq for differential transcript usage analysis, and annotation tables and visualizations of results

##
# Libraries & Functions----
##

#DRIMSeq is used for differential transcript usage analysis.
library(DRIMSeq)
#Transcript annotations retrieved from biomaRt
library(biomaRt)
#ggplot2, reshape2, and used in plots of transcript proportions
library(ggplot2)
library(reshape2)

#Source custom functions

#This functions modifies the source code for DRIMSeq's plotProportions function to allow for user input to specific attributes
source("Functions/5_plot_proportions_DRIMmod.R")
#If creating a reference table of all results with multiple comparisons (not done in this example script), can use following:
source("Functions/5_multiple_comparison_table_dtu.R")

#Data required for this script:
##   > transcript-level counts & sample information in DRIMSeq format (output of script 3_Data_Exploration)
load("Output/tx_ct_formatted.Rdata")


##
# Define Key Variables ----
##

#Output directory for saving files
out_dir = paste0(getwd(),"/Output/")

#Filtering at the gene level ensures that the observed transcript ratios have some minimal reliability. Using recommendations from "Swimming Downstream" paper by Love et al.
#Recommendation for gene expression is all samples, "n_large"
#Recommendation for transcript expression is smallest replicates in a group, "n_small"
#Also define min_{gene/feature}_expr as the minimum sum count of a gene or transcript across all samples. Using suggested value of 10
#Also define min_prop as minimum proportion of a transcript in at least n_small number of samples to be considered. Using suggested value of 0.1
n_large <- 8
n_small <- 2
min_gene_expr <- 10
min_feature_expr <- 10
min_feature_prop <- 0.1

#padj thresh
padj_th <- 0.1

##
# DTU Analysis ----
## 

#Define model matrix
#See discussion of model matrix & design in script 4_DGE_DESeq2
#Using combined factor design for model matrix
mm.full <- model.matrix(~ -1 + group, samples)

#Construct DRIMSeq dataset object
drimds <- dmDSdata(counts=tx_cts, samples=samples)

#Filter for desired DTU stringency: filter values described & defined above, following standard recommendations
drimds <- dmFilter(drimds,
                    min_samps_gene_expr = n_large, 
                    min_samps_feature_expr = n_small,  
                    min_gene_expr = min_gene_expr,
                    min_feature_expr = min_feature_expr,
                    min_samps_feature_prop = n_small, 
                    min_feature_prop = min_feature_prop)



#Perform "Precision" analysis
#This step is quite slow. Depending on number of genes and samples, complexity of design matrix, data within the R workspace, and system capabilities, could take 45 min - 7 hours.

#To test for an estimate of how long the full step will take, start with the first 200 genes:
drimds_sub <- drimds[1:200,]
system.time(drimds_sub <- dmPrecision(drimds_sub, design = mm.full))

#Took 89.866 seconds
#Total time expected to be
(89.866/60)*(length(drimds)/200)
#66 minutes expected

#Save files here, because next step is very slow, so it is a good idea to be able to return to before the command
save(drimds, mm.full, file=paste0(out_dir, "drim_prep_files.Rdata"))
#Clean what is not needed from environment to improve speed of next step
rm(list = setdiff(ls(), c("drimds", "mm.full")))

system.time(drimds <- dmPrecision(drimds, design = mm.full))


#Calculate Fit
#Takes up to a few minutes
drimds <- dmFit(drimds, design = mm.full, verbose = 1)

#Locate coefficients for groups of interest
#Using experiment specific names of "control" and "asd" for groups of interest. 
coef.control <- ifelse(grepl("CON", colnames(design(drimds))), 1, 0)
coef.asd <- ifelse(grepl("ASD", colnames(design(drimds))), 1, 0)

#Define contrast coefficients with appropriate weighting. 
coef_asd.control <- coef.asd/sum(coef.asd) - coef.control/sum(coef.control)

#Tests
test_asd.control <- dmTest(drimds, contrast = coef_asd.control, one_way = FALSE, verbose = TRUE)

res_drim_asd.control <- results(test_asd.control)

sig_drim_asd.control <- res_drim_asd.control[which(res_drim_asd.control$adj_pvalue < padj_th),]
#For convenience of downstream analysis, will use gene_id column as rownames instead
rownames(sig_drim_asd.control) <- sig_drim_asd.control$gene_id
sig_drim_asd.control$gene_id <- NULL

#Also included is a function for making a results table of all genes for a dataset with multiple comparisons similar to the one for DGE in script 4. Inputs are the same, but instead of choice for log2FoldChange is the choice between including p-adj for each comparison
#See function 5_multiple_comparison_table

#Save fit because it took a long time to make
save(drimds, file = paste0(out_dir, "drimds.Rdata"))
#Save test file because it has counts used in precision. Might not need if I use my modified function 
save(test_asd.control, file = paste0(out_dir, "test_drim_asd_control.Rdata"))
#Save significant genes table
save(sig_drim_asd.control, file = paste0(out_dir,"sig_drim_asd_control.Rdata"))

##
# Annotations and plots ----
##

#To prep annotation table with biomaRt, need organism's gene ensembl
org_ens <- "mmusculus_gene_ensembl"
#Other common species are: drerio_gene_ensembl ; hsapiens_gene_ensembl ; dmelanogaster_gene_ensembl
#Can check other available datasets with:
#listDatasets(useMart("ENSEMBL_MART_ENSEMBL"))

ensembl = useDataset(org_ens, mart=useMart("ENSEMBL_MART_ENSEMBL"))

annot_transcript <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'transcript_length', 'external_transcript_name', 'transcript_biotype'),
                   filters = "ensembl_gene_id",
                   values = rownames(sig_drim_asd.control),
                   mart = ensembl)

#Save transcript annotations
save(annot_transcript, file = paste0(out_dir, "annot_transcript.Rdata"))


#DRIMSeq has a function plotProportions, but it doesn't provide much opportunity to customize, so this repository has modified the source code and provided updated version as a function. See script 5_plot_proportions_DRIMmod for full description of parameters. Some examples are below
#For drim_obj can input either drimds after dmFit function, or any result from dmTest, as the information used for transcript counts and precision will be the same.

#Basic plot, similar to DRIMSeq's function:
plotProp_mod(drim_obj = drimds, 
             gene = rownames(sig_drim_asd.control)[2], 
             group.name = "treatment")
#For comparison, this is the original plot
#Main difference is value assigned to group mean annotations, which are calculated in the original function based on group input to design matrix, and not grouping variable supplied to plot like with modified function
plotProportions(x = drimds, 
                gene_id = rownames(sig_drim_asd.control)[2],
                group_variable = "treatment")

#To add annotations to transcripts & if transcript annotation table is provided, can change title to gene symbol instead:
plotProp_mod(drim_obj = test_asd.control, 
             gene = rownames(sig_drim_asd.control)[2], 
             group.name = "treatment",
             tx_annots = annot_transcript,
             main = "symbol")
#Saved as 5_plotProp_Uty_annots

#Can remove gene annotations and points with group means
plotProp_mod(drim_obj = drimds, 
             gene = rownames(sig_drim_asd.control)[2], 
             group.name = "treatment",
             gp.mean = FALSE,
             gene_annot = FALSE,
             main = "Transcript proportions without annotations")

#See function's script to get further information, including on additional arguments for only plotting a subset of samples, colour options, and choosing whether to/ how to order features and samples

##