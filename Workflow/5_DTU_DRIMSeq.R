##
## Differential Transcript Usage Analysis with DRIMSeq ----
##

#Data required for this script:
##   > transcript-level counts & sample information in DRIMSeq format (output of script 3_Data_Exploration)
load("Output/tx_ct_formatted.Rdata")


##
# Libraries & Functions----
##

#DRIMSeq is used for differential transcript usage analysis.
library(DRIMSeq)


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

#Took 91.745 seconds
#Total time expected to be
(91.745/60)*(length(drimds)/200)
#68 minutes expected

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
coef.control <- ifelse(grepl("control", colnames(design(drimds))), 1, 0)
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

#Save fit because it took a long time to make
save(drimds, file = "Output/drimds.Rdata")
#Save test file because it has counts used in precision. Might not need if I use my modified function 
save(test_asd.control, file = "Output/test_drim_asd_control.Rdata")
#Save significant genes table
save(sig_drim_asd.control, file = "Output/sig_drim_asd_control.Rdata")

##
# Annotations and plots ----
##

# In progress still, currently includes rough pasted code from previous scripts

#Can get transcript annotations from biomaRt (pasted from metabolite rescue exploration)
library("biomaRt")
ensembl_dr = useDataset("drerio_gene_ensembl", mart=useMart("ENSEMBL_MART_ENSEMBL"))

annot_tr1 <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 'transcript_length', 'external_transcript_name', 'transcript_biotype'),
                   filters = "ensembl_gene_id",
                   values = all_drim_mr,
                   mart = ensembl_dr)

#Can use transcript annotations and plotProportions (DRIMSeq version, or my modified code) to plot and compare transcript usage

