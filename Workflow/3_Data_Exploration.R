##
## Data Exploration ----
##

#Data required for this script (Outputs from R script "2_Data_Import"):
##   > transcript to gene mapping,
load("Output/txmap_mm.Rdata")
##   > Table of sample information,
load("Output/sample_info.Rdata")
##   > tximport object of scaled transcript-level quantifications, 
load("Output/tx_txi_stpm.Rdata")
##   > tximport object of gene-level quantifications
load("Output/gene_txi.Rdata")


##
# Libraries & Functions----
##

#DESeq2 is used to combine technical replicates, and store gene, sample, and count information in a convenient object
library(DESeq2)
#ggplot2 is used to plot graphics for data exploration
library(ggplot2)
#sva is used for defining batch variables
library(sva)

#Source import functions
source("Functions/PCA_from_dds.R")


##
# Define Key Variables ----
##

#Output directory for saving files
out_dir = paste0(getwd(),"/Output/")

#Pre-filtering minimum counts of gene expression: At least {min.samps} samples will have at least {min.counts} counts for the gene
#These filters are applied to the gene-level counts. The DRIMSeq package for differential transcript usage has more complex filtering that will be done in script 5
min.samps <- 2
min.counts <- 5


##
# Format Data for Downstream Analyses ----
##

#Import gene-level counts to DESeq
#Initializing dataset with no design specified. This will be updated prior to analysis
gene_dds <- DESeqDataSetFromTximport(txi = gene_txi, colData = samples, design = ~1)
#Pre-filter for very low count genes
keep <- rowSums(counts(gene_dds) >= min.counts) >= min.samps
gene_dds <- gene_dds[keep,]

#This data has technical replicates. Before combining the data, preview the distribution of samples to verify that technical replicates are nearly identical
#Using "sample" col in samples dataframe for colour (same colour represents technical replicate, and should cluster together)
#Using "treatment" col in samples dataframe for shape
PCs_plot_dds(gene_dds, color.var = "sample", shape.var = "treatment")

#Technical replicates do present as nearly identical with no outliers, so the samples can be combined
#Use collapseReplicates function to combine runs of same sample:
gene_dds <- collapseReplicates(gene_dds, groupby =  gene_dds$sample, run = gene_dds$run)

#Check collapse did as expected
colnames(gene_dds)
#colnames are from the "groupby" argument above ("sample")
colData(gene_dds)
#Looks good; will remove unnecessary "run" col in colData, since the information is more accurately contained in the new "runsCollapsed" column
gene_dds$run <- NULL

#Factor "group" column of coldata for use in design 
gene_dds$group <- factor(gene_dds$group)
#Reference level is not necessary to set for this analysis, but can relevel to a control sample
gene_dds$group <- relevel(gene_dds$group, "Striatum_control_1")

#Combine technical replicates of transcript-level counts with same method
tx_dds <- DESeqDataSetFromTximport(txi = tx_txi_stpm, colData = samples, design = ~1)
tx_dds <- collapseReplicates(tx_dds, groupby =  tx_dds$sample, run = tx_dds$run)
tx_dds$run <- NULL
tx_dds$group <- factor(gene_dds$group)
tx_dds$group <- relevel(tx_dds$group, "Striatum_control_1")

#Use counts & sample information with combined replicates to produce table in format used for DRIMSeq differential transcript usage
tx_cts <- counts(tx_dds)
samples <- as.data.frame(colData(tx_dds))

#Remove tx version code suffix from transcript id in rowname
rownames(tx_cts) <- gsub("\\.[0-9]*$", "", rownames(tx_cts))
#Only keep transcripts that are mapped to genes
keep <- which(rownames(tx_cts) %in% txmap_mm$TXNAME)
tx_cts <- tx_cts[keep,]
#Set up transcript mapping in same order as counts table
txmap_mm <- txmap_mm[match(rownames(tx_cts),txmap_mm$TXNAME),]

#Set up dataframe of counts as required by DRIMSeq
#DRIMSeq for DTU requires first 2 columns of count matrix as "gene_id" and "feature_id"
tx_cts <- data.frame(gene_id=txmap_mm$GENEID,
                          feature_id=txmap_mm$TXNAME,
                          tx_cts)
#Sample dataframe for DRIMSeq requires a column labelled "sample_id" with same sample names as used in columns of count matrix
colnames(samples)[1] <- "sample_id"

## Save files
save(gene_dds, file=paste0(out_dir, "gene_dds.Rdata"))
save(tx_cts, samples, file=paste0(out_dir, "tx_ct_colRepsamp.Rdata"))


## 