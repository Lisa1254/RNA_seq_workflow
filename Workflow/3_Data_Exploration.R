##
## Data Exploration ----
##

#Data required for this script (Outputs from R script "2_Data_Import"):
##   > transcript to gene mapping
load("Output/txmap_mm.Rdata")
##   > Table of sample information
load("Output/sample_info.Rdata")
##   > tximport object of scaled transcript-level quantifications (NOTE: this is a large file, and is not included with saved output files) 
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
#Complex heatmap is used to draw sample distance heatmap, with circlize being used for color function in the heatmap
library(ComplexHeatmap)
library(circlize)

#Source import functions
#Function to plot samples by principle components from DESeq data set, or to return PC loadings. 
source("Functions/PCA_from_dds.R")
#Functions to modify and plot count data according to batch control with Surrogate Variable Analysis
source("Functions/SV_data_mods.R")
#Draw heatmap of relationships between samples
source("Functions/samp_dist_heatmap.R")


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

#Format gene-level counts to DESeq
#Initializing dataset with no design specified. This will be updated prior to analysis
gene_dds <- DESeqDataSetFromTximport(txi = gene_txi, colData = samples, design = ~1)
#Pre-filter for very low count genes
keep <- rowSums(counts(gene_dds) >= min.counts) >= min.samps
gene_dds <- gene_dds[keep,]

#This data has technical replicates. Before combining the data, preview the distribution of samples to verify that technical replicates are nearly identical
#Using "sample" col in samples dataframe for colour (same colour represents technical replicate, and should cluster together)
#Using "treatment" col in samples dataframe for shape
PCs_plot_dds(gene_dds, color.var = "sample", shape.var = "treatment")

#Figure saved as 3_PC_tech_reps

#Technical replicates do present as nearly identical with no outliers, so the samples can be combined

## If your data has no technical replicates to combine, ##
## skip to line 79 for continuing data format ##


#Use collapseReplicates function to combine runs of same sample:
gene_dds <- collapseReplicates(gene_dds, groupby =  gene_dds$sample, run = gene_dds$run)

#Check collapse did as expected
colnames(gene_dds)
#colnames are from the "groupby" argument above ("sample")
colData(gene_dds)
#Looks good; will remove unnecessary "run" col in colData, since the information is more accurately contained in the new "runsCollapsed" column
gene_dds$run <- NULL

## Return to script here if no technical replicates ##

#Factor "group" column of coldata for use in design 
gene_dds$group <- factor(gene_dds$group)
#Reference level is not necessary to set for this analysis, but can relevel to a control sample
gene_dds$group <- relevel(gene_dds$group, "CON_1")

#Format transcript level counts using DESeq2
tx_dds <- DESeqDataSetFromTximport(txi = tx_txi_stpm, colData = samples, design = ~1)

#If no technical replicates to combine, skip to line 95 ##

#Combine technical replicates of transcript-level counts with same method
tx_dds <- collapseReplicates(tx_dds, groupby =  tx_dds$sample, run = tx_dds$run)
tx_dds$run <- NULL

## Return here if not combining technical replicates ##
tx_dds$group <- factor(tx_dds$group)
tx_dds$group <- relevel(tx_dds$group, "CON_1")

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
save(tx_cts, samples, file=paste0(out_dir, "tx_ct_formatted.Rdata"))


## 
# Batch Correction ----
##

#If your data has known batches (e.g. multiple technicians, data collection on different dates), you can account for them in the design formula (e.g. design = ~ group + batch)
#If batch influences are expected from the experimental design or detected in the exploratory PCA plots, but not specifically identified, the sva package can estimate surrogate variables to be used in the design

#Define full model matrix of intended design
#In this case group variable is being used rather than treatment because of the nested design - "Control" and "ASD" samples are not fully independent, so relationship between multiple inputs needs to be accounted for
#Since there are more than 2 levels to the group factor, the full model matrix is defined without intercept by using -1
mm.full <- model.matrix(~ -1 + group, colData(gene_dds))
#Define null model matrix
mm.null <- model.matrix(~ 1, colData(gene_dds))

#SVA-Seq method
#Uses normalized counts matrix:
gene_dds <- estimateSizeFactors(gene_dds)
norm.cts <- counts(gene_dds, normalized=TRUE) 

#Construct svaseq fit
#As the original paper used a sequence metric PC analysis to include PC1 and PC2 as batch in the design, the number of surrogate variables defined here will be 2 (n.sv=2).
#Remove the n.sv argument to let the function determine the number of surrogate variables.
fit.sva <- svaseq(norm.cts, mod=mm.full, mod0=mm.null, n.sv = 2)

#Integrate surrogate variables to data set
gene_dds$SV1 <- fit.sva$sv[,1]
gene_dds$SV2 <- fit.sva$sv[,2]

#Update design matrix
design(gene_dds) <- ~ -1 + SV1 + SV2 + group

# Use provided SV plot to check if samples cluster according to a known variable
#If there is an expected batch to emerge with the surrogate variables, use the color or shape argument to verify clustering.
#Clustering is not expected to occur by the experimental treatment variable
#This data set does not have supplied information for technical variation that would be expected to cluster, so this plot is mostly for demonstrating its use
SVs_plot(gene_dds, color = "group", shape = "treatment")


#Check clustering of experimental treatment on data that has been adjusted with the surrogate variables ("cleaned" for batch)

#Use provided function to adjust count data for surrogate variables
adj_cts <- rem_SVs(norm.counts = norm.cts, mm.full, fit.sva)

#Use provided PCs_plot function to visualize adjusted data
PCs_plot(cts = adj_cts, sample.data = colData(gene_dds), color.var = "group", shape.var = "treatment", main = "PCA plot from SVA adjusted data")
#Saved as 3_PC_sva_adj

#And for comparison to the unadjusted count data:
PCs_plot_dds(gene_dds, color.var = "group", shape.var = "treatment", main = "PCA from VST of raw count data")
#The differences in these plots show that the SVs are acting to enhance the connection between samples that used the same bacterial donor. The high percentage of variance explained by a single principle component in the adjusted data also demonstrates why adjusted counts should only be considered in data exploration and not for use in analysis.

## Save files
#gene_dds will be used in script 4_DGE_DESeq2, and is being saved as gene_dds_sva to distinguish from the DESeq data set object saved with replicates collapsed, but no surrogate variable information.
save(gene_dds, file=paste0(out_dir, "gene_dds_sva.Rdata"))



##
# Exploratory Plots ----
##

## Add some additional useful & common plots for exploration
## Section in progress - still need to streamline code & generalize further where possible
##Update comments

#Library Size
#Show variation in library sizes, and impact of normalization
#Set to display two plots side by side
par(mfrow = c(1, 2))
#Counts log-transformed, and not normalized
boxplot(counts(gene_dds)+1, xlab="", log = "y", las = 2, col="lightblue", cex.axis = 0.75)
abline(h=median(counts(gene_dds)),col="blue")
title("Boxplots of log counts (unnormalised)", cex.main = 0.75)
#Counts log-transformed, and normalized with DESeq2's RLE method
boxplot(norm.cts+1, xlab="", log = "y", las = 2, col="darkred", cex.axis = 0.75)
abline(h=median(norm.cts),col="blue")
title("Boxplots of log counts (normalized: RLE)", cex.main=0.75)
#Reset plotting to default parameters
par(mfrow = c(1,1))

#Saved as 3_box_lib_size

#Heatmap

#Using provided function to visualize sample distances
sample_dist_heatmap(gene_dds)
#Saved as 3_samp_dist
#Samples from the same donor individual are similar to each other, control samples are similar to each other, but there is variation in the ASD samples. This distribution could indicate that there are few specific differences between the experimental groups, and may also reflect the subsampling from the original data

#Possibly add scatterplot

##
# Test vidger package ----
##

#Might delete, in progress to see
BiocManager::install("vidger")
library(vidger)

vsScatterPlot(x="CON", y="ASD", data=gene_dds, d.factor="treatment", type="deseq")
vsScatterMatrix(data=gene_dds, d.factor="group", type="deseq")

#Actually really like this. Super easy, but not a lot of options to modify or colour

##
