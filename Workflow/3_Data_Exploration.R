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

#Source import functions
#Function to plot samples by principle components from DESeq data set, or to return PC loadings. 
source("Functions/PCA_from_dds.R")
#Functions to modify and plot count data according to batch control with Surrogate Variable Analysis
source("Functions/SV_data_mods.R")


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
gene_dds$group <- relevel(gene_dds$group, "Striatum_control_1")

#Format transcript level counts using DESeq2
tx_dds <- DESeqDataSetFromTximport(txi = tx_txi_stpm, colData = samples, design = ~1)

#If no technical replicates to combine, skip to line 95 ##

#Combine technical replicates of transcript-level counts with same method
tx_dds <- collapseReplicates(tx_dds, groupby =  tx_dds$sample, run = tx_dds$run)
tx_dds$run <- NULL

## Return here if not combining technical replicates ##
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
#And for comparison to the unadjusted count data:
PCs_plot_dds(gene_dds, color.var = "group", shape.var = "treatment", main = "PCA from VST of raw count data")
#The differences in these plots show that the SVs are acting to enhance the connection between samples that used the same bacterial donor. The high percentage of variance explained by a single principle component in the adjusted data also demonstrates why adjusted counts should only be considered in data exploration and not for use in analysis.

## Save files
#gene_dds will be used in script 4_DGE_DESeq2, and is being saved as gene_dds_sva to distinguish from the DESeq data set object saved with replicates collapsed, but no surrogate variable information.
save(gene_dds, file=paste0(out_dir, "gene_dds_sva.Rdata"))



##
# Exploratory Plots ----
##

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




## Add some additional useful & common plots for exploration
## Suggestions include: 
##   > Scatterplot matrix
##   > Heatmaps
##   > See section 4.2 of Love's RNA-Seq workflow for examples of Mean SD plot, & 4.6 for MDS plot.


##
# Test of BigPint package ----
##
#Testing out package bigPint for visualizations
#I don't think I like this package, will probably delete.
library(bigPint)
#See full manual and help at 
#https://lindsayrutter.github.io/bigPint

#Format counts
bp_data <- data.frame(ID=rownames(gene_dds), assay(gene_dds))
#First column of count dataframe must be named ID and include unique identifiers
#Remainder of columns represent samples, and must be named according to the Perl expression ^[a-zA-Z0-9]+\\.[0-9]+
#that is, an alpha-numeric code for treatment type (e.g. ASD or CON) followed by "." as delimiter, and number for replicate
# e.g. ASD.1
#Going to only compare ASD to CON for the purposes of this example, and neglect nested design with donor. Might come back to change.
colnames(bp_data)
bp_samp_names <- colData(gene_dds)$treatment

c <- 1
a <- 1
while ((c<5) & (a<5)) {
  for (name in 1:length(bp_samp_names)) {
    if (bp_samp_names[name] == "Control") {
      bp_samp_names[name] <- paste("CON", c, sep = ".")
      c <- c + 1
    } else {
      bp_samp_names[name] <- paste("ASD", a, sep = ".")
      a <- a + 1
    }
  }
}

colnames(bp_data)[-1] <- bp_samp_names

coldata = data.frame(row.names = colnames(bp_data)[-1], treatment = unlist(lapply(colnames(bp_data)[-1], function (x) unlist(strsplit(x, "[.]"))[1])))

dds = DESeqDataSetFromMatrix(countData = bp_data[,-1], colData = coldata,
                             design = ~ treatment)
dds <- DESeq(dds)

uTreat = unique(unlist(lapply(colnames(bp_data)[-1], function (x) unlist(strsplit(
  x, "[.]"))[1])))
sub_metrics <- list()

for (i in 1:(length(uTreat)-1)){
  for (j in (i+1):length(uTreat)){
    res <- results(dds, contrast=c("treatment", uTreat[i], uTreat[j]))
    metrics = as.data.frame(res@listData)
    metrics = cbind(ID = res@rownames, metrics)
    metrics$ID = as.character(metrics$ID)
    metrics <- metrics[order(metrics$padj), ]
    sub_metrics[[paste0(uTreat[i], "_", uTreat[j])]] <- metrics
  }
}

for (df in seq_len(length(sub_metrics))){
  whichPadj = which(colnames(sub_metrics[[df]])=="pvalue")
  colnames(sub_metrics[[df]])[whichPadj] = "PValue"
  whichPadj = which(colnames(sub_metrics[[df]])=="padj")
  colnames(sub_metrics[[df]])[whichPadj] = "FDR"
  whichPadj = which(colnames(sub_metrics[[df]])=="log2FoldChange")
  colnames(sub_metrics[[df]])[whichPadj] = "logFC"
}

str(sub_metrics, strict.width = "wrap")
names(sub_metrics)

library(dplyr)
tenSigGenes <- sub_metrics[["CON_ASD"]] %>% select(ID) %>%
  filter(row_number() <= 10)
tenSigGenes <- tenSigGenes[,1]
bp_data[,-1] <- log(bp_data[,-1] + 1)
ret <- plotLitre(data=bp_data, geneList = tenSigGenes, saveFile = FALSE)
names(ret)

ret[["CON_ASD_ENSMUSG00000094065"]]

ret2 <- plotSM(bp_data, sub_metrics, pointColor = "pink",
              saveFile = FALSE)
names(ret2)




##