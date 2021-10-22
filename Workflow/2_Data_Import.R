##
## Data Import ----
##

#Data required for this script:
##   > salmon quantifications, at gene and transcript level identification
##   > absolute path to quantifications
##   > dataframe of Transcript to Gene mappings OR GTF file corresponding to organism
##   > dataframe of sample information, including groups of interest such as batch and condition

##
# Libraries & Functions----
##

#tximport is used to organize .sf  files from Salmon into count information
library(tximport)
#GenomicFeatures is used to translate the gtf file into a transcript to gene mapping
library(GenomicFeatures)
#biomaRt is used for collecting gene annotations that will be useful downstream
library(biomaRt)

#Source import functions
source("Functions/tximport_salmon.R")
source("Functions/gtf_2_txmap.R")

##
# Define Key Variables ----
##

#By defining these variables at the beginning of the script, the rest can be followed through with minimal modifications

#Using location of files for Sharon et al's study as example
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109827
#Parent directory
dir <-"/Users/lisa/Documents/Bioinformatics/6999/Data/"
#Subdirectory with Salmon count files. Only count files intended for import should be in this folder
src <- "Sharon/Salmon_quants/"
#Location of gtf file within parent directory {dir}
gtf <- "Mus_musculus.GRCm39.104.gtf"

#Output directory
out_dir = "/Users/lisa/Documents/Bioinformatics/6999/Results_gen/"

#For RNA-Seq data acquired from NCBI's SRA, there will be a metadata file that can be downloaded as a csv to use for sample information
#Location of metadata for sample information
metadata <- "Sharon/Sharon_SraRunTable.txt"

#If preparing annotation tables through biomaRt, specify organism's gene ensembl
org_ens <- "mmusculus_gene_ensembl"
#Other common species are: drerio_gene_ensembl ; hsapiens_gene_ensembl ; dmelanogaster_gene_ensembl
#Can check other available datasets with:
#listDatasets(useMart("ENSEMBL_MART_ENSEMBL"))

##
# Basic Import ----
##

#Set up sample information

# Use metatdata table downloaded from SRA
samples <- read.csv(paste0(dir, metadata))
#Many columns contain information not required downstream, such as the organism, or sequencing platform.
#Check which cols have unique information to retain
apply(samples,2,unique)

#Subset table for desired info
samples <- samples[,c(1,24,20,29,26)]
#Define row and column names to something simple and descriptive
rownames(samples) <- samples$Run
colnames(samples) <- c("run", "sample", "treatment", "tissue", "group")

#Clean content of cells, simplifying descriptions used
samples[grep("Control",samples$treatment),3] <- "Control"
samples[grep("ASD",samples$treatment),3] <- "ASD"
samples[grep("Striatum",samples$tissue),4] <- "Striatum"
samples[grep("cortex",samples$tissue),4] <- "PFC"
samples$group <- gsub(" ", "_", samples$group)
samples$group <- gsub("PrefrontalCortex", "PFC", samples$group)

#Preview table
View(samples)

#Construct gene to transcript mapping
#Be sure to include accurate organism
txmap_mm <- gtf2txmap(dir=dir, file=gtf, organism = "Mus musculus")

#Import transcripts
#Using default of scaledTPM for counts from abundance, as recommended for downstream differential transcript usage analysis
tx_txi_stpm <- salmon2tx(paste0(dir, src), samples$run)

#Import counts at gene-level without scaling
gene_txi <- salmon2gene(paste0(dir, src), samples$run, tx2gene = txmap_mm)


#If you want to save any files, use or modify the following:
#To save Rdata files:
save(samples, file = paste0(out_dir, "sample_info.Rdata"))
save(txmap_mm, file=paste0(out_dir, "txmap_mm.Rdata"))
save(tx_txi_stpm, file=paste0(out_dir, "tx_txi_stpm.Rdata"))
save(gene_txi, file=paste0(out_dir, "gene_txi.Rdata"))

#To save csv of charts:
write.csv(samples, file=paste0(out_dir, "sample_info.csv"))


##
# Other Useful Tables to Prep ----
##

#Identify ensembl being used
ensembl <- useDataset(org_ens, mart=useMart("ENSEMBL_MART_ENSEMBL"))

#Prepare table of GO annotations associated with genes in dataset
#This will be used in downstream topGO gene ontology analysis
go_annotations <- getBM(attributes=c('go_id', 'ensembl_gene_id', 
                                          'external_gene_name', 'namespace_1003'),
                             filters='ensembl_gene_id',
                             values=rownames(gene_txi$counts), mart=ensembl)

#Set up table of gene symbol and description mapping for all genes in set
#This will be used to help interpret results of differential gene expression and differential transcript usage analyses, as well as facilitate comparisons between different datasets
gene_symbols <- getBM(attributes=c('ensembl_gene_id', 
                                             'external_gene_name', 'description'),
                                filters='ensembl_gene_id',
                                values=rownames(gene_txi$counts), mart=ensembl)

#Save files:
save(go_annotation, file = paste0(out_dir, "go_annotations.Rdata"))
save(gene_symbols, file = paste0(out_dir, "gene_symbols.Rdata"))