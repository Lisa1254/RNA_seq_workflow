# RNA_seq_workflow

This collection of scripts and functions describes an RNA-Seq pipeline from fastq files through analysis.

## Intended final product:  
1. Example script for using Salmon to quantify fastq files  
2. Script to import and set up data in preparation for analysis  
3. Script for exploring data and quality control, and formatting data in preparation of downstream analysis  
4. Differential gene expression analysis using DESeq2
5. Differential transcript usage analysis using DRIMSeq  
6. Gene ontology analysis using topGO  
7. Example script for comparing multiple datasets, or multiple attributes in a single dataset

## TO DO:  
Functions - verify that functions from specific packages are identified using double colon format (package::function)  
Functions - switch to internal variables assigned with "=" instead of "<-" as per convention  
All scripts - read through for consistency and clarity  
  


## Data  

For illustration and following along with the example scripts, I've used a selection of samples from the publicly available data from:  

Sharon G, Cruz NJ, Kang DW, Gandal MJ et al. Human Gut Microbiota from Autism Spectrum Disorder Promote Behavioral Symptoms in Mice. Cell 2019 May 30;177(6):1600-1618.e17. PMID: 31150625  
  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109827
  
This data was selected because the samples include multiple runs on the same source RNA data to demonstrate combining technical replicates, and nested design of multiple factors. Original data includes multiple tissue types (Prefrontal cortex & striatum), but only striatum samples are included here.
  
Ensembl ftp site was used for current genome build (GRCm39):  
Transcriptome:  
http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz  
GTF:  
http://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz  
