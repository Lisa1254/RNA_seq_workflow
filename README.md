# RNA_seq_workflow

This collection of scripts and functions describes an RNA-Seq pipeline from fastq files through analysis.

Intended final product:  
1. Example script for using Salmon to quantify fastq files  
2. Script to import and set up data in preparation for analysis  
3. Script for exploring data and quality control  
4. Differential gene expression analysis using DESeq2
5. Differential transcript usage analysis using DRIMSeq  
6. Gene ontology analysis using topGO  
7. Example script for comparing multiple datasets

For illustration and following along with the example scripts, I've used a selection of samples from the publicly available data from:  

Sharon G, Cruz NJ, Kang DW, Gandal MJ et al. Human Gut Microbiota from Autism Spectrum Disorder Promote Behavioral Symptoms in Mice. Cell 2019 May 30;177(6):1600-1618.e17. PMID: 31150625  

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109827