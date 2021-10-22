#!/bin/bash

## This script uses Salmon to quantify transcripts from paired-end reads
## Time request for 17 pairs (index previously done) on graham.computecanada took ~40 min, so adjust time request according to your system & files
## Add account & user before running. 
## As set up, have this script in a parent directory that has the following: 
## >"Read_data" folder with fastq files,
## >"Transcript_index" as an empty folder for index output,
## >"Salmon_quants" as an empty folder for Salmon outputs,
## >fasta of transcriptome for indexing
## run with command: sbatch 1_Salmon.sh


#SBATCH --account=def-{add-account}
#SBATCH --time=0-05:00:00 ## days-hours:minutes:seconds
#SBATCH -c 8
#SBATCH --mem=32000 # requested memory (in MB)
#SBATCH --mail-type=ALL
#SBATCH --mail-user={add-user-email}


#Load environment for Salmon application
#Can verify modules to load prior to running script with command
# module spider salmon/1.4.0

module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 salmon/1.4.0

#Create transcript index
#Replace {fasta_file} with file location
#Example, Danio_rerio.GRCz11.cdna.all.fa.gz
#Can browse ftp site of Ensembl for species of interest
#https://useast.ensembl.org/info/data/ftp/index.html
#Copy link and use wget prior to running this script to transfer file
#wget http://ftp.ensembl.org/pub/release-104/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz

salmon index -t {fasta_file} -i Transcript_index


#Quantification
#This assumes files are gzipped fastq files (end with .fastq.gz)
#Filenames for paired reads were identical, except R1/R2
#Filenames used in example ended in "_R1_001.fastq.gz" and _R2_001.fastq.gz
#Modify to suit your data's naming convention

#Flags chosen & description from the documentation:
#https://salmon.readthedocs.io/en/latest/
# > "-l A" indicates that the library type should be automatically detected
# > "validateMappings" can improve both the sensitivity and specificity of 
##mapping and, as a result, can improve quantification accuracy
# > "seqBias" learn and correct for sequence-specific biases in the input data. 
##Specifically, this model will attempt to correct for random hexamer priming 
##bias, which results in the preferential sequencing of fragments starting with 
##certain nucleotide motifs.
# > "gcBias" will attempt to correct for biases in how likely a sequence is to 
##be observed based on its internal GC content.
## When both --seqBias and --gcBias are enabled, Salmon will learn a 
##conditional fragment-GC bias model

for each in Read_data/*fastq*
do
 file=`echo $each | sed 's/Read_data\///'`
 if [[ $file == *"R1"* ]]
 then
  basename=`echo $file | sed 's/_R1_001\.fastq\.gz//'`
  echo "Aligning pair: ${basename}_R1 and ${basename}_R2"
  salmon quant -i Transcript_index -l A \
        -1 Read_data/${basename}_R1_001.fastq.gz \
        -2 Read_data/${basename}_R2_001.fastq.gz \
        -p 8 --validateMappings \
        --seqBias --gcBias -o Salmon_quants/${basename}_quant
 fi
done

