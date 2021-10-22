##
## Data Exploration ----
##

#Data required for this script:
##   > salmon quantifications, at gene and transcript level identification
##   > absolute path to quantifications

##
# Libraries & Functions----
##

#DESeq2 is used to combine technical replicates, and store gene, sample, and count information in a convenient object
library(DESeq2)
#ggplot2 is used to plot graphics for data exploration
library(ggplot2)

#Source import functions
source()


##
# Define Key Variables ----
##

