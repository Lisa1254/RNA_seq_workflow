#Construct dataframe of transcript and gene mappings

#Wrapper of inputs for GenomicFeatures to create AnnotationDb object and select columns of gene and transcript ids
#GenomicFeatures library should be loaded to properly use select command

#Default behaviour to return dataframe of TXNAME and GENEID
#Set return_db=TRUE to return full Genomic Features database and select keys and mappings manually

gtf2txmap <- function(dir, file, organism, return_db = FALSE) {
  gtf_file <- file.path(dir, file)
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, organism = organism)
  if (return_db == TRUE) {
    return(txdb)
  } else if (return_db == FALSE) {
    txdf <- select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")
    return(txdf)
  }
}
