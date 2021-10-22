# Wrapper functions for using tximport to import salmon counts

#ignoreTxVersion has been set to TRUE
#dir = filepath to salmon quantification files
#samp_names = vector of names in same order as sf files provided
#If samp_names is not provided, file name will be used for names, with ".sf" removed
#tx2gene is dataframe with columns TXNAME and GENEID, as made with GenomicFeatures package

#Import counts at gene level ----
#Import done without scaling. Can be scaled later in analysis
salmon2gene <- function(dir, samp_names, tx2gene) {
  lsf <- list.files(dir)
  files <- file.path(dir, lsf)
  if (missing(samp_names)) {
    names(files) <- gsub("\\.sf", "", lsf)
  } else {
    names(files) <- samp_names
  }
  if (all(file.exists(files))) {
    txi_tx <- tximport::tximport(files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
    return(txi_tx)
  } else {
    cat("All files do not exist. Check filepath.")
  }
}

#Import counts at transcript level ----
#Scaled-TPM recommended for DTU analysis
salmon2tx <- function(dir, samp_names, cts = "scaledTPM") {
  lsf <- list.files(dir)
  files <- file.path(dir, lsf)
  if (missing(samp_names)) {
    names(files) <- gsub("\\.sf", "", lsf)
  } else {
    names(files) <- samp_names
  }
  if (all(file.exists(files))) {
    txi_tx <- tximport::tximport(files, type="salmon", countsFromAbundance = cts, txOut = TRUE, ignoreTxVersion = TRUE)
    return(txi_tx)
  } else {
    cat("All files do not exist. Check filepath.")
  }
}
