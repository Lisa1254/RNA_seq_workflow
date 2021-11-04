# Function to visualize gene counts of significant genes from differential gene expression analysis



#Uses ggplot to show count differences of specified genes between groups of interest
#Inputs are dds object for count data and sample information, genes of interest for plotting, annotation frame of gene ensembl id to symbol mapping (as made with biomaRt), and factor variable from sample information for colour of points
#dds object should have already been through DESeq() function, or estimateSizeFactors() in order to use normalized counts in plot
#Default binwidth is set to 1/30 for use in ggplot2's geom_dotplot, but this value may need adjustment depending on counts, and number of genes in plot.

count.plot <- function(dds, sig.genes, annot_frame, col.var, binwidth = 1/30) {
  symbols <- annot_frame[match(sig.genes, annot_frame$ensembl_gene_id),"external_gene_name"]
  count.data <- counts(dds, normalized=TRUE)[sig.genes,]
  sample.data <- as.data.frame(colData(dds))
  col.var.n <- which(colnames(sample.data) == col.var)
  count.frame <- data.frame(gene=rep(sig.genes, ncol(count.data)),
                            symbol=rep(symbols, ncol(count.data)),
                            counts=as.vector(count.data),
                            sample=rep(colnames(count.data), 
                                       each=length(sig.genes)),
                            group=rep(sample.data[,col.var.n], 
                                      each=length(sig.genes)))
  count.frame$symbol <- factor(count.frame$symbol, levels = count.frame[1:length(sig.genes),2])
  
  ggplot(count.frame, mapping = aes(x=symbol, y=counts, fill=group, color=group)) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = binwidth) +
    scale_y_log10() +
    coord_flip() +
    labs(title = "Gene counts by group") +
    labs(y = "Log10 normalized counts")

}

