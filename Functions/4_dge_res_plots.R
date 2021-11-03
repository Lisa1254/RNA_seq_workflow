# Function to visualize significant genes from differential gene expression analysis


#Gene counts ----

#count differences of specified genes between groups of interest

count.plot <- function(count.data, sample.data, annot_frame, col.var) {
  
}

count.table <- function(count.data, sample.data, annot_frame) {
  genes <- annot_frame[,1]
  symbols <- annot_frame[,2]
  count.frame <- data.frame(gene=rep(genes, ncol(count.data)),
                            symbol=rep(symbols, ncol(count.data)),
                            counts=as.vector(count.data[genes,]),
                            sample=rep(colnames(count.data), 
                                       each=length(genes)),
                            group=rep(sample.data[,2], 
                                      each=length(genes)))
  count.frame$symbol <- factor(count.frame$symbol, levels = count.frame[1:length(genes),2])
  return(count.frame)
}
annot_sub <- gene_symbol_all[which(gene_symbol_all$ensembl_gene_id %in% resc_genes_update_ens),]
samples <- as.data.frame(colData(gene_dds_mr))
gene_dds_mr <- estimateSizeFactors(gene_dds_mr)
counts <- counts(gene_dds_mr, normalized=TRUE)

prev_resc_table <- count.table(counts, samples, annot_sub)

ggplot(prev_resc_table, mapping = aes(x=symbol, y=counts, fill=group, color=group)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1/30) +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Gene counts for previously identified rescued genes") +
  labs(y = "Log10 normalized counts")

## Volcano plot ----

# This is also a good place to add a volcano plot
ggplot(as.data.frame(res_asd.control), mapping=aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point()