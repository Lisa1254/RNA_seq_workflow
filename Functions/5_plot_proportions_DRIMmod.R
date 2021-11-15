## Paste & mod of source code for plotProportions (barplot) in the DRIMSeq package

#Libraries are DRIMSeq for input objects, ggplot2 for figure, reshape for formatting data
#Required inputs are:
#   > drim_obj = DRIMSeq Fit object
#   > gene for plotting
#   > group.name for sample information variable to be used as colour
#Optional inputs are:
#   > sample.sub = vector of sample names to include in plot
#   > tx_annots = annotation table of transcript information as made with biomaRt to provide description and size annotation information to transcript labels
#   > order_features, default = TRUE, to order transcripts in plot by decreasing highest mean proportion
#   > order_samples, default = "default", will order samples within each transcript by same group as identified for colours. Alternatives include "none" for using same order as samples in DRIMSeq count table, or any other character string that matches a factor variable of sample information
#   > group_colors, default is to use 24 colours provided sequentially, repeated if necessary, or can provide
#   > main, default is to use gene ID. Alternatives include providing character vector of choice, or, if tx_annots is provided, can specify "symbol" to use gene symbol

plotProp_mod <- function(drim_obj, gene, sample.sub, group.name, tx_annots, order_features = TRUE, order_samples = "default", group_colors = NULL, main = "gene", gp.mean = TRUE, gene_annot = TRUE) {
  
  #Get count table for gene of interest
  count.table <- counts(drim_obj[gene,])
  feat_names <- count.table[,2]
  
  if (!missing(tx_annots)) {
    feat_descrip <- tx_annots[match(feat_names, tx_annots$ensembl_transcript_id),"transcript_biotype"]
    feat_size <- tx_annots[match(feat_names, tx_annots$ensembl_transcript_id),"transcript_length"]
    feat_names_annot <- paste0(feat_names, ": ", feat_descrip, " (", feat_size, " bp)")
    rownames(count.table) <- feat_names_annot
  } else {
    rownames(count.table) <- feat_names
  }
  
  
  
  count.table <- count.table[,-c(1,2)]
  #Identify group factor
  group <- samples(drim_obj)[,group.name]
  
  #Subset by samples of interest
  if(!missing(sample.sub)){
    count.table <- count.table[,sample.sub]
    group <- group[sample.sub]
    group <- factor(as.character(group))
  }
  
  ## Calculate observed proportions
  proportions <- prop.table(as.matrix(count.table), 2)
  proportions[proportions == "NaN"] <- NA
  
  ## Add gene annotations for subtitle
  if (gene_annot) {
    #Define genewise precision
    gwp <- genewise_precision(drim_obj)[which(genewise_precision(drim_obj)$gene_id == gene),2]
    gwp <- round(gwp, 3)
    mn_ct <- round(mean(apply(count.table,2,sum)))
    sub_string <- paste0("Genewise precision: ", gwp, "  Mean expression: ", mn_ct)
  }

  
  ## Set up table
  prop_samp <- data.frame(feature_id = rownames(proportions), proportions, 
                          stringsAsFactors = FALSE)
  
  
  ## Order transcripts by decreasing proportion
  # This step takes the median expression per transcript per group factor, then sorts by order of max proportion in the transcript in any group
  if(order_features){
    oo <- order(apply(aggregate(t(prop_samp[, -1]), 
                                by = list(group = group), median)[, -1], 2, max), decreasing = TRUE)
    feature_levels <- rownames(prop_samp)[oo]  
  }else{
    feature_levels <- rownames(count.table)
  }
  
  ## Order samples by group
  if(order_samples == "default"){
    o <- order(group)
    sample_levels <- colnames(count.table)[o]
  }else if (order_samples == "none"){
    sample_levels <- colnames(count.table)
  } else {
    order.group <- samples(drim_obj)[,order_samples]
    o <- order(order.group)
    sample_levels <- colnames(count.table)[o]
  }
  
  ## Melt prop_samp
  prop_samp <- reshape2::melt(prop_samp, id.vars = "feature_id", 
                    variable.name = "sample_id", value.name = "proportion", 
                    factorsAsStrings = FALSE)
  
  #Add group variable
  prop_samp$group <- rep(group, each = nrow(count.table))
  
  # Add group-wise mean for each transcript
  if (gp.mean) {
    df_means <- data.frame(feat_id = feat_names)
    for (var in seq(1,length(levels(group)))) {
      mean_col <- apply(proportions[,which(group == levels(group)[var])], 1, mean)
      df_means <- cbind(df_means, mean_col)
      colnames(df_means)[var+1] <- levels(group)[var]
    }
    
    prop_vec <- vector()
    for (each in seq(1:length(prop_samp$sample_id))) {
      col.ind <- which(colnames(df_means) == prop_samp$group[each])
      row.ind <- which(rownames(df_means) == prop_samp$feature_id[each])
      prop_vec <- c(prop_vec, df_means[row.ind, col.ind])
    }
    
    prop_samp$gp_mean <- prop_vec
  }

  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = feature_levels)
  
  prop_samp$sample_id <- factor(prop_samp$sample_id, levels = sample_levels)
  
  ## Prepare colors for groups
  if(is.null(group_colors)){
    clrs <- c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1" ,
              "blueviolet", "firebrick2", "deepskyblue",  "orchid2", 
              "chartreuse3", "gold", "slateblue1", "tomato" , "blue", 
              "magenta", "green3", "yellow", "purple3", "red" ,
              "darkslategray1", "lightpink1", "lightgreen", "khaki1", 
              "plum3", "salmon")
    nc <- length(clrs)
    n <- nlevels(group)
    if (n > nc) {
      clrs <- rep(clrs, ceiling(n/nc))
    }
    group_colors <- clrs[1:n]
  }
  names(group_colors) <- levels(group)
  
  
  #Adding line to use gene name for main if main not provided
  if (main == "gene") {
    main <- gene
  } else if (main == "symbol") {
    main <- tx_annots[which(tx_annots$ensembl_gene_id == gene), "external_gene_name"][1]
  }
  
  #Plot
  prop_plot <- ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", 
                                          group = "sample_id", fill = "group"), 
             stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.text=element_text(size=12), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16), 
          legend.position = "right", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    ggtitle(main) +
    scale_fill_manual(name = group.name, values = group_colors, 
                      breaks = names(group_colors)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    xlab("Features") +
    ylab("Proportions")
    
  
  if (gp.mean) {
    prop_plot <- prop_plot + 
      geom_point(data = prop_samp, 
                 mapping = aes_string(x="feature_id", y = "gp_mean", 
                                      group = "sample_id", fill = "group"), 
                 position = position_dodge(width = .9),
                 shape = 23, size = 3,
                 color = "black")
  }
  
  if (gene_annot) {
    prop_plot + labs(subtitle = sub_string)
  } else {
    prop_plot
  }
  
}


