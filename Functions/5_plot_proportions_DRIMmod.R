## Paste & mod of source code for plotProportions (barplot) in the DRIMSeq package

plotProp_mod <- function(drim_obj, gene, sample.sub, group.name, order_features = TRUE, order_samples = TRUE, group_colors = NULL, main = NULL) {
  
  #Get count table for gene of interest
  count.table <- counts(drim_obj[gene,])
  rownames(count.table) <- count.table[,2]
  count.table <- count.table[,-c(1,2)]
  #Identify group factor
  group <- samples(drim_obj)[,group.name]
  
  #Subset by samples of interest
  if(!missing(sample.sub)){
    count.table <- count.table[,sample.sub]
    group <- group[sample.sub]
  }
  
  ## Calculate observed proportions
  proportions <- prop.table(as.matrix(count.table), 2)
  proportions[proportions == "NaN"] <- NA
  
  prop_samp <- data.frame(feature_id = rownames(proportions), proportions, 
                          stringsAsFactors = FALSE)
  
  prop_fit <- NULL
  
  #Skipping this part of fit_full for now, this is the copied code
  #if(!is.null(fit_full))
  #prop_fit <- data.frame(feature_id = rownames(fit_full), fit_full, 
  #stringsAsFactors = FALSE)
  
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
  if(order_samples){
    o <- order(group)
    sample_levels <- colnames(count.table)[o]
  }else{
    sample_levels <- colnames(count.table)
  }
  
  ## Melt prop_samp
  #Note: melt is a function of library(reshape2)
  prop_samp <- melt(prop_samp, id.vars = "feature_id", 
                    variable.name = "sample_id", value.name = "proportion", 
                    factorsAsStrings = FALSE)
  
  prop_samp$feature_id <- factor(prop_samp$feature_id, levels = feature_levels)
  prop_samp$group <- rep(group, each = nrow(count.table))
  prop_samp$sample_id <- factor(prop_samp$sample_id, levels = sample_levels)
  
  #Skipping steps for md (additional sample information)
  #Skipping steps for prop_fit (from fit_full skipped above)
  
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
  
  ## Prepare colors for features
  # Only seems to be in boxplot2, so I will skip this step
  # If I return, will need to fix the colorb function that is internal to the package, as I fixed above
  #if(is.null(feature_colors)){
  #feature_colors <- colorb(nrow(count.table))
  #}
  #names(feature_colors) <- rownames(count.table)
  
  #Adding line to use gene name for main if main not provided
  if (is.null(main)) {
    main <- gene
  }
  
  #Plot
  ggplot() +
    geom_bar(data = prop_samp, aes_string(x = "feature_id", y = "proportion", 
                                          group = "sample_id", fill = "group"), 
             stat = "identity", position = position_dodge(width = 0.9)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.text=element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16), 
          legend.position = "right", 
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 14)) +
    ggtitle(main) +
    scale_fill_manual(name = "Groups", values = group_colors, 
                      breaks = names(group_colors)) +
    xlab("Features") +
    ylab("Proportions")
  
}

