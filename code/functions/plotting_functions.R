#Ridge plot
plot_ridge.fx <- function(mydf, var){  
  
  ridgeline.plot <- ggplot() +
    geom_density_ridges(data = mydf, aes(x = eval(parse(text = var)), y = immune_cluster, fill = immune_cluster),
                        scale = 2) + 
    myplot +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(size = 30, color = "black"),
          axis.title = element_blank(), 
          plot.title = element_text(size = 30, hjust = 0.5),
          legend.position = "none") + #coord_flip() + #flip axes
    scale_fill_manual(values = cluster_col) + scale_y_discrete(expand = c(0.1,0))
  
  return(ridgeline.plot) 
}

# stats for ridgeplots for Fig3B
hallmark_IC_stats_ridge <- function(proteinmat, metadata, hallmark){
  
  protein_mat_pathway <- as.matrix(proteinmat[rownames(proteinmat) %in% Hs.H[[hallmark]], rownames(metadata)])
  print(dim(protein_mat_pathway))
  pathway_mean <- colMeans(protein_mat_pathway)
  
  metadata$pathway_mean <- pathway_mean
  
  print(pairwise.t.test(metadata$pathway_mean, metadata$immune_cluster, p.adjust = "bonferroni",
                        pool.sd = FALSE, paired = FALSE))
  
  ridge_plot <- plot_ridge.fx(metadata, "pathway_mean")
  return(ridge_plot)
}


#density plot for tpm values
densplot.fx <- function(dat, pheno){
  groups_df <- as.data.frame(t(dat))  
  matchsamples <- match(rownames(groups_df), rownames(pheno))
  groups_df$group <- NA
  groups_df$group <- pheno$group[matchsamples] 
  print(table(groups_df$group))  
  groups_df_m <- reshape2::melt(groups_df)    
  densplot <- ggplot(data = groups_df_m) + geom_density(aes(value, group = group, color = group)) + 
    myplot + myaxis  
  
  return(densplot)  
}



# Stacked barplots for cancer subtypes
subgroupcount_IC.fx <- function(metadata, tumour){
  #Subset to tumour
  mytumour <- metadata[metadata$cohort == tumour,]
  message(tumour)
  # Remove ND subtypes
  mytumour <- mytumour[!grepl("ND", mytumour$tumour_subtype),]
  # Create tables   
  subtype_tab <- as.data.frame(table(mytumour$tumour_subtype), stringsAsFactors = F)
  # Remove ND subtypes
  subtype_tab <- subtype_tab[!grepl("ND", subtype_tab$Var1),]
  
  tumour_tab <- mytumour %>% group_by(tumour_subtype,immune_cluster) %>% summarise(n = n()) %>% 
    mutate(freq = n / sum(n))  
  
  #Count plot  
  subtype_tab$Var1 <- factor(subtype_tab$Var1, levels = subtype_tab$Var1[order(subtype_tab$Freq, decreasing = T)])
  
  subtype_count_plot <- ggplot(subtype_tab,aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", width = 0.8) + myaxis + myplot +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20, hjust = 0.5)) + 
    labs(y = "Frequency")
  
  return(subtype_count_plot)
}


subgroupfreq_IC.fx <- function(metadata, tumour){
  #Subset to tumour
  mytumour <- metadata[metadata$cohort == tumour,]
  message(tumour)
  # Remove ND subtypes
  mytumour <- mytumour[!grepl("ND", mytumour$tumour_subtype),]
  # Create tables   
  subtype_tab <- as.data.frame(table(mytumour$tumour_subtype), stringsAsFactors = F)
  
  tumour_tab <- mytumour %>% group_by(tumour_subtype,immune_cluster) %>% summarise(n = n()) %>% 
    mutate(freq = n / sum(n))
  print(subtype_tab)
  #Fisher test for each subtype and immune cluster
  mysubtypes <- subtype_tab$Var1
  myk <- c("Pediatric Inflamed", "Myeloid Predominant", "Immune Neutral", "Immune Excluded")
  # if subtype is the subtype of interest == 0, if not 1. This is just to get a table with easier interpretation
  # ie. The test shows the odds of samples in the subtype of interest being in the cluster of interest (first cell in the contingency table)
  for( type in mysubtypes){
    message(type)
    mytumour$subtypegroup <- NA
    mytumour$subtypegroup[ mytumour$tumour_subtype ==  type] <- 0
    mytumour$subtypegroup[ mytumour$tumour_subtype !=  type] <- 1
    
    for(k in myk){
      message(k)
      mytumour$immunegroup <- NA
      mytumour$immunegroup[ mytumour$immune_cluster ==  k] <- 0
      mytumour$immunegroup[ mytumour$immune_cluster !=  k] <- 1    
      mytab <- table(mytumour$subtypegroup, mytumour$immunegroup , 
                     dnn = c("tumour subtype", "immune cluster"))
      print(mytab)
      print(fisher.test(mytab, alternative = "greater"))
    }
  }    
  
  #Freq plot
  tumour_tab$tumour_subtype <- factor(tumour_tab$tumour_subtype, levels = subtype_tab$Var1[order(subtype_tab$Freq, decreasing = T)])
  
  freq_plot <- ggplot(tumour_tab, aes(fill=immune_cluster, y=freq, x= tumour_subtype)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = cluster_col) +
    theme(axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 20, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 20, color = "black")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(size = 20, hjust = 0.5)) +
    theme(legend.position = "none") + labs(y = "Fraction")
  
  return(freq_plot)
  
}



# Volcano plot for DEG
volcano_DEG_plot <- function(genetable, contrast, foldchange, padj){
  res <- genetable[ genetable$contrast == contrast,]      
  res$threshold <- NA
  res$threshold[res$log2FoldChange > foldchange & res$padj < padj] <- "Up-regulated"
  res$threshold[res$log2FoldChange < -(foldchange) & res$padj < padj] <- "Down-regulated"
  res$threshold[is.na(res$threshold)] <- "not significant"    
  res_upreg <- res[ res$threshold == "Up-regulated",]    
  res_upreg <- res_upreg[order(res_upreg$log2FoldChange, decreasing = T),]    
  
  res_downreg <- res[ res$threshold == "Down-regulated",]    
  res_downreg <- res_downreg[order(res_downreg$log2FoldChange, decreasing = F),]  
  
  if(nrow(res_upreg) < 10){
    res$genelabels[res$Gene %in% res_upreg$Gene] <- "UP"
  }
  
  if(nrow(res_downreg) < 10){
    res$genelabels[res$Gene %in% res_downreg$Gene] <- "DOWN"    
  }    
  
  if(nrow(res_upreg) >= 10){
    res$genelabels[res$Gene %in% res_upreg$Gene[1:10]] <- "UP"   
  }
  
  if(nrow(res_downreg) >= 10){
    res$genelabels[res$Gene %in% res_downreg$Gene[1:10]] <- "DOWN"      
  }
  
  myvolcano_plot <- ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(aes(color = threshold), size=2.5) +
    scale_colour_manual(values = c("Down-regulated"= "blue", "Up-regulated"="red",  "not significant"= "black")) +
    geom_text_repel(label = ifelse(res$genelabels == "UP", as.character(res$Gene),""), 
                    size = 8, box.padding = 1, max.overlaps = Inf, direction = "y", nudge_x = 0.9, 
                    vjust = 0.9, min.segment.length = 0) + 
    geom_text_repel(label = ifelse(res$genelabels == "DOWN", as.character(res$Gene),""), 
                    size = 8, box.padding = 1, max.overlaps = Inf, direction = "y", nudge_x = -1.2, 
                    vjust = 0.5, min.segment.length = 0) +     
    myplot + myaxis + 
    theme(axis.text.x = element_text(size = 50, angle = 0, hjust = 0.5),
          axis.title = element_text(size = 50), axis.text.y = element_text(size = 50),
          plot.title = element_text(size = 50, hjust = 0.5),
          legend.position = "none") + 
    labs(x = "Fold change (Log2)" ,y = "p-value (-Log10)") 
  
  return(myvolcano_plot)    
}


#to align plots (from stackoverflow)
# Function to align plots (from stackoverflow) 
align_plots1 <- function (...) {
  pl <- list(...)
  stopifnot(do.call(all, lapply(pl, inherits, "gg")))
  gl <- lapply(pl, ggplotGrob)
  bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
  combined <- Reduce(bind2, gl[-1], gl[[1]])
  wl <- lapply(gl, "[[", "widths")
  combined$widths <- do.call(grid::unit.pmax, wl)
  grid::grid.newpage()
  grid::grid.draw(combined)
}


#circle plots for repertoire analysis

circlepack.reads.fx <- function(inputfile, sample_id){
  #inputfile is: all_clustered_IGH_12rm
  sample_df <- inputfile[inputfile$sample_id == sample_id,]
  
  octamer_tab <- as.data.frame(table(sample_df$octamer), stringsAsFactors = F)
  # Clusters Should have only octamers > 3 sequences 
  cluster_tab <- octamer_tab[octamer_tab$Freq > 2,]
  # if there is a cluster, group cdr3s and use sum of reads
  
  if(nrow(cluster_tab) != 0){    
    cluster_tab$octamerreads <- NA
    for(j in 1:nrow(cluster_tab)){
      myoctamer <- cluster_tab$Var1[j]
      octamerreads <- sum(sample_df$cloneCount[which(sample_df$octamer == myoctamer)])
      cluster_tab$octamerreads[j] <- octamerreads
    }
    #clean up cluster_tab
    cluster_tab <- cluster_tab[, c("Var1", "octamerreads")]
    colnames(cluster_tab) <- c("name", "size")
  }
  
  # Make edge df
  octamer_cdr3 <- sample_df[, c("octamer", "nSeqCDR3")]
  colnames(octamer_cdr3) <- c("from", "to")
  #if no clusters, edge df is octamer_cdr3 df and all octamers will be replaced by sample_id
  if(nrow(cluster_tab) == 0) {
    myedge <- octamer_cdr3
    myedge[,"from"] <- unique(sample_df$sample_id)
  }
  
  # if there is a cluster, make a list of sample_id and clusters
  if(nrow(cluster_tab) != 0){
    myclusters <- cbind.data.frame(NA,cluster_tab$name, stringsAsFactors = F)
    colnames(myclusters) <- c("from", "to")
    myclusters$from <- unique(sample_df$sample_id)
    # replace octamers < 3 sequences with sample_id
    octamer_cdr3$from[!octamer_cdr3$from %in% cluster_tab$name] <- unique(sample_df$sample_id)
    myedge <- rbind(myclusters, octamer_cdr3)   
  }
  
  #vertix df  
  
  # get cdr3 and reads    
  cdr3_freq <- sample_df[, c("nSeqCDR3", "cloneCount")]
  colnames(cdr3_freq) <- c("name", "size")
  # get sample frequency        
  sample_tab <- as.data.frame(table(sample_df$sample_id), stringsAsFactors = F)
  colnames(sample_tab) <- c("name", "size")
  
  #bind all and cleanup
  myvertex <- rbind(sample_tab, cdr3_freq)
  
  # if there is a cluster, include it to vertex
  if(nrow(cluster_tab) != 0){
    myvertex <- rbind(myvertex, cluster_tab)
  }
  
  # first row is NA so remove it
  myvertex <- myvertex[!is.na(myvertex$name),]
  myvertex$size <- as.numeric(myvertex$size)
  
  # Make a type variable for colors
  myvertex$type <- NA
  myvertex$type[myvertex$name %in% cdr3_freq$name] <- "CDR3"
  # If there is cluster, name it as octamer to color
  if(nrow(cluster_tab) != 0){
    myvertex$type[myvertex$name %in% cluster_tab$name] <- "octamer"
  }
  
  #Make graph and plot    
  mygraph <- graph_from_data_frame(myedge, vertices = myvertex)  
  colpal <- c("CDR3" = "light grey","octamer" = "blue")
  alphapal <- c("CDR3" = 1,"octamer" = 0.4)
  
  circlep <- ggraph(mygraph, layout = 'circlepack', weight = size) + 
    geom_node_circle(aes(fill = type, alpha = type)) +
    theme_void() + scale_fill_manual(values = colpal ,na.value="transparent") +
    scale_alpha_manual(values = alphapal)
  
  pdf(file = paste0(plotpath,sample_id,"IGHreads_circles.pdf"),
      width = 15, 
      height = 15,
      useDingbats = FALSE)
  print(circlep)
  dev.off()
}

# Scatter plot with linear regression comparing estimated div from RNA-seq and observed div from CapTCR-seq
TCRcap_rnaplot.fx <- function(dta, met_rna, met_cap, labels, x,y){
  
  lmreg <- lm(log10(eval(as.name(met_rna))) ~ log10(eval(as.name(met_cap))), data = dta)
  print(broom::glance(lmreg))
  reg <- broom::glance(lmreg)
  cap_rna_plot <- ggplot(data = dta, aes(y = eval(as.name(met_rna)), x = eval(as.name(met_cap)))) + 
    geom_point(size = 5, aes(color = group)) +
    geom_smooth(method = "lm", se = FALSE) + myplot +
    scale_y_continuous(trans = "log10") + 
    scale_x_continuous(trans = "log10") +
    annotation_logticks(sides = "bl") +
    theme(legend.position = "none", legend.title = element_blank()) +
    theme(axis.title = element_text(size = 45),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 45, color = "black"),
          axis.text.y = element_text(size = 45, color = "black"),
          plot.title = element_text(hjust = 0.5, size = 45)) +
    labs(y = paste0("Estimated ", labels, "\n(RNAseq)"), x = paste0(labels, "\n(CapTCR-seq)")) +
    annotate("text", x=x, y=y, label= paste("R^2 = ", round(reg$r.squared,2)), size = 15) 
  
  return(cap_rna_plot)
}

