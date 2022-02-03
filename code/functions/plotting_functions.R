#Ridge plot
plot_ridge.fx <- function(mydf, var){  
  
  ridgeline.plot <- ggplot() +
    geom_density_ridges(data = mydf, aes(x = eval(parse(text = var)), y = immune_cluster, fill = immune_cluster)) + 
    myplot +
    theme(axis.line = element_line(color = "black"),
          axis.text = element_text(size = 30, color = "black"),
          axis.title = element_blank(), 
          axis.text.x = element_blank(),
          plot.title = element_text(size = 30, hjust = 0.5),
          legend.position = "none") + coord_flip() + #flip axes
    
    scale_fill_manual(values = cluster_col) 
  
  return(ridgeline.plot) 
}

# ridgeplot for proteomics Fig3B
hallmark_IC_stats_ridge <- function(proteinmat, metadata, hallmark){
  
  Hs.H <- read.gmt("/data/h.all.v7.1.symbols.gmt")
  protein_mat_pathway <- as.matrix(proteinmat[rownames(proteinmat) %in% Hs.H[[hallmark]], rownames(metadata)])
  print(dim(protein_mat_pathway))
  pathway_mean <- colMeans(protein_mat_pathway)
  
  metadata$pathway_mean <- pathway_mean
  
  print(pairwise.t.test(metadata$pathway_mean, metadata$immune_cluster, p.adjust = "bonferroni"))
  
  ridge_plot <- plot_ridge.fx(metadata, "pathway_mean")
  return(ridge_plot)
}


# Stacked barplots for cancer subtypes
subgroup_IC.fx <- function(metadata, tumour, color){
  mytumour <- metadata[metadata$cohort == tumour,]
  tumour_tab <- as.data.frame(table(mytumour$tumour_subtype, mytumour$immune_cluster),
                              stringsAsFactors = T)

  if(length(which(tumour_tab$Var2 == "Pediatric inflamed")) == 0){
    tumour_tab <- rbind(tumour_tab, NA)
    tumour_tab[ is.na(tumour_tab$Var1),1:3] <- list(tumour_tab$Var1[1], "C1", 0)
  }
  
  
  myplot <- ggplot(tumour_tab, aes(fill=Var1, y=Freq, x=Var2)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_brewer(palette = color) + 
    theme(axis.title.y = element_text(size = 65),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(size = 50, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 50, color = "black")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(size = 50, hjust = 0.5)) +
    theme(legend.position = "bottom", legend.direction="vertical",
          legend.text = element_text(size = 50),
          legend.key.height= unit(2, 'cm'),
          legend.key.width= unit(2, 'cm'),
          legend.title = element_blank()) + labs(y = "counts")
  return(myplot)
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


#KM plot for survival
KM_plot <- function(metadata, fit_model, col, title, ylab, legendlabs){
  kmplot <- ggsurvplot(fit_model, data = metadata,
                       palette = as.vector(col),conf.int=FALSE, 
                       xlim = c(0,5000),break.x.by = 1000,
                       #Risk table
                       risk.table = TRUE,
                       # pvalue
                       pval = TRUE, pval.size = 10, pval.coord = c(200, 0.1),
                       # legend
                       legend.title="", font.legend = 25, legend.labs = legendlabs, legend = c(0.75, 0.9),
                       # fonts
                       font.main = 30, font.x = 30, font.y = 30, font.tickslab = 30, 
                       # titles
                       title = title, xlab = "Time (days)", ylab = ylab)
  
  kmplot$table <- ggrisktable(fit_model, data = metadata, 
                              color = "strata", palette = as.vector(col),
                              fontsize = 10, risk.table.title = "",
                              xlim = c(0,5000),break.time.by = 1000,
                              y.text = TRUE, ylab = "",  xlab = "",legend.labs = legendlabs,
                              tables.theme = theme_cleantable(), font.tickslab = 20)
  
  return(kmplot)
}
