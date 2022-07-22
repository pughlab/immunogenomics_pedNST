# heatmap for groups, TARGET, ICGC, CBTN

group_hm.fx <- function(group_vector_ordered){
  
  group_mat <- t(as.matrix(group_vector_ordered))
  rownames(group_mat) <- "Group"
  group_hm = Heatmap(group_mat,
                     #titles and names 
                     name = "Group",
                     show_row_names = TRUE,
                     show_column_names = FALSE,    
                     #clusters
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     #aesthestics
                     col = group_col,
                     column_names_gp = gpar(fontsize = 15),
                     height = unit(1, "cm"),
                     row_names_gp = gpar(fontsize = 20),
                     row_title = NULL)
  return(group_hm)
}

# CRI iATLAS clusters

cri_hm.fx <- function(class_vector_ordered){
  
  class_mat <- t(as.matrix(class_vector_ordered))
  rownames(class_mat) <- "CRI-iAtlas cluster"
  class_hm = Heatmap(class_mat,
                     #titles and names
                     name = "CRI-iAtlas cluster",
                     show_row_names = TRUE,
                     show_column_names = TRUE,    
                     #clusters
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     #aesthestics
                     col = cri_col,
                     column_names_gp = gpar(fontsize = 15),
                     height = unit(1, "cm"),
                     row_names_gp = gpar(fontsize = 20))
  return(class_hm)
}

# immune clusters
class_hm.fx <- function(class_vector_ordered){
  
  class_mat <- t(as.matrix(class_vector_ordered))
  rownames(class_mat) <- "Cluster"
  class_hm = Heatmap(class_mat,
                     #titles and names
                     name = "Immune cluster",
                     show_row_names = TRUE,
                     show_column_names = TRUE,    
                     #clusters
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     #aesthestics
                     col = cluster_col,
                     column_names_gp = gpar(fontsize = 15),
                     height = unit(1, "cm"),
                     row_names_gp = gpar(fontsize = 20))
  return(class_hm)
}

# Cohort
cohorts_hm.fx <- function(cohorts_mat){
  
  cohort_hm = Heatmap(mycohorts,
                      #titles and names
                      name = "Cohort",
                      show_row_names = TRUE,
                      show_column_names = FALSE,    
                      #clusters
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      #aesthestics
                      col = cohort_col,
                      column_names_gp = gpar(fontsize = 15),
                      height = unit(1, "cm"),
                      row_names_gp = gpar(fontsize = 20))
  return(cohort_hm)   
}


# genesets for immune celltypes
cells_hm.fx <- function(cells_mat){
  col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
  cells_hm = Heatmap(cells_mat,
                     #titles and names   
                     name = "Cell-types z score",   
                     show_row_names = TRUE,
                     show_column_names = FALSE,  
                     col = col_fun,
                     #clusters and orders  
                     cluster_columns = TRUE,
                     cluster_rows = TRUE,
                     show_column_dend = TRUE,
                     #aesthestics
                     column_names_gp = gpar(fontsize = 25),
                     row_names_gp = gpar(fontsize = 25),
                     height = unit(nrow(cells_mat), "cm"),
                     column_title_gp = gpar(fontsize = 25),
                     column_title = NULL,
                     row_title = NULL,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param = list(col_fun = col_fun, 
                                                 at = c(-3,0,3),
                                                 labels = c("<-3", "0", ">3"),
                                                 title = "Cell-type\nscore")
  )
  return(cells_hm)   
}

