###############
# Figure 6
###############

source("/code/functions/dependencies.R")
source("/code/functions/ggplot2_theme.R")
source("/code/functions/color_schemes.R")
source("/code/functions/Heatmap_functions.R")
source("/code/functions/Oncoprint_functions.R")
source("/code/functions/plotting_functions.R")

datapath <- "/data/"
plotpath <- "/results/"

###############
# Figure 6A
###############

load(file = paste0(datapath, "immunereg_genmat.RData"))
load(file = paste0(datapath,"metadata_IC.RData"))

# median scaled gene expression for each cluster
median_mat <- matrix(nrow = 59, ncol = 4,
                     dimnames = list(rownames(immunereg_genmat), 
                                     c("Pediatric Inflamed", "Myeloid Predominant", 
                                        "Immune Neutral", "Immune Excluded")))
#make a matrix
immunereg_genmat <- as.matrix(immunereg_genmat)
#scale
immunereg_genmat_z <- t(scale(t(immunereg_genmat)))
# median of each gene in z-score in each immune cluster
for( g in rownames(mygen_mat)){
  median_mat[g,1] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Pediatric Inflamed"]])
  median_mat[g,2] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Myeloid Predominant"]])    
  median_mat[g,3] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Immune Neutral"]])
  median_mat[g,4] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Immune Excluded"]])
}

#Heatmap
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
fig6a = Heatmap(t(median_mat),
                #titles and names   
                name = "Median z score",   
                show_row_names = TRUE,
                show_column_names = TRUE,     
                #clusters and orders  
                cluster_columns = TRUE,
                cluster_rows = TRUE,
                show_column_dend = TRUE,
                show_row_dend = TRUE,    
                #aesthestics
                col = col_fun,
                column_names_gp = gpar(fontsize = 25),
                column_names_rot = 45,
                show_heatmap_legend = TRUE,
                row_names_gp = gpar(fontsize = 25),
                width = unit(nrow(median_mat), "cm"),
                height = unit(ncol(median_mat)*2, "cm"),
                column_title_gp = gpar(fontsize = 25),
                column_title = NULL,
                row_title = NULL,
                heatmap_legend_param = list(col_fun = col_fun, 
                                            at = c(-1, 0, 1), 
                                            title = "Median\nz-score"))

pdf(paste0(plotpath, "Fig6_A.pdf"),
    width = 40, height = 20)
draw(fig6a)
dev.off()


###############
# Figure 6B
###############

###############
# Figure 6C
###############

###############
# Figure 6D
###############

###############
# Figure 6E
###############


###############
# Compile in one file
###############

setwd(plotpath)

plotflow:::mergePDF(
  in.file = list.files(file.path(plotpath), pattern = "Fig6_", full.names = TRUE),
  file="Figure6.pdf"
)
do.call(file.remove, list(list.files(plotpath, pattern = "Fig6_", full.names = TRUE)))
