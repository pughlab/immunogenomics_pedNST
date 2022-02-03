###############
# Figure 3
###############

library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(tidyr)

# library(dplyr)
# library(reshape2)
# library(ggridges)
# library(ggbeeswarm)
# library(ggsignif)
# library(cowplot)
# library(gridExtra)
# library(gtable)
# library(grid)
# library(ComplexHeatmap)
# library(circlize)

source("/code/functions/ggplot2_theme.R")
source("/code/functions/color_schemes.R")
source("/code/functions/Heatmap_functions.R")
source("/code/functions/plotting_functions.R")

datapath <- "/data/"

###############
# Figure 3A
###############

load(file = paste0(datapath,"pathwaytable_qusage.RData"))

fc_tab <- pathwaytable_qusage[,c("pathway.name", "contrast","log.fold.change")]
fc_mat <- spread(fc_tab, key = contrast, value = log.fold.change)
rownames(fc_mat) <- fc_mat$pathway.name
fc_mat$pathway.name <- NULL
fc_mat <- as.matrix(fc_mat)
fc_mat_pathways <- fc_mat[rownames(fdr_mat)[rowSums(fdr_mat > 0.1) <= 2],]
#clean up row and colnames
rownames(fc_mat_pathways) <- gsub("HALLMARK_", "", rownames(fc_mat_pathways))
rownames(fc_mat_pathways) <- gsub("_", " ", rownames(fc_mat_pathways))
colnames(fc_mat_pathways) <- gsub("-Others", "", colnames(fc_mat_pathways))
colnames(fc_mat_pathways) <- gsub("[.]", " ", colnames(fc_mat_pathways))

col_fun= colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
pathways_hm = Heatmap(fc_mat_pathways,
                      #titles and names   
                      name = "Log2_FC",   
                      show_row_names = TRUE,
                      show_column_names = TRUE,     
                      #clusters and orders  
                      cluster_columns = TRUE,
                      cluster_rows = TRUE,
                      show_column_dend = TRUE,
                      #aesthestics
                      col = col_fun,
                      column_names_gp = gpar(fontsize = 10),
                      row_names_gp = gpar(fontsize = 10),
                      height = unit(29, "cm"), width = unit(4, "cm"),
                      column_title_gp = gpar(fontsize = 10),
                      column_title = NULL,
                      row_title = NULL,
                      column_names_rot = 45,
                      show_heatmap_legend = FALSE
)

pdf("/results/Fig3A.pdf",
    width = 10, height = 45)
draw(pathways_hm)
decorate_heatmap_body("Log2_FC", {
  grid.rect(x = 0, y = 1, just = c("left", "top"),
            width = 1,height = 0.345, gp = gpar(col = "black", lwd = 5))
})
dev.off()

col_fun= colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
qusage_lgd = Legend(col_fun = col_fun, title = "Fold change\n(Log2)",
                    at = c(-0.6, 0, 0.6), legend_height = unit(2, "cm"))
pdf(file = "/results/Fig3A_legend.pdf",
    width = 8, 
    height = 8,
    useDingbats = FALSE)
draw(qusage_lgd)
dev.off() 

###############
# Figure 3B
###############
