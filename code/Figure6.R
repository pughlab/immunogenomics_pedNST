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
for( g in rownames(immunereg_genmat)){
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
                show_heatmap_legend = FALSE,
                row_names_gp = gpar(fontsize = 25),
                width = unit(nrow(median_mat), "cm"),
                height = unit(ncol(median_mat)*2, "cm"),
                column_title_gp = gpar(fontsize = 25),
                column_title = NULL,
                row_title = NULL)


col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, 
             at = c(-1,0,1), 
             title = "Median\nz-score", title_gp = gpar(fontsize = 15),
             grid_height = unit(2, "cm"),
             labels_gp = gpar(fontsize = 15))

pdf(paste0(plotpath, "Fig6_A.pdf"),
    width = 40, height = 20)
draw(fig6a,annotation_legend_side =  "bottom", 
     annotation_legend_list = list(lgd))
dev.off()

###############
# Figure 6B
###############

load(file = paste0(datapath, "metadata_inflamed.RData"))
load(file = paste0(datapath, "gsea_Tcellgroups_norm.RData"))

# list the included T-cell signatures 
Tcellclusters <- c('CD4.c01(Tn)','CD4.c03(ADSL+ Tn)','CD4.c04(IL7R- Tn)', 'CD8.c01(Tn)', #Naive
                   'CD4.c21(ISG+ Treg)','CD4.c22(ISG+ Th)','CD8.c15(ISG+ CD8+ T)', #ISG
                   'CD4.c23(NME1+CCR4- T)','CD4.c24(NME1+CCR4+ T)','CD8.c17(NME1+ T)',#NME1+
                   'CD4.c07(TIMP1+ Tm)','CD4.c08(CREM+ Tm)','CD4.c09(CCL5+ Tm)','CD4.c10(CAPG+ Tm)', #Tm
                   'CD4.c11(CAPG+CREM- Tm)','CD4.c06(AREG+ Tm)', 'CD8.c02(IL7R+ Tm)','CD8.c04(ZNF683+CXCR6- Tm)','CD8.c10(ZNF683+CXCR6+ Trm)',#Tm
                   'CD4.c12(GZMK+ Tem)', 'CD4.c13(Temra)', #Tem
                   'CD8.c05(GZMK+ early Tem)','CD8.c06(GZMK+ Tem)','CD8.c07(Temra)',#Tem
                   'CD4.c02(CXCR5+ pre-Tfh)','CD4.c16(IL21+ Tfh)','CD4.c17(IFNG+ Tfh/Th1)', 'CD4.c14(CCR6+ Th17)','CD4.c15(IL26+ Th17)', #Th
                   'CD8.c13(OXPHOS- Tex)','CD8.c14(TCF7+ Tex)','CD8.c11(GZMK+ Tex)','CD8.c12(terminal Tex)',#Tex
                   'CD4.c18(TNFRSF9- Treg)','CD4.c19(S1PR1+ Treg)','CD4.c20(TNFRSF9+ Treg)',#Tregs
                   'CD8.c08(KIR+EOMES+ NK-like)','CD8.c09(KIR+TXK+ NK-like)','CD8.c16(Tc17)',
                   'CD4.c05(TNF+ T)')

# order by t-cell groups and cohort
metadata_inflamed <- metadata_inflamed[order(metadata_inflamed$Tcellgroups,
                                             metadata_inflamed$cohort),]

# TG cluster membership for inflamed samples
myTcluster <- as.character(metadata_inflamed$Tcellgroups)
names(myTcluster) <- rownames(metadata_inflamed)
class_mat <- t(as.matrix(myTcluster))
rownames(class_mat) <- "Cluster"

# cohort track
mycohort <- metadata_inflamed$cohort
names(mycohort) <- metadata_inflamed$sample_id
mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"
cohorts_hm <- cohorts_hm.fx(mycohorts)

cells_mat <- gsea_Tcellgroups_norm[, rownames(metadata_inflamed)]

#Group signatures
Tcellsgroups <- c(rep("Naive", 4), rep("ISG+", 3), rep("NME1+", 3), rep("Tm", 9), rep("Tem", 5),
                  rep("Tfh/h", 5),rep("Tex", 4),rep("Treg", 3), rep("Others", 4))
names(Tcellsgroups) <- Tcellclusters

#order groups for heatmap
Tcellsgroups <- factor(Tcellsgroups, levels = c("Naive","NME1+","ISG+", "Tm", "Tem","Treg","Tfh/h","Tex","Others"))

# T-cell group heatmap annotation
Tha <- HeatmapAnnotation(`T-cell group` = anno_block(labels = c("TG1", "TG2","TG3", "TG4", "TG5"),
                                                     labels_gp = gpar(fontsize = 15),
                                                     show_name = FALSE, height = unit(1,"cm")))

#Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Tcellscells_hm <- Heatmap(cells_mat,
                          #titles and names   
                          name = "Cell-types z-score",   
                          show_row_names = TRUE,
                          show_column_names = TRUE,  
                          col = col_fun,
                          #clusters and orders  
                          cluster_columns = FALSE,
                          cluster_rows = FALSE,
                          show_column_dend = FALSE,
                          #aesthestics
                          row_names_gp = gpar(fontsize = 20),
                          row_names_max_width = unit(10,"cm"),
                          column_title_gp = gpar(fontsize = 25),
                          row_title_gp = gpar(fontsize = 25),
                          column_title = "Pediatric inflamed (n = 90)",
                          show_heatmap_legend = TRUE,
                          row_title_rot = 90,
                          column_split = myTcluster,                       
                          row_split = Tcellsgroups, 
                          cluster_row_slices = FALSE,
                          top_annotation = Tha,
                          heatmap_legend_param = list(col_fun = col_fun, 
                                                      at = c(-2, 0, 2),
                                                      labels = c("<-2", "0", ">2"),
                                                      title = "Cell-type\nz-score"))

pdf(file = paste0(plotpath, "Fig6_B.pdf"),
    width = 12, height = 18)
Tcellscells_hm %v% cohorts_hm
dev.off()


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
