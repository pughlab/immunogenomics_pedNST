###############
# Figure 5
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

load(file = paste0(datapath, "geneset_Tcells_norm.RData"))
load(file = paste0(datapath, "vars_inflamed.RData"))

#order
vars_inflamed <- vars_inflamed[order(vars_inflamed$Tcellgroups, vars_inflamed$cohort),]

myTcluster <- as.character(vars_inflamed$Tcellgroups)
names(myTcluster) <- rownames(vars_inflamed)
class_mat <- t(as.matrix(myTcluster))
rownames(class_mat) <- "Cluster"

mycohort <- vars_inflamed$cohort
names(mycohort) <- vars_inflamed$sample_id

mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"
cohorts_hm <- cohorts_hm.fx(mycohorts)

cells_mat <- geneset_Tcells_norm[,rownames(vars_inflamed)]

#Group signatures
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

Tcellsgroups <- c(rep("Naive", 4), rep("ISG+", 3), rep("NME1+", 3), rep("Tm", 9), rep("Tem", 5),
                  rep("Tfh/h", 5),rep("Tex", 4),rep("Treg", 3), rep("Others", 4))
names(Tcellsgroups) <- Tcellclusters

# get summary
summary(as.vector(cells_mat))

#order groups for heatmap
Tcellsgroups <- factor(Tcellsgroups, levels = c("Naive","NME1+","ISG+", "Tm", "Tem","Treg","Tfh/h","Tex","Others"))

# annotation for TGs
Tha <- HeatmapAnnotation(`T-cell group` = anno_block(labels = c("TG1", "TG2","TG3", "TG4", "TG5"),
                                                     labels_gp = gpar(fontsize = 15),
                                                     show_name = T, height = unit(1,"cm")))

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
                          column_title = "Pediatric Inflamed (n = 90)",
                          row_title_rot = 0,
                          column_split = myTcluster,                       
                          row_split = Tcellsgroups, 
                          cluster_row_slices = FALSE,
                          # annotations
                          top_annotation = Tha,
                          # legends
                          show_heatmap_legend = TRUE,
                          heatmap_legend_param = list(
                                                      at = c(-2,0,2),
                                                      labels = c("<-2", "0", ">2"),
                                                      title = "Cell-type\nz-score")
                          )

pdf(file = paste0(plotpath, "CC_Tcells_inflamed.pdf"),
    width = 12, height = 18)
Tcellscells_hm %v% cohorts_hm
dev.off()


###############
# Figure 6B
###############



###############
# Compile in one file
###############

# setwd(plotpath)
# 
# plotflow:::mergePDF(
#   in.file = list.files(file.path(plotpath), pattern = "Fig5_", full.names = TRUE),
#   file="Figure5.pdf"
# )
# 
# do.call(file.remove, list(list.files(plotpath, pattern = "Fig5_", full.names = TRUE)))
