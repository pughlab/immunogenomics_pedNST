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

load(file = paste0(datapath, "gsea_Tcellgroups_norm.RData"))
load(file = paste0(datapath, "metadata_inflamed.RData"))

#order
metadata_inflamed <- metadata_inflamed[order(metadata_inflamed$Tcellgroups, metadata_inflamed$cohort),]

myTcluster <- as.character(metadata_inflamed$Tcellgroups)
names(myTcluster) <- rownames(metadata_inflamed)
class_mat <- t(as.matrix(myTcluster))
rownames(class_mat) <- "Cluster"

mycohort <- metadata_inflamed$cohort
names(mycohort) <- metadata_inflamed$sample_id

mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"
cohorts_hm <- cohorts_hm.fx(mycohorts)

cells_mat <- gsea_Tcellgroups_norm[,rownames(metadata_inflamed)]

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
                                                     labels_gp = gpar(fontsize = 10),
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
                          row_names_gp = gpar(fontsize = 10),
                          row_names_max_width = unit(10,"cm"),
                          column_title_gp = gpar(fontsize = 15),
                          row_title_gp = gpar(fontsize = 15),
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

pdf(file = paste0(plotpath, "Fig6_A.pdf"),
    width = 12, height = 18)
Tcellscells_hm %v% cohorts_hm
dev.off()

###############
# Figure 6B
###############

load(file = paste0(datapath, "metadata_TRB.RData"))
load(file = paste0(datapath, "metadata_inflamed.RData"))

# Div plot
TG_TRB_inflamed <- merge(metadata_TRB, metadata_inflamed[, c("sample_id", "Tcellgroups")], by = "sample_id")

pairwise.t.test(log10(TG_TRB_inflamed$estimated_Shannon), 
                TG_TRB_inflamed$Tcellgroups,
                p.adjust = "bonferroni", pool.sd = F)

TG_div_plot <- ggplot(data = TG_TRB_inflamed,
                      aes(x = Tcellgroups, y = estimated_Shannon)) + 
  geom_beeswarm(cex = 3,aes(color = cohort), size = 5) + 
  geom_boxplot(width = 0.5, outlier.colour = NA, fill = NA) + 
  myaxis + myplot +
  scale_color_manual(values = cohort_col) +
  scale_y_continuous(trans = "log10") + annotation_logticks(sides = "l") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        plot.title = element_text(size = 30, hjust = 0.5)) + 
  geom_signif(y_position = c(4.3, 4.5, 4.7), 
              xmin = c(1,2,2), xmax = c(2,4,5), annotation = c("**", "**", "**"),
              textsize = 15, vjust = 0.5) +
  labs(y = "Estimated\nShannon diversity") + ggtitle(~underline("Pediatric Inflamed (n = 82)"))


pdf(file = paste0(plotpath,"Fig6_B1.pdf"),
    width = 8, 
    height = 10,
    useDingbats = FALSE)
TG_div_plot
dev.off()

# TMB
load(file = paste0(datapath,"pedNST_TMB.RData"))
load(file = paste0(datapath,"metadata_IC.RData"))

ped_tmb_IC <- merge(ped_tmb, metadata_IC[,c("sample_id", "immune_cluster")], by = "sample_id")
TG_TMB_inflamed <- merge(ped_tmb_IC, metadata_inflamed[, c("sample_id", "Tcellgroups")], by = "sample_id")

dim(TG_TMB_inflamed)

TG_tmb_plot <- ggplot(data = TG_TMB_inflamed,
                      aes(x = Tcellgroups, y = mutpermb)) + 
  geom_beeswarm(cex = 1.8,aes(color = cohort), size = 5) + 
  geom_boxplot(width = 0.5, outlier.colour = NA, fill = NA) + 
  myaxis + myplot +
  scale_color_manual(values = cohort_col) +
  scale_y_continuous(trans = "log10") + annotation_logticks(sides = "l") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        plot.title = element_text(size = 30, hjust = 0.5)) + 
  geom_signif(comparisons = list(c("TG1", "TG2")), y_position = 1, map_signif_level=TRUE,
              textsize = 15, test = "wilcox.test", vjust = 0.5) +
  labs(y = "SNV + Indel / Mb") + ggtitle(~underline("Pediatric Inflamed (n = 72)"))

pdf(file = paste0(plotpath,"Fig6_B2.pdf"),
    width = 8, 
    height = 10,
    useDingbats = FALSE)
TG_tmb_plot
dev.off()

###############
# Figure 6C
###############
load(file = paste0(datapath, "metadata_inflamed.RData"))
load(file = paste0(datapath, "geneset_cc_normalized.RData"))

geneset_cc_norm_t <- t(geneset_cc_norm)
metadata_inflamed_genesets <- cbind(metadata_inflamed, geneset_cc_norm_t[ vars_inflamed$sample_id,])

vars_inflamed_genesets$Tcellgroups <- as.factor(vars_inflamed_genesets$Tcellgroups)

# bases plots
dc <- celltype_baseplot.fx(vars_inflamed_genesets, "DC", "Dendritic cells")
tcells <- celltype_baseplot.fx(vars_inflamed_genesets, "T_cells", "T cells")
bcells <- celltype_baseplot.fx(vars_inflamed_genesets, "B_cells", "B cells")
nkcells <- celltype_baseplot.fx(vars_inflamed_genesets, "NK_cells", "NK cells")
granuls <- celltype_baseplot.fx(vars_inflamed_genesets, "Granulocytes", "Granulocytes")
mono <- celltype_baseplot.fx(vars_inflamed_genesets, "Monocytes", "Monocytes")

# for each celltype, do pairwise_comparison and save as DF
DC_sig <- pairwise_comparisons(vars_inflamed_genesets, Tcellgroups, DC, p.adjust.method = "bonferroni") %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1)
DC_sig <- signif_column(DC_sig, p.value)

T_sig <- pairwise_comparisons(vars_inflamed_genesets, Tcellgroups, T_cells, p.adjust.method = "bonferroni") %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1)
T_sig <- signif_column(T_sig, p.value)

B_sig <- pairwise_comparisons(vars_inflamed_genesets, Tcellgroups, B_cells, p.adjust.method = "bonferroni") %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1)
B_sig <- signif_column(B_sig, p.value)

NK_sig <- pairwise_comparisons(vars_inflamed_genesets, Tcellgroups, NK_cells, p.adjust.method = "bonferroni") %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1)
NK_sig <- signif_column(NK_sig, p.value)

gran_sig <- pairwise_comparisons(vars_inflamed_genesets, Tcellgroups, Granulocytes, p.adjust.method = "bonferroni") %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1)
gran_sig <- signif_column(gran_sig, p.value)

mono_sig <- pairwise_comparisons(vars_inflamed_genesets, Tcellgroups, Monocytes, p.adjust.method = "bonferroni") %>%
  dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
  dplyr::arrange(group1)
mono_sig <- signif_column(mono_sig, p.value)

# add to base plots
dc <- dc + geom_signif(comparisons= DC_sig$groups, map_signif_level = TRUE, textsize = 7, vjust = 0.5, 
                       annotations= DC_sig$significance, na.rm = TRUE, step_increase = 0.1, y_position = 1.2, test = NULL, tip_length = 0.01)

nkcells <- nkcells + geom_signif(comparisons= NK_sig$groups, map_signif_level = TRUE, textsize = 7, vjust = 0.5, annotations= NK_sig$significance, 
                                 na.rm = TRUE, y_position =c(4,0,0,6,0,4.5,5,0,5.5,0), test = NULL, tip_length = 0.01) + scale_y_continuous(expand = c(0.1,0.1))

tcells <- tcells + geom_signif(comparisons= T_sig$groups, map_signif_level = TRUE, textsize = 7, vjust = 0.5, annotations= T_sig$significance, 
                               na.rm = TRUE, y_position = 6.6, test = NULL, tip_length = 0.01) + scale_y_continuous(expand = c(0.1,0.1))

bcells <- bcells + geom_signif(comparisons= B_sig$groups, map_signif_level = TRUE, textsize = 7, vjust = 0.5, annotations= B_sig$significance, 
                               na.rm = TRUE, y_position = 4.3, test = NULL, tip_length = 0.01) + scale_y_continuous(expand = c(0.1,0.1))

mono <- mono + geom_signif(comparisons= mono_sig$groups, map_signif_level = TRUE, textsize = 7, vjust = 0.5, annotations= mono_sig$significance, 
                           na.rm = TRUE, y_position = c(0,0,0,7.5,0,6,7,5.5,6.5,0), test = NULL, tip_length = 0.01) + scale_y_continuous(expand = c(0.1,0.1))

granuls <- granuls + geom_signif(comparisons= gran_sig$groups, map_signif_level = TRUE, textsize = 7, vjust = 0.5, annotations= gran_sig$significance, 
                                 na.rm = TRUE, test = NULL, y_position = c(6.1,0,0,0,0,6.6,7,0,0,0), tip_length = 0.01) #+ scale_y_continuous(expand = c(0.1,0.1))


# bind all
all <- cowplot::plot_grid(tcells + theme(axis.text.x = element_blank()), 
                          mono + theme(axis.text.x = element_blank()), 
                          bcells + theme(axis.text.x = element_blank()), 
                          granuls + theme(axis.text.x = element_blank()),
                          nkcells, 
                          dc, 
                          ncol = 2, nrow = 3, align = "vh")

#save with ggarrange
pdf(file = paste0(plotpath,"Fig6_C.pdf"),
    width = 10, height = 12, useDingbats = FALSE)
grid.draw(ggarrange(plots=list(
  tcells + theme(axis.text.x = element_blank())+ ggtitle(~underline("Lymphoid cells")) ,
  mono + theme(axis.text.x = element_blank())+ ggtitle(~underline("Myeloid cells")) ,
  bcells + theme(axis.text.x = element_blank()) ,
  granuls + theme(axis.text.x = element_blank()) ,
  nkcells,
  dc), align = "hv"))
dev.off()

###############
# Figure 6D
###############

load(file = paste0(datapath, "metadata_inflamed.RData"))
load(file = paste0(datapath, "tpm_inflamed_selectgenes.RData"))

# heatmap median of select genes different in TG2 vs TG5
mygenes <- c("HLA-A", "HLA-B", "HLA-C", "FOS","FOSB", "JUN", "JUNB", "ATF3",
             "CXCL8", "CXCL1", "IL6","CD86", "CD69", "TLR2", "TLR4",
             "CD22", "CD19", "BTLA", "FCER2", "FCRL2")
# bind tpms and metadata
metadata_inflamed_genes <- cbind(metadata_inflamed, tpm_inflamed_selectgenes[ metadata_inflamed$sample_id, mygenes])
# order
metadata_inflamed_genes <- metadata_inflamed_genes[order(metadata_inflamed_genes$Tcellgroups, metadata_inflamed_genes$cohort),]

# make a matrix sample x gene
genmat <- matrix( ncol = nrow(vars_inflamed_genes), nrow = length(mygenes))
rownames(genmat) <- mygenes
colnames(genmat) <- rownames(metadata_inflamed_genes)

for(g in 1:nrow(genmat)){
# for each gene, get median for each cancer
  gen <- rownames(genmat)[g]
  mymed <- lapply(unique(vars_inflamed_genes$cohort), function(x){
    median(vars_inflamed_genes[[gen]][ vars_inflamed_genes$cohort == x])})
  names(mymed) <- unique(vars_inflamed_genes$cohort)
  
# center gene tpm for each sample using median of the corresponding cancer
  myscaledgene <- vector(mode = "numeric", length = ncol(genmat))
  names(myscaledgene) <- rownames(vars_inflamed_genes)
  for(i in 1:nrow(vars_inflamed_genes)){
    myscaledgene[i] <- vars_inflamed_genes[[gen]][i]-unlist(mymed[vars_inflamed_genes$cohort[i]])  
  }
  # fill the matrix with centered gene tpm
  genmat[g,] <- myscaledgene[colnames(genmat)]
}

# heatmap median of gene tpms centered with median of each cancer type
medmat <- matrix( ncol = 2, nrow = length(mygenes))
rownames(medmat) <- mygenes
colnames(medmat) <- c("TG2", "TG5")

TG2s <- vars_inflamed_genes$sample_id[ vars_inflamed_genes$Tcellgroups == "TG2"]
TG5s <- vars_inflamed_genes$sample_id[ vars_inflamed_genes$Tcellgroups == "TG5"]

for(g in 1:nrow(medmat)){
  gen <- rownames(medmat)[g]
  medmat[g, "TG2"] <- median(genmat[gen, TG2s], na.rm = T)
  medmat[g, "TG5"] <- median(genmat[gen, TG5s], na.rm = T)   
  }

genegroup <- c(rep("TG5", 15), rep("TG2",5))
names(genegroup) <- rownames(medmat)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
med_hm <- Heatmap(t(medmat),
                  #titles and names   
                  name = "Median\nnormalized\nexpression",   
                  show_row_names = TRUE,
                  show_column_names = TRUE,  
                  col = col_fun,
                  #clusters and orders  
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  show_column_dend = FALSE,
                  #aesthestics
                  column_names_gp = gpar(fontsize = 20),
                  row_names_gp = gpar(fontsize = 20),
                  width = unit(19,"cm"),
                  height = unit(2, "cm"),
                  column_title_gp = gpar(fontsize = 20),
                  row_title_gp = gpar(fontsize = 20),
                  column_title = " ",
                  row_title_rot = 90,
                  column_split = factor(genegroup, levels = c("TG5", "TG2")),  
                  column_names_rot = 45,
                  column_gap = unit(0.5, "cm"),
                  cluster_column_slices = FALSE,
                  # legends
                  show_heatmap_legend = TRUE,
                  heatmap_legend_param = list(col_fun = col_fun, 
                                              at = c(-1.5,0,1.5),
                                              labels = c("<-1.5", "0", ">1.5"),
                                              title = "Median\nnormalized\nexpression")
                  )

pdf(paste0(plotpath, "Fig_6D.pdf"),
    width = 10, height = 10)
draw(med_hm)
dev.off()

###############
# Figure 6E
###############

load(file = paste0(datapath, "vars_myeloid.RData"))
load(file = paste0(datapath, "myeloid_geneset_norm.RData"))

# cohort
vars_myeloid <- vars_myeloid[order(vars_myeloid$Myeloidgroups, vars_myeloid$cohort),]
myMcluster <- as.character(vars_myeloid$Myeloidgroups)
names(myMcluster) <- rownames(vars_myeloid)
class_mat <- t(as.matrix(myMcluster))
rownames(class_mat) <- "Cluster"
mycohort <- vars_myeloid$cohort
names(mycohort) <- rownames(vars_myeloid)

mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"
cohorts_hm <- cohorts_hm.fx(mycohorts)

#Order
myeloid_cells_mat <- myeloid_geneset_norm[,rownames(vars_myeloid)]

#Some reordering for heatmap
myeloid_cells_mat <- myeloid_cells_mat[c('hM01_Mast.TPSAB1',
                                         'hM04_cDC1.BATF3','hM03_cDC2.CD1C','hM02_pDC.LILRA4',
                                         'hM06_Mono.CD16',
                                         'hM09_Macro.PLTP','hM10_Macro.IL1B',
                                         'hM12_TAM.C1QC','hM13_TAM.SPP1'),]

summary(as.vector(myeloid_cells_mat))

#Group signatures
Myeloidgroups <- c("Mast", rep("DC", 3), "Mono", rep("Mac", 2), rep("TAM", 2))
names(Myeloidgroups) <- rownames(myeloid_cells_mat)

Myeloidgroups <- factor(Myeloidgroups, levels = c("Mast","DC","Mono","Mac","TAM"))

# annotation

Mha <- HeatmapAnnotation(`Myeloid group` = anno_block(labels = c("MG1", "MG2","MG3", "MG4", "MG5"),
                                                      labels_gp = gpar(fontsize = 15),
                                                      show_name = T, height = unit(1,"cm")))                                         
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Myeloid_hm <- Heatmap(myeloid_cells_mat,
                      #titles and names   
                      name = "Cell-types z-score",   
                      show_row_names = TRUE,
                      show_column_names = TRUE,  
                      col = col_fun,
                      #clusters and orders  
                      cluster_columns = FALSE,
                      cluster_rows = FALSE,
                      show_column_dend = FALSE,
                      # annotations
                      top_annotation = Mha,
                      #aesthestics
                      column_names_gp = gpar(fontsize = 5),
                      row_names_gp = gpar(fontsize = 10),
                      row_names_max_width = unit(8,"cm"),                      
                      height = unit(9, "cm"),
                      column_title_gp = gpar(fontsize = 10),
                      row_title_gp = gpar(fontsize = 10),
                      row_title_rot = 0,
                      column_split = myMcluster, 
                      column_title = "Myeloid Predominant (n = 279)",                  
                      row_split = Myeloidgroups, 
                      cluster_row_slices = FALSE,
                      #legend
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(col_fun = col_fun, 
                                                  at = c(-2,0,2),
                                                  labels = c("<-2", "0", ">2"),
                                                  title = "Cell-type\nz-score")
)

pdf(file = paste0(plotpath, "Fig6_E.pdf"),
    width = 10, height = 10)
Myeloid_hm %v% cohorts_hm
dev.off()

###############
# Figure 6F
###############

load(file = paste0(datapath, "TME_clusters/vars_myeloid.RData"))
load(file = paste0(datapath, "/TME_clusters/microglia_geneset_norm.RData"))

C2samples_microglia <- metadata_IC$sample_id[metadata_IC$immune_cluster == "Myeloid Predominant" &
                                               metadata_IC$cohort != "NBL"]

tpms_microglia <- tpms_myeloid[,C2samples_microglia]
vars_microglia <- vars_myeloid[C2samples_microglia,]


vars_microglia <- vars_microglia[order(vars_microglia$Myeloidgroups, vars_microglia$cohort),]

myMcluster <- as.character(vars_microglia$Myeloidgroups)
names(myMcluster) <- rownames(vars_microglia)
class_mat <- t(as.matrix(myMcluster))
rownames(class_mat) <- "Cluster"

mycohort <- vars_microglia$cohort
names(mycohort) <- rownames(vars_microglia)

mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"
cohorts_hm <- cohorts_hm.fx(mycohorts)

#Order
microglia_cells_mat <- microglia_geneset_norm[,rownames(vars_microglia)]


Mha <- HeatmapAnnotation(`Myeloid group` = anno_block(labels = c("MG1", "MG2","MG3", "MG4", "MG5"),
                                                      labels_gp = gpar(fontsize = 10),
                                                      show_name = T, height = unit(1,"cm")))

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
Microglia_hm <- Heatmap(microglia_cells_mat,
                        #titles and names   
                        name = "Cell-types z score",   
                        show_row_names = TRUE,
                        show_column_names = TRUE,  
                        col = col_fun,
                        #clusters and orders  
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        show_column_dend = FALSE,
                        #aesthestics
                        column_names_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 20),
                        height = unit(2, "cm"),
                        column_title_gp = gpar(fontsize = 20),
                        row_title_gp = gpar(fontsize = 20),
                        show_heatmap_legend = FALSE,
                        row_title = "Microglia",   
                        row_title_rot = 90,
                        column_split = myMcluster, 
                        column_title = "Myeloid-driven - pedCNS (n = 234)",                  
                        cluster_row_slices = FALSE,
                        top_annotation = Mha)

pdf(file = paste0(plotpath, "CC_microgliacells.pdf"),
    width = 8, height = 10)
Microglia_hm %v% cohorts_hm
dev.off()

summary(as.vector(microglia_cells_mat))





###############
# Figure 6G
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
