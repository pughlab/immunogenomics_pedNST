###############
# Figure 3
###############

source("/code/functions/dependencies.R")
source("/code/functions/ggplot2_theme.R")
source("/code/functions/color_schemes.R")
source("/code/functions/plotting_functions.R")

datapath <- "/data/"
plotpath <- "/results/"

###############
# Figure 3A
###############

load(paste0(datapath, "pathwaytable_qusage.RData"))

# get fold change and fdr for pathways
fc_tab <- pathwaytable_qusage[,c("pathway.name", "contrast","log.fold.change")]
fdr_tab <- pathwaytable_qusage[,c("pathway.name", "contrast","FDR")]

fc_mat <- spread(fc_tab, key = contrast, value = log.fold.change)
fdr_mat <- spread(fdr_tab, key = contrast, value = FDR)

rownames(fc_mat) <- fc_mat$pathway.name
fc_mat$pathway.name <- NULL
fc_mat <- as.matrix(fc_mat)

rownames(fdr_mat) <- fdr_mat$pathway.name
fdr_mat$pathway.name <- NULL
fdr_mat <- as.matrix(fdr_mat)

#Remove those pathways with fdr > 0.1 in all clusters
fc_mat_pathways <- fc_mat[rownames(fdr_mat)[rowSums(fdr_mat > 0.1) <= 2],]

#clean up row and colnames
rownames(fc_mat_pathways) <- gsub("HALLMARK_", "", rownames(fc_mat_pathways))
rownames(fc_mat_pathways) <- gsub("_", " ", rownames(fc_mat_pathways))
colnames(fc_mat_pathways) <- gsub("-Others", "", colnames(fc_mat_pathways))
colnames(fc_mat_pathways) <- gsub("[.]", " ", colnames(fc_mat_pathways))

#transpose
fc_mat_pathways_t <- t(fc_mat_pathways)

#Heatmap
ht_opt$HEATMAP_LEGEND_PADDING = unit(1,"cm")

col_fun= colorRamp2(c(-0.6, 0, 0.6), c("blue", "white", "red"))
pathways_hm = Heatmap(fc_mat_pathways_t,
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
                      column_names_gp = gpar(fontsize = 20),
                      row_names_gp = gpar(fontsize = 20),
                      height = unit(8, "cm"), width = unit(44, "cm"),
                      column_title = NULL,
                      row_title = NULL,
                      column_names_rot = 45,
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(col_fun = col_fun,
                                                  at = c(-0.6, 0, 0.6),
                                                  legend.height = unit(4, "cm"),
                                                  title_gp = gpar(fontsize = 18),
                                                  title = "Fold change (Log2)",
                                                  labels_gp = gpar(fontsize = 18)))

pdf(paste0(plotpath,"Fig3_A.pdf"), width = 30, height = 10)
draw(pathways_hm)
decorate_heatmap_body("Log2_FC", {
  grid.rect(x = 0, y = 1, just = c("left", "top"),
            width = 0.45,height = 1, gp = gpar(col = "black", lwd = 5))
})
dev.off()


###############
# Figure 3B
###############

load(file = paste0(datapath,"metadata_IC.RData"))
load(file = paste0(datapath, "proteomics.RData"))
Hs.H <- read.gmt(paste0(datapath, "h.all.v7.1.symbols.gmt"))  

# subset to those with matched RNAseq in pedNST
metadata_IC_proteom <- metadata_IC[metadata_IC$sample_id %in% colnames(proteomics),]

#as factor
metadata_IC_proteom$immune_cluster <- factor(metadata_IC_proteom$immune_cluster,
                                             levels =c("Immune Desert", "Immune Neutral",
                                                       "Myeloid Predominant", "Pediatric Inflamed"))
#asterisks added manually based on the test
coag <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_COAGULATION")
label.df <- data.frame(x = c(1.2,1.2), y = c(2.3,1.3))
coag <- coag + 
  geom_text(data = label.df,aes(x = x, y = y), label = c("*", "***"), size = 15) +
  labs(title = "Coagulation\n(n = 47/138)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1))

emt <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
emt <- emt + geom_text(data = label.df,aes(x = x, y = y), label = c("***", "**"), size = 15) +
  labs(title = "EMT\n(n = 27/200)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1), limits = c(-1.5,1.5))

angio <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_ANGIOGENESIS")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
angio <- angio + geom_text(data = label.df,aes(x = x, y = y), label = c("***", "*"), size = 15) +
  labs(title = "Angiogenesis\n(n = 7/36)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1))

ifng <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_INTERFERON_GAMMA_RESPONSE")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
ifng <- ifng + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(title = "IFN-gamma\n(n = 100/200)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1))

ifna <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_INTERFERON_ALPHA_RESPONSE")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
ifna <- ifna + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(title = "IFN-alpha\n(n = 56/97)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1))

inflam <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_INFLAMMATORY_RESPONSE")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
inflam <- inflam + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(title = "Inflammatory\n(n = 54/200)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1), limits = c(-1.5,1.5))

il6 <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_IL6_JAK_STAT3_SIGNALING")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
il6 <- il6 + geom_text(data = label.df, aes(x = x, y = y), label = c("***", "**"), size = 15) +
  labs(title = "IL6/JAK/STAT3\n(n = 23/87)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1), labels = scales::label_number(accuracy = 1), limits = c(-1.5,1.5))

complement <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_COMPLEMENT")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
complement <- complement + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(title = "Complement\n(n = 123/200)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1), limits = c(-1.5,1.5))

tnfa <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
label.df <- data.frame(x = c(1.2,1.2), y = c(1.3,2.3))
tnfa <- tnfa + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(title = "TNFa\n(n = 51/200)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1), limits = c(-1.5,1.5))

allograft <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_ALLOGRAFT_REJECTION")
label.df <- data.frame(x = c(1,1), y = c(1.3,2.3))
allograft <- allograft + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(title = "Allograft rejection\n(n = 72/200)") + 
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks = c(-1,0,1),labels = scales::label_number(accuracy = 1))

fig3b <- plot_grid(
  ifng, ifna, allograft, il6, inflam, tnfa,complement, emt, angio, coag, 
  rel_widths = c(2.5,1,1,1,1,1,1,1,1,1,1), 
  align = "h", ncol = 11)

y.grob <- textGrob("Average protein z-score", gp=gpar(fontsize=30), rot=90)

pdf(paste0(plotpath,"Fig3_B.pdf"), width = 30, height = 15, useDingbats = FALSE)
fig3b
dev.off()

###############
# Figure 3C
###############

load(file = paste0(datapath, "DESeq2_genetable.RData"))

# remove pvalue = NA
genetable <- genetable[ !is.na(genetable$pvalue),]

C1_volcano <- volcano_DEG_plot(genetable, "C1-Others", 1.5, 0.1) + 
  scale_x_continuous(expand = c(0.4,0)) + 
  scale_y_continuous(expand = c(0.1,0)) +
  ggtitle(expression(~underline("Pediatric Inflamed"))) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype="dashed", color = "black")

C2_volcano <- volcano_DEG_plot(genetable, "C2-Others", 1.5, 0.1) + 
  scale_x_continuous(expand = c(0.5,0)) +
  scale_y_continuous(expand = c(0.2,0)) +
  ggtitle(expression(~underline("Myeloid Predominant"))) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype="dashed", color = "black")

C3_volcano <- volcano_DEG_plot(genetable, "C3-Others", 1.5, 0.1) + 
  scale_x_continuous(expand = c(0.3,0)) +
  scale_y_continuous(expand = c(0.2,0)) +
  ggtitle(expression(~underline("Immune Neutral"))) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype="dashed", color = "black")

C4_volcano <- volcano_DEG_plot(genetable, "C4-Others", 1.5, 0.1) + 
  scale_x_continuous(expand = c(0.3,0)) +
  scale_y_continuous(expand = c(0.1,0)) +
  ggtitle(expression(~underline("Immune Desert"))) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype="dashed", color = "black")

fig3c <- plot_grid(C1_volcano, C2_volcano, C3_volcano, C4_volcano,
                   ncol = 4, nrow =1, align = "h")

pdf(paste0(plotpath,"Fig3_C.pdf"), width = 50, height = 15, useDingbats = FALSE)
fig3c
dev.off()

###############
# Figure 3D
###############

load(file = paste0(datapath, "immunereg_genmat.RData"))
load(file = paste0(datapath,"metadata_IC.RData"))

# median scaled gene expression for each cluster
median_mat <- matrix(nrow = 59, ncol = 4,
                     dimnames = list(rownames(immunereg_genmat), 
                                     c("Pediatric Inflamed", "Myeloid Predominant", 
                                       "Immune Neutral", "Immune Desert")))
#make a matrix
immunereg_genmat <- as.matrix(immunereg_genmat)
#scale
immunereg_genmat_z <- t(scale(t(immunereg_genmat)))

# median of each gene in z-score in each immune cluster
for( g in rownames(immunereg_genmat)){
  median_mat[g,1] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Pediatric Inflamed"]])
  median_mat[g,2] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Myeloid Predominant"]])    
  median_mat[g,3] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Immune Neutral"]])
  median_mat[g,4] <- median(immunereg_genmat_z[g, metadata_IC$sample_id[metadata_IC$immune_cluster == "Immune Desert"]])
}

#Heatmap
ht_opt$HEATMAP_LEGEND_PADDING = unit(1,"cm")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
fig3d = Heatmap(t(median_mat),
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
                column_names_gp = gpar(fontsize = 20),
                column_names_rot = 45,
                row_names_gp = gpar(fontsize = 20),
                width = unit(nrow(median_mat), "cm"),
                height = unit(ncol(median_mat)*2, "cm"),
                column_title = NULL,
                row_title = NULL,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(col_fun = col_fun,
                                            at = c(-1, 0, 1),
                                            legend.height = unit(4, "cm"),
                                            title_gp = gpar(fontsize = 18),
                                            title = "Median z-score",
                                            labels_gp = gpar(fontsize = 18)))

pdf(paste0(plotpath, "Fig3_D.pdf"), width = 30, height = 10)
draw(fig3d)
dev.off()

###############
# Compile in one file
###############

setwd(plotpath)

plotflow:::mergePDF(
  in.file = list.files(file.path(plotpath), pattern = "Fig3_", full.names = TRUE),
  file="Figure3.pdf"
)

do.call(file.remove, list(list.files(plotpath, pattern = "Fig3_", full.names = TRUE)))





