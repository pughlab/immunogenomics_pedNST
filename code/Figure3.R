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
fdr_tab <- pathwaytable_qusage[,c("pathway.name", "contrast","FDR")]

fc_mat <- spread(fc_tab, key = contrast, value = log.fold.change)
fdr_mat <- spread(fdr_tab, key = contrast, value = FDR)

rownames(fc_mat) <- fc_mat$pathway.name
fc_mat$pathway.name <- NULL
fc_mat <- as.matrix(fc_mat)

rownames(fdr_mat) <- fdr_mat$pathway.name
fdr_mat$pathway.name <- NULL
fdr_mat <- as.matrix(fdr_mat)

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

load(file = paste0(datapath,"metadata_IC.RData"))
load(file = paste0(datapath, "proteomics.RData"))

metadata_IC_proteom <- metadata_IC[metadata_IC$sample_id %in% colnames(proteomics),]

coag <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_COAGULATION")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
coag <- coag + 
  geom_text(data = label.df,aes(x = x, y = y), label = c("*", "***"), size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm"), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))

emt <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
emt <- emt + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))

angio <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_ANGIOGENESIS")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
angio <- angio + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))

ifng <- hallmark_IC_stats_ridge(proteomics, metadata_IC_proteom, "HALLMARK_INTERFERON_GAMMA_RESPONSE")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
ifng <- ifng + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))

ifna <- hallmark_IC_stats_ridge(proteomics,metadata_IC_proteom, "HALLMARK_INTERFERON_ALPHA_RESPONSE")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
ifna <- ifna + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))

inflam <- hallmark_IC_stats_ridge(proteomics,metadata_IC_proteom, "HALLMARK_INFLAMMATORY_RESPONSE")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
inflam <- inflam + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-1.1,1.3))

il6 <- hallmark_IC_stats_ridge(proteomics,metadata_IC_proteom, "HALLMARK_IL6_JAK_STAT3_SIGNALING")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
il6 <- il6 + geom_text(data = label.df, aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-1.1,1.3))

complement <- hallmark_IC_stats_ridge(proteomics,metadata_IC_proteom, "HALLMARK_COMPLEMENT")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
complement <- complement + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-1.1,1.3))

tnfa <- hallmark_IC_stats_ridge(proteomics,metadata_IC_proteom, "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
tnfa <- tnfa + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-1.1,1.3))

allograft <- hallmark_IC_stats_ridge(proteomics,metadata_IC_proteom, "HALLMARK_ALLOGRAFT_REJECTION")
label.df <- data.frame(x = c(1.2,1.2), y = c(3.3,4.3))
allograft <- allograft + geom_text(data = label.df,aes(x = x, y = y), label = "***", size = 15) +
  labs(x = "Average protein z score") +
  theme(plot.margin = unit(c(1,0,0,0), "cm")) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1))

proteom_title <- ggdraw() + draw_label(~underline("Pathway protein levels (CBTN n = 141)"), size = 30)

myplot <- plot_grid(proteom_title,
                    ifng, ifna, allograft, il6, inflam, tnfa,complement, emt, angio, coag, 
                    rel_heights = c(1/4,1,1,1,1,1,1,1,1,1,2), 
                    align = "v", nrow = 11)
y.grob <- textGrob("Average protein z-score", gp=gpar(fontsize=30), rot=90)

pdf("/results/Fig3B.pdf",
    width = 15, height = 30, useDingbats = FALSE)
grid.arrange(arrangeGrob(myplot, left = y.grob))
dev.off()






