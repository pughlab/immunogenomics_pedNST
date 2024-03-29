###############
# Figure 4
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
# Figure 4A-C
###############

load(file = paste0(datapath, "pedNST_TMB.RData"))

#recode 0 to 0.001 for plotting in log scale
ped_tmb$snvpermb[ped_tmb$snvpermb == 0] <- 0.001
ped_tmb$mutpermb[ped_tmb$mutpermb == 0] <- 0.001
ped_tmb$indelpermb[ped_tmb$indelpermb == 0] <- 0.001

# Fig4A
fig4a <- ggplot(data = ped_tmb, aes(x = immune_cluster, y = snvpermb)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 0.7, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  myplot + myaxis +
  theme(axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30)) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  scale_y_continuous(trans = "log10", breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = scales::label_number(accuracy = 0.01)) + annotation_logticks(sides = "l") +
  labs(y = "SNV / Mb") +
  ggtitle(expression(~underline("pedNST (n = 763)")))

pdf(file = paste0(plotpath,"Fig4_A.pdf"), width = 10, height = 12, useDingbats = FALSE)
fig4a
dev.off()

# Fig4B
fig4b <- ggplot(data = ped_tmb, aes(x = immune_cluster, y = mutpermb)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 0.7, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  myplot + myaxis +
  theme(axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30)) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  scale_y_continuous(trans = "log10", breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = scales::label_number(accuracy = 0.01)) + annotation_logticks(sides = "l") +
  labs(y = "SNV + InDel / Mb") +
  ggtitle(expression(~underline("pedNST (n = 763)")))

pdf(file = paste0(plotpath,"Fig4_B.pdf"), width = 10, height = 12, useDingbats = FALSE)
fig4b
dev.off()

# Fig4C
hgg <- ped_tmb[ped_tmb$cohort == "pedHGG",]
fig4c <- ggplot(data = hgg, aes(x = immune_cluster, y = mutpermb)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 3, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  myplot + myaxis +
  theme(axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30)) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  scale_y_continuous(trans = "log10", breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = scales::label_number(accuracy = 0.01)) + annotation_logticks(sides = "l") +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Desert")), 
              map_signif_level=TRUE, textsize = 20, test = "wilcox.test", vjust = 0.5, y_position = 3) +
  geom_signif(comparisons = list(c("Immune Neutral", "Immune Desert")), 
              map_signif_level=TRUE, textsize = 20, test = "wilcox.test", vjust = 0.5, y_position = 2.6) +
  labs(y = "SNV + InDel / Mb") +
  ggtitle(expression(~underline("pedHGG (n = 63)")))

pdf(file = paste0(plotpath,"Fig4_C.pdf"), width = 10, height = 12, useDingbats = FALSE)
fig4c
dev.off()

###############
# Figure 4D-E
###############

load(file = paste0(datapath, "pedNST_strongpeptides.RData"))

#Fig4D
# column binders shows number of strong binding peptides by Mupexi with the BA rank of <= 0.5%
fig4d <- ggplot(data = metadata_SB, aes(x = immune_cluster, y = binders)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 1, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  myplot + myaxis +
  theme(axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30)) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::label_number(accuracy = 1)) + annotation_logticks(sides = "l") +
  labs(y = paste0("Strong binding peptides")) +
  ggtitle(expression(~underline("pedNST (n = 459)")))

pdf(file = paste0(plotpath,"Fig4_D.pdf"), width = 10,height = 12, useDingbats = FALSE)
fig4d
dev.off()

#Fig4E
load(file = paste0(datapath, "pedNST_strongpeptides.RData"))

hgg <- metadata_SB[ metadata_SB$cohort == "pedHGG",]

# make color transparent if number of data point are <=2
mytab <- table(hgg$immune_cluster) 
cluster_col[ names(cluster_col) %in% names(mytab)[ mytab <= 2] ] <- "transparent"

fig4e <- ggplot(data = hgg, aes(x = immune_cluster, y = binders)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 4, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5, aes(color = immune_cluster)) +
  myplot + myaxis +
  theme(axis.title.y = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30)) +
  theme(plot.title = element_text(size = 30, hjust = 0.5), legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::label_number(accuracy = 0.01)) + annotation_logticks(sides = "l") +
  labs(y = paste0("Strong binding peptides")) +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Neutral")), y_position = 3.1,
              map_signif_level=TRUE, textsize = 20, test = "wilcox.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Desert")), y_position = 3.5,
              map_signif_level=TRUE, textsize = 20, test = "wilcox.test", vjust = 0.5) +
  ggtitle(expression(~underline("pedHGG (n = 35)")))  

pdf(file = paste0(plotpath,"Fig4_E.pdf"), width = 10, height = 12, useDingbats = FALSE)
fig4e
dev.off()

###############
# Figure 4F
###############

load(file = paste0(datapath, "pedNST_TMB.RData"))
load(file = paste0(datapath, "pedNST_mutatedgenes.RData"))
load(file = paste0(datapath, "oncogenic_pathways_genes.RData"))

ped_tmb_IC_pathway <- ped_tmb

mypaths <- c("wnt", "notch", "rtk", "pi3k", "cellcycle", "mmr", "myc", "hippo", "tp53", "nrf2", "tgfb")

# add a column for samples with alterations in each pathway
for( i in mypaths){
  message(i)
  # subset genmat to genes in oncogenic pathway
  genmat_pathway <- genmat[,colnames(genmat)[colnames(genmat) %in% oncopath_list[[i]]]] 
  message("number of altered genes:")
  print(ncol(genmat_pathway))
  print(colnames(genmat_pathway))
  # total number of mutations
  mut_sum <- rowSums(genmat_pathway, na.rm = T)
  message("total number of mutations:")
  print(sum(mut_sum))
  # count number of samples with at least one mutations. recode any non zero value to 1.
  mut_sum[mut_sum > 0] <- 1 
  message("total number of samples with at least one mutation:")
  print(sum(mut_sum))
  # append number of genes in each pathway to the colname
  pathway_name <- paste0(i, "_", ncol(genmat_pathway))
  ped_tmb_IC_pathway[[pathway_name]] <- mut_sum  
}

# make a clusters x pathways matrix to hold fraction of samples with alterations in oncogenic pathways in immune clusters
pathway_mat <- matrix(nrow = 4, ncol = 11)
rownames(pathway_mat) <- c("Pediatric Inflamed", "Myeloid Predominant", "Immune Neutral", "Immune Desert")
colnames(pathway_mat) <- c("wnt_48", "notch_51", "rtk_71", "pi3k_26", "cellcycle_10", "mmr_5",
                           "myc_10", "hippo_32", "tp53_5", "nrf2_3", "tgfb_7")

# for each cluster, add fraction of samples with mutations in each pathway
for(k in 1:nrow(pathway_mat)){
  mycluster <- rownames(pathway_mat)[k]
  tmp <- ped_tmb_IC_pathway[ ped_tmb_IC_pathway$immune_cluster == mycluster,]
  for(p in 1:ncol(pathway_mat)){
    mypathway <- colnames(pathway_mat)[p]
    freqtab <- as.data.frame(table(tmp[[mypathway]]), stringsAsFactors = F)
    freqtab$perc <- freqtab$Freq/sum(freqtab$Freq)   
    pathway_mat[k,p] <- freqtab$perc[2]
  }}    

# for each pathway, add total number of samples with alterations across pedNST
pathwaycount <- vector()
for(i in 1:ncol(pathway_mat)){
  mypathway <- colnames(pathway_mat)[i]
  mypathwaycount <- sum(ped_tmb_IC_pathway[[mypathway]] == 1)
  pathwaycount <- c(pathwaycount, mypathwaycount)
}
names(pathwaycount) <- colnames(pathway_mat)

#sort based on total n altered
pathwaycount <- sort(pathwaycount, decreasing = T)

# Heatmap annotation for total altered samples for oncogenic pathways
ha = HeatmapAnnotation(`n altered` = anno_barplot(pathwaycount,
                                                  border = FALSE, height = unit(3,"cm"),
                                                  axis_param=list(gp = gpar(fontsize=20), labels_rot = 0)),
                       annotation_name_side = "left", annotation_name_rot = 0,
                       annotation_name_gp = gpar(fontsize = 20))

# make a cohort x pathway matrix
cohort_mat <- matrix(nrow = 12, ncol = 11)
rownames(cohort_mat) <- unique(ped_tmb_IC_pathway$cohort)
colnames(cohort_mat) <- c("wnt_48", "notch_51", "rtk_71", "pi3k_26", "cellcycle_10", "mmr_5",
                          "myc_10", "hippo_32", "tp53_5", "nrf2_3", "tgfb_7")

# for each pathway, add number of samples with alterations in each pathway
for(i in 1:ncol(cohort_mat)){
  mypathway <- colnames(cohort_mat)[i]
  tmp <- ped_tmb_IC_pathway[ped_tmb_IC_pathway[[mypathway]] == 1,]
  mytab <- as.data.frame(table(tmp$cohort), stringsAsFactors = F)
  cohort_mat[ match(mytab$Var1, rownames(cohort_mat)),i] <- mytab$Freq
}
# If NA convert to 0
cohort_mat[is.na(cohort_mat)] <- 0

#convert to fractions of altered samples in cancer types
cohort_mat_frac <- t(apply(cohort_mat,2, function(x){x/sum(x)}))

#color pallette same as before
mycol <- c("ETMR" = "#76afa9",
           "MNG" = "#a6cee3", 
           "MB" =  "#1f78b4",
           "SCHW" =  "#678ba5",
           "ATRT" =  "#33a02c",
           "EPN" = "#e31a1c",
           "pedHGG" = "#fdbf6f",
           "pedLGG" = "#8763ae",
           "NFB" = "#cab2d6",
           "CPH" = "#cccc7a",
           "CP" = "#a9a9a9",
           "NBL" = "#b2df8a")

# order based on total n altered
cohort_mat_frac <- cohort_mat_frac[names(pathwaycount),]

# Heatmap annotation for fraction of altered samples in cancer types
ha1 = HeatmapAnnotation(`Tumour type` = anno_barplot(cohort_mat_frac, 
                                                     gp = gpar(fill = mycol[ colnames(cohort_mat_frac)]),
                                                     border = FALSE, height = unit(3,"cm"),
                                                     axis_param = list(gp = gpar(fontsize=20), labels_rot = 0)),
                        annotation_name_side = "left", annotation_name_rot = 0,
                        annotation_name_gp = gpar(fontsize = 20))

# cleanup pathway labels
col_labels = structure(c("RTK (71)", "PI3K (26)", "Wnt (48)", "Notch (51)", "HIPPO (32)",
                         "TP53 (5)", "Cell cycle (10)", "MYC (10)", "MMR (5)",
                         "TGFb (7)", "NRF2 (3)"),
                       names = names(pathwaycount))

#order with cohort_mat
pathway_mat <- pathway_mat[,rownames(cohort_mat_frac)]

# heatmap for fraction of mutated samples in each immune cluster (z score)
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
my_hm = Heatmap(scale(pathway_mat), #scale
                #titles and names
                name = "% mutated samples (z-score)",
                col = col_fun, 
                show_row_names = TRUE,
                show_column_names = TRUE,    
                #clusters
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                #aesthestics
                column_title_gp = gpar(fontsize = 20),
                column_names_gp = gpar(fontsize = 20),
                column_names_rot = 45,
                height = unit(4, "cm"),
                width = unit(11,"cm"),
                row_names_gp = gpar(fontsize = 20),
                row_names_side = "left", 
                column_labels = col_labels,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(col_fun = col_fun, 
                                            title = "Mutated samples\n(z-score)", 
                                            labels_gp = gpar(fontsize = 15), 
                                            title_gp = gpar(fontsize = 15)))

lgd_cohort = Legend(labels = names(mycol), title = "", nrow = 2, legend_gp = gpar(fill = mycol))

pdf(file = paste0(plotpath,"Fig4_F.pdf"), width = 12, height = 10, useDingbats = FALSE)
draw(ha %v% ha1 %v% my_hm, annotation_legend_list = lgd_cohort, annotation_legend_side = "bottom")
dev.off() 

###############
# Figure 4G
###############

load(file = paste0(datapath, "common_snv_fusion.RData"))
load(file = paste0(datapath, "pedNST_TMB.RData"))

mostfreq <- rowSums(!is.na(oncomat))
topgenesalters <- names(sort(mostfreq, decreasing = T)[1:15])
oncomat <- oncomat[topgenesalters,]

C1 <- IC_oncoprint("Pediatric Inflamed", ped_tmb, oncomat, "TRUE")
C2 <- IC_oncoprint("Myeloid Predominant", ped_tmb, oncomat, "FALSE")
C3 <- IC_oncoprint("Immune Neutral", ped_tmb, oncomat, "FALSE")
C4 <- IC_oncoprint("Immune Desert", ped_tmb, oncomat, "FALSE")

#Add total percentage mutations as row annotation on the left
percmutations <- rowSums(!is.na(oncomat))/ncol(oncomat)*100

right_ha = rowAnnotation(`% mutated` = anno_barplot(percmutations, border = FALSE, 
                                                    axis_param=list(gp = gpar(fontsize=35), labels_rot = 0)),
                         width = unit(8,"cm"), annotation_name_side = "top", annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 35))

#legend
lgd = Legend(labels = names(col), title = "", 
             grid_height = unit(1, "cm"), grid_width = unit(1, "cm"),
             legend_gp = gpar(fill = col), labels_gp = gpar(fontsize = 30))

pdf(file = paste0(plotpath,"Fig4_G.pdf"), width = 40,  height = 17, useDingbats = FALSE)
draw(C1 + C2 + C3 + C4 + right_ha, gap = unit(c(0.3,0.3,0.3,1),"cm"), annotation_legend_side =  "right", annotation_legend_list = list(lgd))
dev.off()  

###############
# Compile in one file
###############

setwd(plotpath)

plotflow:::mergePDF(
  in.file = list.files(file.path(plotpath), pattern = "Fig4_", full.names = TRUE),
  file="Figure4.pdf"
)

do.call(file.remove, list(list.files(plotpath, pattern = "Fig4_", full.names = TRUE)))

