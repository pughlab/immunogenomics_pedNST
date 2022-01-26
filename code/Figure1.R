###############
# Figure 1
###############

library(dplyr)
library(reshape2)
library(ggbeeswarm)
library(gridExtra)
library(gtable)
library(grid)

source("/code/functions/ggplot2_theme.R")
source("/code/functions/color_schemes.R")
source("/code/functions/Heatmap_functions.R")

datapath <- "/data/"

###############
# Figure 1A
###############

load(paste0(datapath,"estimate_ped_pdx.RData"))

ped <- estimate_ped_pdx[ estimate_ped_pdx$group != "TCGA",]

tab <- as.data.frame(table(ped$cohort), stringsAsFactors = F)
tab <- tab[order(tab$Freq, decreasing = T),]
ped$cohort <- factor(ped$cohort, levels = tab$Var1)
ped$group <- factor(ped$group, levels = c("CBTN", "ICGC", "TARGET", "PDX (ITCC)"))

fig1A <- ggplot(data = ped) + geom_bar(aes(y = cohort, fill = group)) + myaxis + myplot +
theme(axis.title = element_blank(), axis.text.x = element_text(size = 25, angle = 0, hjust = 0.5), 
      legend.position = c(0.8,0.9), legend.title = element_blank()) + scale_fill_manual(values = group_col)

pdf(file = "/results/Fig1A.pdf",
    width = 10, height = 10, useDingbats = FALSE)
fig1A
dev.off()

###############
# Figure 1B
###############

load(paste0(datapath,"estimate_ped_pdx.RData"))
#Create two dummy variables to split adult and ped
emptyvar <- as.data.frame(matrix(ncol = 7, nrow = 2))
colnames(emptyvar) <- colnames(estimate_ped_pdx)

emptyvar$group <- as.character(emptyvar$group)
emptyvar$cohort <- as.character(emptyvar$cohort)

emptyvar[1,"aliquot_id"] <- "empty1"
emptyvar[2,"aliquot_id"] <- "empty2"
emptyvar[1,"sample_id"] <- "empty1"
emptyvar[2,"sample_id"] <- "empty2"
emptyvar[1,"ImmuneScore"] <- 3500
emptyvar[2,"ImmuneScore"] <- 3500
emptyvar[1,"cohort"] <- "EMPTY1"
emptyvar[2,"cohort"] <- "EMPTY2"

estimate_ped_pdx <- rbind(estimate_ped_pdx,emptyvar)
#calculate immune reads
estimate_ped_pdx$percread <- 8.0947988*exp(estimate_ped_pdx$ImmuneScore*0.0006267)

immune.cohorts <- cbind(NA, unique(estimate_ped_pdx$cohort))
colnames(immune.cohorts) <- c("group","cohort")
immune.cohorts <- as.data.frame(immune.cohorts)

immune.cohorts$group <- as.character(immune.cohorts$group)

adults <- c("PRAD", "LGG", "OV", "SKCM", "COAD", "GBM", "LUAD")
peds <- c("PDX","ETMR", "MB", "ATRT", "EPN", "pedHGG", "CP", "NBL", "pedLGG", "CPH", "MNG", "SCHW", "NFB")

immune.cohorts[immune.cohorts$cohort %in% adults, 1] <- "Adult"
immune.cohorts[immune.cohorts$cohort %in% peds, 1] <- "Pediatric"

immune.cohorts[immune.cohorts$cohort == "EMPTY1",1] <- "Pediatric"
immune.cohorts[immune.cohorts$cohort == "EMPTY2",1] <- "Pediatric"

for(i in 1:nrow(immune.cohorts)){
  immune.cohorts$median_immunereads[i]<-median(estimate_ped_pdx$percread[estimate_ped_pdx$cohort == immune.cohorts$cohort[i]])}

tmp <- immune.cohorts[which(immune.cohorts$group == "Pediatric"),]
tmp1 <- immune.cohorts[which(immune.cohorts$group == "Adult"),]
immune.cohorts <- rbind(tmp,tmp1)

immune.cohorts$cohort <- factor(immune.cohorts$cohort, levels = c("PDX","ETMR", "MB", "ATRT", "EPN", "pedHGG", "CP", 
                                                                  "NBL", "pedLGG", "CPH",  "MNG", "SCHW", "NFB", 
                                                                  "EMPTY1","EMPTY2", "PRAD", "LGG", "OV", "SKCM", 
                                                                  "COAD", "GBM", "LUAD"))                                          

immune.cohorts <- immune.cohorts[order(immune.cohorts$cohort),]


disease.width <- (nrow(estimate_ped_pdx)/nrow(immune.cohorts)) 
sorted.estimate_ped_pdx <- estimate_ped_pdx[0,]

start = 0
for(i in 1:(nrow(immune.cohorts))){
  tmp <- estimate_ped_pdx[estimate_ped_pdx$cohort == immune.cohorts$cohort[i],]
  tmp <- tmp[order(tmp$percread),] 
  #create range of x values to squeeze dots into equal widths of the plot for each Disease regardless of the number of samples
  div <- disease.width/nrow(tmp)
  #If there is only one sample, put the dot in the middle of the alloted space
  if(dim(tmp)[1]==1)
  {
    tmp$Xpos<-start+(disease.width/2)
  } else tmp$Xpos<-seq(from = start, to = start+disease.width, by = div)[-1]
  
  sorted.estimate_ped_pdx<-rbind(sorted.estimate_ped_pdx,tmp)  
  immune.cohorts$Median.start[i] <- tmp$Xpos[1]
  immune.cohorts$Median.stop[i] <- tmp$Xpos[nrow(tmp)]
  immune.cohorts$N[i]<-nrow(tmp)
  
  start <- start+disease.width+30
  
}

immune.cohorts$medianloc <- immune.cohorts$Median.start+((immune.cohorts$Median.stop-immune.cohorts$Median.start)/2)

sorted.estimate_ped_pdx$cohort <- factor(sorted.estimate_ped_pdx$cohort, levels = levels(immune.cohorts$cohort))     

#Color dummy variables white
rmEMPTY <- rep("black",22)
rmEMPTY[14:15] <- "white"

immune.cohorts$color_crossbar <- NA

immune.cohorts$color_crossbar[immune.cohorts$cohort == "EMPTY1"] <- "white"
immune.cohorts$color_crossbar[immune.cohorts$cohort == "EMPTY2"] <- "white"

immune.cohorts$color_crossbar[is.na(immune.cohorts$color_crossbar)] <- "black"

fig1B <- ggplot() + 
  geom_point(data = sorted.estimate_ped_pdx, aes(x = Xpos ,y = percread, color = cohort), 
             size = 7, shape = 20) +
  geom_crossbar(data = immune.cohorts, 
                aes(x = medianloc,y = median_immunereads, color = color_crossbar,
                    ymin = median_immunereads,ymax = median_immunereads), width = disease.width) +
  myaxis + myplot +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = rmEMPTY)) +
  scale_color_manual(values = c(cohort_col, "white" = "white", "black" = "black"), guide = FALSE) +
  scale_x_continuous(breaks = seq((disease.width)/2,max(sorted.estimate_ped_pdx$Xpos),
                                  disease.width+30), labels = immune.cohorts$cohort, expand = c(0,20)) + 
  scale_y_continuous(breaks = seq(0, 70, by = 10)) + 
  labs(y = "% Immune Reads") 

pdf(file = "/results/Fig1B.pdf",
    width = 20, height = 8, useDingbats = FALSE)
print(fig1B)
dev.off()

###############
# Figure 1C
###############

load(file = paste0(datapath,"metadata_IC.RData"))
load( file = paste0(datapath, "geneset_cc_normalized.RData"))

#order based on immune clusters
cluster_cohort <- metadata_IC[order(metadata_IC$immune_cluster, metadata_IC$cohort),]

#Clusters
mycluster <- as.character(cluster_cohort$immune_cluster)
names(mycluster) <- rownames(cluster_cohort)
cluster_hm <- class_hm.fx(mycluster)

#Cohort
mycohort <- cluster_cohort$cohort
names(mycohort) <- rownames(cluster_cohort)
mycohorts <- t(as.matrix(mycohort))
rownames(mycohorts) <- "Cohort"
cohorts_hm <- cohorts_hm.fx(mycohorts)

cells_mat <- geneset_cc_norm[,rownames(cluster_cohort)]

# Heatmap
cells_hm <- cells_hm.fx(cells_mat)

#annotation
annotation_order <- c("Pediatric inflammed", "Myeloid-driven", "Pediatric cold", "Immune excluded")

cluster_ha = HeatmapAnnotation(clusters = anno_mark(at = c(50, 235, 566, 844), labels_rot = 0,
                                                    labels = annotation_order, side = "top",
                                                    labels_gp = gpar(fontsize = 20), 
                                                    link_height = unit(0.5, "cm")))

pdf("/results/Fig1C.pdf",width = 18, height = 10)
draw(cluster_ha%v% cluster_hm %v% cells_hm %v% cohorts_hm, gap = unit(0.5, "cm"))
dev.off()

#Legends
lgd_cohort = Legend(labels = names(cohort_col)[2:13], title = "Cohort",
                    legend_gp = gpar(fill = cohort_col[2:13]))

pdf(file = "/results/1C_legend_cohort.pdf",
    width = 8,  height = 8, useDingbats = FALSE)
draw(lgd_cohort)
dev.off() 

col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
celltype_lgd = Legend(col_fun = col_fun, 
                      at = c(-3,0,3), labels = c("<3", "0", ">3"), title = "Cell-type\nscore")

pdf(file = "/results/1C_legend_heatmap.pdf",
    width = 8, height = 8, useDingbats = FALSE)
draw(celltype_lgd)
dev.off()  

###############
# Figure 1D
###############

load(file = paste0(datapath,"metadata_IC.RData"))

tab <- as.data.frame(table(metadata_IC$cohort), stringsAsFactors = F)
tab <- tab[order(tab$Freq, decreasing = F),]

cancer_IC_mat <- matrix(nrow = 12, ncol = 4,
                        dimnames = list(tab$Var1, 
                                        c("Pediatric inflamed", "Myeloid-driven", "Pediatric cold", "Immune excluded")))

for(i in 1:nrow(cancer_IC_mat)){
  mycancer <- metadata_IC[ metadata_IC$cohort == rownames(cancer_IC_mat)[i],]    
  freq_tab <- as.data.frame(table(mycancer$immune_cluster), stringsAsFactors = F)
  freq_tab$perc <- freq_tab$Freq/sum(freq_tab$Freq)    
  cancer_IC_mat[i, freq_tab$Var1] <- freq_tab$perc *100
}

cancer_IC_mat[is.na(cancer_IC_mat)] <- 0

col_fun= colorRamp2(c(0, 100), c("white", "red"))
hm1D = Heatmap(cancer_IC_mat,
                    #titles and names   
                    name = "% cancer",   
                    show_row_names = TRUE,
                    show_column_names = TRUE,     
                    #clusters and orders  
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    show_column_dend = TRUE,
                    #aesthestics
                    row_names_side = "left",
                    col = col_fun,
                    column_names_rot = 45,
                    column_names_gp = gpar(fontsize = 20),
                    row_names_gp = gpar(fontsize = 20),
                    height = unit(nrow(cancer_IC_mat), "cm"),
                    width = unit(ncol(cancer_IC_mat), "cm"),
                    column_title_gp = gpar(fontsize = 20),
                    column_title = NULL,
                    row_title = NULL,
                    show_heatmap_legend = FALSE)

ha = rowAnnotation(`cohort size` = anno_barplot(tab$Freq, bar_width = 1,
                                                gp = gpar(col = "white", fill = "#4d4d4d"), 
                                                border = FALSE,
                                                axis_param = list(at = c(0, 100,200,300), labels_rot = 45),
                                                width = unit(4, "cm")), show_annotation_name = FALSE)

pdf("/results/Fig1D.pdf", width = 10, height = 10)
draw(hm1D + ha)
dev.off()

#Legend
col_fun = colorRamp2(c(0, 100), c("white", "red"))
cancer_lgd = Legend(col_fun = col_fun, title = "% cancer")

pdf(file = "/results/1D_legend.pdf",
    width = 8, height = 8,useDingbats = FALSE)
draw(cancer_lgd)
dev.off() 

###############
# Figure 1E
###############

load(file = paste0(datapath,"metadata_IC.RData"))

tab <- as.data.frame(table(metadata_IC$CRI_cluster), stringsAsFactors = F)
tab <- tab[order(tab$Freq, decreasing = F),]

cri_IC_mat <- matrix(nrow = 6, ncol = 4,
                     dimnames = list(tab$Var1,c("Pediatric inflamed", "Myeloid-driven", "Pediatric cold", "Immune excluded")))

for(i in 1:nrow(cri_IC_mat)){
  mycancer <- metadata_IC[ metadata_IC$CRI_cluster == rownames(cri_IC_mat)[i],]    
  freq_tab <- as.data.frame(table(mycancer$immune_cluster), stringsAsFactors = F)
  freq_tab$perc <- freq_tab$Freq/sum(freq_tab$Freq)    
  cri_IC_mat[i, freq_tab$Var1] <- freq_tab$perc *100
}

cri_IC_mat[is.na(cri_IC_mat)] <- 0

col_fun= colorRamp2(c(0, 100), c("white", "red"))
hm1E = Heatmap(cri_IC_mat,
                 #titles and names   
                 name = "% CRI-iAtlas cluster",   
                 show_row_names = TRUE,
                 show_column_names = TRUE,     
                 #clusters and orders  
                 cluster_columns = FALSE,
                 cluster_rows = FALSE,
                 show_column_dend = TRUE,
                 #aesthestics
                 row_names_side = "left",
                 col = col_fun,
                 column_names_rot = 45,                 
                 column_names_gp = gpar(fontsize = 20),
                 row_names_gp = gpar(fontsize = 20),
                 height = unit(nrow(cri_IC_mat), "cm"),
                 width = unit(ncol(cri_IC_mat), "cm"),
                 column_title_gp = gpar(fontsize = 20),
                 column_title = NULL,
                 row_title = NULL,
                 show_heatmap_legend = FALSE)

ha = rowAnnotation(
  `cohort size` = anno_barplot(tab$Freq, bar_width = 1, 
                               gp = gpar(col = "white", fill = "#4d4d4d"), 
                               border = FALSE,
                               axis_param = list(at = c(0, 100,250,500), labels_rot = 45),
                               width = unit(4, "cm")), 
  show_annotation_name = FALSE)

pdf("/results/Fig1E.pdf",width = 10, height = 10)
draw(hm1E + ha)
dev.off()

# Legend
col_fun = colorRamp2(c(0, 100), c("white", "red"))
cri_lgd = Legend(col_fun = col_fun, 
                 title = "% CRI-iAtlas\ncluster")

pdf(file = "/results/1E_legend.pdf",
    width = 8, height = 8, useDingbats = FALSE)
draw(cri_lgd)
dev.off() 


###############
# Figure 1F
###############


###############
# Figure 1G
###############

load(file = paste0(datapath,"metadata_IC.RData"))

lggp <- subgroup_IC.fx(metadata_IC, "pedLGG", "Purples")
hggp <- subgroup_IC.fx(metadata_IC, "pedHGG", "Oranges")
atrtp <- subgroup_IC.fx(metadata_IC, "ATRT", "Greens")
epnp <- subgroup_IC.fx(metadata_IC, "EPN", "Reds")
mbp <- subgroup_IC.fx(metadata_IC, "MB", "Blues")
nblp <- subgroup_IC.fx(metadata_IC, "NBL", "GnBu")

p1 <- plot_grid(atrtp+ ggtitle(expression(~underline("ATRT"))),
                epnp + ggtitle(expression(~underline("EPN"))),
                hggp+ ggtitle(expression(~underline("pedHGG"))),
                nblp + ggtitle(expression(~underline("NBL"))), 
                mbp+ ggtitle(expression(~underline("MB"))), 
                lggp+ ggtitle(expression(~underline("pedLGG"))),
                nrow=1 ,
                align="h")
pdf(paste0(plotpath, "Fig1G.pdf"),
    width = 50, height = 15, useDingbats = FALSE)
p1
dev.off()

