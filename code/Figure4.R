###############
# Figure 4
###############

source("/code/functions/ggplot2_theme.R")
source("/code/functions/color_schemes.R")
source("/code/functions/Heatmap_functions.R")
source("/code/functions/plotting_functions.R")

datapath <- "/data/"

###############
# Figure 4A
###############


###############
# Figure 4B
###############

###############
# Figure 4C
###############

###############
# Figure 4D
###############

###############
# Figure 4E
###############


###############
# Figure 4F
###############


###############
# Figure 4G
###############

load(file = paste0(datapath, "common_snv_fusion.RData"))
laod(file = paste0(datapath, "pedNST_TMB.RData"))

mostfreq <- rowSums(!is.na(oncomat))
topgenesalters <- names(sort(mostfreq, decreasing = T)[1:15])
oncomat <- oncomat[topgenesalters,]

C1 <- IC_oncoprint("Pediatric Inflamed", ped_tmb, oncomat, "TRUE")
C2 <- IC_oncoprint("Myeloid Predominant", ped_tmb, oncomat, "FALSE")
C3 <- IC_oncoprint("Immune Neutral", ped_tmb, oncomat, "FALSE")
C4 <- IC_oncoprint("Immune Excluded", ped_tmb, oncomat, "FALSE")

#Add total percentage mutations as row annotation on the left
percmutations <- rowSums(!is.na(oncomat))/ncol(oncomat)*100

right_ha = rowAnnotation(`% mutated` = anno_barplot(percmutations, border = FALSE, 
                                                    axis_param=list(gp = gpar(fontsize=50), labels_rot = 0)),
                         width = unit(8,"cm"), annotation_name_side = "top", annotation_name_rot = 0,
                         annotation_name_gp = gpar(fontsize = 50))


#legend
lgd = Legend(labels = names(col), title = "Alteration",
             legend_gp = gpar(fill = col))

pdf(file = "/results/Fig4_G.pdf",
    width = 40, 
    height = 17,
    useDingbats = FALSE)
draw(C1 + C2 + C3 + C4 + right_ha, 
     gap = unit(c(0.3,0.3,0.3,1),"cm"),
     annotation_legend_side =  "right", 
     annotation_legend_list = list(lgd))
dev.off()  





###############
# Compile in one file
###############

setwd("/results")

plotflow:::mergePDF(
  in.file = list.files(file.path("/results"), pattern = "Fig3_", full.names = TRUE),
  file="Figure3.pdf"
)

do.call(file.remove, list(list.files("/results/", pattern = "Fig3_", full.names = TRUE)))

