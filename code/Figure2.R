###############
# Figure 2
###############

library(ggplot2)
library(survival)
library(survminer)

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
# Figure 2A
###############



###############
# Figure 2B
###############



pdf(file = "/results/Fig1B.pdf",
    width = 20, height = 8, useDingbats = FALSE)
print(fig1B)
dev.off()

###############
# Figure 2C
###############


###############
# Figure 2D
###############

load(file = paste0(datapath,"metadata_IC.RData"))
metadata_IC$vital_status <- as.numeric(as.character(metadata_IC$vital_status))
cluster_names <- c("Pediatric inflamed", "Myeloid driven", "Pediatric cold", "Immune excluded")
sfit <- survfit(Surv(days_to_death, vital_status)~ immune_cluster, data=metadata_IC)
IC_KM_OS <- KM_plot(metadata_IC, sfit, cluster_col, "Immune cluster", "Overall survival", cluster_names)

pdf("/results/Fig2D.pdf",
    width = 10, height = 12, onefile = F)
IC_KM_OS
dev.off()

###############
# Figure 2E
###############

load(file = paste0(datapath,"metadata_IC.RData"))
metadata_IC$recurrence <- as.numeric(as.character(metadata_IC$recurrence))
sfit <- survfit(Surv(days_to_progress, recurrence)~ immune_cluster, data=metadata_IC)
IC_KM_PFS <- KM_plot(metadata_IC, sfit, cluster_col, "Immune clusters", "Progression-free survival", cluster_names)

pdf("/results/Fig2E.pdf",
    width = 10, height = 12, onefile = F)
IC_KM_PFS
dev.off()


