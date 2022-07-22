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


###############
# Figure 5A
###############




###############
# Compile in one file
###############

setwd("/results")

plotflow:::mergePDF(
  in.file = list.files(file.path("/results"), pattern = "Fig4_", full.names = TRUE),
  file="Figure4.pdf"
)

do.call(file.remove, list(list.files("/results/", pattern = "Fig4_", full.names = TRUE)))
