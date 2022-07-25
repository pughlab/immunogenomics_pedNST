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




###############
# Compile in one file
###############

setwd(plotpath)

plotflow:::mergePDF(
  in.file = list.files(file.path(plotpath), pattern = "Fig6_", full.names = TRUE),
  file="Figure6.pdf"
)
do.call(file.remove, list(list.files(plotpath, pattern = "Fig6_", full.names = TRUE)))
