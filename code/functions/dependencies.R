
mypackages <- c("dplyr", "reshape2", 
                "ggplot2", "ggridges", "ggbeeswarm", "ggsignif",
                "grid", "gridExtra", "gtable", "cowplot", 
                "ComplexHeatmap", "circlize",
                "plotflow")

invisible(lapply(mypackages, library, character.only = TRUE))
