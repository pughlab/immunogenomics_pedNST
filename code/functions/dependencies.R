
mypackages <- c("dplyr", "reshape2", "tidyr",
                "ggplot2", "ggridges", "ggbeeswarm", "ggsignif", "ggrepel",
                "grid", "gridExtra", "gtable", "cowplot", 
                "ComplexHeatmap", "circlize",
                "survival", "survminer",
                "plotflow", "qusage")

invisible(lapply(mypackages, library, character.only = TRUE))
