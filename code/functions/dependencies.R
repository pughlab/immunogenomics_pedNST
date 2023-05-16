
mypackages <- c("dplyr", "reshape2", "tidyr",
                "ggplot2", "ggridges", "ggbeeswarm", "ggsignif", "ggrepel",
                "grid", "gridExtra", "gtable", "cowplot", "patchwork",
                "ComplexHeatmap", "circlize",
                "survival", "survminer", "egg",
                "plotflow", "qusage", "broom")

suppressPackageStartupMessages(invisible(lapply(mypackages, library, character.only = TRUE)))
