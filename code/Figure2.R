###############
# Figure 2
###############

source("/code/functions/dependencies.R")
source("/code/functions/color_schemes.R")
source("/code/functions/survival_functions.R")

datapath <- "/data/"

###############
# Figure 2A
###############

load(file = paste0(datapath,"metadata_IC.RData"))

metadata_gender <- metadata_IC[!is.na(metadata_IC$gender),]
metadata_gender <- metadata_gender[metadata_gender$gender != "Unknown",]

freqtab <- metadata_gender %>% group_by(immune_cluster,gender) %>%
  summarise(n = n()) %>% mutate(freq = n / sum(n)) 

fig2a <- stacked_plots(freqtab, "immune_cluster", "gender") + labs(y = "Frequency")

pdf("/results/Fig2_A.pdf",
    width = 12, height = 10, onefile = F)
fig2a
dev.off()

###############
# Figure 2B
###############

load(file = paste0(datapath,"metadata_IC.RData"))

metadata_race <- metadata_IC[!is.na(metadata_IC$race),]
metadata_race <- metadata_race[metadata_race$race != "Unknown",]

freqtab <- metadata_race %>% group_by(immune_cluster,race) %>%
  summarise(n = n()) %>% mutate(freq = n / sum(n)) 

fig2b <- stacked_plots(freqtab, "immune_cluster", "race") + labs(y = "Frequency")

pdf("/results/Fig2_B.pdf",
    width = 12, height = 10, onefile = F)
fig2b
dev.off()

###############
# Figure 2C
###############

load(file = paste0(datapath,"metadata_IC.RData"))

metadata_age <- metadata_IC[!is.na(metadata_IC$age_at_diagnosis),]

fig2c <- ggplot(data = metadata_age, aes(x = immune_cluster, y = age_at_diagnosis)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 0.7) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  theme(axis.title.y = element_text(size = 45),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 45,angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 45, color = "black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(size = 45)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  labs(y = "Age at diagnosis (years)")

pdf("/results/Fig2_C.pdf",
    width = 12, height = 10, onefile = F)
fig2c
dev.off()

###############
# Figure 2D
###############

load(file = paste0(datapath,"metadata_IC.RData"))

metadata_IC$vital_status <- as.numeric(as.character(metadata_IC$vital_status))

cluster_names <- c("Pediatric Inflamed", "Myeloid Predominant", "Immune Neutral", "Immune Excluded")

sfit <- survfit(Surv(days_to_death, vital_status) ~ immune_cluster, data=metadata_IC)
fig2d <- KM_plot(metadata_IC, sfit, cluster_col, "", "Overall survival", cluster_names)

pdf("/results/Fig2_D.pdf",
    width = 10, height = 12, onefile = F)
fig2d
dev.off()

###############
# Figure 2E
###############

load(file = paste0(datapath,"metadata_IC.RData"))

metadata_IC$recurrence <- as.numeric(as.character(metadata_IC$recurrence))

sfit <- survfit(Surv(days_to_progress, recurrence)~ immune_cluster, data=metadata_IC)
fig2e <- KM_plot(metadata_IC, sfit, cluster_col, "", "Progression-free survival", cluster_names)

pdf("/results/Fig2_E.pdf",
    width = 10, height = 12, onefile = F)
fig2e
dev.off()

###############
# Compile in one file
###############

setwd("/results")

plotflow:::mergePDF(
  in.file = list.files(file.path("/results/"), pattern = "Fig2_", full.names = TRUE),
  file="Figure2.pdf"
)

do.call(file.remove, list(list.files("/results/", pattern = "Fig2_", full.names = TRUE)))


