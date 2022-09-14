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
plotpath <- "/results/"

###############
# Figure 5A
###############

load(file = paste0(datapath, "TRB_CapTCR_RNAseq.RData"))

shan_trb <- TCRcap_rnaplot.fx(TRB_CapTCR_RNAseq, "estimated_Shannon_RNAseq", "observed_Shannon_TCRCap", "Shannon diversity", 40,3)

#add title
fig5a <- shan_trb + ggtitle(~underline("TRB diversity inference")) + 
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) 

pdf(file = paste0(plotpath,"Fig5_A.pdf"),
    width = 10, height = 10, useDingbats = FALSE)
fig5a
dev.off()

###############
# Figure 5B
###############

load(file = paste0(datapath, "metadata_TRB.RData"))

# linear regression
regression_trb <- lm(log10(estimated_Shannon) ~ log10(Reads), data = metadata_TRB)
summary(regression_trb)

#add residuals
metadata_TRB$residuals <- residuals(regression_trb)

# define outliers as more than abs(sd) for residuals
metadata_TRB$outliergroup <- "none"
metadata_TRB$outliergroup[ metadata_TRB$residuals <= -2*sd(metadata_TRB$residuals)] <- "Outlier_down"
metadata_TRB$outliergroup[ metadata_TRB$residuals >= 2*sd(metadata_TRB$residuals)] <- "Outlier_up"

fig5b <- ggplot(data = metadata_TRB, aes(y = estimated_Shannon, x = Reads, label = cohort)) + 
  geom_point(aes(color = outliergroup), size = 5) +
  scale_color_manual(values = c("Outlier_down" = "#74add1", "Outlier_up" = "#d73027", "none" = "black")) +
  geom_smooth(method = "lm", se = FALSE) + myplot +
  scale_y_continuous(trans = "log10") + scale_x_continuous(trans = "log10") + annotation_logticks(sides = "bl") +
  geom_text_repel(data = subset(metadata_TRB, outliergroup == "Outlier_down"), segment.size = 0.2,
                  box.padding = 0.5, direction = "x", hjust = 1,nudge_x  = 1, min.segment.length = 0, size = 6) +
  geom_text_repel(data = subset(metadata_TRB, outliergroup == "Outlier_up"), segment.size = 0.2,
                  box.padding = 0.5,  direction = "x", hjust = 1, nudge_y  = 1, min.segment.length = 0, size = 6) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 45),
        axis.title = element_text(size = 35),
        plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 35, color = "black"),
        axis.text.y = element_text(size = 35, color = "black")) +
  labs(y = "Estimated Shannon diversity", x = "TCRb reads") +
  ggtitle(~underline("pedNST (n = 361)"))

pdf(file = paste0(plotpath,"Fig5_B.pdf"),
    width = 10, 
    height = 10,
    useDingbats = FALSE)
fig5b
dev.off()

###############
# Figure 5C
###############

load(file = paste0(datapath, "metadata_TRB.RData"))

# nbl
nbl <- metadata_TRB[ metadata_TRB$cohort == "NBL",]

trbplot_nbl <- ggplot(data = nbl, aes(x = immune_cluster, y = estimated_Shannon)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 4, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5, aes(color = immune_cluster)) +
  theme(axis.title.y = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 40,angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 40, color = "black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(size = 40, hjust = 0.5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::label_number(accuracy = 1)) + annotation_logticks(sides = "l") +
  labs(y = paste0("Estimated\nShannon diversity")) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Excluded")), y_position = 5.5,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Neutral")), y_position = 4.5,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Neutral")), y_position = 4,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Excluded")), y_position = 5,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Immune Neutral", "Immune Excluded")), y_position = 4,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  ggtitle(expression(~underline("NBL (n = 116)"))) 

# cns
cns <- metadata_TRB[ metadata_TRB$cohort != "NBL",]

# make color transparent if number of data point are <=2
mytab <- table(cns$immune_cluster) 

clustcol <- cluster_col
clustcol[ names(clustcol) %in% names(mytab)[ mytab <= 2] ] <- "transparent"

trbplot_cns <- ggplot(data = cns, aes(x = immune_cluster, y = estimated_Shannon)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 2, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 40,angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 40, color = "black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.title = element_text(size = 40, hjust = 0.5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = clustcol) +
  scale_y_continuous(trans = "log10", 
                     labels = scales::label_number(accuracy = 1)) + annotation_logticks(sides = "l") +
  labs(y = paste0("Estimated\nShannon diversity")) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Myeloid Predominant")), y_position = 4.3,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Neutral")), y_position = 5.5,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Neutral")), y_position = 4.8,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  ggtitle(expression(~underline("pedCNS (n = 245)"))) 

fig5c <- plot_grid(trbplot_nbl, trbplot_cns, nrow = 1, align = "h", ncol = 2)

pdf(file = paste0(plotpath,"Fig5_C.pdf"), width = 14, height = 10, useDingbats = FALSE)
fig5c
dev.off()

###############
# Figure 5D
###############


###############
# Figure 5E
###############

load(file = paste0(datapath, "metadata_IGH.RData"))

alligs <- c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "IGHM")

mytab <- metadata_igh[,c("immune_cluster","TotalIsotypes", alligs)]

# remove rare IGHs
mytab$IGHD <- NULL
mytab$IGHE <- NULL

# replace NA with 0
mytab[is.na(mytab)] <- 0

#convert to fractions
mytab[,3:9] <- mytab[,3:9]/mytab$TotalIsotypes

#melt and order
mytab_melted <- melt(mytab[,c(1,3:9)])

mytab_melted$variable <- factor(mytab_melted$variable,
                                levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM"))

pairwise.t.test(mytab$IGHG1, mytab$immune_cluster, "none", paired=FALSE, pool.sd=TRUE)
pairwise.t.test(mytab$IGHG3, mytab$immune_cluster, "none", paired=FALSE, pool.sd=TRUE)

fig5e <- ggplot(data = mytab_melted, aes(x = variable, y = value, fill = immune_cluster)) + 
  geom_boxplot(outlier.color = NA, position = "dodge") + 
  myplot + myaxis + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 45),
        legend.title = element_blank(), axis.text.y = element_text(size = 45), 
        axis.text.x = element_text(size = 45),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 45),
        plot.margin = margin(1,0.5,1.5,0.5, "cm")) +
  scale_fill_manual(values = cluster_col) +
  geom_signif(annotation="*",y_position= 0.55, xmin= 0.9, xmax=1.1, textsize = 10, vjust = 0.5) + # for IGHG1
  geom_signif(annotation="*",y_position= 0.65, xmin= 0.9, xmax=1.3, textsize = 10, vjust = 0.5) +# for IGHG1
  geom_signif(annotation="*",y_position= 0.3, xmin= 2.9, xmax=3.3, textsize = 10, vjust = 0.5) +# for IGHG3
  ggtitle(~underline("PedNST (n = 742)")) + 
  labs(y = "Isotype fraction") #+ coord_cartesian(ylim = c(0,0.75))

#plot legend in two rows
lgd <- get_legend(fig5e + theme(legend.position = "bottom") +  guides(fill=guide_legend(nrow=2, byrow=TRUE))) 

pdf(file = paste0(plotpath,"Fig5_E.pdf"), width = 12, height = 10, useDingbats = FALSE)
plot_grid(fig5e, lgd, nrow=2, ncol =1, rel_heights = c(3,0.4))
dev.off()

###############
# Figure 5F
###############

load(file = paste0(datapath, "metadata_IGrep.RData"))

#nbl
nbl <- metadata_igrep[ metadata_igrep$cohort == "NBL",]
giniplot_nbl <- ggplot(data = nbl, aes(x = immune_cluster, y = gini)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 3, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  theme(axis.title.y = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 40,angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 40, color = "black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(1,1,1,2.5, "cm"),
        plot.title = element_text(size = 40, hjust = 0.5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  labs(y = paste0("gini index (Ig)")) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Neutral")), y_position = 1.1,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Excluded")), y_position = 1.2,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Myeloid Predominant", "Immune Excluded")), y_position = 1,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  ggtitle(expression(~underline("NBL (n = 113)"))) +
  scale_y_continuous(expand = c(0.1, 0))

#cns
cns <- metadata_igrep[ metadata_igrep$cohort != "NBL",]

giniplot_cns <- ggplot(data = cns, aes(x = immune_cluster, y = gini)) + 
  geom_beeswarm(color = "grey", size = 5, cex = 2, alpha = 0.7, shape = 16) + 
  geom_boxplot(outlier.shape = NA, fill = NA, lwd = 1.5,aes(color = immune_cluster)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 40,angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 40, color = "black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(1,2.5,1,1, "cm"),
        plot.title = element_text(size = 40, hjust = 0.5)) +
  theme(legend.position = "none") +
  scale_color_manual(values = cluster_col) +
  labs(y = paste0("gini index (Ig)")) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Myeloid Predominant")), y_position = 1,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Neutral")), y_position = 1.1,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  geom_signif(comparisons = list(c("Pediatric Inflamed", "Immune Excluded")), y_position = 1.2,
              map_signif_level=TRUE, textsize = 12, test = "t.test", vjust = 0.5) +
  ggtitle(expression(~underline("pedCNS (n = 248)"))) +
  scale_y_continuous(expand = c(0.1, 0))

pdf(file = paste0(plotpath,"Fig5_F.pdf"),
    width = 14, 
    height = 10,
    useDingbats = FALSE)
plot_grid(giniplot_nbl, giniplot_cns, nrow = 1, align = "h", ncol = 2)
dev.off()

###############
# Figure 5G
###############


###############
# Compile in one file
###############

setwd(plotpath)

plotflow:::mergePDF(
  in.file = list.files(file.path(plotpath), pattern = "Fig5_", full.names = TRUE),
  file="Figure5.pdf"
)

do.call(file.remove, list(list.files(plotpath, pattern = "Fig5_", full.names = TRUE)))
