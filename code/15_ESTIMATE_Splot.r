
library(dplyr)
library(reshape2)
library(ggbeeswarm)
library(gridExtra)
library(gtable)
library(grid)

source("~/OneDrive - UHN/R_src/ggplot2_theme.R")
source("~/OneDrive - UHN/R_src/color_schemes.R")

manifestpath <- "/Users/anabbi/OneDrive - UHN/Documents/IPD2/Manifests/"
datapath <- "/Users/anabbi/OneDrive - UHN/Documents/IPD2/Data/"
plotpath <- "/Users/anabbi/OneDrive - UHN/Documents/IPD2/Plots/"

load(paste0(datapath,"ESTIMATE/estimate_manifest_primary_clean_final.RData"))

load(paste0(datapath,"ESTIMATE/PDX_estimate.RData"))

head(PDX_estimate)

PDX_estimate_itcc <- read.table(paste0(datapath, "ESTIMATE/ESTIMATE_output/RNAseq_itcc-P4_TPM_values_210726_estimateOutput.txt"),
                               sep = "\t", stringsAsFactors = F, header = T)

PDX_estimate_itcc <- PDX_estimate_itcc[grep(paste(c("EP", "MB", "HG", "NB", "RT"),collapse="|"),PDX_estimate_itcc$NAME),]

PDX_estimate_itcc <- PDX_estimate_itcc[grep(paste(c("PP", "PT", "PR", "PM"),collapse="|"),PDX_estimate_itcc$NAME),]

dim(PDX_estimate_itcc)

#remove
PDX_estimate_itcc <- PDX_estimate_itcc[!PDX_estimate_itcc$NAME %in% c("ITCC.P4_s07_HG0068_PT01_F01_R01_A03", "ITCC.P4_s15_NB0537_PT03_F01_R01_03"), ]

colnames(PDX_estimate_itcc)[1] <- "SampleID"

PDX_estimate <- rbind(PDX_estimate, PDX_estimate_itcc)

dim(PDX_estimate)

mycol <- c("aliquot_id", "StromalScore", "ImmuneScore","ESTIMATEScore", "sample_id", "group", "cohort")

estimate_df <- estimate_manifest_primary_clean[,mycol]

summary(PDX_estimate$ImmuneScore)

PDX_estimate$group <- "ICGC"
PDX_estimate$cohort <- "PDX"
colnames(PDX_estimate)[1] <- "aliquot_id"
PDX_estimate$sample_id <- PDX_estimate$aliquot_id

head(estimate_df)

estimate_ped_pdx <- rbind(estimate_df,PDX_estimate)

head(estimate_ped_pdx)

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

immune.cohorts

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

(immune.cohorts)

immune.cohorts$medianloc <- immune.cohorts$Median.start+((immune.cohorts$Median.stop-immune.cohorts$Median.start)/2)

sorted.estimate_ped_pdx$cohort <- factor(sorted.estimate_ped_pdx$cohort, levels = levels(immune.cohorts$cohort))     

head(sorted.estimate_ped_pdx)

rmEMPTY <- rep("black",22)
rmEMPTY[14:15] <- "white"

immune.cohorts$color_crossbar <- NA

immune.cohorts$color_crossbar[immune.cohorts$cohort == "EMPTY1"] <- "white"
immune.cohorts$color_crossbar[immune.cohorts$cohort == "EMPTY2"] <- "white"

immune.cohorts$color_crossbar[is.na(immune.cohorts$color_crossbar)] <- "black"

cohort_col

immuneplot <- ggplot() + 
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

immuneplot

pdf(file = paste0(plotpath,"Immunereads_Splot.pdf"),
        width = 20, height = 8, useDingbats = FALSE)
print(immuneplot)
dev.off()




