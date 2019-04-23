#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyr)

setwd("Z:/Paula/Mothur")
#Assign varibles for the paths of the data to import
shamSD <- "stability.opti_mcc.0.03.subsample.0.03.Sham_DMM_SD.metastats.txt"
shamDMM <- "stability.opti_mcc.0.03.subsample.0.03.Sham_DMM.metastats.txt" #Phylo table
SDDMM <- "stability.opti_mcc.0.03.subsample.0.03.DMM_SD_DMM.metastats.txt"
taxfile <- "stability.cons.taxonomy" #Phylo tax
metadata <- "rat.group.design" #metadata

shamSD <- read.table(shamSD, sep = "\t", header = T) #Read the OTU table into R
shamDMM <- read.table(shamDMM, sep = "\t", header = T) #Read the OTU table into R
SDDMM <- read.table(SDDMM, sep = "\t", header = T) #Read the OTU table into R

taxonomy <- read.table(taxfile, sep = "\t", header = T) #Read the taxonomy file into R
taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("kingdom", "phylum", "class", "family", "order", "genus", "species"), sep = ";")
str(taxonomy)

meta <- read.table(file = metadata, sep = '\t', header = TRUE)
str(meta)

C1 <- merge(shamSD, taxonomy, by.x="OTU", by.y="OTU")
C2 <- merge(shamDMM, taxonomy, by.x="OTU", by.y="OTU")
C3 <- merge(SDDMM, taxonomy, by.x="OTU", by.y="OTU")

C1$q.value <- p.adjust(C1$p.value, method = "fdr", n =10756)
C2$q.value <- p.adjust(C2$p.value, method = "fdr", n =10756)
C3$q.value <- p.adjust(C3$p.value, method = "fdr", n =10756)

qplot(C1$p.value, geom="histogram")
qplot(C1$q.value, geom="histogram") 
qplot(C2$p.value, geom="histogram")
qplot(C2$q.value, geom="histogram") 
qplot(C3$p.value, geom="histogram")
qplot(C3$q.value, geom="histogram") 

subC1 <- C1[ which(C1$q.value<1
                   & C1$Weigh.Mean>0.001), ]
subC2 <- C2[ which(C2$q.value<1
                   & C2$Weigh.Mean>0.001), ]
subC3 <- C3[ which(C3$q.value<1
                   & C3$Weigh.Mean>0.001), ]

subC1 <- subC1[order(subC1$q.value),]
subC2 <- subC2[order(subC2$q.value),]
subC3 <- subC3[order(subC3$q.value),]

subC1$b1 <- subC1$mean.group1>subC1$mean.group2
subC1$gr <- ifelse((subC1$b1), subC1$mean.group1,-1*subC1$mean.group2)
subC2$b1 <- subC2$mean.group1>subC2$mean.group2
subC2$gr <- ifelse((subC2$b1), subC2$mean.group1,-1*subC2$mean.group2)
subC3$b1 <- subC3$mean.group1>subC3$mean.group2
subC3$gr <- ifelse((subC3$b1), subC3$mean.group1,-1*subC3$mean.group2)

subC1 <- subC1[order(subC1$gr),]
subC2 <- subC2[order(subC2$gr),]
subC3 <- subC3[order(subC3$gr),]

# Fitting Labels 
par(las=2) # make label text perpendicular to axis
par(mar=c(5,16,5,3)) # increase y-axis margin.

#counts <- table(mtcars$gear)
barplot(subC2$gr, main="Sham (+) & DMM (-) Diferential OTUs", horiz=TRUE, names.arg=subC2$genus, cex.names=0.8, col=c(rep("salmon",length(which(subC2$gr<0))),rep("cornflowerblue",length(which(subC2$gr>0)))))

#counts <- table(mtcars$gear)
barplot(subC1$gr, main="Sham (+) & DMM-SD (-) Diferential OTUs", horiz=TRUE, names.arg=subC1$genus, cex.names=0.8, col=c(rep("mediumseagreen",length(which(subC1$gr<0))),rep("cornflowerblue",length(which(subC1$gr>0)))))

#counts <- table(mtcars$gear)
barplot(subC3$gr, main="DMM-SD (+) & DMM (-) Diferential OTUs", horiz=TRUE, names.arg=subC3$genus, cex.names=0.8 , col=c(rep("salmon",length(which(subC3$gr<0))),rep("mediumseagreen",length(which(subC3$gr>0)))))

