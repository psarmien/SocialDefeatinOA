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
sharedfile <- "stability.opti_mcc.0.03.subsample.shared" #Phylo table
taxfile <- "stability.cons.taxonomy" #Phylo tax
metadata <- "rat.group.design" #metadata
mdsfile <- "nmds.txt" #nmds for sorting

#first the otu table file
otu_subsample <- read.table(sharedfile, sep = "\t", header = T) #Read the OTU table into R
#setting rownames to sampleID and removing that column from array contents
row.names(otu_subsample) <- otu_subsample$Group
otu_subsample <- otu_subsample[,-c(1:3)]

taxonomy <- read.table(taxfile, sep = "\t", header = T) #Read the taxonomy file into R
taxonomy <- separate(data = taxonomy, col = Taxonomy, into = c("kingdom", "phylum", "class", "family", "order", "genus", "species"), sep = ";")
str(taxonomy)

meta <- read.table(file = metadata, sep = '\t', header = TRUE)
str(meta)

#then the nmds file to merge with the meta data
nmds <- read.table(file = mdsfile, sep =" ", header = TRUE)
metanmds <- merge(nmds, meta, by.x="group", by.y="Group")

#If you want different taxonomic level, find and replace the taxonomic level listed here
rm(otu.summary)
otu.summary <- prop.table(as.matrix(otu_subsample), 1) 
otu_abund <- colSums(otu.summary)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]

#top 30 most abundant genera
num_genera <- 300 # enter the number of genera you want

melt_otu <- melt(otu.summary_sorted[,c(1:num_genera)])
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")
tail(melt_otu)

#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu <- merge(metanmds, melt_otu, by.x = "group", by.y = "Sample")
meta_otu_tax <- merge(meta_otu, taxonomy)
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$Treatment),]
meta_otu_tax$Treatment<- factor(meta_otu_tax$Treatment, levels=unique(as.character(meta_otu_tax$Treatment)) )
cols <- c(1,2,8,9,10,11,12,13,14,15,16)
meta_otu_tax<- meta_otu_tax[which(meta_otu_tax$family=="Bifidobacteriales(100)" |meta_otu_tax$family=="Lactobacillales(100)"),cols ]

ggplot(meta_otu_tax, aes(x =Treatment, y = Abundance, fill=Treatment)) + 
  geom_bar(stat = "identity") +
  facet_grid(meta_otu_tax$Treatment, scales = "free", space = "free") +
  facet_wrap(meta_otu_tax$family)+
  theme(axis.title.x = element_blank()) + 
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Family Abundance (top ", num_genera, " genera)")) +
  ggtitle("Family Composition by treatment and family") 
ggsave("FamilySpecies.png", width = 10, height = 4.5)
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$Treatment, decreasing = FALSE),]
a<- meta_otu_tax[which(meta_otu_tax$family=="Lactobacillales(100)"), ]
b<- meta_otu_tax[which(meta_otu_tax$family=="Bifidobacteriales(100)"),]
ggplot(b, aes(Treatment, Abundance, fill=Treatment)) + 
  geom_boxplot()+
  coord_cartesian(ylim=c(0, 0.001))
  #facet_wrap(meta_otu_tax$family)
  
