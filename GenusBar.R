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
names(taxonomy) <- c("OTU", "Size","kindom","phylum","class","family","order","genus")
taxonomy<-separate(taxonomy,phylum,into=c("phylum","PNumber"),sep=" ")
taxonomy<-separate(taxonomy,family,into=c("family","FNumber"),sep=" ")
taxonomy<-separate(taxonomy,genus,into=c("genus","GNumber"),sep=" ")
str(taxonomy)

meta <- read.table(file = metadata, sep = '\t', header = TRUE)
str(meta)

#then the nmds file to merge with the meta data
nmds <- read.table(file = mdsfile, sep =" ", header = TRUE)
metanmds <- merge(nmds, meta, by.x="group", by.y="Group")

my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black", "violetred2", "darkslategray1", "gold"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
rm(otu.summary)
otu.summary <- prop.table(as.matrix(otu_subsample), 1) 
otu_abund <- colSums(otu.summary)
otu.summary <- rbind(otu_abund, otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]

#top 30 most abundant genera
num_genera <- 30 # enter the number of genera you want

melt_otu <- melt(otu.summary_sorted[,c(1:num_genera)])
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")
tail(melt_otu)


#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu <- merge(metanmds, melt_otu, by.x = "group", by.y = "Sample")
meta_otu_tax <- merge(meta_otu, taxonomy)
str(meta_otu_tax)
summary(meta_otu_tax$group)
#sorting based on MDS1 from negative to positive (NMDS axis 1)
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$MDS1),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$kingdom),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$phylum),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$class),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$family),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$order),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$Treatment),]
#ordering samples based on NMDS axis 1
meta_otu_tax$group <- factor(meta_otu_tax$group, levels=unique(as.character(meta_otu_tax$group)) )

#MAKE A GRAPH! Plot individuals not group means
ggplot(meta_otu_tax, aes(x = group, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  facet_grid(meta_otu_tax$Treatment, scales = "free", space = "free") +
  #facet_wrap(~Treatment)+
  scale_fill_manual(values=my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Genus Composition by mice sorted by NMDS Axis 1") 
ggsave("GenusBarPlotNMDS1_Phylo.png", width = 10, height = 4)

meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$MDS1),]
ggplot(meta_otu_tax, aes(x = group, y = Abundance, fill = OTU)) + 
  geom_bar(stat = "identity") +
  facet_grid(meta_otu_tax$Treatment, scales = "free", space = "free") +
  #facet_wrap(~Treatment)+
  scale_fill_manual(values=my_colors, label=meta_otu_tax$genus) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Genus Composition by mice sorted by NMDS Axis 1") 
ggsave("GenusBarPlotNMDS1_Phylo2.png", width = 10, height = 4)

#top 300 most abundant genera
num_genera <- 300 # enter the number of genera you want

melt_otu <- melt(otu.summary_sorted[,c(1:num_genera)])
colnames(melt_otu) <- c("Sample", "OTU", "Abundance")
tail(melt_otu)


#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu <- merge(metanmds, melt_otu, by.x = "group", by.y = "Sample")
meta_otu_tax <- merge(meta_otu, taxonomy)
str(meta_otu_tax)
summary(meta_otu_tax$group)
#sorting based on MDS1 from negative to positive (NMDS axis 1)
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$MDS1),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$class),]
meta_otu_tax <- meta_otu_tax[order(meta_otu_tax$Treatment),]
#ordering samples based on NMDS axis 1
meta_otu_tax$group <- factor(meta_otu_tax$group, levels=unique(as.character(meta_otu_tax$group)) )

#MAKE A GRAPH! Plot individuals not group means
ggplot(meta_otu_tax, aes(x = group, y = Abundance, fill = family)) + 
  geom_bar(stat = "identity") +
  facet_grid(meta_otu_tax$Treatment, scales = "free", space = "free") +
  #facet_wrap(~Treatment)+
  #scale_fill_manual(values=my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Family Composition by mice sorted by NMDS Axis 1") 
ggsave("GenusBarPlotNMDS1_Phylo_Family2.png", width = 10, height = 4.5)

#MAKE A GRAPH! Plot individuals not group means
ggplot(meta_otu_tax, aes(x = group, y = Abundance, fill =phylum)) + 
  geom_bar(stat = "identity") +
  facet_grid(meta_otu_tax$Treatment, scales = "free", space = "free") +
  #facet_wrap(~Treatment)+
  scale_fill_manual(values=my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Phylum Composition by mice sorted by NMDS Axis 1") 
ggsave("GenusBarPlotNMDS1_Phylo_Phylum.png", width = 10, height = 4)

# 
# 
# #making graph of next 16 most abundant genera
# num_genera <- 20 # enter the number of genera you want
# 
# melt_otu32 <- melt(otu.summary_sorted[,c(21:40)])
# colnames(melt_otu32) <- c("Sample", "OTU", "Abundance")
# tail(melt_otu32)
# 
# 
# #Putting it all together: merge melt_otu, metadata, taxonomy tables
# meta_otu32 <- merge(metanmds, melt_otu32, by.x = "group", by.y = "Sample")
# meta_otu32_tax <- merge(meta_otu32, taxonomy)
# str(meta_otu32_tax)
# #sorting based on MDS1 from negative to positive (NMDS axis 1)
# meta_otu32_tax <- meta_otu32_tax[order(meta_otu32_tax$MDS1),]
# #ordering samples based on NMDS axis 1
# meta_otu32_tax$group <- factor(meta_otu32_tax$group, levels=unique(as.character(meta_otu32_tax$group)) )
# 
# ggplot(meta_otu32_tax, aes(x = group, y = Abundance, fill = genus)) + 
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = my_colors) +
#   # Remove x axis title
#   theme(axis.title.x = element_blank()) + 
#   ylim(c(0,1)) +
#   guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
#   theme(legend.text=element_text(size=8)) +
#   #theme(legend.position="bottom") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ylab(paste0("Relative Abundance (second ", num_genera, " top genera)")) +
#   ggtitle("Genus Composition by Bull sorted by NMDS Axis 1") 
# ggsave("GenusBarPlot_17to32_NMDS1_Phylo.png", width = 10, height = 4)
# 
# 
