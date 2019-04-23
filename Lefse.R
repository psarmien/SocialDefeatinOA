#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyr)
library(readxl)

setwd("Z:/Paula/Mothur")

# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   require(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }

lefse_results <- read_excel("lefse_results.xlsx",sheet = "Sheet1", col_types = c("text", "numeric", "text", "numeric", "numeric","blank", "numeric", "text", "text", "text", "text", "text", "text"))
lefse_results<-separate(lefse_results,Genus,into=c("Genus","Number"),sep="_")

# DMM_lefse <- subset(lefse_results, Class...3=="DMM")
# DMM_lefse <- DMM_lefse[1:5,]
# 
# Sham_lefse <- subset(lefse_results, Class...3=="Sham")
# Sham_lefse <- Sham_lefse[1:7,]
# Sham_lefse <- Sham_lefse[-2,]
# Sham_lefse <- Sham_lefse[-5,]
# 
# DMMSD_lefse <- subset(lefse_results, Class...3=="DMM_SD")
# DMMSD_lefse <- DMMSD_lefse[1:6,]
# DMMSD_lefse <- DMMSD_lefse[-4,]
# 
# 
# p1 <- ggplot(Sham_lefse, aes(x = Genus, y = LDA))+
#   geom_bar(stat = "identity", fill = "#a6cee3")+
#   xlab("Genus")+
#   ylab("LDA Score")+
#   ggtitle("Sham")+
#   geom_text(aes(label=Genus), vjust=1, angle=90, size=3)+
#   guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
#   ylim(0,5.5) +
#   theme(legend.position="none", axis.text.x=element_blank())
# p2 <- ggplot(DMM_lefse, aes(x = Genus, y = LDA))+
#   geom_bar(stat = "identity", fill = "#1f78b4")+
#   xlab("Genus")+
#   ylab("LDA Score")+
#   ggtitle("DMM")+
#   geom_text(aes(label=Genus), vjust=1, angle=90, size=3)+
#   guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
#   ylim(0,5.5) +
#   theme(legend.position="none", axis.text.x=element_blank())
# p3 <- ggplot(DMMSD_lefse, aes(x = Genus, y = LDA))+
#   geom_bar(stat = "identity", fill = "#FF6666")+
#   xlab("Genus")+
#   ylab("LDA Score")+
#   ggtitle("DMM_DS")+
#   geom_text(aes(label=Genus), vjust=1, angle=90, size=3)+
#   guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
#   ylim(0,5.5) +
#   theme(legend.position="none", axis.text.x=element_blank())
# 
# multiplot(p1,p2,p3, cols=3)


DMM_lefse <- subset(lefse_results, Class...3=="DMM")
DMM_lefse <- DMM_lefse[1:15,]

Sham_lefse <- subset(lefse_results, Class...3=="Sham")
Sham_lefse <- Sham_lefse[1:15,]

DMMSD_lefse <- subset(lefse_results, Class...3=="DMM SD")
DMMSD_lefse <- DMMSD_lefse[1:15,]


totleft <- rbind(Sham_lefse,DMM_lefse,DMMSD_lefse)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,16,5,3)) # increase y-axis margin.
barplot(totleft$LDA, main="LEFSE Differential OTUs", horiz=TRUE, names.arg=totleft$Genus, cex.names=0.8 , col=c(rep("cornflowerblue",length(Sham_lefse$OTU)),rep("salmon",length(DMM_lefse$OTU)),rep("mediumseagreen",length(DMMSD_lefse$OTU))))
legend("bottomleft", c("Sham", "DMM", "DMM-SD"), col=c("cornflowerblue","salmon","mediumseagreen"), lwd=10)

                    