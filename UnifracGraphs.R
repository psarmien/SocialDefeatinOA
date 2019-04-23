#if (!requireNamespace("BiocManager", quietly = TRUE))
#+     install.packages("BiocManager")
#BiocManager::install("phyloseq", version = "3.8")

library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(plotly)
library(qgraph)
library(phyloseq)

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

setwd("Z:/Paula/Mothur")

###### DATA #####
unifrac <-"stability.1.unweighted.ave.dist" #Distance matrix from unweighted unifrac
metadata <- "rat.group.design" #metadata

M <- import_mothur_dist(unifrac)
unifrac <- metaMDS(M, distance = M, k = 3, trymax=100)
meta <- read.table(file = metadata, sep = '\t', header = TRUE)

#Calculation of the irdination stress

unifrac$stress

#merging MDS and metadata

nmds <-as.data.frame(unifrac$points)    
nmds$group <- rownames(nmds)
metanmds <- merge(meta, nmds, by.x = 'Group', by.y = 'group')
metanmds$Treatment <- factor(metanmds$Treatment)
str(metanmds)

ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Treatment)) +
  labs(x='PC1', y= 'PC2', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("unif_unw.png", height = 5, width = 7)
ggplot(metanmds, aes(x=MDS1, y=MDS3)) + geom_point(aes(color=Treatment)) +
  labs(x='PC1', y= 'PC3', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("unif_unw2.png", height = 5, width = 7)

#Adding in a 95% confidence ellipse
# this generates a dataframe containing the group centroids
metanmds <- metanmds[,1:5]
cols <- c(3,5)
NMDS.mean <- aggregate(metanmds[,cols], list(group=metanmds$Treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean , by.x = 'Treatment', by.y='design')
nmds$group == metanmds$Group  # this is a problem, we need these two dataframes in the same order...
metanmds <- metanmds[match(nmds$group,metanmds$Group),] # ok I think this works as long as $group isnt a factor...
nmds$group == metanmds$Group  # hurray!
colms <- c(1,3)
taxis <- unifrac$points
taxis <- taxis[,colms]
ord <- ordiellipse(taxis, metanmds$Treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')

# this little loop generates a dataframe containing the ellipse data for each group
df_ell <- data.frame()
for (d in levels(metanmds$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds$Treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),NMDS=d))
}
colnames(df_ell) <- c('MDS1', 'MDS3', 'design') # just making it so our column names are consistent

# now we are adding metadata to the ellipse dataframe
# probably an easier way to do this but oh well...
str(df_ell)
df_ell$design <- factor(df_ell$design)
ggplot(metanmds, aes(x=MDS1, y=MDS3)) +
  geom_point(aes(color=Treatment)) + 
  geom_segment(aes(x=MDS1, y=MDS3, xend=groupX, yend=groupY, color=Treatment), size = .3) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS3, group=design, fill=design), alpha = 0.25)+ 
  labs(x='PC1', y= 'PC3', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("Ufrifrac_unwellip3.png", height = 3, width = 8.5)

#Adding in a 95% confidence ellipse
# this generates a dataframe containing the group centroids
metanmds <- metanmds[,1:5]
cols <- c(3,4)
NMDS.mean <- aggregate(metanmds[,cols], list(group=metanmds$Treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean , by.x = 'Treatment', by.y='design')
nmds$group == metanmds$Group  # this is a problem, we need these two dataframes in the same order...
metanmds <- metanmds[match(nmds$group,metanmds$Group),] # ok I think this works as long as $group isnt a factor...
nmds$group == metanmds$Group  # hurray!

ord <- ordiellipse(unifrac, metanmds$Treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')

# this little loop generates a dataframe containing the ellipse data for each group
df_ell <- data.frame()
for (d in levels(metanmds$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds$Treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),NMDS=d))
}
colnames(df_ell) <- c('MDS1', 'MDS2', 'design') # just making it so our column names are consistent

# now we are adding metadata to the ellipse dataframe
# probably an easier way to do this but oh well...
str(df_ell)
df_ell$design <- factor(df_ell$design)
ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Treatment)) + 
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, color=Treatment), size = .3) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS2, group=design, fill=design), alpha = 0.25)+
  labs(x='PC1', y= 'PC2', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("Ufrifrac_unwellip2.png", height = 3, width = 8.5)




unifrac <-"stability.tre1.weighted.ave.dist" #Distance matrix from unweighted unifrac
M <- import_mothur_dist(unifrac)
unifrac <- metaMDS(M, distance = M, k = 3, trymax=100)

#Calculation of the irdination stress

unifrac$stress

#merging MDS and metadata

nmds <-as.data.frame(unifrac$points)    
nmds$group <- rownames(nmds)
metanmds <- merge(meta, nmds, by.x = 'Group', by.y = 'group')
metanmds$Treatment <- factor(metanmds$Treatment)
str(metanmds)

ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Treatment)) +
  labs(x='PC1', y= 'PC2', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("unif_w.png", height = 5, width = 7)
ggplot(metanmds, aes(x=MDS1, y=MDS3)) + geom_point(aes(color=Treatment)) +
  labs(x='PC1', y= 'PC3', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("unif_w2.png", height = 5, width = 7)

#Adding in a 95% confidence ellipse
# this generates a dataframe containing the group centroids
metanmds <- metanmds[,1:5]
cols <- c(3,5)
NMDS.mean <- aggregate(metanmds[,cols], list(group=metanmds$Treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean , by.x = 'Treatment', by.y='design')
nmds$group == metanmds$Group  # this is a problem, we need these two dataframes in the same order...
metanmds <- metanmds[match(nmds$group,metanmds$Group),] # ok I think this works as long as $group isnt a factor...
nmds$group == metanmds$Group  # hurray!
taxis <- unifrac$points
taxis <- taxis[,colms]
ord <- ordiellipse(taxis, metanmds$Treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')

# this little loop generates a dataframe containing the ellipse data for each group
df_ell <- data.frame()
for (d in levels(metanmds$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds$Treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),NMDS=d))
}
colnames(df_ell) <- c('MDS1', 'MDS3', 'design') # just making it so our column names are consistent

# now we are adding metadata to the ellipse dataframe
# probably an easier way to do this but oh well...
str(df_ell)
df_ell$design <- factor(df_ell$design)
ggplot(metanmds, aes(x=MDS1, y=MDS3)) +
  geom_point(aes(color=Treatment)) + 
  geom_segment(aes(x=MDS1, y=MDS3, xend=groupX, yend=groupY, color=Treatment), size = .3) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS3, group=design, fill=design), alpha = 0.25)+ 
  labs(x='PC1', y= 'PC3', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("Ufrifrac_wellip3.png", height = 3, width = 8.5)

#Adding in a 95% confidence ellipse
# this generates a dataframe containing the group centroids
metanmds <- metanmds[,1:5]
cols <- c(3,4)
NMDS.mean <- aggregate(metanmds[,cols], list(group=metanmds$Treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean , by.x = 'Treatment', by.y='design')
nmds$group == metanmds$Group  # this is a problem, we need these two dataframes in the same order...
metanmds <- metanmds[match(nmds$group,metanmds$Group),] # ok I think this works as long as $group isnt a factor...
nmds$group == metanmds$Group  # hurray!

ord <- ordiellipse(unifrac, metanmds$Treatment, label = TRUE, conf = .95, kind = 'se', draw = 'none')

# this little loop generates a dataframe containing the ellipse data for each group
df_ell <- data.frame()
for (d in levels(metanmds$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds$Treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),NMDS=d))
}
colnames(df_ell) <- c('MDS1', 'MDS2', 'design') # just making it so our column names are consistent

# now we are adding metadata to the ellipse dataframe
# probably an easier way to do this but oh well...
str(df_ell)
df_ell$design <- factor(df_ell$design)
ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Treatment)) + 
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, color=Treatment), size = .3) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS2, group=design, fill=design), alpha = 0.25)+
  labs(x='PC1', y= 'PC2', caption = paste('Ordination stress: ', round(unifrac$stress, digits = 2)))
ggsave("Ufrifrac_wellip2.png", height = 3, width = 8.5)



