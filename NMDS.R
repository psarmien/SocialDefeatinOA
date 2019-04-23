library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("Z:/Paula/Mothur")

##### FUNCTIONS #####

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method, permutations = 9999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

###### DATA #####
otu_table <-"stability.opti_mcc.0.03.subsample.shared" #rarefied OTU table
metadata <- "rat.group.design" #metadata

#Read in OTU table
otu_subsample <- read.table(otu_table, header = TRUE)

#Remove sequence ID and leave just sample names
#otu_subsample <- separate(data = otu_subsample, col = Group, into = c("sequencing_id", "Group"), sep = "_")
#otu_subsample <- otu_subsample[-c(1,2),-2] #remove unneeded rows and columns

#Stores the sample name info as the rownames of the dataframe rather
rownames(otu_subsample) <- otu_subsample$Group

#Read in metadata
meta <- read.table(file = metadata, sep = '\t', header = TRUE)

#Makes sure that the meta table and the otu table have the same samples
meta <- meta[meta$Group %in% rownames(otu_subsample),] 
otu_subsample <- otu_subsample[rownames(otu_subsample) %in% meta$Group,]

# removes extra info that mothur includes in their OTU tables and outlier points
otu_subsample <- otu_subsample[,-c(1:3)]  

##################################################
##################################################

# this calculates the distance matrix using Bray-Curtis distances with vegan 
dist.matr.bray <- vegdist(otu_subsample, method = 'bray')

# this is vegan's function to make an NMDS ordination using k=2 dimensions

mds <- metaMDS(dist.matr.bray, k = 5,trymax = 1000, autotransform = FALSE)

#Calculation of the irdination stress

mds$stress

#merging MDS and metadata

nmds <-as.data.frame(mds$points)    
nmds$group <- rownames(nmds)
metanmds <- merge(meta, nmds, by.x = 'Group', by.y = 'group')
metanmds$Treatment <- factor(metanmds$Treatment)
str(metanmds)
#write.csv(metanmds)

###### PLOTS #######

#General plots with basic facets

ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Treatment))


#Plot for manuscript
ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Treatment)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(mds$stress, digits = 2)))
ggsave("nmds_all.png", height = 5, width = 7)
ordiellipse(mds, metanmds$Treatment, display = "sites", kind = "sd")
ordiellipse(sol, MyMeta$amt, display = "sites", kind = "sd", label = T)
#Save MDS data for later
write.table(nmds, file="nmds.txt")


#Adding in a 95% confidence ellipse

############

# this generates a dataframe containing the group centroids
NMDS.mean <- aggregate(metanmds[,3:4], list(group=metanmds$Treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean , by.x = 'Treatment', by.y='design')

str(metanmds)
#metanmds$day <- factor(metanmds$day)


nmds$group == metanmds$Group  # this is a problem, we need these two dataframes in the same order...
#metanmds$group <- as.character(metanmds$group)
metanmds <- metanmds[match(nmds$group,metanmds$Group),] # ok I think this works as long as $group isnt a factor...
nmds$group == metanmds$Group  # hurray!

ord <- ordiellipse(mds, metanmds$Treatment, label = TRUE, conf = .99, kind = 'se', draw = 'none')

# this little loop generates a dataframe containing the ellipse data for each group

df_ell <- data.frame()
for (d in levels(metanmds$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(metanmds[metanmds$Treatment == d,],
                                                   veganCovEllipse(ord[[d]]$cov, ord[[d]]$center, ord[[d]]$scale))),NMDS=d))
}

colnames(df_ell) <- c('MDS1', 'MDS2', 'design') # just making it so our column names are consistent

# now we are adding metadata to the ellipse dataframe
# probably an easier way to do this but oh well...

#meta_sub <- meta[,-1]
#meta_sub2 <- unique(meta_sub)
#df_ell2 <- merge(df_ell, meta_sub2, by.x = 'design', by.y = 'day_location_treatment')
str(df_ell)
df_ell$design <- factor(df_ell$design)



ggplot(metanmds, aes(x=MDS1, y=MDS2)) + geom_point(aes(color=Treatment)) +
  labs(x='Axis 1', y= 'Axis 2', caption = paste('Ordination stress: ', round(mds$stress, digits = 2))) +
  geom_path(data = df_ell, aes(x=MDS1, y=MDS2, group=design, color=design))


ggplot(metanmds, aes(x=groupX, y=groupY)) +
  geom_point(aes(color=Treatment)) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS2, group=design, fill=design), alpha = 0.25) + 
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(mds$stress, digits = 2))) 

ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=Treatment)) + 
  geom_segment(aes(x=MDS1, y=MDS2, xend=groupX, yend=groupY, color=Treatment), size = .3) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS2, group=design, fill=design), alpha = 0.25) 
ggsave("ellipses2.png", height = 3, width = 8.5)


#Adding in a 95% confidence ellipse

############

# this generates a dataframe containing the group centroids
metanmds <- metanmds[,1:7]
cols <- c(3,5)
NMDS.mean <- aggregate(metanmds[,cols], list(group=metanmds$Treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean , by.x = 'Treatment', by.y='design')

str(metanmds)
#metanmds$day <- factor(metanmds$day)


nmds$group == metanmds$Group  # this is a problem, we need these two dataframes in the same order...
#metanmds$group <- as.character(metanmds$group)
metanmds <- metanmds[match(nmds$group,metanmds$Group),] # ok I think this works as long as $group isnt a factor...
nmds$group == metanmds$Group  # hurray!
colms<-c(1,3)
taxis <- mds$points
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

#meta_sub <- meta[,-1]
#meta_sub2 <- unique(meta_sub)
#df_ell2 <- merge(df_ell, meta_sub2, by.x = 'design', by.y = 'day_location_treatment')
str(df_ell)
df_ell$design <- factor(df_ell$design)

ggplot(metanmds, aes(x=MDS1, y=MDS3)) + geom_point(aes(color=Treatment)) +
  labs(x='NMDS1', y= 'NMDS3', caption = paste('Ordination stress: ', round(mds$stress, digits = 2))) +
  geom_path(data = df_ell, aes(x=MDS1, y=MDS3, group=design, color=design))


ggplot(metanmds, aes(x=groupX, y=groupY)) +
  geom_point(aes(color=Treatment)) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS3, group=design, fill=design), alpha = 0.25) + 
  labs(x='NMDS 1', y= 'NMDS 3', caption = paste('Ordination stress: ', round(mds$stress, digits = 2))) 

ggplot(metanmds, aes(x=MDS1, y=MDS3)) +
  geom_point(aes(color=Treatment)) + 
  geom_segment(aes(x=MDS1, y=MDS3, xend=groupX, yend=groupY, color=Treatment), size = .3) + 
  geom_polygon(data = df_ell, aes(x=MDS1, y=MDS3, group=design, fill=design), alpha = 0.25) 
ggsave("ellipses3.png", height = 3, width = 8.5)
