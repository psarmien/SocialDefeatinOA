library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("Z:/Paula/Mothur")

metadatarats <- read.table(file = "rat.group.design", sep = "\t", header = TRUE)

V1rfact <- read_tsv(file = "stability.opti_mcc.groups.rarefaction")%>%
  select(-contains("lci-"), -contains("hci-")) %>%
  gather(-numsampled, key=sample, value=coverage) %>%
  mutate(sample=str_replace_all(sample, pattern="0.03-", replacement="")) %>%
  drop_na()
#V1rfact <- separate(data = V4rfact, col = sample, into = c("sequencing_id", "sample"), sep = "_") #removing sequencing ID
#V1rfact <- V4rfact[,-2]

V1meta_rare <- metadatarats%>%
  #sample_n(20) %>%
  merge(., V1rfact, by.x= "Group", by.y = "sample")

#graph rarefaction plot with vertical line where subsampling cutoff 
ggplot(V1meta_rare, aes(x=numsampled, y=coverage, group=Group, color=Treatment)) +
  geom_line()+
  geom_vline(xintercept=65008) +
  coord_cartesian(ylim=c(0,2000)) +
  labs(x="Number of Sequences Sampled per Subject",
       y="Number of OTUs per Subject") +
  theme_classic()
ggsave("Rarefaction_all.png")