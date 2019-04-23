library(readr)
library(stringr)
library(ggplot2)
library(tidyr)

setwd("Z:/Paula/Mothur")
alpha_div <- read.table(file = "stability.opti_mcc.groups.summary", sep = "\t", header = T)
meta <- read.table(file = "rat.group.design", sep = "\t", header = TRUE)

alpha_div_merge <- merge(meta, alpha_div, by.x = "Group", by.y = "group")

# #checking boxplots
# qplot(Treatment, chao, geom = "boxplot", colour = Treatment, data = alpha_div_merge, size = I(0.3))
# qplot(Treatment, simpsoneven, geom = "boxplot", colour = Treatment, data = alpha_div_merge, size = I(0.3))
# qplot(Treatment, shannoneven, geom = "boxplot", colour = Treatment, data = alpha_div_merge, size = I(0.3))
# qplot(Treatment, chao, geom = "boxplot", colour = Treatment, data = alpha_div_merge, size = I(0.3))
# qplot(Treatment, simpson, geom = "boxplot", colour = Treatment, data = alpha_div_merge, size = I(0.3))
# qplot(Treatment, shannon, geom = "boxplot", colour = Treatment, data = alpha_div_merge, size = I(0.3))

#Get figures for manuscript.
chao <- ggplot(alpha_div_merge, aes(Treatment, chao)) + 
  geom_boxplot(aes(color = Treatment)) + 
ggsave("chao.png")

shannon <- ggplot(alpha_div_merge, aes(Treatment, shannon)) + 
  geom_boxplot(aes(color = Treatment)) + 
ggsave("shannon.png")

shannon <- ggplot(alpha_div_merge, aes(Treatment, simpson)) + 
  geom_boxplot(aes(color = Treatment)) + 
ggsave("simpson.png")

ad_metrics <- c("chao", "shannon", "simpson")
  for(m in ad_metrics){
    print(m)
    aov_temp <- aov(get(m) ~ Treatment, data = alpha_div_merge)
    summary(aov_temp)
    anova_summary <- as.data.frame(summary(aov_temp)[[1]])
    write.table(anova_summary, file = paste0("anova_", m, ".txt"), sep = "\t", quote = FALSE)
    if (summary(aov_temp)[[1]][["Pr(>F)"]][[1]] < 0.05){
      tukey_out <- TukeyHSD(aov_temp)
      tukey_out_df <- as.data.frame(tukey_out$Treatment)
      tukey_out_df$ad_metric <- m
      if (exists("tukey_summary")) {
        tukey_summary <- rbind(tukey_summary, tukey_out_df)
      } else {
        tukey_summary <- tukey_out_df
      }
    }
  }

tukey_summary$q.value <- p.adjust(tukey_summary$`p adj`, method = "BH")
write.table(tukey_summary, file = "tukey_summary.txt", sep = "\t", quote = FALSE)
