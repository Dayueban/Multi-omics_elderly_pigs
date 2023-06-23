# represent the ratio of Firmicutes and Bacteroidetes
rm(list = ls())

# load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(cowplot)

# load dataset
# feces
a <- read.csv("stool_metag/Feces_Firmicutes_vs_Bacteroidetes_ratio.csv", header = T, check.names = F)

# saliva
#b <- read.csv("saliva_metag/Saliva_Firmicutes_vs_Bacteroidetes_ratio.csv", header = T, check.names = F)

# combine two dataset based on row
#ab <- rbind(a, b)

# add group information
ab <- a
ab <- ab %>%
  mutate(log2 = log2(ab$`Firmicutes/Bacteroidetes`))
ab$group <- NA
ab$group[grep("^Ad", ab$ID)] <- "Ad"
ab$group[grep("^Ma", ab$ID)] <- "Ma"
ab$group[grep("^Oa", ab$ID)] <- "Oa"

# add bodysite information
#ab$bodysite <- NA
#ab$bodysite[grep("*Fe", ab$ID)] <- "Feces"
#ab$bodysite[grep("*Sa", ab$ID)] <- "Saliva"


my_comparisons_a <- list(c("Ad","Ma"), c("Ad", "Oa"), c("Ma", "Oa"))
p1 <- ggplot(ab, aes(x = group,y = log2, fill = group)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color = "white", size = 0.75) +
  #facet_grid(~bodysite) +
  
  #Change the order of items in the legend
  theme_classic() +
  scale_x_discrete(limits=c("Ad", "Ma", "Oa")) +
  stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Log2 Firmucutes/Bacteroidetes ratio") +
  xlab("") +
  labs(title = "Feces") +
        #strip.background = element_rect(fill=c("#C2C4C5")), # 分面的背景色设置
        #strip.text = element_text(size = 15,colour = "black", face = "bold")) +  # 分面的字体设置
  scale_color_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  scale_fill_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E")) + 
  geom_jitter(position=position_jitter(0.2))
p1
ggsave("Plots/Feces_Firmucutes_vs_Bacteroidetes_ratio.png", plot = p1, width = 3, height = 6)


