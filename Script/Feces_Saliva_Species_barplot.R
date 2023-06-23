rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(cowplot)

a <- read.csv("stool_metag/table.s10.tran.group.csv", header = T, check.names = F)

# 将数据进行长宽矩阵转换

a_t <- melt(a)

names(a_t) <- c("Species", "Group", "Abundance")
a_t$Group <- factor(a_t$Group, levels = c("Ad", "Ma", "Oa"))
a_t$Species <- factor(a_t$Species, levels = c("Others","s__Lactobacillus reuteri","s__Corynebacterium xerosis","s__Clostridium sp. CAG:226","s__Clostridium sp. CAG:138",
                                            "s__Lactobacillus amylovorus","s__Lactobacillus johnsonii","s__Bacteroides fragilis","s__Firmicutes bacterium CAG:110","s__Prevotella sp. CAG:279",
                                            "s__Ruminococcus flavefaciens"))

# plot
p <- ggplot(a_t, aes(x = Group, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"),
        #legend.position = "bottom",
        axis.ticks=element_line(size=1.2),
        panel.background = element_rect(colour = "black", size = 1)) +
  #guides(fill = guide_legend(ncol = 2)) +
  ylab("Relative Abundance") +
  xlab("") +
  labs(title = "Feces") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + # 以百分数进行Y轴的展示
  
  scale_fill_manual(values = c("#b0b9b9","#ffe327","#2f4e87","#f69001","#f0eedf",
                                        "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
                                        "#262a35"))
#p  
#ggsave("Plots/Feces_Plylum_percent_barplot.png", plot = p, width = 4, height = 6, dpi = 300)  
    

b <- read.csv("saliva_metag/table.s10.tran.group.csv", header = T, check.names = F)

# 将数据进行长宽矩阵转换

b_t <- melt(b)

names(b_t) <- c("Species", "Group", "Abundance")
b_t$Group <- factor(b_t$Group, levels = c("Ad", "Ma", "Oa"))
b_t$Species <- factor(b_t$Species, levels = c("Others","s__Chryseobacterium taklimakanense","s__Moraxella pluranimalium","s__[Haemophilus] parasuis","s__Flavobacterium ummariense",
                                            "s__Streptococcus suis","s__Moraxella porci","s__Chryseobacterium gleum","s__Rothia nasimurium","s__Bergeyella zoohelcum",
                                            "s__Neisseria dentiae"))

# plot
p1 <- ggplot(b_t, aes(x = Group, y = Abundance, fill = Species)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"),
        #legend.position = "bottom",
        axis.ticks=element_line(size=1.2),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Relative Abundance") +
  xlab("") +
  labs(title = "Saliva") + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + # 以百分数进行Y轴的展示
  
  scale_fill_manual(values = c("#b0b9b9","#ffe327","#2f4e87","#f69001","#f0eedf",
                                        "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c","#262a35")) 
                                        # 拼图，参考：https://mp.weixin.qq.com/s/qycXFTSt0VtHmm7z3NQcpg
p_total <- p / p1 +
  plot_annotation(tag_levels = c("A","C"))

ggsave("Plots/Feces_saliva_barplot_Species.png", plot = p_total, width = 5, height = 10, dpi = 300)                                    
