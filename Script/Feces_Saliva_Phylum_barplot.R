rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(cowplot)

a <- read.csv("stool_metag/Feces.Unigenes.group.relative.p.csv", header = T, check.names = F)

# 将数据进行长宽矩阵转换

a_t <- melt(a)

names(a_t) <- c("Phylum", "Group", "Abundance")
a_t$Group <- factor(a_t$Group, levels = c("Ad", "Ma", "Oa"))
a_t$Phylum <- factor(a_t$Phylum, levels = c("Others","Ascomycota","Planctomycetes","Tenericutes","Lentisphaerae",
                                            "Chlamydiae","Synergistetes","Fibrobacteres","Verrucomicrobia","Fusobacteria",
                                            "Spirochaetes","Euryarchaeota","Proteobacteria","Actinobacteria","Bacteroidetes",
                                            "Firmicutes"))

# plot
p <- ggplot(a_t, aes(x = Group, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"),
        axis.ticks=element_line(size=1.2),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Relative Abundance") +
  xlab("") +
  labs(title = "Feces") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + # 以百分数进行Y轴的展示
  
  scale_fill_manual(values = c("#b0b9b9","#ffe327","#2f4e87","#f69001","#f0eedf",
                               "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
                               "#262a35","#c5942e","#C696C4","#006b31","#fddd00",
                               "#d80b13")) 
p  
#ggsave("Plots/Feces_Plylum_percent_barplot.png", plot = p, width = 4, height = 6, dpi = 300)  


## for saliva
b <- read.csv("saliva_metag/Saliva.Unigenes.group.relative.p.csv", header = T, check.names = F)

# 将数据进行长宽矩阵转换

b_t <- melt(b)

names(b_t) <- c("Phylum", "Group", "Abundance")
b_t$Group <- factor(b_t$Group, levels = c("Ad", "Ma", "Oa"))
b_t$Phylum <- factor(b_t$Phylum, levels = c("Others","Chlamydiae","Candidatus Gracilibacteria","Basidiomycota","Mucoromycota",
                                            "Elusimicrobia","Chloroflexi","Euryarchaeota","Spirochaetes","Synergistetes",
                                            "Candidatus Saccharibacteria","Fusobacteria","Firmicutes","Actinobacteria","Proteobacteria",
                                            "Bacteroidetes"))

# plot
p1 <- ggplot(b_t, aes(x = Group, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  theme_classic() + 
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16, color = "black", face = "bold"),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width=unit(0.5,"cm"),
        axis.ticks=element_line(size=1.2),
        panel.background = element_rect(colour = "black", size = 1)) +
  ylab("Relative Abundance") +
  xlab("") +
  labs(title = "Saliva") + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + # 以百分数进行Y轴的展示
  
  scale_fill_manual(values = c("#b0b9b9","#ffe327","#2f4e87","#f69001","#f0eedf",
                               "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
                               "#262a35","#c5942e","#C696C4","#006b31","#fddd00",
                               "#d80b13")) 
# 拼图，参考：https://mp.weixin.qq.com/s/qycXFTSt0VtHmm7z3NQcpg
p_total <- p + p1 +
  plot_annotation(tag_levels = "A")

ggsave("Plots/Feces_saliva_barplot_phylum.png", plot = p_total, width = 8, height = 6, dpi = 300)















