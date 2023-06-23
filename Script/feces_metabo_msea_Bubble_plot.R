# dotplot or barplot of enrichment analysis for metabolomics
rm(list = ls())
library(argparser)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggtext)
library(magrittr)
library(reshape)
library(psych)

# 导入数据
a <- read.csv("Feces_metabo/Differential_Ma_Ad_pathway_kegg.csv", header = T, check.names = F)
b <- read.csv("Feces_metabo/Differential_Oa_Ad_pathway_kegg.csv", header = T, check.names = F)
c <- read.csv("Feces_metabo/Differential_Oa_Ma_pathway_kegg.csv", header = T, check.names = F)

# management of data a, b, and c
a$Enrichment_ratio <- a$hits / a$expected
a <- a[order(a$Enrichment_ratio, decreasing = T), ]

b$Enrichment_ratio <- b$hits / b$expected
b <- b[order(b$Enrichment_ratio, decreasing = T), ]

c$Enrichment_ratio <- c$hits / c$expected
c <- c[order(c$Enrichment_ratio, decreasing = T), ]



# top 20 metabolic pathways were kept for subsequent analysis in a and b dataset
a2 <- a[, c(1,7,8)]
b2 <- b[, c(1,7,8)]
c2 <- c[, c(1,7,8)] # c因为本身的pathway条目就比较少，因此全部纳入分析

# add group information
a2$Group <- rep("Ma_vs_Ad", dim(a)[1])
b2$Group <- rep("Oa_vs_Ad", dim(b)[1])
c2$Group <- rep("Oa_vs_Ma", dim(c)[1])

# 将3个表进行合并
abc30 <- rbind(a2, b2, c2)
dim(abc30)

# plot
xx = ggplot(abc30, aes(x = Group, y = Pathway)) + 
  geom_point(aes(size = Enrichment_ratio, fill = FDR), alpha = 0.75, shape = 21) + 
  #scale_size_continuous(limits = c(1, 130), range = c(1,13), breaks = c(1,5,10,13)) + 
  labs( x= "", y = "KEGG Enrichment Pathways", size = "Enrichment ratio", fill = "FDR")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.95, hjust = 0.95), 
        axis.text.y = element_text(colour = "black", size = 10), 
        axis.title = element_text(color = "black", face = "bold", size = 13),
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_continuous(type = "viridis")
xx
ggsave("Feces_metabo/Feces_metabolic_KEGG_pathway_AdMaOa_20230103.png", plot = xx, width = 6, height = 8, dpi = 300)



