rm(list = ls())

library(ggplot2)
library(vegan)
library(ape)
library(plyr)
library(ggExtra)

# (1)加载数据
# 此方法适合用于已经计算好主成分的数据
a <- read.csv("stool_metag/PCoA.csv", header = T, row.names = 1, check.names = F)
names(a)[1:2] <- c("PCA1", "PCA2")

# dim(a)
# [1] 45  2
# 分组标签按照所需要的顺序
a$group <- "NA"
a$group[grep("^Ad",rownames(a))] <- "Ad"
a$group[grep("^Ma",rownames(a))] <- "Ma"
a$group[grep("^Oa",rownames(a))] <- "Oa"
a$group <- factor(a$group, levels = c("Ad", "Ma", "Oa"))

### 出图
pp <- ggplot(a, aes(PCA1, PCA2))+ 
  geom_point(aes(color = group), size = 2.5, alpha = 0.8)+######设置处理group，点的大小
  scale_shape_manual(aes(values = group)) + 
  scale_fill_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))+ ##底色设置
  scale_color_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))+ ##颜色设置
  xlab(paste("PCoA1 ( ",22.68,"%"," )", sep = "")) + ##坐标轴设置
  ylab(paste("PCoA2 ( ",16.54,"%"," )", sep = ""))+##坐标轴设置
  stat_ellipse(aes(fill = group),geom = "polygon",level = 0.95,alpha = 0.1)+
  #geom_polygon(data = border, aes(fill = group), alpha = 0.1, show.legend = FALSE) +
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color = "black",size =9))
pp

# 添加组间ANOSIM检验的结果，以表格展示
#anosim_all <- data.frame("Age"=c("Ad Vs. Ma","Ad Vs. Oa","Ma Vs. Oa"),
#                         "R"=c(0.42, 0.20, 0.22),
#                         "Pvalue"=c(0.001,0.001,0.002),stringsAsFactors = FALSE)
#anosim_all_p <- ggtexttable(anosim_all,rows = NULL, theme=ttheme("minimal", base_size = 6))

# 图形拼接
#pp1 <- pp + annotation_custom(ggplotGrob(anosim_all_p),xmin = 0.1, ymin = 0.06, xmax = 0.2)

#添加边界密度图
pp2 <- ggMarginal(pp,type="density",size=3,margins="both",groupColour=TRUE,groupFill=TRUE)
ggsave("Plots/Feces_Pcoa_plot.png", plot = pp2, width = 6, height = 6, dpi = 300)


# for saliva metagenomic

rm(list = ls())

library(ggplot2)
library(vegan)
library(ape)
library(plyr)
library(ggExtra)

# (1)加载数据
# 此方法适合用于已经计算好主成分的数据
a <- read.csv("saliva_metag/PCoA.csv", header = T, row.names = 1, check.names = F)
names(a)[1:2] <- c("PCA1", "PCA2")

# dim(a)
# [1] 45  2
# 分组标签按照所需要的顺序
a$group <- "NA"
a$group[grep("^Ad",rownames(a))] <- "Ad"
a$group[grep("^Ma",rownames(a))] <- "Ma"
a$group[grep("^Oa",rownames(a))] <- "Oa"
a$group <- factor(a$group, levels = c("Ad", "Ma", "Oa"))

### 出图
pp <- ggplot(a, aes(PCA1, PCA2))+ 
  geom_point(aes(color = group), size = 2.5, alpha = 0.8)+######设置处理group，点的大小
  scale_shape_manual(aes(values = group)) + 
  scale_fill_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))+ ##底色设置
  scale_color_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))+ ##颜色设置
  xlab(paste("PCoA1 ( ",36.62,"%"," )", sep = "")) + ##坐标轴设置
  ylab(paste("PCoA2 ( ",23.80,"%"," )", sep = ""))+##坐标轴设置
  stat_ellipse(aes(fill = group),geom = "polygon",level = 0.95,alpha = 0.1)+
  #geom_polygon(data = border, aes(fill = group), alpha = 0.1, show.legend = FALSE) +
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text = element_text(color = "black",size =9))
pp

# 添加组间ANOSIM检验的结果，以表格展示
#anosim_all <- data.frame("Age"=c("Ad Vs. Ma","Ad Vs. Oa","Ma Vs. Oa"),
#                         "R"=c(0.42, 0.20, 0.22),
#                         "Pvalue"=c(0.001,0.001,0.002),stringsAsFactors = FALSE)
#anosim_all_p <- ggtexttable(anosim_all,rows = NULL, theme=ttheme("minimal", base_size = 6))

# 图形拼接
#pp1 <- pp + annotation_custom(ggplotGrob(anosim_all_p),xmin = 0.1, ymin = 0.06, xmax = 0.2)

#添加边界密度图
pp2 <- ggMarginal(pp,type="density",size=3,margins="both",groupColour=TRUE,groupFill=TRUE)
ggsave("Plots/Saliva_Pcoa_plot.png", plot = pp2, width = 6, height = 6, dpi = 300)




