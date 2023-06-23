# for PCA analsis with Lipid data
rm(list = ls())

# load packages
library(ade4)
library(MASS)
library(ggplot2)
library(ggpubr)
library(devtools)
library("FactoMineR")
library("factoextra")
library(ggbiplot)
library(vegan)
library(readxl)
library(ggExtra)

# load data
a <- read.csv("stool_metag/Unigenes.absolute.level3.csv", header = T, row.names = 1, check.names = F)
#a_t <- data.frame(t(a))
# (2)数据处理
sum1 <- apply(a[,1:ncol(a)], 1, sum)
sum2 <- apply(a[,1:ncol(a)], 1, function(x)length(x[x == 0]))
sumall <- sum(sum1)
num_sample <- ncol(a)-1

## 
a$ratio1 <- (sum1 / sumall) * 100
a$ratio1 <- paste(sprintf("%.2f", a$ratio1), "%", sep = "")
a$ratio2 <- (sum2 / num_sample)
a$ratio2 <- paste(sprintf("%.2f", a$ratio2), "%", sep = "")

## to exclude species which with abundance less than 0.01% and present 
## in at least 20% of the study cohort
keep <- intersect(rownames(a)[which(a$ratio1 > 0.05)],
                  rownames(a)[which(a$ratio2 < 0.8)])

## 保留keep变量里面的那些species
a2 <- a[rownames(a) %in% keep, ]
a2$ratio1 <- NULL
a2$ratio2 <- NULL

a2_rel <- sweep(a2, 2, colSums(a2), "/")


# PCA analysis
a.bc <- vegdist(t(a2_rel), na.rm = TRUE)  

a.bc.ma <- as.matrix(a.bc)
a.bc.df <- as.data.frame(a.bc.ma)
a.pca <- PCA(a.bc.df,graph = FALSE)

#' add the group information according to the rownames
groups <- rep("NA", ncol(a2))
groups[grep("^Ad",names(a2))] <- "Ad"
groups[grep("^Ma",names(a2))] <- "Ma"
groups[grep("^Oa",names(a2))] <- "Oa"

# 只取前面两个主成分
b <- as.data.frame(a.pca$ind$coord[,1:2])
b_PC <- a.pca$eig[,2][1:2]
# comp 1   comp 2 
# 54.72515 38.49288
b$group <- groups
# 分组标签按照所需要的顺序
b$group <- factor(b$group, levels = c("Ad", "Ma", "Oa"))
names(b)[1:2] <- c("PCA1", "PCA2")

### 出图
pp <- ggplot(b, aes(PCA1, PCA2))+ 
  geom_point(aes(color = group), size = 2.5, alpha = 0.8)+######设置处理group，点的大小
  scale_shape_manual(aes(values = group)) + 
  scale_fill_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))+ ##底色设置
  scale_color_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))+ ##颜色设置
  xlab(paste("PCoA1 ( ",round(b_PC[1], digits = 2),"%"," )", sep = "")) + ##坐标轴设置
  ylab(paste("PCoA2 ( ",round(b_PC[2], digits = 2),"%"," )", sep = ""))+##坐标轴设置
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
ggsave("Serum_lipid/Lipid_PCA_plot.png", plot = pp2, width = 6, height = 6, dpi = 300)


# 另外一种作图方式
library(doBy)
se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]
stat <- summaryBy(PCA1 + PCA2 ~ group, b, FUN = c(mean, sd, se, conf_95))

# (3)添加各组样本与其质心点的连线
# 各组样本与质心的连线
# 对于质心点，本示例直接使用各组的均值，已经在上一步中计算得到
plot_data <- merge(b, stat, by = 'group')

# (4) 作图1
p <- ggplot(plot_data) +
  geom_segment(aes(x = PCA1.mean, y = PCA2.mean, xend = PCA1, yend = PCA2, color = group), show.legend = FALSE) +
  geom_point(aes(x = PCA1, y = PCA2, color = group)) +
  scale_color_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'),
        panel.grid = element_blank(), legend.key = element_blank(),
        axis.text = element_text(color = 'black'), axis.ticks = element_line(color = 'black')) +
  labs(x = paste("PC1 ( ",round(b_PC[1], digits = 2),"%"," )", sep = ""), 
       y = paste("PC2 ( ",round(b_PC[2], digits = 2),"%"," )", sep = ""))

# 或者在质心点位置添加样本分组标签
p1 <- p +
  geom_point(data = stat, aes(x = PCA1.mean, y = PCA2.mean, color = group), shape = 22, fill = 'white', size = 10, show.legend = FALSE) +
  geom_text(data = stat, aes(x = PCA1.mean, y = PCA2.mean, color = group, label = group), size = 4, show.legend = FALSE) +
  stat_ellipse(aes(x = PCA1, y = PCA2, color = group),linetype = 2, level = 0.95, show.legend = FALSE) +
  scale_color_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  geom_vline(xintercept = 0, lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = 0, lty=4,col="black",lwd=0.5) +
  theme(legend.position = 'top',
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        axis.text = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10, color = "black"))

# 添加组间ANOSIM检验的结果，以表格展示
#anosim_all <- data.frame("Treatment"=c("Con Vs. AB","Con Vs. FMT","AB Vs. FMT","FMT Vs.DC"),
#                         "R"=c(0.892, 0.548, 0.012,0.26),
#                         "Pvalue"=c(0.015,0.006,0.314,0.07),stringsAsFactors = FALSE)

#anosim_all_p <- ggtexttable(anosim_all,rows = NULL, theme=ttheme("minimal", 
#                                                                 base_size = 12))

# 图形拼接
#pp <- p1 + annotation_custom(ggplotGrob(anosim_all_p),xmin = 3, ymin = 2.2, xmax = 5)

ggsave("stool_metag/stool_level3_PCA_plot.png", plot = p1, width = 6, height = 6, dpi = 600)








