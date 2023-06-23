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
library(stringr)
library(qpcR)
library(VennDiagram)

# feces和saliva样本KO genes在三组间差异的结果取交集发现有14个是重合的
# 所以绘制下venn图
a <- read.csv("Output_data/Feces_level3_tidy_kruskal_test.csv",
              row.names = 1, header = T, check.names = F)
b <- read.csv("Output_data/Saliva_level3_tidy_kruskal_test20230303.csv",
              row.names = 1, header = T, check.names = F)

a1 <- rownames(a)[which(a$fdr.kw < 0.05)]
b1 <- rownames(b)[which(b$fdr.kw < 0.05)]

# 绘制韦恩图
#求四个向量的交集
keep_KO <- Reduce(intersect,  list(v1 = a1,
                                   v2 = b1))
#[1] "NEG691"  "POS2085" "POS248"
par(mar=c(1,1,1,1))
venn.diagram(
  x = list(
    a1, b1
  ),
  category.names = c("Feces" , "Saliva"),
  filename = "Plots/venn_level3_intersect_feces_saliva.png",
  output = TRUE ,
  imagetype="png" ,
  height = 1200 , 
  width = 1200 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#e20612","#00b0eb"),
  fill = c(alpha("#e20612",0.8), alpha("#00b0eb",0.8)),
  cex = 1.2,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 0.8,
  #cat.default.pos = "inter",
  cat.pos = c(-30, 30),
  margin = 0.1,
  cat.dist = c(0.05, 0.05),
  cat.fontfamily = "sans",
  #cat.col = c("#7292AA","#90CBD3","#EAA58E"), # 标签的颜色
  cat.col = c("black", "black"),
  cat.fontface = "bold",
  #rotation = 1 #绘制两个韦恩图的时候该参数不能用，否则报错
  #rotation = 1
)


# 分别对feces和saliva交集的26个Ko map取出来作柱状图
feces_ko13 <- a[keep_KO, 1:45]
feces_ko13 <- as.data.frame(t(feces_ko13))
feces_ko13 <- sweep(feces_ko13, 1, rowSums(feces_ko13), "/")
saliva_ko13 <- b[keep_KO, 1:45]
saliva_ko13 <- as.data.frame(t(saliva_ko13))
saliva_ko13 <- sweep(saliva_ko13, 1, rowSums(saliva_ko13), "/")


# feces ploting
library(reshape2)
feces_ko13.2 <- cbind(rownames(feces_ko13), feces_ko13)
names(feces_ko13.2)[1] <- "ID"
## reshape dataframe, from long sheet to wild sheet
feces_ko13.3 <- melt(feces_ko13.2,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")

# add a new cloumn as the group information
feces_ko13.3$group <- NA
feces_ko13.3$group[grep("^Ad", feces_ko13.3$ID)] <- "Ad"
feces_ko13.3$group[grep("^Ma", feces_ko13.3$ID)] <- "Ma"
feces_ko13.3$group[grep("^Oa", feces_ko13.3$ID)] <- "Oa"
feces_ko13.3$group <- factor(feces_ko13.3$group, levels = c("Ad","Ma","Oa"))

# 计算mean,sd,se的方法，感觉更简便，用的是dplyr包里面的函数，加上管道操作更直观
library(dplyr)
library(plyr)
feces_ko13.35 <- feces_ko13.3 %>%
  group_by(Taxa, group) %>% # 根据两个指标进行分组
  dplyr::summarise(value_mean=mean(Relative_abundance),sd=sd(Relative_abundance),
                   se=sd(Relative_abundance)/sqrt(n())) # summarise函数进行计算

# draw barplot a
library(ggplot2)
library(ggbreak)
#par(mar=c(2,2,2,2))
p <- ggplot(data=feces_ko13.35,aes(x = Taxa,y=value_mean,fill=group)) + 
  geom_bar(stat = "identity", color="black",position = position_dodge()) + 
  geom_errorbar(aes(ymin=value_mean - se, ymax = value_mean + se), 
                width=0.2,position = position_dodge(0.9)) + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=90,color="black",size=8, vjust = 0.95, hjust = 0.95),
        axis.text.y=element_text(size = 10, colour = "black"),
        axis.title.y=element_text(size = 14,color = "black"),
        axis.ticks.y = element_line(colour = 'transparent'),
        axis.line.y = element_blank(),
        legend.position = "top",
        plot.title = element_text(size = 20)) + 
  guides(fill=guide_legend(title=NULL)) + 
  xlab("") + 
  labs(y = "Relative abundance") +
  scale_fill_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E"),
                    labels = c("Ad", "Ma","Oa")) + 
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_y_break(breaks = c(0.0001, 0.0002)) + 
  theme(legend.key.height=unit(0.6,"cm"),
        legend.key.width=unit(0.6,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=20),
        plot.margin = unit(c(1,0.5,1,0.5),"cm"))
#p <- p + scale_y_continuous(trans = "sqrt", expand = c(0,0))
p <- p + scale_y_continuous(trans = "reverse", expand = c(0,0)) ## reverse the 'Percent' axis using trans = "reverse"
##p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Feces_Saliva_intersect_in_feces_level3.png",plot = p, width = 4, height = 6, dpi = 300)


# saliva ploting
library(reshape2)
saliva_ko13.2 <- cbind(rownames(saliva_ko13), saliva_ko13)
names(saliva_ko13.2)[1] <- "ID"
## reshape dataframe, from long sheet to wild sheet
saliva_ko13.3 <- melt(saliva_ko13.2,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")

# add a new cloumn as the group information
saliva_ko13.3$group <- NA
saliva_ko13.3$group[grep("^Ad", saliva_ko13.3$ID)] <- "Ad"
saliva_ko13.3$group[grep("^Ma", saliva_ko13.3$ID)] <- "Ma"
saliva_ko13.3$group[grep("^Oa", saliva_ko13.3$ID)] <- "Oa"
saliva_ko13.3$group <- factor(saliva_ko13.3$group, levels = c("Ad","Ma","Oa"))

# 计算mean,sd,se的方法，感觉更简便，用的是dplyr包里面的函数，加上管道操作更直观
library(dplyr)
library(plyr)
saliva_ko13.35 <- saliva_ko13.3 %>%
  group_by(Taxa, group) %>% # 根据两个指标进行分组
  dplyr::summarise(value_mean=mean(Relative_abundance),sd=sd(Relative_abundance),
                   se=sd(Relative_abundance)/sqrt(n())) # summarise函数进行计算

# draw barplot a
library(ggplot2)
library(ggbreak)
#par(mar=c(2,2,2,2))
p <- ggplot(data=saliva_ko13.35,aes(x = Taxa,y=value_mean,fill=group)) + 
  geom_bar(stat = "identity", color="black",position = position_dodge()) + 
  geom_errorbar(aes(ymin=value_mean - se, ymax = value_mean + se), 
                width=0.2,position = position_dodge(0.9)) + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=90,color="black",size=8, vjust = 0.95, hjust = 0.95),
        axis.text.y=element_text(size = 10, colour = "black"),
        axis.title.y=element_text(size = 14,color = "black"),
        axis.ticks.y = element_line(colour = 'transparent'),
        axis.line.y = element_blank(),
        legend.position = "top",
        plot.title = element_text(size = 20)) + 
  guides(fill=guide_legend(title=NULL)) + 
  xlab("") + 
  labs(y = "Relative abundance (sqrt)") +
  scale_fill_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E"),
                    labels = c("Ad", "Ma","Oa")) + 
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_y_break(breaks = c(0.0001, 0.0002)) + 
  theme(legend.key.height=unit(0.6,"cm"),
        legend.key.width=unit(0.6,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=20),
        plot.margin = unit(c(1,0.5,1,0.5),"cm"))
#p <- p + scale_y_continuous(trans = "sqrt", expand = c(0,0))
p <- p + scale_y_continuous(expand = c(0,0))
##p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Feces_Saliva_intersect_in_saliva_level3.png",plot = p, width = 4, height = 6, dpi = 300)




