rm(list = ls())
options(stringsAsFactors = FALSE)

# stool
alpha_df <- read.csv("saliva_metag/Unigenes.absolute.s.csv",header = T,sep = ',', row.names = 1)
alpha_df <- as.data.frame(t(alpha_df))

# 如果输入数据不是整数的话，estimateR函数会报错
alpha_df2 <- sapply(alpha_df, as.integer)
alpha_df2 <- as.data.frame(alpha_df2)
rownames(alpha_df2) <- rownames(alpha_df)  

library(ggplot2)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(vegan)
library(picante)
library("RColorBrewer")
display.brewer.all()

# calculation for alpha diversity based on diversity function in R package vegan
shannon_index <- diversity(alpha_df2, index = "shannon", base = exp(1))
shannon_diversity <- exp(1)^shannon_index


#定义函数:https://mp.weixin.qq.com/s/1HEEoWQi8mb1FotZep_Z7Q
library(picante)       #picante 包加载时默认同时加载 vegan

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

#现在直接使用定义好的命令 alpha()，一步得到多种 Alpha 多样性指数
#加载 OTU 丰度表和进化树文件
#otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#otu <- t(otu)
#tree <- read.tree('otu_tree.tre')

#不包含谱系多样性，无需指定进化树；Shannon 公式的 log 底数我们使用 2
alpha_all <- alpha(alpha_df2 + 0.00001, base = 2)
#包含谱系多样性时，指定进化树文件；Shannon 公式的 log 底数我们使用 2
#alpha_all <- alpha(otu, tree, base = 2)
write.csv(alpha_all, "Output_data/Saliva_alpha_stat.csv")

#alpha_df <- alpha_df[-8,]
#httr::set_config( httr::config( ssl_verifypeer = 0L ) )
#violin plot
alpha_all$group <- "NA"
#rownames(alpha_all) <- as.character(rownames(alpha_all))
alpha_all$group[grep("^Ad",rownames(alpha_all))] <- "Ad"
alpha_all$group[grep("^Ma",rownames(alpha_all))] <- "Ma"
alpha_all$group[grep("^Oa",rownames(alpha_all))] <- "Oa"

# do alpha analysis without group Severe
#alpha_all <- alpha_all[5:nrow(alpha_all), ]

my_comparisons_a <- list(c("Ad","Ma"), c("Ad", "Oa"), c("Ma", "Oa"))
#my_comparisons_a <- list(c("Health", "Mild"))
#stat_t <- t_test(data = alpha_all, Chao1~group)
#stat_t <- add_significance(stat_t, 'p')
#stat_t.test <-  add_xy_position(stat_t, x = 'group', dodge = 0.8)
p1 <- ggplot(alpha_all, aes(x=group,y=Chao1,fill=group)) +
  geom_boxplot(position = position_dodge(1))
#Change the order of items in the legend
p1 <- p1 + theme_classic()
p1 <- p1 +  scale_x_discrete(limits=c("Ad", "Ma", "Oa")) +
  stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "top")
p1 <-  p1 + scale_color_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  scale_fill_manual(values = c("#00B9C0", "#CA8BDB", "#C29C5E"))
p1 <- p1 + geom_jitter(position=position_jitter(0.2))
p1




##Shannon index
p2 <- ggplot(alpha_all, aes(x=group,y=Shannon,fill=group)) +
  geom_boxplot(position = position_dodge(1))
#Change the order of items in the legend
p2 <- p2 + theme_classic()
p2 <- p2 +  scale_x_discrete(limits=c("Ad", "Ma", "Oa"))+
  stat_compare_means(comparisons = my_comparisons_a, method = "t.test",label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "top")
p2 <- p2 + scale_color_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  scale_fill_manual(values= c("#00B9C0", "#CA8BDB", "#C29C5E"))
p2 <- p2 + geom_jitter(position=position_jitter(0.2))
p2
#Arrange on one page
p_total <- ggarrange(p1,p2,labels = "A")
ggsave("Plots/Saliva_alpha_compare20221107.png", plot = p_total, width = 6, height = 6, dpi = 300)



## 将feces和saliva的alpha多样性结果放在一个表里面，通过分面的方式一起展示两个部位的
a <- read.csv("Output_data/Feces_Saliva_alpha_stat_merge.csv", header = T, row.names = 1,
              check.names = F)
#dim(a)
a$group <- "NA"
#rownames(alpha_all) <- as.character(rownames(alpha_all))
a$group[grep("^Ad",rownames(a))] <- "Ad"
a$group[grep("^Ma",rownames(a))] <- "Ma"
a$group[grep("^Oa",rownames(a))] <- "Oa"

# do alpha analysis without group Severe
#alpha_all <- alpha_all[5:nrow(alpha_all), ]

my_comparisons_a <- list(c("Ad","Ma"), c("Ad", "Oa"), c("Ma", "Oa"))
#my_comparisons_a <- list(c("Health", "Mild"))
#stat_t <- t_test(data = alpha_all, Chao1~group)
#stat_t <- add_significance(stat_t, 'p')
#stat_t.test <-  add_xy_position(stat_t, x = 'group', dodge = 0.8)
p1 <- ggplot(a, aes(x = group,y = Richness, color = group)) +
  geom_rect(data = a, aes(xmin = 0,xmax = 4,
                          ymin = -Inf,ymax = Inf, fill = Bodysite), alpha=0.03) +
  geom_boxplot(position = position_dodge(1), color = "black") +
  facet_grid(~Bodysite) +

#Change the order of items in the legend
  theme_classic() +
  scale_x_discrete(limits=c("Ad", "Ma", "Oa")) +
  stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none",
        strip.background = element_rect(fill=c("#C2C4C5")), # 分面的背景色设置
        strip.text = element_text(size = 15,colour = "black", face = "bold")) +  # 分面的字体设置
  scale_color_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  scale_fill_manual(values = c("#E9F5FA", "#E2EFE1")) + 
  geom_jitter(position=position_jitter(0.2))
p1




##Shannon index
p2 <- ggplot(a, aes(x = group,y = Shannon, color = group)) +
  geom_rect(data = a, aes(xmin = 0,xmax = 4,
                          ymin = -Inf,ymax = Inf, fill = Bodysite), alpha=0.03) +
  geom_boxplot(position = position_dodge(1), color = "black") +
  facet_grid(~Bodysite) +
  
  #Change the order of items in the legend
  theme_classic() +
  scale_x_discrete(limits=c("Ad", "Ma", "Oa")) +
  stat_compare_means(comparisons = my_comparisons_a,method = "t.test", label = "p")+
  #stat_compare_means(label.y = 3500)+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.position = "none",
        strip.background = element_rect(fill=c("#C2C4C5")), # 分面的背景色设置
        strip.text = element_text(size = 15,colour = "black", face = "bold")) +  # 分面的字体设置
  scale_color_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E")) +
  scale_fill_manual(values = c("#E9F5FA", "#E2EFE1")) + 
  geom_jitter(position=position_jitter(0.2))
p2
#Arrange on one page
p_total <- ggarrange(p1,p2,labels = c("A","B"))
ggsave("Plots/Saliva_Feces_alpha_compare20221109-2.png", plot = p_total, width = 10, height = 6, dpi = 300)


