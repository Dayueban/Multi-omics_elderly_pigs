# (1) dataset cleaning with dataset at species
rm(list = ls())
library(readxl)
library(ggplot2)
library(stringr)
library(dplyr)
library(stringr)

#######################################################################
#######################################################################
#######################################################################
# (1) 读取数据, 种水平的相对丰度表
a <- read.csv("stool_metag/Feces_Unigenes.absolute.s.csv", header = T, row.names = 1, check.names = F)

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





# (2) kruskal  t-test (读取完数据后直接到这一步)
nr <- nrow(a2)
Pvalue_kw <- c(rep(0,nr))
#log2_FC <- c(rep(0,nr))
#FC <- c(rep(0,nr))
group_list <- rep(c("Ad", "Ma", "Oa"), c(15,15,15))  
for (i in 1:nr){
  #Pvalue_NC_HC[i] <- wilcox.test(as.numeric(a[i,grep("^NC", names(a))]), as.numeric(a[i,grep("^HC", names(a))]))$p.value
  #Pvalue_NC_HE[i] <- wilcox.test(as.numeric(a[i,grep("^NC", names(a))]), as.numeric(a[i,grep("^HE", names(a))]))$p.value
  #Pvalue_HC_HE[i] <- wilcox.test(as.numeric(a[i,grep("^HC", names(a))]), as.numeric(a[i,grep("^HE", names(a))]))$p.value
  #Pvalue[i] <- wilcox.test(as.numeric(a3[i,]) ~ group_list, data = a3)$p.value
  Pvalue_kw[i] <- kruskal.test(as.numeric(a2[i,]) ~ group_list, data = a)$p.value
  #FC[i] <- (mean(as.numeric(a3[i,1:6]))+0.001)/ (mean(as.numeric(a3[i,7:13]))+0.001)
  #log2_FC[i] <- log2((mean(as.numeric(a3[i,1:6]))+0.001)/ (mean(as.numeric(a3[i,7:13]))+0.001))
}

#Pvalue <- as.numeric(Pvalue)  #在操作的过程中发现秩和检验的包输出来不是列表的形式，这里变成列表
#fdr.NC_HC <- p.adjust(Pvalue_NC_HC,method="fdr",length(Pvalue_NC_HC))    #p.adjust就是计算FDR的包，这个可要记得了
#fdr.NC_HE <- p.adjust(Pvalue_NC_HE,method="fdr",length(Pvalue_NC_HE))
fdr.kw <- p.adjust(Pvalue_kw, method="BH", length(Pvalue_kw))
#a3_res <- cbind(a3,Pvalue,fdr.w, FC,log2_FC)
a3 <- cbind(a2, Pvalue_kw, fdr.kw)
write.csv(a3, file="Output_data/Feces_species_tidy_kruskal_test20230228.csv")

# 1) 选择差异最大的前30个进行展示
a4 <- a3[order(a3$fdr.kw, decreasing = FALSE), ] # 根据fdr值进行排序
Taxo_df1 <- a4[1:30, 1:45]
rownames(Taxo_df1) <- str_split(rownames(Taxo_df1), ";", simplify = TRUE)[,2] # 去掉属信息，只保留种
Taxo_df2 <- Taxo_df1
Taxo_df3 <- as.data.frame(t(Taxo_df2))

library(reshape2)
Taxo_df4 <- cbind(rownames(Taxo_df3), Taxo_df3)
names(Taxo_df4)[1] <- "ID"
## reshape dataframe, from long sheet to wild sheet
Taxo_df5 <- melt(Taxo_df4,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")

# add a new cloumn as the group information
Taxo_df5$group <- NA
Taxo_df5$group[grep("^Ad", Taxo_df5$ID)] <- "Ad"
Taxo_df5$group[grep("^Ma", Taxo_df5$ID)] <- "Ma"
Taxo_df5$group[grep("^Oa", Taxo_df5$ID)] <- "Oa"
Taxo_df5$group <- factor(Taxo_df5$group, levels = c("Ad","Ma","Oa"))

# 计算mean,sd,se的方法，感觉更简便，用的是dplyr包里面的函数，加上管道操作更直观
library(dplyr)
Taxo_df55 <- Taxo_df5 %>%
  group_by(Taxa, group) %>% # 根据两个指标进行分组
  summarise(value_mean=mean(Relative_abundance),sd=sd(Relative_abundance),
            se=sd(Relative_abundance)/sqrt(n())) # summarise函数进行计算


# draw barplot a
library(ggplot2)
library(ggbreak)
#par(mar=c(2,2,2,2))
p <- ggplot(data=Taxo_df55,aes(x = Taxa,y=value_mean,fill=group)) + 
  geom_bar(stat = "identity", color="black",position = position_dodge()) + 
  geom_errorbar(aes(ymin=value_mean - se, ymax = value_mean + se), 
                width=0.2,position = position_dodge(0.9)) + 
  theme_classic() + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=90,color="black",size=12,vjust = 0.5,hjust=0.95),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.y=element_text(size = 16,color = "black"),
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
p <- p + scale_y_continuous(trans = "sqrt", expand = c(0,0))
##p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Feces_species_tidy_top30_kw20230228.png",plot = p, width = 8, height = 12, dpi = 300)



###################################################################
###################################################################
###################################################################
# (3) wilcoxon t-test
nr <- nrow(a2)
Pvalue_ad_ma <- c(rep(0,nr))
Pvalue_ad_oa <- c(rep(0,nr))
Pvalue_ma_oa <- c(rep(0,nr))
#log2_FC <- c(rep(0,nr))
#FC <- c(rep(0,nr))
#group_list <- rep(c("Ctrl","RSV"), c(6,7))  
for (i in 1:nr){
  Pvalue_ad_ma[i] <- wilcox.test(as.numeric(a2[i,1:15]), as.numeric(a2[i,16:30]))$p.value
  Pvalue_ad_oa[i] <- wilcox.test(as.numeric(a2[i,1:15]), as.numeric(a2[i,31:45]))$p.value
  Pvalue_ma_oa[i] <- wilcox.test(as.numeric(a2[i,16:30]), as.numeric(a2[i,31:45]))$p.value
  #Pvalue[i] <- wilcox.test(as.numeric(a3[i,]) ~ group_list, data = a3)$p.value
  #Pvalue[i] <- kruskal.test(as.numeric(a[i,]) ~ group_list, data = a)$p.value
  #FC[i] <- (mean(as.numeric(a3[i,1:6]))+0.001)/ (mean(as.numeric(a3[i,7:13]))+0.001)
  #log2_FC[i] <- log2((mean(as.numeric(a3[i,1:6]))+0.001)/ (mean(as.numeric(a3[i,7:13]))+0.001))
}

#Pvalue <- as.numeric(Pvalue)  #在操作的过程中发现秩和检验的包输出来不是列表的形式，这里变成列表
fdr.ad_ma <- p.adjust(Pvalue_ad_ma,method="BH",length(Pvalue_ad_ma))    #p.adjust就是计算FDR的包，这个可要记得了
fdr.ad_oa <- p.adjust(Pvalue_ad_oa,method="BH",length(Pvalue_ad_oa))
fdr.ma_oa <- p.adjust(Pvalue_ma_oa,method="BH",length(Pvalue_ma_oa))
#a3_res <- cbind(a3,Pvalue,fdr.w, FC,log2_FC)
a3 <- cbind(a2,Pvalue_ad_ma,Pvalue_ad_oa, Pvalue_ma_oa, fdr.ad_ma,
            fdr.ad_oa, fdr.ma_oa)
write.csv(a3, file="Output_data/Feces_species_tidy_wilcoxn_test_BH.csv")

# plot01 Ad and Ma
a4 <- a3[order(a3$fdr.ad_ma, decreasing = FALSE), ] # Ad vs Ma
Taxo_df1 <- a4[1:30, 1:30] # 选择差异最显著的前20个进行展示
rownames(Taxo_df1) <- str_split(rownames(Taxo_df1), ";", simplify = TRUE)[,2] # 去掉属信息，只保留种
Taxo_df2 <- Taxo_df1
Taxo_df3 <- as.data.frame(t(Taxo_df2))

library(reshape2)
Taxo_df4 <- cbind(rownames(Taxo_df3), Taxo_df3)
names(Taxo_df4)[1] <- "ID"
## reshape dataframe, from long sheet to wild sheet
Taxo_df5 <- melt(Taxo_df4,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")

# add a new cloumn as the group information
Taxo_df5$group <- NA
Taxo_df5$group[grep("^Ad", Taxo_df5$ID)] <- "Ad"
Taxo_df5$group[grep("^Ma", Taxo_df5$ID)] <- "Ma"
#Taxo_df5$group[grep("^Oa", Taxo_df5$ID)] <- "Oa"
Taxo_df5$group <- factor(Taxo_df5$group, levels = c("Ad","Ma"))

# 计算mean,sd,se的方法，感觉更简便，用的是dplyr包里面的函数，加上管道操作更直观
library(dplyr)
Taxo_df55 <- Taxo_df5 %>%
  group_by(Taxa, group) %>% # 根据两个指标进行分组
  summarise(value_mean=mean(Relative_abundance),sd=sd(Relative_abundance),
            se=sd(Relative_abundance)/sqrt(n())) # summarise函数进行计算


# draw barplot a
library(ggplot2)
library(ggbreak)
#par(mar=c(2,2,2,2))
p <- ggplot(data=Taxo_df55,aes(x = Taxa,y=value_mean,fill=group)) + 
  geom_bar(stat = "identity", color="black",position = position_dodge()) + 
  geom_errorbar(aes(ymin=value_mean - se, ymax = value_mean + se), 
                width=0.2,position = position_dodge(0.9)) + 
  theme_classic() + 
  #coord_flip() + 
  theme(axis.text.x=element_text(angle=60,color="black",size=10,vjust = 0.95,hjust=0.95),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.y=element_text(size = 16,color = "black"),
        legend.position = "top",
        plot.title = element_text(size = 20)) + 
  guides(fill=guide_legend(title=NULL)) + 
  xlab("") + 
  labs(y = "Relative abundance (sqrt)") +
  #scale_fill_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E"),
  #                  labels = c("Ad", "Ma","Oa")) + 
  scale_fill_manual(values=c("#00B9C0", "#CA8BDB"),
                    labels = c("Ad", "Ma")) +
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_y_break(breaks = c(0.0001, 0.0002)) + 
  theme(legend.key.height=unit(0.6,"cm"),
        legend.key.width=unit(0.6,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=20),
        plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "sqrt", expand = c(0,0))
##p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Feces_species_tidy_Ad_Ma_top30_wilcoxn20230228.png",plot = p, width = 12, height = 8, dpi = 300)

# plot02 Ad and Oa
a4 <- a3[order(a3$fdr.ad_oa, decreasing = FALSE), ] # Ad vs Oa
Taxo_df1 <- a4[1:30, c(1:15, 31:45)] # 选择差异最显著的前20个进行展示
rownames(Taxo_df1) <- str_split(rownames(Taxo_df1), ";", simplify = TRUE)[,2] # 去掉属信息，只保留种
Taxo_df2 <- Taxo_df1
Taxo_df3 <- as.data.frame(t(Taxo_df2))

library(reshape2)
Taxo_df4 <- cbind(rownames(Taxo_df3), Taxo_df3)
names(Taxo_df4)[1] <- "ID"
## reshape dataframe, from long sheet to wild sheet
Taxo_df5 <- melt(Taxo_df4,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")

# add a new cloumn as the group information
Taxo_df5$group <- NA
Taxo_df5$group[grep("^Ad", Taxo_df5$ID)] <- "Ad"
#Taxo_df5$group[grep("^Ma", Taxo_df5$ID)] <- "Ma"
Taxo_df5$group[grep("^Oa", Taxo_df5$ID)] <- "Oa"
Taxo_df5$group <- factor(Taxo_df5$group, levels = c("Ad","Oa"))

# 计算mean,sd,se的方法，感觉更简便，用的是dplyr包里面的函数，加上管道操作更直观
library(dplyr)
Taxo_df55 <- Taxo_df5 %>%
  group_by(Taxa, group) %>% # 根据两个指标进行分组
  summarise(value_mean=mean(Relative_abundance),sd=sd(Relative_abundance),
            se=sd(Relative_abundance)/sqrt(n())) # summarise函数进行计算


# draw barplot a
library(ggplot2)
library(ggbreak)
#par(mar=c(2,2,2,2))
p <- ggplot(data=Taxo_df55,aes(x = Taxa,y=value_mean,fill=group)) + 
  geom_bar(stat = "identity", color="black",position = position_dodge()) + 
  geom_errorbar(aes(ymin=value_mean - se, ymax = value_mean + se), 
                width=0.2,position = position_dodge(0.9)) + 
  theme_classic() + 
  #coord_flip() + 
  theme(axis.text.x=element_text(angle=60,color="black",size=10,vjust = 0.95,hjust=0.95),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.y=element_text(size = 16,color = "black"),
        legend.position = "top",
        plot.title = element_text(size = 20)) + 
  guides(fill=guide_legend(title=NULL)) + 
  xlab("") + 
  labs(y = "Relative abundance (sqrt)") +
  #scale_fill_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E"),
  #                  labels = c("Ad", "Ma","Oa")) + 
  scale_fill_manual(values=c("#00B9C0", "#C29C5E"),
                    labels = c("Ad", "Oa")) +
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_y_break(breaks = c(0.0001, 0.0002)) + 
  theme(legend.key.height=unit(0.6,"cm"),
        legend.key.width=unit(0.6,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=20),
        plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "sqrt", expand = c(0, 0))
##p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Feces_species_tidy_Ad_Oa_top30_wilcoxn20230228.png",plot = p, width = 12, height = 8, dpi = 300)

# plot03 Ma and Oa
a4 <- a3[order(a3$fdr.ma_oa, decreasing = FALSE), ] # Ma vs Oa
Taxo_df1 <- a4[1:30, 16:45] # 选择差异最显著的前20个进行展示
rownames(Taxo_df1) <- str_split(rownames(Taxo_df1), ";", simplify = TRUE)[,2] # 去掉属信息，只保留种
Taxo_df2 <- Taxo_df1
Taxo_df3 <- as.data.frame(t(Taxo_df2))

library(reshape2)
Taxo_df4 <- cbind(rownames(Taxo_df3), Taxo_df3)
names(Taxo_df4)[1] <- "ID"
## reshape dataframe, from long sheet to wild sheet
Taxo_df5 <- melt(Taxo_df4,id.vars = "ID",variable.name = "Taxa", value.name = "Relative_abundance")

# add a new cloumn as the group information
Taxo_df5$group <- NA
#Taxo_df5$group[grep("^Ad", Taxo_df5$ID)] <- "Ad"
Taxo_df5$group[grep("^Ma", Taxo_df5$ID)] <- "Ma"
Taxo_df5$group[grep("^Oa", Taxo_df5$ID)] <- "Oa"
Taxo_df5$group <- factor(Taxo_df5$group, levels = c("Ma","Oa"))

# 计算mean,sd,se的方法，感觉更简便，用的是dplyr包里面的函数，加上管道操作更直观
library(dplyr)
Taxo_df55 <- Taxo_df5 %>%
  group_by(Taxa, group) %>% # 根据两个指标进行分组
  summarise(value_mean=mean(Relative_abundance),sd=sd(Relative_abundance),
            se=sd(Relative_abundance)/sqrt(n())) # summarise函数进行计算


# draw barplot a
library(ggplot2)
library(ggbreak)
#par(mar=c(2,2,2,2))
p <- ggplot(data=Taxo_df55,aes(x = Taxa,y=value_mean,fill=group)) + 
  geom_bar(stat = "identity", color="black",position = position_dodge()) + 
  geom_errorbar(aes(ymin=value_mean - se, ymax = value_mean + se), 
                width=0.2,position = position_dodge(0.9)) + 
  theme_classic() + 
  #coord_flip() + 
  theme(axis.text.x=element_text(angle=60,color="black",size=10,vjust = 0.95,hjust=0.95),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.y=element_text(size = 16,color = "black"),
        legend.position = "top",
        plot.title = element_text(size = 20)) + 
  guides(fill=guide_legend(title=NULL)) + 
  xlab("") + 
  labs(y = "Relative abundance (sqrt)") +
  #scale_fill_manual(values=c("#00B9C0", "#CA8BDB", "#C29C5E"),
  #                  labels = c("Ad", "Ma","Oa")) + 
  scale_fill_manual(values=c("#CA8BDB", "#C29C5E"),
                    labels = c("Ma", "Oa")) +
  #scale_y_continuous(expand = c(0, 0)) +
  #scale_y_break(breaks = c(0.0001, 0.0002)) + 
  theme(legend.key.height=unit(0.6,"cm"),
        legend.key.width=unit(0.6,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=12),
        legend.title=element_text(size=20),
        plot.margin = unit(c(1,0.5,1,0.5),"cm"))
p <- p + scale_y_continuous(trans = "sqrt", expand = c(0, 0))
##p <- p + geom_jitter() ##添加数据点
ggsave("Plots/Feces_species_tidy_Ma_Oa_top30_wilcoxn20230228.png",plot = p, width = 12, height = 8, dpi = 300)

