rm(list = ls())
#temp_nm2 <- "filtered"
#library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (WGCNA) ### -	Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations
library(stringr) ### Mainly Dataframe/Character manipulation
library(tidyverse)
library(circlize)

# load working datasheet
options (stringsAsFactors = FALSE)
#work_df <- read.csv("Ileum_OTU.csv",header = T,row.names = 1,check.names = F)
a <- read.csv("Output_data/Feces_species_tidy_wilcoxn_test_BH.csv", header = T,
              row.names = 1, check.names = F)
a0.05 <- a[which(a$fdr.ad_ma < 0.05), 1:30]
a0.05t <- as.data.frame(t(a0.05))

dat <- a0.05t


# 计算模块后的代谢物与16S_ASV之间的相关性以及作图
metabo <- read.csv("Serum_metabo/Differential_Ad_Ma.csv", header = T,
                    row.names = 1, check.names = F)
metabo_t <- as.data.frame(t(metabo))

# 要保证两组数据的个体ID是一样的
dat2 <- dat[c(1:10, 18:22, 24:27, 29),]
library(plyr)
library(psych)
library(reshape2)
corr_df <- corr.test(dat2, metabo_t, method = "spearman", adjust = "BH", alpha = .05)
corr_df_cor <- corr_df$r
corr_df_p <- corr_df$p.adj

#corr_df_cor <- corr_df_cor[abs(corr_df_cor) > 0.6 | corr_df_p < 0.05 ]

#' jiont the r and p matrix
asso_df <- rbind(corr_df_cor,corr_df_p)
#' transpose the data frame
asso_df <- as.data.frame(t(asso_df))
# rename the column variables
colnames(asso_df) <- paste(colnames(asso_df),rep(c("cor","P_adj"),c(62,62)),sep = "_")
#' re-order the asso_df according to the adjusted p value
#asso_df <- asso_df[order(asso_df$P_adj,decreasing = FALSE),]


## use ggplot2
# Reset rownames
corr_df_cor <- data.frame(row=rownames(corr_df_cor),corr_df_cor,check.names = F) # create a column called "row" 
rownames(corr_df_cor) <- NULL
corr_df_p <- data.frame(row=rownames(corr_df_p),corr_df_p,check.names = F) # create a column called "row" 
rownames(corr_df_p) <- NULL
# Melt
nbar.m <- melt(corr_df_cor)
nbap.m <- melt(corr_df_p)
# Classify (you can classify differently for nbar and for nbap also)         
nbar.m$value2<-cut(nbar.m$value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend
nbap.m$value2<-cut(nbap.m$value,breaks=c(-Inf, 0.001, 0.01, 0.05),label=c("***", "**", "*")) 
nbar.m<-cbind.data.frame(nbar.m,nbap.m$value,nbap.m$value2) # adding the p value and its cut to the first dataset of R coefficients
names(nbar.m)[5]<-paste("valuep") # change the column names of the dataframe 
names(nbar.m)[6]<-paste("signif.")


#' part one: chord diagram plot
# 添加一列数值1，分别计算相同代谢物的个数，用于后续排序
df3 <- nbar.m
#write.csv(df3, file = "metabo_assoc_16S.csv") # 用作文章的Supplementary Table
# remove column which the |correlation coefficients| less than 0.3
df4 <- df3[-which(abs(df3$value) <= 0.7) ,]
df4 <- df4[which(df4$valuep < 0.01), ]
df4 <- df4[,1:3]
names(df4) <- c("Species", "metabo", "Coeff")
df4$value <- rep(1, nrow(df4))
# 根据列名taxo将相同细菌凑一起，然后根据value计算个数

data_total <- df4 %>%
  group_by(Species) %>% 
  transmute(Total=sum(value))
# 根据列名metabo将相同代谢物凑一起，然后根据value计算个数
data_total2 <- df4 %>%
  group_by(metabo) %>% 
  transmute(Total2=sum(value))
# 将data_total中的Total列合并到df4数据框中
df4 <- cbind(df4, data_total$Total, data_total2$Total2)
names(df4)[c(ncol(df4)-1, ncol(df4))] <- c("Total", "Total2")

# change the first column by paste metabo names with numbers in column "Total"
df4$Species <- paste(df4$Species, "\t(",df4$Total, ")",sep = "")
df4$metabo <- paste(df4$metabo, "\t(", df4$Total2,")", sep = "")
# sort df4 based on the numbers of links
df4 <- df4[order(df4$Total2, decreasing = T),]
# remove column 'value' and 'Total'
df4$value <- NULL
df4$Total <- NULL
df4$Total2 <- NULL
## set colours for segments
#df4 <- df4[-which(abs(df4$Corr_score) < 5), ] # only |Corr_score| > 3 were used for exhibition
circos.clear()
# 第一种图是将links根据正负相关显示的
#pdf(file = paste(temp_nm2,"_otus_metabo_circlize.pdf",sep = ""), height = 10, width = 10)
tiff(file = "Faeces_microbe_serum_metabolites_circlize202321.tiff", height = 2400, width = 2400,
     res = 300)
#circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
#Put horizontally or vertically symmetric
# gap.after是为了将Otu和代谢物之间分一个gap
circos.par(start.degree = 0, gap.degree = 1, points.overflow.warning = FALSE,
           gap.after = c(rep(1.5, length(unique(df4$Species))-1), 5, 
                         rep(1.5, length(unique(df4$metabo))-1), 5),
           circle.margin = c(0.2,0.3,0.2,0.2)) #controls the margins on the left, right, bottom and top sides of the circle
#par(mar = rep(-2, 4))
#par(mar = c(-8,-8,0,0))

# 给代谢物上色，而Otu全部给予灰色
grid.col <- setNames(c(topo.colors(length(unique(df4$Species))),rep("#BEBEBE", length(unique(df4$metabo)))), 
                     c(unique(df4$Species), unique(df4$metabo)))
# now, plot the image with rotated labels
chordDiagram(df4,  
             preAllocateTracks = 1, 
             annotationTrack = "grid",
             annotationTrackHeight = c(0.03, 0.1),
             grid.col = grid.col,
             col = ifelse(df4$Coeff > 0, "#F76F72", "blue"), # link为正负相关可以用两种颜色表示
             link.sort = TRUE, link.decreasing = FALSE)
#directional = 1, 
#direction.type = c("diffHeight", "arrows")
#link.arr.type = "big.arrow")

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.4), cex = 0.6)
}, bg.border = NA)
title(paste("Metabolites assocaited with" ,"Microbial species", sep = " "))
legend("topleft", pch = 15, bty = "n",col = c("#F76F72", "blue"), title = "Connection type",
       legend = c("Positive","Negative"), cex = 1.2)
dev.off()

# 图形后期通过PS修正一下

# 另外一种作图方法
nbar.m <- nbar.m %>%
  mutate(value_abs = abs(value))
# 将物种名进行裁剪
nbar.m$row <- str_split(nbar.m$row, ";", simplify = TRUE)[,2]
# ab(value) > 0.7的才进行展示
nbar.m0.7 <- nbar.m[which(nbar.m$value_abs > 0.7 | nbar.m$valuep < 0.01), ]
# 结果储存
write.csv(nbar.m0.7,"assoc_results/Feces_microbial_assoc_metabolites.csv")
pb <- ggplot(nbar.m0.7, aes(row, variable)) +
  #geom_tile(aes(fill=value),colour="white") +
  geom_point(aes(size = value_abs, fill = value), shape = 21) +
  #scale_fill_brewer(palette = "RdYlGn",name="Correlation")# RColorBrewer package
  scale_fill_gradient2(low="#4da0a0", high="#9b3a74", guide="colorbar",name="Correlation") +
  theme_classic() +
  theme(axis.text.x=element_text(face="bold",angle=90,color="black",vjust = 0.5,hjust = 0.95,size=16),
        axis.text.y=element_text(face = "bold",size = 16, color = "black"),
        axis.title=element_text(size = 28))+
  labs(title="Metabolites associated with Microbial species") + 
  xlab("Microbial species") +
  ylab("Metabolites") +
  scale_size(range = c(6,12)) +
  theme(plot.title = element_text(size = 24,color = "black", hjust = 0.5)) +
  theme(axis.line.x=element_line(colour = "black"),axis.line.y=element_line(colour = "black"),
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(lineheight=0.8,face="bold",size=14),
        legend.title=element_text(size=18))
# Adding the significance level stars using geom_text 
pp<- pb +
  geom_text(aes(label=signif.),size=6,color = "white",na.rm=TRUE,nudge_y = -0.15)
ggsave("Microbial_species_assoc_metabo_Ma_Ad2.png", plot = pp, width = 24,height = 24, dpi = 300)

# 将1-Palmitoyllysophosphatidylcholine抽出来，和与其显著相关的菌做相关性
# 散点图
dat3 <- dat2
names(dat3) <- str_split(names(dat2), ";", simplify = TRUE)[,2]
dat3 <- dat3[, names(dat3) %in% c("s__Bacteroides fragilis",
                                  "s__Bacteroides plebeius CAG:211",
                                  "s__Enterococcus cecorum",
                                  "s__Oscillibacter sp. CAG:241",
                                  "s__Oscillibacter sp. CAG:241_62_21",
                                  "s__Firmicutes bacterium CAG:83",
                                  "s__Fusobacterium mortiferum",
                                  "s__Fusobacterium necrophorum")]
metabo_t.s <- metabo_t[, names(metabo_t) %in% "1-Palmitoyllysophosphatidylcholine"]
dddf <- cbind(dat3, metabo_t.s)
names(dddf)[ncol(dddf)] <- "Palmitoyllysophosphatidylcholine"

dddf_rel <- sweep(dddf, 2, colSums(dddf), '/')

#Min-Max Normalization
#trans_max_min <- function(x){
#  (x-min(x))/(max(x)-min(x))
#}
library(ggpubr)
#Two_variable_one[,2:5] <- as.data.frame(apply(Two_variable_one[,2:5],2,trans_max_min))
#dddf[,1:ncol(dddf)] <- as.data.frame(apply(dddf[,1:ncol(dddf)],2,trans_max_min))
myplots <- list()  # new empty list
for (i in 1:8) {
  p1 <- ggscatter(dddf_rel,x= "Palmitoyllysophosphatidylcholine",y=colnames(dddf_rel)[i],
                  color = "black",fill="black", shape = 21, size = 1.5, # Points color, shape and size
                  add.params = list(color = "red", fill = "grey66",size = 2), # Customize reg. line
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.sep = "\n",size=6),
                  xlab = "1-Palmitoyllysophosphatidylcholine", ylab = colnames(dddf_rel)[i])
  p1 <- ggpar(p1,
              font.tickslab = c(10,"bold", "black"),font.x = c(12, "bold"),
              font.y = c(12, "bold"))
  print(i)
  print(p1)
  myplots[[i]] <- p1  # add each plot into plot list
}
#Arrange on one page
p_total <- ggarrange(plotlist = myplots,labels = c("B","C","D","E","F",
                                                   "G","H","I"),
                     ncol = 2,nrow = 4)
p_total
ggsave("assoc_results/Feces_microbial_KO_metabolites_selected.png",width = 8,height = 14)


# 将Pro-hyp抽出来，和与其显著相关的菌做相关性
# 散点图
dat3 <- dat2
names(dat3) <- str_split(names(dat2), ";", simplify = TRUE)[,2]
dat3 <- dat3[, names(dat3) %in% c("s__Clostridium botulinum",
                                  "s__Clostridium celatum",
                                  "s__Clostridium disporicum")]
metabo_t.s <- metabo_t[, names(metabo_t) %in% "Pro-hyp"]
dddf <- cbind(dat3, metabo_t.s)
names(dddf)[ncol(dddf)] <- "Pro_hyp"

dddf_rel <- sweep(dddf, 2, colSums(dddf), '/') #计算相对丰度

#Min-Max Normalization
#trans_max_min <- function(x){
#  (x-min(x))/(max(x)-min(x))
#}
library(ggpubr)
#Two_variable_one[,2:5] <- as.data.frame(apply(Two_variable_one[,2:5],2,trans_max_min))
#dddf[,1:ncol(dddf)] <- as.data.frame(apply(dddf[,1:ncol(dddf)],2,trans_max_min))
myplots <- list()  # new empty list
for (i in 1:3) {
  p1 <- ggscatter(dddf_rel,x= "Pro_hyp",y=colnames(dddf_rel)[i],
                  color = "black",fill="black", shape = 21, size = 1.5, # Points color, shape and size
                  add.params = list(color = "red", fill = "grey66",size = 2), # Customize reg. line
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  cor.coeff.args = list(method = "spearman", label.sep = "\n",size=6),
                  xlab = "Pro-hyp", ylab = colnames(dddf_rel)[i])
  p1 <- ggpar(p1,
              font.tickslab = c(10,"bold", "black"),font.x = c(12, "bold"),
              font.y = c(12, "bold"))
  print(i)
  print(p1)
  myplots[[i]] <- p1  # add each plot into plot list
}
#Arrange on one page
p_total <- ggarrange(plotlist = myplots,labels = c("J","K","L"),ncol = 2,nrow = 2)
p_total
ggsave("assoc_results/Feces_microbial_KO_metabolites_selected2.png",width = 8,height = 8)
