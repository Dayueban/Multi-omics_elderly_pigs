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

# load working datasheet
options (stringsAsFactors = FALSE)
#work_df <- read.csv("Ileum_OTU.csv",header = T,row.names = 1,check.names = F)
a <- read.csv("Output_data/Saliva_species_tidy_wilcoxn_test_BH.csv", header = T,
              row.names = 1, check.names = F)
a0.05 <- a[which(a$fdr.ma_oa < 0.05), c(16:45)]
a0.05t <- as.data.frame(t(a0.05))

dat <- a0.05t


# 计算模块后的代谢物与16S_ASV之间的相关性以及作图
metabo <- read.csv("Serum_metabo/Differential_Ma_Oa.csv", header = T,
                   row.names = 1, check.names = F)
metabo_t <- as.data.frame(t(metabo))

# 要保证两组数据的个体ID是一样的
dat2 <- dat[c(3:7,9:12,14,16:25),]
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
colnames(asso_df) <- paste(colnames(asso_df),rep(c("cor","P_adj"),c(8,8)),sep = "_")
#' re-order the asso_df according to the adjusted p value
#asso_df <- asso_df[order(asso_df$P_adj,decreasing = FALSE),]
# write.csv(asso_df, "M3C_metabo_transc_asso.csv",row.names = TRUE)

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

# 另外一种作图方法
nbar.m <- nbar.m %>%
  mutate(value_abs = abs(value))
# 将物种名进行裁剪
nbar.m$row <- str_split(nbar.m$row, ";", simplify = TRUE)[,2]
# ab(value) > 0.6的才进行展示
nbar.m0.6 <- nbar.m[which(nbar.m$value_abs > 0.6), ]
write.csv(nbar.m0.6, "assoc_results/Saliva_microbial_assoc_metabolites_Ma_Oa.csv")
pb <- ggplot(nbar.m0.6, aes(row, variable)) +
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
ggsave("Microbial_species_assoc_metabo_Oa_Ma.png", plot = pp, width = 22,height = 14, dpi = 300)



