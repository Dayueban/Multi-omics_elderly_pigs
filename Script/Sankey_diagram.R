#Alluvial diagram
####################################
library(ggplot2)
library(ggalluvial)
connect_table <- zotu_occ_core[, c(1,2,3,10)]
connect_table[85,4] <- "Unknown"
connect_table <-na.omit(connect_table)
connect_table$occ <- connect_table$occ/822317
colnames(connect_table)[4] <- "zotuGenus"
library(dplyr)
test <- left_join(connect_table, HOMD_1522[, c(1,7)]) 
test=test[order(test$occ,decreasing = T),]
test$otu = factor(test$otu, levels = unique(test$otu))
levels(test$otu)
test$zOTU = factor(test$zOTU, levels = c("Zotu1",  "Zotu2",  "Zotu3",  "Zotu4",  "Zotu5",  "Zotu6",  "Zotu7",  "Zotu8",  "Zotu15", "Zotu9", "Zotu10", "Zotu11", "Zotu13", "Zotu14","Zotu21", "Other" ))
levels(test$zOTU)
test$Genus = factor(test$Genus, levels = unique(test$Genus))
levels(test$Genus)
write.table(test, file = "connecttable.txt",row.names=FALSE,col.names=TRUE, sep="\t")

library(viridis)
library(readxl)
connecttable <- read_excel("connecttable.xlsx")
connecttable=connecttable[order(connecttable$occ,decreasing = T),]

colnames(connecttable)[1] <- "otu"
connecttable$otu = factor(connecttable$otu, levels = rev(unique(connecttable$otu)))
levels(connecttable$otu)
connecttable$Genus = factor(connecttable$Genus, levels = c( "999_999", "Rothia", "Leptotrichia" , "Gemella" , "Campylobacter"  
                                                            ,"Aggregatibacter", "Streptococcus",   "Veillonella" ,    "Fusobacterium"  , "Prevotella"  
                                                            ,"Neisseria" ,      "Haemophilus"   ))
levels(connecttable$Genus)
connecttable$zOTU = factor(connecttable$zOTU, levels = rev(unique(connecttable$zOTU)))
levels(connecttable$zOTU)
connecttable$zotuGenus = factor(connecttable$zotuGenus, levels = c("999_999","Haemophilus" ,    "Neisseria",       "Prevotella" ,     "Fusobacterium"  ,"Veillonella","Streptococcus" ,  "Aggregatibacter", "Campylobacter",   "Gemella","Leptotrichia","Rothia" ))
levels(connecttable$zotuGenus)
View(connecttable)    

#Fig7A
alluvial_white <-ggplot(data = connecttable,
                        aes(axis1 = Genus, axis2 = otu, 
                            axis3 = zOTU,
                            axis4 =zotuGenus,
                            y = occ)) +
  scale_x_discrete(limits = c("Genus"
                              ,"otu"
                              , "zOTU"
                              , "zotuGenus"
  ), expand = c(.1, .1)) +
  xlab("Demographic") +
  scale_fill_viridis(discrete=T)+
  geom_alluvium(aes(fill = zOTU)) +
  scale_y_continuous(expand=c(0.01,0.01))+
  geom_stratum(color = "grey") +
  #  geom_flow(stat = "alluvium", lode.guidance = "rightleft", 
  #            color = "darkgray") + 
  geom_text(stat = "stratum", size=3, aes(label = after_stat(stratum))) +
  theme(
    #text=element_text(family="sans", size=10, face = "bold"),
    axis.title = element_text(family="sans", size=10, face = "bold"),
    axis.text.y = element_text(family="sans", size=8, face = "bold"),
    axis.text.x = element_text(family="sans", size=8),
    #axis.line = element_line(colour = "black", size = 0.5),
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    panel.background = element_blank(),
    legend.title = element_text(size = 8, face = "bold")
  )
alluvial_white