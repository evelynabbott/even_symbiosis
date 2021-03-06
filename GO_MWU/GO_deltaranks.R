#GO delta ranks
library(ggrepel)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)
library(tidyverse)
library(cowplot)
rm(list=ls())
setwd("~/Dropbox/even/GO_MWU")

#combine molecular function, cellular component, and biological processes GO_MWU files
MF=read.table("MWU_MF_GO_C.csv",header=T)
CC=read.table("MWU_CC_GO_C.csv",header=T)
BP=read.table("MWU_BP_GO_C.csv",header=T)
data1=rbind(MF,CC,BP)

MF=read.table("MWU_MF_GO_D.csv",header=T)
CC=read.table("MWU_CC_GO_D.csv",header=T)
BP=read.table("MWU_BP_GO_D.csv",header=T)
data2=rbind(MF,CC,BP)

#also do this for the resampled counts

goods=intersect(data1$term,data2$term)
#goods=unique(as.character(c(data1$term[data1$p.adj<=0.1],data2$term[data2$p.adj<=0.1])))
length(goods)

data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]

# all overlapping GO terms
ress=merge(data1,data2,by="term")
#assign colors based on genes

ress$color = ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.x,ignore.case = T),"red",
                    ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.x,ignore.case = T),"blue",
                           ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.y,ignore.case = T),"red",
                                  ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.y,ignore.case = T),"blue","black"))))

ress$Function <- ifelse(ress$color == "red", "translation",
                        ifelse(ress$color == "blue", "photosynthesis", "other"))


#photosynthesis/translation
ggplot(ress, aes(x = delta.rank.x, y = delta.rank.y, color=Function))+
  #geom_point(size = 4, aes(alpha=Function))+ geom_text_repel(aes(label=name.x), size = 3, point.padding = 1)+
  geom_point(size = 4, aes(alpha=Function))+
  scale_color_manual(breaks=c("translation","photosynthesis","other"), values = c("firebrick1","dodgerblue2","grey80"))+
  scale_alpha_manual(breaks = c("translation","photosynthesis","other"), values = c(1, 1,.05))+
  scale_x_continuous(name="C baseline")+ 
  scale_y_continuous(name="D baseline")+
  #coord_cartesian(xlim=c(-2000,2000),ylim = c(-2000,2000))+
  geom_hline(yintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept = 0, linetype = "dashed", color="gray70")

#subset and plot alone
ress1 <- ress

# GO terms highly signifcant in any of the two datasets
ress1 = ress %>% filter(p.adj.x<=0.001 | p.adj.y<=0.001)

ggplot(ress1, aes(x = delta.rank.x, y = delta.rank.y, color=Function))+ 
  #geom_point(size = 4)+ #geom_text_repel(aes(label=name.y), size = 2.0, point.padding = 0.9)+
  geom_point(size = 4, aes(alpha=Function))+ geom_text_repel(aes(label=name.y), size = 3.0, point.padding = 1.2)+
  scale_color_manual(breaks=c("translation","photosynthesis", "other"), values = c("firebrick1","dodgerblue2","gray40"))+
  scale_alpha_manual(breaks = c("translation","photosynthesis", "other"), values = c(1, 1, 1))+
  scale_x_continuous(name="Cladocopium")+ 
  scale_y_continuous(name="Durusdinium")+
  coord_cartesian(xlim=c(-2300,2200),ylim = c(-2200,2300))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept = 0, linetype = "dashed", color="gray70")

