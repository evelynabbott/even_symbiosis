#resamp GO DELTA
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

sym="C"

MF=read.table(paste0("MWU_MF_GO_",sym,".csv"),header=T)
CC=read.table(paste0("MWU_CC_GO_",sym,".csv"),header=T)
BP=read.table(paste0("MWU_BP_GO_",sym,".csv"),header=T)
data1=rbind(MF,CC,BP)

targetcount=8000
MF=read.table(paste0("MWU_MF_GO_",sym,"_",targetcount,".csv"),header=T)
CC=read.table(paste0("MWU_CC_GO_",sym,"_",targetcount,".csv"),header=T)
BP=read.table(paste0("MWU_BP_GO_",sym,"_",targetcount,".csv"),header=T)
data2=rbind(MF,CC,BP)

goods=intersect(data1$term,data2$term)
length(goods)
data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]
# all overlapping GO terms
ress=merge(data1,data2,by="term")

ress$color = ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.x,ignore.case = T),"red",
                    ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.x,ignore.case = T),"blue",
                           ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.y,ignore.case = T),"red",
                                  ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.y,ignore.case = T),"blue","black"))))

ress$Function <- ifelse(ress$color == "red", "translation",
                        ifelse(ress$color == "blue", "photosynthesis", "other"))

ress1=ress

#----------------------------------------------
targetcount=15000
MF=read.table(paste0("MWU_MF_GO_",sym,"_",targetcount,".csv"),header=T)
CC=read.table(paste0("MWU_CC_GO_",sym,"_",targetcount,".csv"),header=T)
BP=read.table(paste0("MWU_BP_GO_",sym,"_",targetcount,".csv"),header=T)
data2=rbind(MF,CC,BP)

goods=intersect(data1$term,data2$term)
length(goods)
data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]
# all overlapping GO terms
ress=merge(data1,data2,by="term")

ress$color = ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.x,ignore.case = T),"red",
                    ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.x,ignore.case = T),"blue",
                           ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.y,ignore.case = T),"red",
                                  ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.y,ignore.case = T),"blue","black"))))

ress$Function <- ifelse(ress$color == "red", "translation",
                        ifelse(ress$color == "blue", "photosynthesis", "other"))

ress2=ress
#----------------------------------------------
targetcount=30000
MF=read.table(paste0("MWU_MF_GO_",sym,"_",targetcount,".csv"),header=T)
CC=read.table(paste0("MWU_CC_GO_",sym,"_",targetcount,".csv"),header=T)
BP=read.table(paste0("MWU_BP_GO_",sym,"_",targetcount,".csv"),header=T)
data2=rbind(MF,CC,BP)

goods=intersect(data1$term,data2$term)
length(goods)
data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]
# all overlapping GO terms
ress=merge(data1,data2,by="term")

ress$color = ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.x,ignore.case = T),"red",
                    ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.x,ignore.case = T),"blue",
                           ifelse(grepl("translat|ribosom|rRNA|ribonu",ress$name.y,ignore.case = T),"red",
                                  ifelse(grepl("photosy|thylak|plastid|chlorop",ress$name.y,ignore.case = T),"blue","black"))))

ress$Function <- ifelse(ress$color == "red", "translation",
                        ifelse(ress$color == "blue", "photosynthesis", "other"))

ress3=ress

#photosynthesis/translation
a=ggplot(ress1, aes(x = delta.rank.x, y = delta.rank.y, color=Function))+
  #geom_point(size = 4, aes(alpha=Function))+ geom_text_repel(aes(label=name.x), size = 3, point.padding = 1)+
  geom_point(size = 3, aes(alpha=Function))+
  scale_color_manual(breaks=c("translation","photosynthesis","other"), values = c("firebrick1","dodgerblue2","grey80"))+
  scale_alpha_manual(breaks = c("translation","photosynthesis","other"), values = c(1, 1,.08))+
  scale_x_continuous(name=paste(sym,"all"))+ 
  scale_y_continuous(name=paste(sym,"8k"))+
  #coord_cartesian(xlim=c(-2000,2000),ylim = c(-2000,2000))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept = 0, linetype = "dashed", color="gray70")

b=ggplot(ress1, aes(x = delta.rank.x, y = delta.rank.y, color=Function))+
  #geom_point(size = 4, aes(alpha=Function))+ geom_text_repel(aes(label=name.x), size = 3, point.padding = 1)+
  geom_point(size =3, aes(alpha=Function))+
  scale_color_manual(breaks=c("translation","photosynthesis","other"), values = c("firebrick1","dodgerblue2","grey80"))+
  scale_alpha_manual(breaks = c("translation","photosynthesis","other"), values = c(1, 1,.1))+
  scale_x_continuous(name=paste(sym,"all"))+ 
  scale_y_continuous(name=paste(sym,"15k"))+
  #coord_cartesian(xlim=c(-2000,2000),ylim = c(-2000,2000))+
  theme(legend.position = "none")+
  geom_hline(yintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept = 0, linetype = "dashed", color="gray70")

c=ggplot(ress1, aes(x = delta.rank.x, y = delta.rank.y, color=Function))+
  #geom_point(size = 4, aes(alpha=Function))+ geom_text_repel(aes(label=name.x), size = 3, point.padding = 1)+
  geom_point(size = 3, aes(alpha=Function))+
  scale_color_manual(breaks=c("translation","photosynthesis","other"), values = c("firebrick1","dodgerblue2","grey80"))+
  scale_alpha_manual(breaks = c("translation","photosynthesis","other"), values = c(1, 1,.08))+
  scale_x_continuous(name=paste(sym,"all"))+ 
  scale_y_continuous(name=paste(sym,"30k"))+
  theme(legend.position = "none")+
  #coord_cartesian(xlim=c(-2000,2000),ylim = c(-2000,2000))+
  geom_hline(yintercept = 0, linetype = "dashed", color="gray70")+
  geom_vline(xintercept = 0, linetype = "dashed", color="gray70")

library(gridExtra)

grid.arrange(a, b,c,
             ncol = 3)
