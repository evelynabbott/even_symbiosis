#plot evenness against study, pool, and host genotype

library(ggplot2)
rm(list=ls())
load("~/Dropbox/even/metadata/coldata.Rdata")

#Turn your 'treatment' column into a character vector
rownames(coldata)=coldata$Run
coldata <- coldata[order(coldata$logCDratio),]
coldata$Run <- as.character(coldata$Run)
#Then turn it back into a factor with the levels in the correct order
coldata$Run <- factor(coldata$Run, levels=unique(coldata$Run)) 
#coldata$Run = rownames(coldata)

coldata$percentc=coldata$cladeC/((coldata$cladeC)+(coldata$cladeD))
coldata <- coldata[order(coldata$percentc),]
coldata <- coldata[order(coldata$logCDratio),]

#make nicer labels
coldata$study=coldata$my_title
coldata$study=gsub("L_Barshis_bleachResillience_PRJNA177515","Barshis_2013",coldata$study)
coldata$study=gsub("k1_Palumbi_lab_heat_resilience_PRJNA274410","Seneca_Palumbi_2015",coldata$study)

coldata$treatment=coldata$c_t
coldata$treatment=gsub("c","control",coldata$treat)
coldata$treatment=gsub("h","heated",coldata$treat)

coldata$genotype=coldata$geno


a=ggplot(coldata,aes(x=Run,y=evenness,color=study,group=1))+
  geom_point(size=2)+
  theme_classic()+
  xlab(" ")+
  theme(legend.position= "none")+ #comment out to view legend
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

b=ggplot(coldata,aes(x=Run,y=evenness,color=treatment,group=1))+
  geom_point(size=2)+
  theme_classic()+
  xlab(" ")+
  theme(legend.position= "none")+ #comment out to view legend
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

c=ggplot(coldata,aes(x=Run,y=evenness,color=genotype,group=1))+
  geom_point(size=2)+
  theme_classic()+
  xlab("samples from least to most Cladocopium")+
  theme(legend.position= "none")+ #comment out to view legend
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

library(gridExtra)

grid.arrange(a, b,c,
             nrow = 3)
