#SYMBIONT GENE EXPRESSION
#DESeq2, WGCNA input, GOMWU input
rm(list=ls())
setwd("~/Dropbox/even")

sym="D" #C=cladocopium; D=durusdinium
load(file="~/Dropbox/even/metadata/coldata.Rdata") #metadata
load(paste0("~/Dropbox/even/metadata/clade",sym,"_counts.Rdata")) #counts

#-----------------------------------------------------------------------------------------
#DESeq2
#-----------------------------------------------------------------------------------------
#symbiont gene expression
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)

#remove samples from counts
keeps <- rownames(coldata)
keeps  
counts = counts[,colnames(counts) %in% keeps]

#remove genes with low coverage
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]

#dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*evenness+logCcount+pool+time) #for cladocopium
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*evenness+logDcount+pool+time) #for durusdinium
dds=DESeq(dds)

resultsNames(dds)
rld = vst(dds)
rld.df=assay(rld)
colnames(rld.df) = coldata$Run
vsd = rld.df

library(limma)
vsd=removeBatchEffect(vsd,batch=coldata$time,covariates=coldata$logDcount)

save(coldata,counts,vsd,dds, file=paste0("~/Dropbox/even/metadata/",sym,"_deseq_results.Rdata"))
#-----------------------------------------------------------------------------------------------
#WGCNA input
#-----------------------------------------------------------------------------------------------
datExpr=as.data.frame(t(vsd))
datTraits=coldata
save(datExpr,datTraits,file=paste0("~/Dropbox/even/wgcna/wgcna_input",sym,".Rdata")) #do this twice, once for cladocopium and once for durusdinium

#-----------------------------------------------------------------------------------------------
#GO_MWU input 
#-----------------------------------------------------------------------------------------------
resultsNames(dds)

#evenness
res = results(dds, name ="evenness")
goout=data.frame("gene"=rownames(res),
                 "logP"=-log(res$pvalue,10))


sign=res$log2FoldChange>0
goout=goout %>% 
  mutate(logP=if_else(sign,
                      logP,
                      -logP)) %>%
  write_csv(path=paste0("~/Dropbox/even/GO_MWU/GO_",sym,".csv"))
