#EVENNESS
#coral host gene expression

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
rm(list=ls())
#deseq================================================================================================

ll=load('~/Dropbox/even/metadata/host.Rdata') 

cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut) #how many genes have mean below the cutoff?
counts=cc[means>cut,]

dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*evenness+pool+time)
dds=DESeq(dds)

resultsNames(dds)
rld = vst(dds)
rld.df=assay(rld)
colnames(rld.df) = coldata$Run
vsd = rld.df
save(vsd,dds,coldata,counts,file="host_deseq_results.Rdata")

#KOG lfc output ------------------------------------------------------------------------
#for KOG, we do lfc 
resultsNames(dds)

res = results(dds, name ="evenness")
summary(res)
Keven = data.frame("gene"=rownames(res),
                      "lfc"=res$log2FoldChange)
save(Keven,file="~/Dropbox/even/KOG/input/KOG_even.Rdata")

res = results(dds, name ="c_th.evenness")
summary(res)
Kixn = data.frame("gene"=rownames(res),
                  "lfc"=res$log2FoldChange)
save(Kixn,file="~/Dropbox/even/KOG/input/KOG_even_ct_ixn.Rdata")

res = results(dds, name ="c_t_h_vs_c")
summary(res)
ct = data.frame("gene"=rownames(res),
                "lfc"=res$log2FoldChange)
save(ct,file="~/Dropbox/even/KOG/input/KOG_ct.Rdata")


#GO output -----------------------------------------------------------------------------
resultsNames(dds)
res = results(dds, name ="evenness")
summary(res)

goout=data.frame("gene"=rownames(res),
                 "logP"=-log(res$pvalue,10))
sign=res$log2FoldChange>0
goout=goout %>%
  mutate(logP=if_else(sign,
                      logP,
                      -logP)) %>%
  write_csv(path= "~/Dropbox/even/GO_MWU/input/GO_host.csv")


