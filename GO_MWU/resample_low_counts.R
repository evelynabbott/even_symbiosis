#-----------------------------------------------------------------------------------------------
#Resample counts to see if low count samples skew results
#-----------------------------------------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
#theme_set(theme_cowplot())
library(data.table)
rm(list=ls())

#choose symbiont
sym="D" #C=cladocopium; D=durusdinium
load(file="~/Dropbox/even/metadata/coldata.Rdata") #metadata
load(paste0("~/Dropbox/even/metadata/clade",sym,"_counts.Rdata")) #counts
keeps = row.names(coldata)
counts = counts[,keeps]

#remove genes with low coverage
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]

zooxsums=apply(counts,2,sum)

#choose read count to subsample
#targetcount = 8000
#targetcount = 15000 
targetcount = 30000


cladeC.resampled=c()
for (s in 1:ncol(counts)) {
  probs= counts[,s]/zooxsums[s]
  cts=hist(sample(c(1:nrow(counts)), targetcount,prob=probs,replace=TRUE),breaks=c(0:nrow(counts)),plot=F)$counts
  cladeC.resampled=data.frame(cbind(cladeC.resampled, cts))
}
#reassign the colnames back to counts
colnames(cladeC.resampled) = colnames(counts)
colnames(cladeC.resampled)
rownames(cladeC.resampled) = rownames(counts)

str(cladeC.resampled)
str(counts)

zooxsums=apply(cladeC.resampled,2,sum)

counts <- cladeC.resampled
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*evenness+logCcount+pool+time)
dds=DESeq(dds)

#get variance stabilized counts and save them
#rld = vst(dds)
rld = varianceStabilizingTransformation(dds)
rld.df=assay(rld)
colnames(rld.df) = coldata$Run
vsd = rld.df

res=results(dds,name="evenness")

#save for GO_MWU
goout=data.frame("gene"=rownames(res),
                 "logP"=-log(res$pvalue,10))

sign=res$log2FoldChange>0
goout=goout %>% 
  mutate(logP=if_else(sign,
                      logP,
                      -logP))%>%
  write_csv(path= paste0("~/Dropbox/even/GO_MWU/GO_",sym,"_",targetcount,".csv"))





#======================================================================================
#GO
setwd("~/Dropbox/even/GO_MWU")
input=paste0("GO_",sym,"_",targetcount,".csv") # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="c_and_d_sym_GO_only_BPMFCCrem.tsv" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")
# Calculating stats. It might take ~3 min for  MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)

goDivision="CC" # either MF, or BP, or CC
source("gomwu.functions.R")
# Calculating stats. It might take ~3 min for  MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")
# Calculating stats. It might take ~3 min for  MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=50,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)


#PLOT RESULTS WITH GO_deltaranks.R
