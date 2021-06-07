#coral host gene expression
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(limma)
library(viridis)
library(hexbin)

#deseq================================================================================================

setwd("~/Dropbox/codominant_symbiosis/")

evenness=function(x, padding=0.0001){
  x[x<padding]=padding
  x[(1-x)<padding]=1-padding
  p=x
  q=1-x
  H=-(p*log(p)+q*log(q))
  Hmax=-2*0.5*log(0.5)
  return(H/Hmax)
}

#ll=load('counts_metadata/coldata_redME_evenness.Rdata') 
ll=load('counts_metadata/coldata.Rdata') 
dim(coldata)

# calculating evenness
cc=10^(coldata$logCcount)
dd=10^(coldata$logDcount)
coldata$Cprop=cc/(cc+dd)
coldata$Dprop=1-coldata$Cprop
coldata$evenness=evenness(coldata$Cprop)

# # remove WGCNA outliers
# coldata <- coldata[order(coldata$logCDratio),]
# coldata = coldata[-c(1:10),]

#remove  Barshis outliers 
outliers=c("SRR629129",
"SRR629150",
"SRR629152",
"SRR629133",
"SRR629130",
"SRR629151",
"SRR629134",
"SRR629153",
"SRR629155")
coldata=coldata[!(coldata$Run %in% outliers),]

ll=load('counts_metadata/host_counts')

rownames(coldata)=coldata$Run
keeps = row.names(coldata)
counts = counts[,keeps]
dim(coldata)

# normalizing pool factors
coldata$pool[coldata$pool=="300"]="HV"
coldata$pool[coldata$pool=="400"]="MV"

# table(coldata$my_title)
# plot(coldata$evenness[coldata$my_title=="L_Barshis_bleachResillience_PRJNA177515"])
# plot(coldata$evenness[coldata$my_title=="k1_Palumbi_lab_heat_resilience_PRJNA274410"])

# 
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut) #how many genes have mean below the cutoff?
counts=cc[means>cut,]

coldata$hostCount=apply(cc,2,sum)
save(coldata,file="~/Dropbox/codominant_symbiosis/counts_metadata/coldata_evenness_noOuts.RData")

# fitting DESeq model with covariates ("time" captures both study and timepoint)  
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*evenness+pool+time)
dds=DESeq(dds)

ddsd = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*Dprop+pool+time)
ddsd=DESeq(ddsd)

vsd=assay(vst(dds))
vsd=removeBatchEffect(vsd,batch=coldata$pool,batch2=coldata$time)
pheatmap(cor(vsd))

# Dixon modules
ll=load("moduleAssignment.Rdata")
# extracting red module genes
reds=row.names(geneModuleMembership[moduleColors=="red" & geneModuleMembership$MMred>0.25,])
length(reds)

red.vsd=data.frame(vsd[row.names(vsd) %in% reds,])
dim(red.vsd)
pheatmap(cor(red.vsd))

# computing red module eigengene
library(vegan)
redpc=rda(t(red.vsd)~1)
coldata$redME=scores(redpc,display="sites")[,1]

# plotting red module eigengene vs evenness and D proportion
coldata$treatment=coldata$c_t
coldata$treatment=factor(coldata$treatment,levels=c("h","c"))

save(coldata,file="~/Dropbox/codominant_symbiosis/counts_metadata/coldata_redME_evenness_noOuts.RData")

gg=ggplot(coldata,aes(evenness,redME,color=treatment))+
  geom_point()+
  geom_smooth(method=lm)
gg
pdf("redME_vs_evenness.pdf",height=3,width=3.4)
plot(gg)
dev.off()

gg=ggplot(coldata,aes(Dprop,redME,color=treatment))+
  geom_point()+
  geom_smooth(method=lm)+
  xlab("D proportion")
gg
pdf("redME_vs_Dprop.pdf",height=3,width=3.4)
plot(gg)
dev.off()

library(lme4)
library(lmerTest)

summary(lmer(redME~Dprop+(Dprop|geno),subset(coldata,treatment=="c")))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.1942     0.0809 71.9253  -14.76   <2e-16 ***
#   Dprop         0.1838     0.1171 44.5353    1.57    0.124  

summary(lmer(redME~Dprop+(Dprop|geno),subset(coldata,treatment=="h")))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.2868     0.1910 10.1266   6.738 4.82e-05 ***
#   Dprop        -0.3174     0.2627 32.8201  -1.209    0.235 

summary(lmer(redME~evenness+(evenness|geno),subset(coldata,treatment=="c")))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.05115    0.06704 28.46623 -15.679 1.55e-15 ***
#   evenness    -0.13083    0.21554  3.69027  -0.607    0.579  

summary(lmer(redME~evenness+(evenness|geno),subset(coldata,treatment=="h")))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.2570     0.1640 23.2406   7.665 8.26e-08 ***
#   evenness     -0.7494     0.3168 39.4378  -2.366    0.023 *

#----------Dixon blobs

mct=results(dds,name="c_t_h_vs_c")
me=results(dds,name="evenness")
mei=results(dds,name="c_th.evenness")
mdi=results(ddsd,name="c_th.Dprop")
mctd=results(ddsd,name="c_t_h_vs_c")

#---- collecting results:
# interaction term
ei=mei$log2FoldChange
# change in highly even case (evenness=1):
ec=mei$log2FoldChange+mct$log2FoldChange
# change in uneven case (evenness=0):
uc=mct$log2FoldChange
# change in C-only case (low D):
duc=mctd$log2FoldChange
# effect of D on stress response
di=mdi$log2FoldChange


# plotting response to heat in even case against response in uneven case
gg=ggplot(data.frame(cbind(even=ec,uneven=uc)),aes(uneven, even))+
  geom_hex(bins=c(20,20))+
  xlab("heat response\nwith uneven symbionts")+
  ylab("heat response\nwith even symbionts")+
  scale_fill_viridis(trans="log")+
  coord_equal()+
  geom_abline(intercept=0,slope=1,color="grey80")+
  geom_smooth(method="lm",color="red")+
  ylim(-5,5)+
  geom_vline(xintercept=0,linetype="dotted")+
  geom_hline(yintercept=0,linetype="dotted")
gg
pdf("even_uneven_heat_change.pdf",height=4,width=4)
plot(gg)
dev.off()

# interaction term is negatively correlated with base response ==> evenness *diminishes* response
gg=ggplot(data.frame(cbind(even.ixn=ei,heat.resp=uc)),aes(heat.resp, even.ixn))+
  geom_hex()+
  ylim(-5,5)+
  #  coord_equal()+
  xlab("heat response\nwith uneven symbionts")+
  ylab("change in response\nwith even symbionts")+
  scale_fill_viridis(trans="log")+
  #  geom_abline(intercept=0,slope=1,color="grey80")+
  geom_smooth(method="lm",color="red")+
  geom_vline(xintercept=0,linetype="dotted")+
  geom_hline(yintercept=0,linetype="dotted")
gg
pdf("interaction_vs_heatResp.pdf",height=3,width=4)
plot(gg)
dev.off()


# heat response with all-D vs heat response with all-C
gg=ggplot(data.frame(cbind(D.ixn=di+duc,heat.resp=duc)),aes(heat.resp, D.ixn))+
  geom_hex()+
  ylim(-5,5)+
  #  coord_equal()+
  xlab("heat response\nwith all-C symbionts")+
  ylab("heat response\nwith all-D symbionts")+
  scale_fill_viridis(trans="log")+
   geom_abline(intercept=0,slope=1,color="grey80")+
  geom_smooth(method="lm",color="red")+
  geom_vline(xintercept=0,linetype="dotted")+
  geom_hline(yintercept=0,linetype="dotted")
gg
pdf("heat_allD_vs_allC.pdf",height=3,width=4)
plot(gg)
dev.off()

# change in response with more D
gg=ggplot(data.frame(cbind(D.ixn=di,heat.resp=duc)),aes(heat.resp, D.ixn))+
  geom_hex(bins=c(20,20))+
  xlab("heat response\nwith uneven symbionts")+
  ylab("heat response\nwith even symbionts")+
  scale_fill_viridis(trans="log")+
  coord_equal()+
  geom_abline(intercept=0,slope=1,color="grey80")+
  geom_smooth(method="lm",color="red")+
  ylim(-5,5)+
  geom_vline(xintercept=0,linetype="dotted")+
  geom_hline(yintercept=0,linetype="dotted")
gg
pdf("Dinteraction_vs_heatResp.pdf",height=3,width=4)
plot(gg)
dev.off()


#--------- pvalues
# significance of negative correlation between base response and interaction term
summary(lm(ei~uc))
# (Intercept)  0.083893   0.003571   23.49   <2e-16 ***
#   uc          -0.268352   0.005223  -51.38   <2e-16 ***

summary(lm(di~duc))
# (Intercept) -0.059632   0.002383  -25.02   <2e-16 ***
#   duc         -0.257432   0.003251  -79.19   <2e-16 ***

#-------- evenness~genotype

summary(lm(evenness~geno,coldata))
# Multiple R-squared:  0.7432,	Adjusted R-squared:  0.7085 
# F-statistic: 21.45 on 17 and 126 DF,  p-value: < 2.2e-16

ggplot(coldata,aes(geno,evenness,color=treatment))+
  geom_jitter(width=0.1)+
  xlab("genotype")

ggplot(coldata,aes(geno,Dprop,color=treatment))+
  geom_jitter(width=0.1)+
  xlab("genotype")+
  ylab("D proportion")

plot(logCcount~logDcount, coldata)
