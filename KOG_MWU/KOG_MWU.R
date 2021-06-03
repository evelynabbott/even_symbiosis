#------------------------------------------------------------------------------
#HOST KOG_MWU
#------------------------------------------------------------------------------
rm(list=ls())
library(KOGMWU)

setwd("~/Dropbox/even/KOG")
data(gene2kog)

agene2kog = read.table("amil_gene2kogClass.tab",header=FALSE,sep = "\t")
ggene2kog = read.table("Amil.v2.eggnogWebsite.gene2kog.tsv",header = TRUE,sep = "\t")
load("~/Dropbox/even/KOG/input/KOG_even_ct_ixn.Rdata")
load("~/Dropbox/even/KOG/input/KOG_even.Rdata")
load("~/Dropbox/even/KOG/input/KOG_ct.Rdata")
GSR = read.csv("~/Dropbox/even/KOG/input/clusterAstress_For_MWU.csv")


#even
even = kog.mwu(Keven,agene2kog)
even

#minorct
heated = kog.mwu(ct,agene2kog)
heated

#even ixn
even_ct_ixn = kog.mwu(Kixn, agene2kog)
even_ct_ixn

#GSR stress
GSR = kog.mwu(GSR,ggene2kog)
GSR

ktable=makeDeltaRanksTable(list("GSR"=GSR,"heated"=heated,"even"=even,"even:heat"=even_ct_ixn))
ktable = na.omit(ktable)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation") 

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#?pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

panel.lm=function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                   cex = 1, col.lmline = "red", ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y[ok]~x[ok]), 
           col = col.lmline, ...)
}

pairs(ktable, lower.panel = panel.lm, upper.panel = panel.cor, cex.labels = 1.5)

# plotting individual delta-rank correlations:
par(mfrow=c(3,1))
corrPlot(y="even",x="GSR",ktable)
corrPlot(y="heated",x="GSR",ktable)
corrPlot(y="even:heat",x="GSR",ktable)

