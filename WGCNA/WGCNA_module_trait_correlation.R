#wgcna4_module-correlations.R
#This script calculates and plots relationships between 
#module eigengenes and sample traits. Inputs come from:
#get_variance_stabilized_counts.R and wgcna3b_step-wise_network_construction.R
#code is adapted from examples given here: https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/

rm(list=ls())

setwd("~/Dropbox/even/wgcna")

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the WGCNA package
library(WGCNA)
library(tidyverse)
source('wgcna_functions.R')
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#pick symbiont
sym="D"

# Load the expression and trait data saved in the first part
lnames = load(file = paste0("wgcna_input",sym,".Rdata"));

#The variable lnames contains the names of loaded variables.
lnames #(when you load a .Rdata object, if can have multiple variables assigned, the list of them is saved in this lnames variable, in this case the two objects you loaded are "datExpr" and "datTraits")


# Load network data saved in the second part.
lnames = load(file=paste0(sym,"_networkdata_signed.RData"))
#lnames = load(file = "/home/evelyn/Dropbox/project/bleww_final/Zoox_expr/fightzone/D/wgcna/test/anothertry2/wgcna3b_manual_sft12_minModSize30_cutDpth0.4_signed.Rdata") #D outliers removed
lnames


#match up the rownames in the expression 
rownames(datExpr)[1:10]     #this tells us the rownames as given to WGCNA
rownames(datTraits)[1:10]   #this is for the entire set of datTraits and includes samples that were not put into WGCNA
rownames(MEs) = rownames(datExpr)
#datTraits = datTraits[rownames(datTraits) %in% rownames(datExpr),]  #reduce datTraits to those we have WGCNA results for
#if everything is matched up this should say TRUE
sum(rownames(datTraits) == rownames(datExpr)) == length(rownames(datExpr))
head(datTraits)
dim(datTraits)



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate module eigengenes with color labels
#This will provide the first principal component for expression
#behavior for the genes in each module. On that principal component,
#each sample will have a loading value. These values can be corelated
#with our known sample traits to get an idea what biological mechanism
#the co-reguated genes in a given module might be responding to
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs1 = MEs0[rownames(datTraits),]


########################################
############ SUBSET BY SIZE ############
########################################

#plot all modules
minModSizeToPlot = 0;outputMM=TRUE

# plot only large modules
#minModSizeToPlot = 600;outputMM=FALSE

########################################
########################################
########################################


module.sizes = table(moduleColors)
passing=names(module.sizes)[module.sizes>minModSizeToPlot]
MEs = orderMEs(MEs1[,paste('ME', passing, sep='')])

# 
########################################
########################################
#look at logCcount, CD ratio, minorLR
datTraits = datTraits %>%
  select(evenness, Dprop)

#use the cor() function to get the correlations between the module eigengenes and the trait data
moduleTraitCor = cor(MEs, datTraits, use = "p");
#get p values as well
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#write out module loadings
mEigs=MEs
rownames(mEigs) = rownames(datTraits)
save(mEigs, file='moduleEigengenes_treatment.Rdata')



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

# # Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

rows = rownames(moduleTraitCor)
sub.colors = substr(rownames(moduleTraitCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = module.sizes,
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               cex.lab.y = .8,
               cex.lab.x = .8,
               zlim = c(-1,1),
               font.lab.y = 1,
               main = paste("Module-trait relationships"))

#====================================================================================
#make a heatmap showing sample vs module correlation

library(pheatmap)
library(RColorBrewer)
#Build heatmap of MEs - samples as rows and MEs as columns
MEsm = as.matrix(MEs)
#make rownames in order of most C to least C
datTraits1 <- datTraits[order(datTraits$logCDratio),]
MEsm = MEsm[match(rownames(datTraits1), rownames(MEsm)),]

coul <- coul <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(25))
pheatmap(t(MEsm),scale="row",cluster_cols=FALSE, cluster_rows=F,labels_col = "",color = coul)

#=====================================================================================
#REPLOT THE MODULE EIGENGENE CLUSTERING TREE
#This shows you how the modules are similar to one another.
#Pushing the merging theshold (argument 3 for wgcna3b_step-wise_network_construction.R) will join modules together
#This is like moving a horizontal line up this tree figure, if the line is above a node, modules below that node will be joined into one
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

#plot them
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.4, lty=2, col='red')

