#----------------------------------------------------------------------------------------
#Get number of differentially expressed genes; supplemental fig S2
#---------------------------------------------------------------------------------------
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(gridExtra)

rm(list=ls())
sym="C"
load(paste0("~/Dropbox/even/metadata/",sym,"_deseq_results.Rdata"))
res = results(dds, name ="evenness")
summary(res)
res
C = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
C
nsig = sum(C$sig)
tot = nrow(C)
a=ggplot(data=C, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) +
  theme(legend.position = "none")+
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle="A") +
  xlim(-2,3)+
  ylim(0,30)+
  guides(color=guide_legend(title="FDR<0.1"))
plot(a)
#-----------------------------------------------------------------------------------
sym="D"
load(paste0("~/Dropbox/even/metadata/",sym,"_deseq_results.Rdata"))
res = results(dds, name ="evenness")
summary(res)
res
D = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
D
nsig = sum(D$sig)
tot = nrow(D)
b=ggplot(data=D, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  theme(legend.position = "none")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle="B") +
  xlim(-2, 2)+
  guides(color=guide_legend(title="FDR<0.1"))
plot(b)

#-----------------------------------------------------------------------------------
load("~/Dropbox/even/metadata/host_deseq_results.Rdata")
res = results(dds, name ="evenness")
summary(res)
res
Host = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
Host
nsig = sum(Host$sig)
tot = nrow(Host)
c=ggplot(data=Host, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  theme(legend.position = "none")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle="C") +
  xlim(-2, 2)+
  guides(color=guide_legend(title="FDR<0.1"))
#-----------------------------------------------------------------------------------
sym="C"
load(paste0("~/Dropbox/even/metadata/",sym,"_deseq_results.Rdata"))
res = results(dds, name ="c_t_h_vs_c")
summary(res)
res
C1 = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
C1
nsig = sum(C1$sig)
tot = nrow(C1)
d=ggplot(data=C1, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) +
  theme(legend.position = "none")+
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle="D") +
  xlim(-2, 2)+
  guides(color=guide_legend(title="FDR<0.1"))

#-----------------------------------------------------------------------------------
sym="D"
load(paste0("~/Dropbox/even/metadata/",sym,"_deseq_results.Rdata"))
res = results(dds, name ="c_t_h_vs_c")
summary(res)
res
D1 = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
D1
nsig = sum(D1$sig)
tot = nrow(D1)
e=ggplot(data=D1, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  theme(legend.position = "none")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle="E") +
  xlim(-2, 2)+
  guides(color=guide_legend(title="FDR<0.1"))

#-----------------------------------------------------------------------------------
load("~/Dropbox/even/metadata/host_deseq_results.Rdata")
res = results(dds, name ="c_t_h_vs_c")
summary(res)
res
Host1 = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
Host1
nsig = sum(Host1$sig)
tot = nrow(Host1)
f=ggplot(data=Host1, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  theme(legend.position = "none")+
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle="F") +
  xlim(-2, 2)+
  guides(color=guide_legend(title="FDR<0.1"))



grid.arrange(a,b,c,d,e,f,
             nrow = 2)

