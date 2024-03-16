rm(list=ls())
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tibble)
library(ggplot2)
library(tidyr)
library(grid)
setwd('~/Documents/nosZ_probecapture/testing_readmapping/')

place.rel<-function(tab) {
  
  tab.sp<-split(tab,tab$edge_num)
  tab.rel<-tibble()

  for(i in 1:length(tab.sp)) {
    
    tab.in<-tab.sp[[i]]
    nseq<-nrow(tab.in)
    nseq.rel<-nseq/nrow(tab)
    tab.rel<-rbind(tab.rel,cbind(tab.in[1,],nseq=nseq,nseq.rel=nseq.rel))
  
  }
  return(tab.rel)
}

shotgun<-read.jplace('epa_result_sh.jplace')
#shotgun@phylo$edge.length<-sqrt(shotgun@phylo$edge.length)

probe<-read.jplace('epa_result_pr.jplace')
#probe@phylo$edge.length<-sqrt(probe@phylo$edge.length)

shotgun.best<-get.placements(shotgun,by='best')
probe.best<-get.placements(probe,by='best')

probe.best<-probe.best[probe.best$node != 3283,]
shotgun.best<-shotgun.best[shotgun.best$node != 3283,]

shotgun.best.rel<-place.rel(shotgun.best)
probe.best.rel<-place.rel(probe.best)

#dropping outgroup?
probe.best.rel<-probe.best.rel[probe.best.rel$node != 3283,]
shotgun.best.rel<-shotgun.best.rel[shotgun.best.rel$node != 3283,]

#2392
#1704
p.shot<-ggtree(shotgun,layout='circular') %<+% shotgun.best.rel
p.shot + geom_point2(aes(subset=node %in% node,size=nseq.rel),color='red',alpha=0.4)+geom_rootedge()+scale_size(range=c(1,15),limits=c(0.00001,0.05))+
         geom_cladelabel(node=2337,label='Clade I', offset.text=0.5, hjust = 1, fontsize=8)+
         geom_cladelabel(node=1652,label='Clade II', offset.text=0.5, hjust = 1, fontsize=8)+
         ggtitle(label='nosZ Placement - shotgun metagenonics',subtitle=paste0('Circle size indicates relative abundance (',nrow(shotgun.best),' seqs. total) '))

p.probe<-ggtree(probe,layout='circular') %<+% probe.best.rel
p.probe + geom_point2(aes(subset=node %in% node,size=nseq.rel),color='blue',alpha=0.4)+geom_rootedge()+scale_size(range=c(1,15),limits=c(0.00001,0.05))+geom_treescale(width=0.5)+
  geom_treescale(width=0.5)+
  geom_cladelabel(node=2337,label='Clade I', offset.text=0.5, hjust = 1, fontsize=8)+
  geom_cladelabel(node=1652,label='Clade II', offset.text=0.5, hjust = 1, fontsize=8)+
  ggtitle(label='nosZ Placement - Probe capture',subtitle=paste0('Circle size indicates relative abundance (',nrow(probe.best),' seqs. total) '))

