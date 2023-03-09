
################################
######Double WGCNA epi stroma tutte le cellule
##############################
load("/home/aurora.savino/crossWGCNA/averagedEpiTNBC.RData")
load("/home/aurora.savino/crossWGCNA/averagedEpiHER2.RData")
load("/home/aurora.savino/crossWGCNA/averagedEpiER.RData")
load("/home/aurora.savino/crossWGCNA/averagedStromaTNBC.RData")
load("/home/aurora.savino/crossWGCNA/averagedStromaHER2.RData")
load("/home/aurora.savino/crossWGCNA/averagedStromaER.RData")
load("~/crossWGCNA/TNmerge.RData")

setwd(paste("results/",at,b,ct,pv, sep=""))

library(WGCNA)
library(limma)
library(pheatmap)
library(Hmisc)
averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedStroma<-cbind(averagedStromaTNBC, averagedStromaHER2, averagedStromaER)

rownames(averagedEpi)<-rownames(TNmerge)
rownames(averagedStroma)<-rownames(TNmerge)

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedStroma))==0)

averagedEpi<-averagedEpi[, ind]
averagedStroma<-averagedStroma[, ind]

genes<-intersect(c(which(apply(averagedStroma,1,var)>=quantile(apply(averagedStroma,1,var),0.5))),
                 c(which(apply(averagedEpi,1,var)>=quantile(apply(averagedEpi,1,var),0.5))))
stroma<-averagedStroma[genes,]
epi<-averagedEpi[genes,]


#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "1", sep="_")
rownames(epi)<-paste(rownames(epi), "2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged<-rbind(stroma, epi)

#Functions da A_L.R
source("../../scripts/crossWGCNA_functions_netdiff.R")

degsc<-network(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2")
save(degsc, file="degsc_netdiff.RData")

source("../../scripts/crossWGCNA_functions_selfloops.R")

degsc<-network(data=data_merged, Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2")
save(degsc, file="degsc_selfloops.RData")

