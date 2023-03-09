
################################
######Double WGCNA epi stroma on all cells
##############################
setwd(dir)

#load the pseudo-bulk data
load(paste(dir, "/averagedEpiTNBC.RData", sep=""))
load(paste(dir, "/averagedEpiHER2.RData"))
load(paste(dir, "/averagedEpiER.RData"))
load(paste(dir, "/averagedStromaTNBC.RData"))
load(paste(dir, "/averagedStromaHER2.RData"))
load(paste(dir, "/averagedStromaER.RData"))
load(paste(dir, "/TNmerge.RData"))

setwd(paste("results/",at,b,ct,pv, sep=""))

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

#crossWGCNA functions
source("../../scripts/crossWGCNA_functions_all.R")

degsc<-network(data=data_merged, method="netdiff", Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2")
save(degsc, file="degsc_netdiff.RData")

degsc<-network(data=data_merged, method="selfloop", Adj_type="signed", cortype="pearson", pval="none", thr=0.05, beta=6, comp1="_1", comp2="_2")
save(degsc, file="degsc_selfloops.RData")

