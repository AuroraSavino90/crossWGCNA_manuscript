GSE161529_degs<-function(at,b,ct, dir){
setwd(dir)
  source("scripts/crossWGCNA.R")

#load the pseudo-bulk data
load(paste(dir, "/data/GSE161529/averagedEpiTNBC.RData", sep=""))
load(paste(dir, "/data/GSE161529/averagedEpiHER2.RData", sep=""))
load(paste(dir, "/data/GSE161529/averagedEpiER.RData", sep=""))
load(paste(dir, "/data/GSE161529/averagedStromaTNBC.RData", sep=""))
load(paste(dir, "/data/GSE161529/averagedStromaHER2.RData", sep=""))
load(paste(dir, "/data/GSE161529/averagedStromaER.RData", sep=""))

setwd(paste("results/",at,b,ct, sep=""))

averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER)
averagedStroma<-cbind(averagedStromaTNBC, averagedStromaHER2, averagedStromaER)

###filtering data
ind<-which(colSums(is.na(averagedEpi))==0 & colSums(is.na(averagedStroma))==0)

averagedEpi<-averagedEpi[, ind]
averagedStroma<-averagedStroma[, ind]

#set.seed(563930)
set.seed(5907596)
averagedEpi<-averagedEpi[, sample(1:ncol(averagedEpi), replace = F)]
averagedStroma<-averagedStroma[, sample(1:ncol(averagedStroma), replace = F)]

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
degsc<-crossWGCNA(data=data_merged, method="netdiff", doClusters=F, Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_1", comp2="_2", )
save(degsc, file="degsc_netdiff_rand.RData")

}
