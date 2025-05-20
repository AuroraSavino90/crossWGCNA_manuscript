GSE161529_degs_clust<-function(at,b,ct, dir, d, res){
setwd(dir)
  source("scripts/crossWGCNA.R")

#load the pseudo-bulk data
  load(paste("results/averagedEpiTNBC",d,res,".RData", sep="_"))
  load(paste("results/averagedEpiHER2",d,res,".RData", sep="_"))
  load(paste("results/averagedEpiER",d,res,".RData", sep="_"))
  load(paste("results/averagedStromaTNBC",d,res,".RData", sep="_"))
  load(paste("results/averagedStromaHER2",d,res,".RData", sep="_"))
  load(paste("results/averagedStromaER",d,res,".RData", sep="_"))

  load(file="GSE161529_RAW/TNmerge.RData")


setwd(paste("results/",at,b,ct, sep=""))

averagedEpi<-cbind(averagedEpiTNBC, averagedEpiHER2, averagedEpiER[,-20])
averagedStroma<-cbind(averagedStromaTNBC, averagedStromaHER2, averagedStromaER[,-20])

rownames(averagedEpi)<-rownames(TNmerge)
rownames(averagedStroma)<-rownames(TNmerge)

###filtering data

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
save(degsc, file=paste("degsc_netdiff",d,res,".RData", sep="_"))
return(degsc)

}
