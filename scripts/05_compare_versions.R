
compare_algorithms<-function(at,b,ct,dir){
  setwd(dir)
  require(openxlsx)
  load("data/metadataAll.RData")
  load("data/datasetsAll.RData")
  source("scripts/crossWGCNA.R")
  setwd(paste("results/",at,b,ct, sep=""))

#Datasets
###dataset GSE5847

stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE5847"]][,stromaID]
epi<-datasetsAll[["GSE5847"]][,epiID]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

#merging stroma and epi
rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")
colnames(epi)<-colnames(stroma)

data_merged_GSE5847<-rbind(stroma, epi)


############ GSE10797
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE10797" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE10797"]][,stromaID]
epi<-datasetsAll[["GSE10797"]][,epiID]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]
##############All genes double WGCNA
###use WGCNA with twice the genes
rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE10797<-rbind(stroma, epi)

##############GSE14548
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE14548" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE14548"]][,stromaID]
epi<-datasetsAll[["GSE14548"]][,epiID]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE14548<-rbind(stroma, epi)


#############GSE83591
load("data/GSE83591.RData")
GSE83591_meta<-read.csv("data/GSE83591_metadata.txt", sep="\t")
GSE83591_meta<-t(GSE83591_meta)
colnames(GSE83591_meta)<-GSE83591_meta[1,]
GSE83591_meta<-GSE83591_meta[-1,]

GSE83591_meta2<-read.xlsx("data/GSE83591_Correspondence_LCM.xlsx",1)
GSE83591_meta2$ID<-unlist(strsplit(GSE83591_meta2$GSE83591, "_"))[seq(1,109*9,9)]

#reorder based on annotation file
GSE83591<-GSE83591[,match(GSE83591_meta2$ID, colnames(GSE83591))]
GSE83591_meta<-GSE83591_meta[match(GSE83591_meta2$ID, GSE83591_meta[,1]),]

GSE83591<-GSE83591[,which(GSE83591_meta2$Cy3 %in%c("TE", "TS"))]
GSE83591_meta<-GSE83591_meta[which(GSE83591_meta2$Cy3 %in%c("TE", "TS")),]
GSE83591_meta2<-GSE83591_meta2[which(GSE83591_meta2$Cy3 %in%c("TE", "TS")),]

GSE83591<-GSE83591[,-c(15,34,67)]
GSE83591_meta<-GSE83591_meta[-c(15,34,67),]
GSE83591_meta2<-GSE83591_meta2[-c(15,34,67),]

stroma<-GSE83591[,which(GSE83591_meta2$Cy3=="TS")]
epi<-GSE83591[,which(GSE83591_meta2$Cy3=="TE")]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE83591<-rbind(stroma, epi)

#######GSE68744
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stroma<-datasetsAll[["GSE68744"]][,stromaID]
epi<-datasetsAll[["GSE68744"]][,epiID]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE68744<-rbind(stroma, epi)

#######GSE88715
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE88715" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE88715"]][,stromaID]
epi<-datasetsAll[["GSE88715"]][,epiID]

###filtering data
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE88715<-rbind(stroma, epi)

#####################################################################
#######################compute crossWGCNA networks
#####################################################################

net_GSE5847<-crossWGCNA(data=data_merged_GSE5847, method="netdiff", Adj_type=at, cortype=ct, pval="none", beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE10797<-crossWGCNA(data=data_merged_GSE10797,method="netdiff", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE14548<-crossWGCNA(data=data_merged_GSE14548,method="netdiff", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE83591<-crossWGCNA(data=data_merged_GSE83591,method="netdiff", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE68744<-crossWGCNA(data=data_merged_GSE68744,method="netdiff", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE88715<-crossWGCNA(data=data_merged_GSE88715, method="netdiff",Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)

net_GSE5847_sl<-crossWGCNA(data=data_merged_GSE5847,method="selfloop", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE10797_sl<-crossWGCNA(data=data_merged_GSE10797,method="selfloop",  Adj_type=at, cortype=ct, pval="none", beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE14548_sl<-crossWGCNA(data=data_merged_GSE14548,method="selfloop",  Adj_type=at, cortype=ct, pval="none", beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE83591_sl<-crossWGCNA(data=data_merged_GSE83591, method="selfloop", Adj_type=at, cortype=ct, pval="none", beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE68744_sl<-crossWGCNA(data=data_merged_GSE68744, method="selfloop", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)
net_GSE88715_sl<-crossWGCNA(data=data_merged_GSE88715, method="selfloop", Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2", doClusters=TRUE, ds=2)


prop_same_genes<-function(net){
prop_same_all<-c()
for(i in 1:length(unique(net[[2]][[1]]))){
  clust<-unique(net[[2]][[1]])[i]
  clust_genes<-names(which(net[[2]][[1]]==clust))
  prop_same<-sum(table(gsub("_tis1|_tis2", "", clust_genes))==2)/length(clust_genes)
  prop_same_all<-c(prop_same_all, prop_same)
}
return(prop_same_all)
}


df_new<-data.frame(proportion=c(prop_same_genes(net_GSE5847),
prop_same_genes(net_GSE10797),
prop_same_genes(net_GSE14548),
prop_same_genes(net_GSE83591),
prop_same_genes(net_GSE68744),
prop_same_genes(net_GSE88715)),
dataset=c(rep("GSE5847", length(prop_same_genes(net_GSE5847))),
rep("GSE10797", length(prop_same_genes(net_GSE10797))),
rep("GSE14548", length(prop_same_genes(net_GSE14548))),
rep("GSE83591", length(prop_same_genes(net_GSE83591))),
rep("GSE68744", length(prop_same_genes(net_GSE68744))),
rep("GSE88715", length(prop_same_genes(net_GSE88715)))
),
algorithm=rep("indirect loops", length(c(prop_same_genes(net_GSE5847),
prop_same_genes(net_GSE10797),
prop_same_genes(net_GSE14548),
prop_same_genes(net_GSE83591),
prop_same_genes(net_GSE68744),
prop_same_genes(net_GSE88715)))))

df_sl<-data.frame(proportion=c(prop_same_genes(net_GSE5847_sl),
prop_same_genes(net_GSE10797_sl),
prop_same_genes(net_GSE14548_sl),
prop_same_genes(net_GSE83591_sl),
prop_same_genes(net_GSE68744_sl),
prop_same_genes(net_GSE88715_sl)),
dataset=c(rep("GSE5847", length(prop_same_genes(net_GSE5847_sl))),
rep("GSE10797", length(prop_same_genes(net_GSE10797_sl))),
rep("GSE14548", length(prop_same_genes(net_GSE14548_sl))),
rep("GSE83591", length(prop_same_genes(net_GSE83591_sl))),
rep("GSE68744", length(prop_same_genes(net_GSE68744_sl))),
rep("GSE88715", length(prop_same_genes(net_GSE88715_sl)))
),
algorithm=rep("self loops", length(c(prop_same_genes(net_GSE5847_sl),
prop_same_genes(net_GSE10797_sl),
prop_same_genes(net_GSE14548_sl),
prop_same_genes(net_GSE83591_sl),
prop_same_genes(net_GSE68744_sl),
prop_same_genes(net_GSE88715_sl)))))

df_tot<-rbind.data.frame(df_sl, df_new)

pdf("algorithm_compare.pdf", 10, 5)
ggplot(df_tot, aes(x=dataset, y=proportion, fill=algorithm))+geom_boxplot()+theme_classic()
dev.off()


#######
##Compare adjacencies: one example
######

adj_GSE5847<-Adjacency(data=data_merged_GSE5847, method="netdiff", Adj_type=at, cortype=ct, pval="none", thr=pv, beta=b, comp1="_tis1", comp2="_tis2")
adj_GSE5847_sl<-Adjacency(data=data_merged_GSE5847, method="selfloop", Adj_type=at, cortype=ct, pval="none", thr=pv, beta=b, comp1="_tis1", comp2="_tis2")

##the correlation of a gene with the others between and within tissues is similar
#because in both tissues it is involved in the same pathways

i<-10
corr_tis1<-adj_GSE5847_sl[i,grep("tis1", colnames(adj_GSE5847_sl))][-i]
corr_tis2<-adj_GSE5847_sl[i,grep("tis2", colnames(adj_GSE5847_sl))][-i]
cor(corr_tis1, corr_tis2)
pdf("Example_cor_selfloops.pdf",5,5)
plot(corr_tis1, corr_tis2, pch=19)
dev.off()

cor_tot<-c()
for(i in 1:(nrow(adj_GSE5847_sl)/2)){
  corr_tis1<-adj_GSE5847_sl[i,grep("tis1", colnames(adj_GSE5847_sl))][-i]
corr_tis2<-adj_GSE5847_sl[i,grep("tis2", colnames(adj_GSE5847_sl))][-i]
cor_tot<-c(cor_tot, cor(corr_tis1, corr_tis2))
}

pdf("All_cor_selfloops.pdf",5,5)
hist(cor_tot)
dev.off()



i<-10
corr_tis1<-adj_GSE5847[i,grep("tis1", colnames(adj_GSE5847))][-i]
corr_tis2<-adj_GSE5847[i,grep("tis2", colnames(adj_GSE5847))][-i]
cor(corr_tis1, corr_tis2)
pdf("Example_cor_netdiff.pdf",5,5)
plot(corr_tis1, corr_tis2, pch=19)
dev.off()

cor_tot_alt<-c()
for(i in 1:(nrow(adj_GSE5847)/2)){
  corr_tis1<-adj_GSE5847[i,grep("tis1", colnames(adj_GSE5847))][-i]
corr_tis2<-adj_GSE5847[i,grep("tis2", colnames(adj_GSE5847))][-i]
cor_tot_alt<-c(cor_tot_alt, cor(corr_tis1, corr_tis2))
}

pdf("All_cor_netdiff.pdf",5,5)
hist(cor_tot_alt)
dev.off()
}
