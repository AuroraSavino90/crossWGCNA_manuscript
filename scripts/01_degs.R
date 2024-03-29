LCM_degrees<-function(at,b,ct,m, dir){
setwd(dir)
require(openxlsx)
load("data/metadataAll.RData")
load("data/datasetsAll.RData")
source("scripts/crossWGCNA.R")
setwd(paste("results/",at,b,ct, sep=""))

########################
#Format data as necessary
##########################

###GSE5847
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

stromaMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiMatch<-metadataAll$matching[which(metadataAll$dataset=="GSE5847" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]

inboth<-intersect(stromaMatch, epiMatch)

stromaID<-stromaID[match(inboth, stromaMatch)]
epiID<-epiID[match(inboth, epiMatch)]

stroma<-datasetsAll[["GSE5847"]][,stromaID]
epi<-datasetsAll[["GSE5847"]][,epiID]

###filter genes
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

#merge stroma and epi
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

###filter genes
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

##merge stroma and epi
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

###filter genes
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

### merge stroma and epi
colnames(epi)<-colnames(stroma)
data_merged_GSE14548<-rbind(stroma, epi)


#############GSE83591
load("../../data/GSE83591.RData")
GSE83591_meta<-read.csv("../../data/GSE83591_metadata.txt", sep="\t")
GSE83591_meta<-t(GSE83591_meta)
colnames(GSE83591_meta)<-GSE83591_meta[1,]
GSE83591_meta<-GSE83591_meta[-1,]

GSE83591_meta2<-read.xlsx("../../data/GSE83591_Correspondence_LCM.xlsx",1)
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

###filter genes
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

##merge stroma and epi
colnames(epi)<-colnames(stroma)
data_merged_GSE83591<-rbind(stroma, epi)

#######GSE68744
stromaID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Stroma" & metadataAll$diseaseStatus=="InvasiveBC")]
epiID<-metadataAll$id[which(metadataAll$dataset=="GSE68744" & metadataAll$compartment=="Epi" & metadataAll$diseaseStatus=="InvasiveBC")]
#checked matching

stroma<-datasetsAll[["GSE68744"]][,stromaID]
epi<-datasetsAll[["GSE68744"]][,epiID]

###filter genes
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

###filter genes
genes<-unique(c(which(apply(stroma,1,var)>=quantile(apply(stroma,1,var),0.25)),
                which(apply(epi,1,var)>=quantile(apply(epi,1,var),0.25))))
stroma<-stroma[genes,]
epi<-epi[genes,]

rownames(stroma)<-paste(rownames(stroma), "tis1", sep="_")
rownames(epi)<-paste(rownames(epi), "tis2", sep="_")

colnames(epi)<-colnames(stroma)
data_merged_GSE88715<-rbind(stroma, epi)

#####################################################################
###############compute degrees for each LCM dataset
#####################################################################

Adj_GSE5847<-Adjacency(data=data_merged_GSE5847, method=m, Adj_type=at, cortype=ct, pval="none", beta=b, comp1="_tis1", comp2="_tis2")
degrees_GSE5847<-degrees(A=Adj_GSE5847, comp1="_tis1", comp2="_tis2")
Adj_GSE10797<-Adjacency(data=data_merged_GSE10797, method=m, Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2")
degrees_GSE10797<-degrees(A=Adj_GSE10797, comp1="_tis1", comp2="_tis2")
Adj_GSE14548<-Adjacency(data=data_merged_GSE14548, method=m, Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2")
degrees_GSE14548<-degrees(A=Adj_GSE14548, comp1="_tis1", comp2="_tis2")
Adj_GSE83591<-Adjacency(data=data_merged_GSE83591, method=m, Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2")
degrees_GSE83591<-degrees(A=Adj_GSE83591, comp1="_tis1", comp2="_tis2")
Adj_GSE68744<-Adjacency(data=data_merged_GSE68744, method=m, Adj_type=at, cortype=ct, pval="none",  beta=b, comp1="_tis1", comp2="_tis2")
degrees_GSE68744<-degrees(A=Adj_GSE68744, comp1="_tis1", comp2="_tis2")
Adj_GSE88715<-Adjacency(data=data_merged_GSE88715, method=m, Adj_type=at, cortype=ct, pval="none", beta=b, comp1="_tis1", comp2="_tis2")
degrees_GSE88715<-degrees(A=Adj_GSE88715, comp1="_tis1", comp2="_tis2")

degs<-list(degrees_GSE5847, degrees_GSE10797, degrees_GSE14548,
           degrees_GSE83591, degrees_GSE68744, degrees_GSE88715)

#save the list of degrees for all LCM datasets
save(degs, file=paste("degs_3rd_", m, ".RData", sep=""))
}
