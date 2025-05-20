require(openxlsx)
require(ggplot2)
require(ggpubr)

require(msigdbr)
require(fgsea)
require(data.table)
##Setting the parameters
at<-"signed"#adjacency type
ct<-"spearman"#correlation type
b<-6#beta
m<-"netdiff"

dir<-"path_to_workdir"# directory with "data","scripts" and "results" folders within
setwd(dir)

setwd(paste("results/",at,b,ct, sep=""))
#load degrees computed in 01_degs.R and GSE161529_degs.R
load(paste("degs_3rd_", m, ".RData", sep=""))
load(paste("degsc_", m, ".RData", sep=""))
degsc_orig<-degsc

load(paste("degsc_", m, "_rand.RData", sep=""))

names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")

########STROMA#####################################################
###rank product of the 6 networks
ranked<-list()
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))
  ranked[[j]]<-rank(-coef_gsea)
}

incommon<-Reduce(intersect, lapply(ranked, names))

rankings<-matrix(nrow=length(incommon), ncol=length(degs))
rownames(rankings)<-incommon
for(j in 1:length(degs)){
  rankings[,j]<-ranked[[j]][incommon]
}
rankprod<-apply(rankings,1,prod)
rankprod<-data.frame(gene=incommon, rankprod=rankprod)
rownames(rankprod)<-rankprod[,1]


#comparison of LCM rankproduct with single cell kRatio
coef_gsea_sc<-degsc_orig$kExt1/degsc_orig$kInt1
names(coef_gsea_sc)<-gsub("_1", "", names(coef_gsea_sc))

inboth<-intersect(rownames(rankprod), names(coef_gsea_sc))


rankprodinboth<-rankprod[inboth,2]
names(rankprodinboth)<-inboth
top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]


coef_gsea_sc_rand<-degsc$kExt1/degsc$kInt1
names(coef_gsea_sc_rand)<-gsub("_1", "", names(coef_gsea_sc_rand))

df<-data.frame(sc=c(coef_gsea_sc[top], coef_gsea_sc_rand[top]), dir=factor(c(rep("Original", length(top)), rep("Random", length(bottom)))))

my_comparisons <- list( c("Original", "Random"))
pdf(paste("singlecell_boxplot_stroma_",m,"_rand.pdf",sep=""),3,5)
ggboxplot(df, x="dir", y="sc")+stat_compare_means(comparisons=my_comparisons,label = "p.signif")+ylab("Comm. score single cell")+xlab("")
dev.off()


####################EPI##############################
###rank product of the 6 networks
ranked<-list()
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))
  ranked[[j]]<-rank(-coef_gsea)
}

incommon<-Reduce(intersect, lapply(ranked, names))

rankings<-matrix(nrow=length(incommon), ncol=length(degs))
rownames(rankings)<-incommon
for(j in 1:length(degs)){
  rankings[,j]<-ranked[[j]][incommon]
}
rankprod<-apply(rankings,1,prod)
rankprod<-data.frame(gene=incommon, rankprod=rankprod)
rownames(rankprod)<-rankprod[,1]


coef_gsea_sc<-degsc_orig$kExt2/degsc_orig$kInt2
names(coef_gsea_sc)<-gsub("_2", "", names(coef_gsea_sc))

#comparison of LCM rankproduct with single cell kRatio
inboth<-intersect(rownames(rankprod), names(coef_gsea_sc))

rankprodinboth<-rankprod[inboth,2]
names(rankprodinboth)<-inboth
top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]

coef_gsea_sc_rand<-degsc$kExt2/degsc$kInt2
names(coef_gsea_sc_rand)<-gsub("_2", "", names(coef_gsea_sc_rand))

df<-data.frame(sc=c(coef_gsea_sc[top], coef_gsea_sc_rand[top]), dir=factor(c(rep("Original", length(top)), rep("Random", length(bottom)))))

my_comparisons <- list( c("Original", "Random"))
pdf(paste("singlecell_boxplot_epi_",m, "_rand.pdf", sep=""),3,5)
print(ggboxplot(df, x="dir", y="sc")+stat_compare_means(comparisons=my_comparisons,label = "p.signif")+ylab("Comm. score single cell")+xlab(""))
dev.off()

