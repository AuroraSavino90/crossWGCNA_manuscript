LCM_vs_sc<-function(at,b,ct,pv,m, dir){
  require(openxlsx)
require(ggplot2)
require(ggpubr)

setwd(dir)

setwd(paste("results/",at,b,ct,pv, sep=""))

load(paste("degs_3rd_", m, ".RData", sep=""))
load(paste("degsc_", m, ".RData", sep=""))

names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")

########STROMA

###rank product dei 6 network
ranked<-list()
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))
  ranked[[j]]<-rank(-coef_gsea)
}

incommon<-Reduce(intersect, lapply(ranked, names))

tocompare<-matrix(nrow=length(incommon), ncol=6)
rownames(tocompare)<-incommon
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))
  tocompare[,j]<-coef_gsea[incommon]
}
cor(tocompare)

rankings<-matrix(nrow=length(incommon), ncol=length(degs))
rownames(rankings)<-incommon
for(j in 1:length(degs)){
  rankings[,j]<-ranked[[j]][incommon]
}
rankprod<-apply(rankings,1,prod)
rankprod<-data.frame(gene=incommon, rankprod=rankprod)
rownames(rankprod)<-rankprod[,1]

write.xlsx(rankprod, file=paste("rankprod_stroma_", m, ".xlsx", sep=""), rowNames=T)


coef_gsea_sc<-degsc$kExt1/degsc$kInt1
names(coef_gsea_sc)<-gsub("_1", "", names(coef_gsea_sc))

#confronto con rankprod
inboth<-intersect(rownames(rankprod), names(coef_gsea_sc))
cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth,2]))

df<-data.frame(rankproduct=log2(rankprod[inboth,2]), sc=coef_gsea_sc[inboth])
pdf(paste("singlecell_scatter_stroma_",m, ".pdf", sep=""),5,4)
print(ggplot(df, aes(x=rankproduct,y=sc))+geom_hex(bins=60)+ylab(label = "Comm. score single cell")+xlab(label = "Rank product")+geom_smooth(method="lm")+theme_classic())
dev.off()

write.csv(cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth,2]))[c(3,4)], paste("singlecell_scatter_stroma_",m, ".csv", sep=""))

library(ggpubr)
rankprodinboth<-rankprod[inboth,2]
names(rankprodinboth)<-inboth
top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]
df<-data.frame(sc=c(coef_gsea_sc[top], coef_gsea_sc[bottom]), dir=factor(c(rep("Top", length(top)), rep("Bottom", length(bottom)))))

my_comparisons <- list( c("Top", "Bottom"))
pdf(paste("singlecell_boxplot_stroma_",m,".pdf",sep=""),3,5)
ggboxplot(df, x="dir", y="sc")+stat_compare_means(comparisons=my_comparisons,label = "p.signif")+ylab("Comm. score single cell")+xlab("")
dev.off()

top50LCM<-rownames(rankprod)[(order(rankprod[,2], decreasing=F))][1:547]
top50sc<-names(sort(coef_gsea_sc, decreasing=T))[1:1519]
intersect(top50LCM, top50sc)
write.csv(intersect(top50LCM, top50sc), paste("stroma_", m, "_1decile.csv", sep=""))

####################EPI
###rank product dei 6 network
ranked<-list()
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))
  ranked[[j]]<-rank(-coef_gsea)
}

incommon<-Reduce(intersect, lapply(ranked, names))

tocompare<-matrix(nrow=length(incommon), ncol=6)
rownames(tocompare)<-incommon
for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))
  tocompare[,j]<-coef_gsea[incommon]
}
cor(tocompare)

rankings<-matrix(nrow=length(incommon), ncol=length(degs))
rownames(rankings)<-incommon
for(j in 1:length(degs)){
  rankings[,j]<-ranked[[j]][incommon]
}
rankprod<-apply(rankings,1,prod)
rankprod<-data.frame(gene=incommon, rankprod=rankprod)
rownames(rankprod)<-rankprod[,1]

write.xlsx(rankprod, file=paste("rankprod_epi_",m, ".xlsx", sep=""), rowNames=T)

coef_gsea_sc<-degsc$kExt2/degsc$kInt2
names(coef_gsea_sc)<-gsub("_2", "", names(coef_gsea_sc))

#confronto con rankprod
inboth<-intersect(rownames(rankprod), names(coef_gsea_sc))
cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth,2]))

df<-data.frame(rankproduct=log2(rankprod[inboth,2]), sc=coef_gsea_sc[inboth])
pdf(paste("singlecell_scatter_epi_",m,".pdf", sep=""),5,4)
print(ggplot(df, aes(x=rankproduct,y=sc))+geom_hex(bins=60)+ylab(label = "Comm. score single cell")+xlab(label = "Rank product")+geom_smooth(method="lm")+theme_classic())
dev.off()

write.csv(cor.test(coef_gsea_sc[inboth], log2(rankprod[inboth,2]))[c(3,4)], paste("singlecell_scatter_epi_",m, ".csv", sep=""))

library(ggpubr)
rankprodinboth<-rankprod[inboth,2]
names(rankprodinboth)<-inboth
top<-names(rankprodinboth)[which(rankprodinboth<quantile(rankprodinboth, c(0.1)))]
bottom<-names(rankprodinboth)[which(rankprodinboth>quantile(rankprodinboth, c(0.9)))]
df<-data.frame(sc=c(coef_gsea_sc[top], coef_gsea_sc[bottom]), dir=factor(c(rep("Top", length(top)), rep("Bottom", length(bottom)))))

my_comparisons <- list( c("Top", "Bottom"))
pdf(paste("singlecell_boxplot_epi_",m, ".pdf", sep=""),3,5)
print(ggboxplot(df, x="dir", y="sc")+stat_compare_means(comparisons=my_comparisons,label = "p.signif")+ylab("Comm. score single cell")+xlab(""))
dev.off()


top50LCM<-rownames(rankprod)[(order(rankprod[,2], decreasing=F))][1:547]
top50sc<-names(sort(coef_gsea_sc, decreasing=T))[1:1519]
intersect(top50LCM, top50sc)
write.csv(intersect(top50LCM, top50sc), paste("epi_", m, "_1decile.csv", sep=""))

###GSEA
library(msigdbr)
library(fgsea)
library(data.table)
library(openxlsx)

m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)


coef_gsea<-degsc$kExt1/degsc$kInt1
names(coef_gsea)<-gsub("_1", "", names(coef_gsea))

fgseaRes <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
fgseaRes<-fgseaRes[fgseaRes$padj<0.05,]
if(nrow(fgseaRes)>0){
  write.xlsx(as.data.frame(fgseaRes[,c(1,2,3,6,7)]), paste("fgsea_ratio_stroma_sc_",m,".xlsx", sep=""))
}

coef_gsea<-degsc$kExt2/degsc$kInt2
names(coef_gsea)<-gsub("_2", "", names(coef_gsea))

fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea, scoreType = "pos")
fgseaRes2<-fgseaRes2[fgseaRes2$padj<0.05,]
if(nrow(fgseaRes2)>0){
  write.xlsx(as.data.frame(fgseaRes2[,c(1,2,3,6,7)]), paste("fgsea_ratio_epi_sc_",m,".xlsx", sep=""))
}
}
