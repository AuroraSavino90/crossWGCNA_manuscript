LCM_vs_cocultures<-function(at,b,ct,m, dir){
require(ggplot2)
require(fgsea)
require(data.table)

setwd(dir)

load("data/CAFs_cocultures_small.RData")
setwd(paste("results/",at,b,ct, sep=""))

load(paste("degs_3rd_", m, ".RData", sep=""))

names(degs)<-c("GSE5847", "GSE10797", "GSE14548", "GSE83591",
               "GSE68744", "GSE88715")

###GSEA
names(DE_stroma)<-c(1:length(DE_stroma))
gsea_stroma<-matrix(nrow=length(DE_stroma), ncol=length(degs))

for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
  names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))

  fgseaRes <- fgseaMultilevel(DE_stroma, coef_gsea)
  gsea_stroma[,j]<-unlist(fgseaRes[,2])

}

names(DE_epi)<-c(1:length(DE_epi))
gsea_epi<-matrix(nrow=length(DE_epi), ncol=length(degs))

for(j in 1:length(degs)){
  coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
  names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))

  fgseaRes <- fgseaMultilevel(DE_epi, coef_gsea)
  gsea_epi[,j]<-unlist(fgseaRes[,2])
 if(j==1){
   pdf(paste("epi_gsea_example_", m, ".pdf", sep=""),4,4)
   print(plotEnrichment(DE_epi[[3]], coef_gsea))
   dev.off()
 }
}

set.seed(46956305)
random_repeated_stroma<-c()
for(i in 1:100){
  names(DE_stroma)<-c(1:length(DE_stroma))
  gsea_stroma_rand<-matrix(nrow=length(DE_stroma), ncol=length(degs))

  for(j in 1:length(degs)){
    coef_gsea<-degs[[j]]$kExt1/degs[[j]]$kInt1
    names(coef_gsea)<-gsub("_tis1", "", names(coef_gsea))
    names(coef_gsea)<-sample(names(coef_gsea), length(names(coef_gsea)), replace = F)

    fgseaRes <- fgseaMultilevel(DE_stroma, coef_gsea)
    gsea_stroma_rand[,j]<-unlist(fgseaRes[,2])

  }

  random_repeated_stroma<-c(random_repeated_stroma, sum(gsea_stroma_rand<0.05, na.rm = T))
}

data <- data.frame(
  name=c("observed", "random"),
  value=c(sum(gsea_stroma<0.05, na.rm=T), median(random_repeated_stroma, na.rm=T)),
  sd=c(0, sd(random_repeated_stroma))
)

pdf(paste("stroma_gsea_signif_", m, ".pdf", sep=""), 3,5)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="darkgreen") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()

save(random_repeated_stroma, file=paste("random_repeated_stroma_", m, ".RData", sep=""))
load(file=paste("random_repeated_stroma_", m, ".RData", sep=""))

set.seed(46956305)
random_repeated_epi<-c()
for(i in 1:100){
  names(DE_epi)<-c(1:length(DE_epi))
  gsea_epi_rand<-matrix(nrow=length(DE_epi), ncol=length(degs))

  for(j in 1:length(degs)){
    coef_gsea<-degs[[j]]$kExt2/degs[[j]]$kInt2
    names(coef_gsea)<-gsub("_tis2", "", names(coef_gsea))
    names(coef_gsea)<-sample(names(coef_gsea), length(names(coef_gsea)), replace = F)

    fgseaRes <- fgseaMultilevel(DE_epi, coef_gsea)
    gsea_epi_rand[,j]<-unlist(fgseaRes[,2])

  }

  random_repeated_epi<-c(random_repeated_epi, sum(gsea_epi_rand<0.05, na.rm = T))
}#5

save(random_repeated_epi, file=paste("random_repeated_epi_", m, ".RData", sep=""))
load(file=paste("random_repeated_epi_", m, ".RData", sep=""))

data <- data.frame(
  name=c("observed", "random"),
  value=c(sum(gsea_epi<0.05, na.rm=T), median(random_repeated_epi, na.rm=T)),
  sd=c(0, sd(random_repeated_epi))
)

pdf(paste("epi_gsea_signif_", m, ".pdf", sep=""), 3,5)
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="brown2") +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, size=1.3) + theme_bw()
dev.off()
}
