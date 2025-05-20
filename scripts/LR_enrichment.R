library(fgsea)
setwd(dir)
load(paste("results/",at,b,ct,"/degs_3rd_netdiff.RData", sep=""))
lr<-read.csv("/data/human_lr_pair.txt", sep="\t")
l<-unique(lr$ligand_gene_symbol)
r<-unique(lr$receptor_gene_symbol)

m_list<-list(lr=unique(c(l,r)))

fgseaRes <- list()
fgseaRes2 <- list()
for (j in 1:length(degs)) {
  coef_gsea <- degs[[j]]$kExt1 / degs[[j]]$kInt1
  names(coef_gsea) <- gsub("_tis1", "", names(coef_gsea))
  if(length(which(names(coef_gsea)==""))>0){
  coef_gsea<-coef_gsea[-which(names(coef_gsea)=="")]
  }
  fgseaRes[[j]] <- fgseaMultilevel(m_list, coef_gsea)

  coef_gsea <- degs[[j]]$kExt2 / degs[[j]]$kInt2
  names(coef_gsea) <- gsub("_tis2", "", names(coef_gsea))
  if(length(which(names(coef_gsea)==""))>0){
    coef_gsea<-coef_gsea[-which(names(coef_gsea)=="")]
  }
  fgseaRes2[[j]] <- fgseaMultilevel(m_list, coef_gsea)

}




library(metap)

meta_dn_rev<-function(x){
  istwo <- rep(T, length(x))
  toinvert <- ifelse(sign(unlist(lapply(x, function(x){x[1,6]})))==(-1), T, F)
  missing<-which(is.na(toinvert))
  if(length(missing)==0){
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]})), two = istwo, invert = toinvert))
  } else {
    p<-sumlog(two2one(unlist(lapply(x, function(x){x[1,3]}))[-missing], two = istwo[-missing], invert = toinvert[-missing]))

  }

  return(p[[3]])
}

p1<-meta_dn_rev(fgseaRes)
p2<-meta_dn_rev(fgseaRes2)

save(list(p1,p2), file=paste("results/",at,b,ct,"/LR_enrichment.RData", sep=""))
