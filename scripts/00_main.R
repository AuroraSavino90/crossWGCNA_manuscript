setwd("~/crossWGCNA/")
dir<-"/home/aurora.savino/crossWGCNA/"

##Setting the parameters
at<-"signed"#adjacency type
ct<-"spearman"#correlation type
pv<-0.05#pvaluje
b<-6#beta

dir.create(paste("results/",at,b,ct,pv, sep=""))

#runs crossWGCNA to compute inter-tissue degrees for six LCM datasets
#the parameter m defines the method to compute the adjacency
source("scripts/01_degs.R")
LCM_degrees(at,b,ct,pv,m="netdiff", dir)
LCM_degrees(at,b,ct,pv,m="selfloop", dir)

#GSEA of kRatio-ranked genes for the six LCM datasets
#the parameter m defines the method used to compute the adjacency
source("scripts/02_GSEA.R")
LCM_GSEA(at,b,ct,pv,m="netdiff", dir)
LCM_GSEA(at,b,ct,pv,m="selfloop", dir)

#compute crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset
source("scripts/GSE161529_degs.R")

##test whether genes with high kRatio in LCM are differentially expressed upon cocultures more often than by chance.
source("scripts/03_cocultures.R")
LCM_vs_cocultures(at,b,ct,pv,m="netdiff", dir)
LCM_vs_cocultures(at,b,ct,pv,m="selfloop", dir)

#compare LCM and single-cell kRatios
#the parameter m defines the method used to compute the adjacency
source("scripts/04_LCM_vs_singlecell.R")
LCM_vs_sc(at,b,ct,pv,m="netdiff", dir)
LCM_vs_sc(at,b,ct,pv,m="selfloop", dir)

#compare the two alternative algorithms
#computing the number of identical genes measured in either stroma or epi within the same co-expression module
#and comparing the correlation between genes' adjacencies within a tissue or inter-tissues
source("scripts/05_compare_versions.R")
compare_algorithms(at,b,ct,pv,dir)

#compare the coherence of kratios between different LCM datasets
source("scripts/06_datasets_coherence.R")
LCM_coherence(at,b,ct,pv,dir)

