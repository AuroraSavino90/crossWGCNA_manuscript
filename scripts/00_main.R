dir<-"path_to_workdir"# directory with "data","scripts" and "results" folders within
setwd(dir)

##Setting the parameters
at<-"signed"#adjacency type
ct<-"spearman"#correlation type
b<-6#beta

dir.create(paste("results/",at,b,ct, sep=""))

#runs crossWGCNA to compute inter-tissue degrees for six LCM datasets
#the parameter m defines the method to compute the adjacency
source("scripts/01_degs.R")
LCM_degrees(at,b,ct,m="netdiff", dir)
LCM_degrees(at,b,ct,m="selfloop", dir)

#GSEA of kRatio-ranked genes for the six LCM datasets
#the parameter m defines the method used to compute the adjacency
source("scripts/02_GSEA.R")
LCM_GSEA(at,b,ct,m="netdiff", dir)
LCM_GSEA(at,b,ct,m="selfloop", dir)

#compute crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset
source("scripts/GSE161529_degs.R")
GSE161529_degs(at,b,ct, dir)

##test whether genes with high kRatio in LCM are differentially expressed upon cocultures more often than by chance.
source("scripts/03_cocultures.R")
LCM_vs_cocultures(at,b,ct,m="netdiff", dir)
LCM_vs_cocultures(at,b,ct,m="selfloop", dir)

#compare LCM and single-cell kRatios
#the parameter m defines the method used to compute the adjacency
source("scripts/04_LCM_vs_singlecell.R")
LCM_vs_sc(at,b,ct,m="netdiff", dir)
LCM_vs_sc(at,b,ct,m="selfloop", dir)

#compare the two alternative algorithms
#computing the number of identical genes measured in either stroma or epi within the same co-expression module
#and comparing the correlation between genes' adjacencies within a tissue or inter-tissues
source("scripts/05_compare_versions.R")
compare_algorithms(at,b,ct,dir)

#compare the coherence of kratios between different LCM datasets
source("scripts/06_datasets_coherence.R")
LCM_coherence(at,b,ct,dir)

#apply crossWGCNA to inter-regulatory-layer interactions
source("scripts/07_mrnavsprot.R")
mrna_vs_prot(at,b,ct,m="netdiff", dir)
