#dir<-"path_to_workdir"# directory with "data","scripts" and "results" folders within
setwd(dir)

dir.create(paste("results/",at,b,ct, sep=""))

#compute crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset
source("scripts/GSE161529_degs_clust.R")
sc_15_01<-GSE161529_degs_clust(at,b,ct, dir, d=15, res=0.1)
sc_15_05<-GSE161529_degs_clust(at,b,ct, dir, d=15, res=0.5)
sc_15_1<-GSE161529_degs_clust(at,b,ct, dir, d=15, res=1)
sc_30_01<-GSE161529_degs_clust(at,b,ct, dir, d=30, res=0.1)
sc_30_05<-GSE161529_degs_clust(at,b,ct, dir, d=30, res=0.5)
sc_30_1<-GSE161529_degs_clust(at,b,ct, dir, d=30, res=1)
sc_50_01<-GSE161529_degs_clust(at,b,ct, dir, d=50, res=0.1)
sc_50_05<-GSE161529_degs_clust(at,b,ct, dir, d=50, res=0.5)
sc_50_1<-GSE161529_degs_clust(at,b,ct, dir, d=50, res=1)


kRatioE_15_01<-sc_15_01$kExt1/sc_15_01$kInt1
kRatioS_15_01<-sc_15_01$kExt2/sc_15_01$kInt2
kRatioE_15_05<-sc_15_05$kExt1/sc_15_05$kInt1
kRatioS_15_05<-sc_15_05$kExt2/sc_15_05$kInt2
kRatioE_15_1<-sc_15_1$kExt1/sc_15_1$kInt1
kRatioS_15_1<-sc_15_1$kExt2/sc_15_1$kInt2
kRatioE_30_01<-sc_30_01$kExt1/sc_30_01$kInt1
kRatioS_30_01<-sc_30_01$kExt2/sc_30_01$kInt2
kRatioE_30_05<-sc_30_05$kExt1/sc_30_05$kInt1
kRatioS_30_05<-sc_30_05$kExt2/sc_30_05$kInt2
kRatioE_30_1<-sc_30_1$kExt1/sc_30_1$kInt1
kRatioS_30_1<-sc_30_1$kExt2/sc_30_1$kInt2
kRatioE_50_01<-sc_50_01$kExt1/sc_50_01$kInt1
kRatioS_50_01<-sc_50_01$kExt2/sc_50_01$kInt2
kRatioE_50_05<-sc_50_05$kExt1/sc_50_05$kInt1
kRatioS_50_05<-sc_50_05$kExt2/sc_50_05$kInt2
kRatioE_50_1<-sc_50_1$kExt1/sc_50_1$kInt1
kRatioS_50_1<-sc_50_1$kExt2/sc_50_1$kInt2

genesE<-Reduce(intersect, list(names(kRatioE_15_01), names(kRatioE_15_05), names(kRatioE_15_1),
                               names(kRatioE_30_01), names(kRatioE_30_05), names(kRatioE_30_1),
                               names(kRatioE_50_01), names(kRatioE_50_05), names(kRatioE_50_1)))
genesS<-Reduce(intersect, list(names(kRatioS_15_01), names(kRatioS_15_05), names(kRatioS_15_1),
                               names(kRatioS_30_01), names(kRatioS_30_05), names(kRatioS_30_1),
                               names(kRatioS_50_01), names(kRatioS_50_05), names(kRatioS_50_1)))

matE<-cbind(kRatioE_15_01[genesE], kRatioE_15_05[genesE], kRatioE_15_1[genesE],
            kRatioE_30_01[genesE], kRatioE_30_05[genesE], kRatioE_30_1[genesE],
            kRatioE_50_01[genesE], kRatioE_50_05[genesE], kRatioE_50_1[genesE])
matS<-cbind(kRatioS_15_01[genesS], kRatioS_15_05[genesS], kRatioS_15_1[genesS],
            kRatioS_30_01[genesS], kRatioS_30_05[genesS], kRatioS_30_1[genesS],
            kRatioS_50_01[genesS], kRatioS_50_05[genesS], kRatioS_50_1[genesS])

colnames(matE)<-c("PC 15; res 0.1", "PC 15; res 0.5", "PC 15; res 1",
                  "PC 30; res 0.1", "PC 30; res 0.5", "PC 30; res 1",
                  "PC 50; res 0.1", "PC 50; res 0.5", "PC 50; res 1")
colnames(matS)<-c("PC 15; res 0.1", "PC 15; res 0.5", "PC 15; res 1",
                  "PC 30; res 0.1", "PC 30; res 0.5", "PC 30; res 1",
                  "PC 50; res 0.1", "PC 50; res 0.5", "PC 50; res 1")

cor(matE)
cor(matS)

library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <- c(seq( 0.6, 1, length.out=paletteLength))

graphics.off()
pdf(paste("cor_matE.pdf", sep=""),10,10)
pheatmap(cor(matE),   cellwidth=15, cellheight=15,  keep.dendro=T, color = myColor, breaks = myBreaks)
dev.off()

pdf(paste("cor_matS.pdf", sep=""),10,10)
pheatmap(cor(matS),   cellwidth=15, cellheight=15,  keep.dendro=T, color = myColor, breaks = myBreaks)
dev.off()

summary(cor(matE)[upper.tri(cor(matE))])

pheatmap(cor(matE))
