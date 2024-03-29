LCM_GSEA<-function(at,b,ct,m, dir){
  require(openxlsx)
require(pheatmap)
require(msigdbr)
require(fgsea)
require(data.table)
require(org.Hs.eg.db)

setwd(dir)
source("scripts/crossWGCNA.R")
setwd(paste("results/",at,b,ct, sep=""))

###load degrees obtained from 01_degs.R
load(paste("degs_3rd_", m, ".RData", sep=""))

names(degs) <- c("GSE5847",
                 "GSE10797",
                 "GSE14548",
                 "GSE83591",
                 "GSE68744",
                 "GSE88715")

###GSEA of kRatio-ranked genes (epi or stroma) for each LCM dataset
m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

for (j in 1:length(degs)) {
  coef_gsea <- degs[[j]]$kExt1 / degs[[j]]$kInt1
  names(coef_gsea) <- gsub("_tis1", "", names(coef_gsea))

  fgseaRes <- fgseaMultilevel(m_list, coef_gsea)
  fgseaRes <- fgseaRes[fgseaRes$padj < 0.05, ]
  if (nrow(fgseaRes) > 0) {
    write.xlsx(
      as.data.frame(fgseaRes[, c(1, 2, 3, 6, 7)]),
      paste("fgsea_ratio_stroma_", m,"_", names(degs)[j], ".xlsx", sep = "")
    )
  }

  coef_gsea <- degs[[j]]$kExt2 / degs[[j]]$kInt2
  names(coef_gsea) <- gsub("_tis2", "", names(coef_gsea))

  fgseaRes2 <- fgseaMultilevel(m_list, coef_gsea)
  fgseaRes2 <- fgseaRes2[fgseaRes2$padj < 0.05, ]
  if (nrow(fgseaRes2) > 0) {
    write.xlsx(
      as.data.frame(fgseaRes2[, c(1, 2, 3, 6, 7)]),
      paste("fgsea_ratio_epi_", m,"_", names(degs)[j], ".xlsx", sep = "")
    )
  }
}


m_df = msigdbr(species = "Homo sapiens", category = "C5")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes <- list()
fgseaRes2 <- list()
for (j in 1:length(degs)) {
  coef_gsea <- degs[[j]]$kExt1 / degs[[j]]$kInt1
  names(coef_gsea) <- gsub("_tis1", "", names(coef_gsea))

  fgseaRes[[j]] <- fgseaMultilevel(m_list, coef_gsea)

  coef_gsea <- degs[[j]]$kExt2 / degs[[j]]$kInt2
  names(coef_gsea) <- gsub("_tis2", "", names(coef_gsea))

  fgseaRes2[[j]] <- fgseaMultilevel(m_list, coef_gsea)

}


##############################################################
####### compare GSEA results for all stroma-based kRatios
##############################################################

allpaths <- c()
for (i in 1:length(fgseaRes)) {
  allpaths <-
    union(allpaths, fgseaRes[[i]]$pathway[fgseaRes[[i]]$pval < 0.05])
}

allpaths_mat <-
  matrix(0, nrow = length(allpaths), ncol = length(fgseaRes))
rownames(allpaths_mat) <- allpaths

for (i in 1:length(fgseaRes)) {
  allpaths_mat[fgseaRes[[i]]$pathway[which(fgseaRes[[i]]$pval < 0.05)], i] <-
    1
}

###plot heatmaps considering different thresholds of number of LCM datasets with the same GO term enriched

shared_paths <- allpaths_mat[rowSums(allpaths_mat) > 3,]
allpaths_mat <-
  matrix(0, nrow = length(allpaths), ncol = length(fgseaRes))
rownames(allpaths_mat) <- allpaths
for (i in 1:length(fgseaRes)) {
  allpaths_mat[fgseaRes[[i]]$pathway[which(fgseaRes[[i]]$pval < 0.05)], i] <-
    fgseaRes[[i]]$NES[which(fgseaRes[[i]]$pval < 0.05)]
}

toplot <- allpaths_mat[rowSums(sign(allpaths_mat)) > 2, ]
colnames(toplot) <- names(degs)
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <-
  c(seq(0, max(unlist(toplot), na.rm = T), length.out = floor(paletteLength)))
length(myBreaks) == length(paletteLength) + 1

pdf(paste("summary_GSEA_stroma_", m, "_2.pdf", sep=""), 20, 30)
pheatmap(
  toplot,
  cellwidth = 10,
  cellheight = 10,
  breaks = myBreaks,
  color = myColor,
  keep.dendro = T
)
dev.off()

toplot <- allpaths_mat[rowSums(sign(allpaths_mat)) > 3, ]
colnames(toplot) <- names(degs)
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <-
  c(seq(0, max(unlist(toplot), na.rm = T), length.out = floor(paletteLength)))
length(myBreaks) == length(paletteLength) + 1

pdf(paste("summary_GSEA_stroma_", m, "_3.pdf",sep=""), 20, 10)
pheatmap(
  toplot,
  cellwidth = 10,
  cellheight = 10,
  breaks = myBreaks,
  color = myColor,
  keep.dendro = T
)
dev.off()

##############################################################
####### compare GSEA results for all epi-based kRatios
##############################################################
allpaths <- c()
for (i in 1:length(fgseaRes2)) {
  allpaths <-
    union(allpaths, fgseaRes2[[i]]$pathway[fgseaRes2[[i]]$pval < 0.05])
}

allpaths_mat <-
  matrix(0, nrow = length(allpaths), ncol = length(fgseaRes2))
rownames(allpaths_mat) <- allpaths

for (i in 1:length(fgseaRes2)) {
  allpaths_mat[fgseaRes2[[i]]$pathway[which(fgseaRes2[[i]]$pval < 0.05)], i] <-
    1
}

###plot heatmaps considering different thresholds of number of LCM datasets with the same GO term enriched

shared_paths <- allpaths_mat[rowSums(allpaths_mat) > 3,]
allpaths_mat <-
  matrix(0, nrow = length(allpaths), ncol = length(fgseaRes2))
rownames(allpaths_mat) <- allpaths
for (i in 1:length(fgseaRes2)) {
  allpaths_mat[fgseaRes2[[i]]$pathway[which(fgseaRes2[[i]]$pval < 0.05)], i] <-
    fgseaRes2[[i]]$NES[which(fgseaRes2[[i]]$pval < 0.05)]
}

toplot <- allpaths_mat[rowSums(sign(allpaths_mat)) > 3, ]
colnames(toplot) <- names(degs)

paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <-
  c(seq(0, max(unlist(toplot), na.rm = T), length.out = floor(paletteLength)))
length(myBreaks) == length(paletteLength) + 1

pdf(paste("summary_GSEA_epi_", m, "_3.pdf", sep=""), 20, 35)
pheatmap(
  toplot,
  cellwidth = 10,
  cellheight = 10,
  breaks = myBreaks,
  color = myColor,
  keep.dendro = T
)
dev.off()

toplot <- allpaths_mat[rowSums(sign(allpaths_mat)) > 4, ]
colnames(toplot) <- names(degs)

paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <-
  c(seq(0, max(unlist(toplot), na.rm = T), length.out = floor(paletteLength)))
length(myBreaks) == length(paletteLength) + 1

pdf(paste("summary_GSEA_epi_", m, "_4.pdf", sep=""), 20, 15)
pheatmap(
  toplot,
  cellwidth = 10,
  cellheight = 10,
  breaks = myBreaks,
  color = myColor,
  keep.dendro = T
)
dev.off()
}
