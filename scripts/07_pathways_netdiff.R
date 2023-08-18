net_pathway<-function(at,b,ct,pv,m, dir){
  setwd(dir)
  require(openxlsx)
  load("data/metadataAll.RData")
  load("data/datasetsAll.RData")
  hippo <- read.xlsx("data/Pathways.xlsx", 2)[-1, 1]
  sign_hippo <- read.xlsx("data/Pathways.xlsx", 2)[-1, 3]

  source("scripts/crossWGCNA_functions_all.R")
  setwd(paste("results/",at,b,ct,pv, sep=""))


sign_hippo[which(sign_hippo == "OG")] <- 1
sign_hippo[which(sign_hippo == "TSG")] <- (-1)
sign_hippo[-which(sign_hippo %in% c("1", "-1"))] <- NA
hippo <- hippo[which(sign_hippo %in% c("1", "-1"))]
sign_hippo <- sign_hippo[which(sign_hippo %in% c("1", "-1"))]


  selgenes <- hippo
  sign <- as.numeric(sign_hippo)
  names(sign) <- selgenes


  ######################################
  ##### LCM datasets preprocessing
  ######################################
     #Datasets
    ###dataset GSE5847

    stromaID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE5847" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE5847" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    stromaMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE5847" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE5847" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    inboth <- intersect(stromaMatch, epiMatch)

    stromaID <- stromaID[match(inboth, stromaMatch)]
    epiID <- epiID[match(inboth, epiMatch)]

    stroma <- datasetsAll[["GSE5847"]][, stromaID]
    epi <- datasetsAll[["GSE5847"]][, epiID]

    ###filtering data
    genes <-
      unique(c(which(
        apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
      ),
      which(
        apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
      )))
    stroma <- stroma[genes, ]
    epi <- epi[genes, ]

    #merging stroma and epi
    rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
    rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")
    colnames(epi) <- colnames(stroma)

    data_merged_GSE5847 <- rbind(stroma, epi)


    ############ GSE10797
    stromaID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE10797" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE10797" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    stromaMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE10797" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE10797" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    inboth <- intersect(stromaMatch, epiMatch)

    stromaID <- stromaID[match(inboth, stromaMatch)]
    epiID <- epiID[match(inboth, epiMatch)]

    stroma <- datasetsAll[["GSE10797"]][, stromaID]
    epi <- datasetsAll[["GSE10797"]][, epiID]

    ###filtering data
    genes <-
      unique(c(which(
        apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
      ),
      which(
        apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
      )))
    stroma <- stroma[genes, ]
    epi <- epi[genes, ]
    ##############All genes double WGCNA
    ###use WGCNA with twice the genes
    rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
    rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

    colnames(epi) <- colnames(stroma)
    data_merged_GSE10797 <- rbind(stroma, epi)

    ##############GSE14548
    stromaID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE14548" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE14548" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    stromaMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE14548" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE14548" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    inboth <- intersect(stromaMatch, epiMatch)

    stromaID <- stromaID[match(inboth, stromaMatch)]
    epiID <- epiID[match(inboth, epiMatch)]

    stroma <- datasetsAll[["GSE14548"]][, stromaID]
    epi <- datasetsAll[["GSE14548"]][, epiID]

    ###filtering data
    genes <-
      unique(c(which(
        apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
      ),
      which(
        apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
      )))
    stroma <- stroma[genes, ]
    epi <- epi[genes, ]

    rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
    rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

    colnames(epi) <- colnames(stroma)
    data_merged_GSE14548 <- rbind(stroma, epi)


    #############GSE83591
    load("data/GSE83591.RData")
    GSE83591_meta <- read.csv("data/GSE83591_metadata.txt", sep = "\t")
    GSE83591_meta <- t(GSE83591_meta)
    colnames(GSE83591_meta) <- GSE83591_meta[1, ]
    GSE83591_meta <- GSE83591_meta[-1, ]

    GSE83591_meta2 <- read.xlsx("data/GSE83591_Correspondence_LCM.xlsx", 1)
    GSE83591_meta2$ID <-
      unlist(strsplit(GSE83591_meta2$GSE83591, "_"))[seq(1, 109 * 9, 9)]

    #reorder based on annotation file
    GSE83591 <- GSE83591[, match(GSE83591_meta2$ID, colnames(GSE83591))]
    GSE83591_meta <-
      GSE83591_meta[match(GSE83591_meta2$ID, GSE83591_meta[, 1]), ]

    GSE83591 <- GSE83591[, which(GSE83591_meta2$Cy3 %in% c("TE", "TS"))]
    GSE83591_meta <-
      GSE83591_meta[which(GSE83591_meta2$Cy3 %in% c("TE", "TS")), ]
    GSE83591_meta2 <-
      GSE83591_meta2[which(GSE83591_meta2$Cy3 %in% c("TE", "TS")), ]

    GSE83591 <- GSE83591[, -c(15, 34, 67)]
    GSE83591_meta <- GSE83591_meta[-c(15, 34, 67), ]
    GSE83591_meta2 <- GSE83591_meta2[-c(15, 34, 67), ]

    stroma <- GSE83591[, which(GSE83591_meta2$Cy3 == "TS")]
    epi <- GSE83591[, which(GSE83591_meta2$Cy3 == "TE")]

    ###filtering data
    genes <-
      unique(c(which(
        apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
      ),
      which(
        apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
      )))
    stroma <- stroma[genes, ]
    epi <- epi[genes, ]

    rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
    rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

    colnames(epi) <- colnames(stroma)
    data_merged_GSE83591 <- rbind(stroma, epi)

    #######GSE68744
    stromaID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE68744" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE68744" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    stroma <- datasetsAll[["GSE68744"]][, stromaID]
    epi <- datasetsAll[["GSE68744"]][, epiID]

    ###filtering data
    genes <-
      unique(c(which(
        apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
      ),
      which(
        apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
      )))
    stroma <- stroma[genes, ]
    epi <- epi[genes, ]

    rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
    rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")

    colnames(epi) <- colnames(stroma)
    data_merged_GSE68744 <- rbind(stroma, epi)

    #######GSE88715
    stromaID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE88715" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiID <-
      metadataAll$id[which(
        metadataAll$dataset == "GSE88715" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    stromaMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE88715" &
          metadataAll$compartment == "Stroma" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]
    epiMatch <-
      metadataAll$matching[which(
        metadataAll$dataset == "GSE88715" &
          metadataAll$compartment == "Epi" &
          metadataAll$diseaseStatus == "InvasiveBC"
      )]

    inboth <- intersect(stromaMatch, epiMatch)

    stromaID <- stromaID[match(inboth, stromaMatch)]
    epiID <- epiID[match(inboth, epiMatch)]

    stroma <- datasetsAll[["GSE88715"]][, stromaID]
    epi <- datasetsAll[["GSE88715"]][, epiID]

    ###filtering data
    genes <-
      unique(c(which(
        apply(stroma, 1, var) >= quantile(apply(stroma, 1, var), 0.25)
      ),
      which(
        apply(epi, 1, var) >= quantile(apply(epi, 1, var), 0.25)
      )))
    stroma <- stroma[genes, ]
    epi <- epi[genes, ]

    rownames(stroma) <- paste(rownames(stroma), "tis1", sep = "_")
    rownames(epi) <- paste(rownames(epi), "tis2", sep = "_")


    colnames(epi) <- colnames(stroma)
    data_merged_GSE88715 <- rbind(stroma, epi)


    ######################################################
    ######### compute degrees
    ########################################################

     degs_GSE5847 <-
      network(
        selgenes = selgenes,
        data = data_merged_GSE5847,
        method=m,
        Adj_type = "keep sign",
        cortype = "spearman",
        pval = "none",
        thr = 0.05,
        beta = 6,
        comp1 = "_tis1$",
        comp2 = "_tis2$",
        sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE5847)[grep("_tis2$", rownames(data_merged_GSE5847))])],
        compartment_sel = "comp2"
      )
    degs_GSE10797 <-
      network(
        selgenes = selgenes,
        data = data_merged_GSE10797,
        method=m,
        Adj_type = "keep sign",
        cortype = "spearman",
        pval = "none",
        thr = 0.05,
        beta = 6,
        comp1 = "_tis1$",
        comp2 = "_tis2$",
        sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE10797)[grep("_tis2$", rownames(data_merged_GSE10797))])],
        compartment_sel = "comp2"
      )
    degs_GSE14548 <-
      network(
        selgenes = selgenes,
        data = data_merged_GSE14548,
        method=m,
        Adj_type = "keep sign",
        cortype = "spearman",
        pval = "none",
        thr = 0.05,
        beta = 6,
        comp1 = "_tis1$",
        comp2 = "_tis2$",
        sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE14548)[grep("_tis2$", rownames(data_merged_GSE14548))])],
        compartment_sel = "comp2"
      )
    degs_GSE83591 <-
      network(
        selgenes = selgenes,
        data = data_merged_GSE83591,
        method=m,
        Adj_type = "keep sign",
        cortype = "spearman",
        pval = "none",
        thr = 0.05,
        beta = 6,
        comp1 = "_tis1$",
        comp2 = "_tis2$",
        sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE83591)[grep("_tis2$", rownames(data_merged_GSE83591))])],
        compartment_sel = "comp2"
      )
    degs_GSE68744 <-
      network(
        selgenes = selgenes,
        data = data_merged_GSE68744,
        method=m,
        Adj_type = "keep sign",
        cortype = "spearman",
        pval = "none",
        thr = 0.05,
        beta = 6,
        comp1 = "_tis1$",
        comp2 = "_tis2$",
        sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE68744)[grep("_tis2$", rownames(data_merged_GSE68744))])],
        compartment_sel = "comp2"
      )
    degs_GSE88715 <-
      network(
        selgenes = selgenes,
        data = data_merged_GSE88715,
        method=m,
        Adj_type = "keep sign",
        cortype = "spearman",
        pval = "none",
        thr = 0.05,
        beta = 6,
        comp1 = "_tis1$",
        comp2 = "_tis2$",
        sign_list = sign[gsub("_tis2$", "", rownames(data_merged_GSE88715)[grep("_tis2$", rownames(data_merged_GSE88715))])],
        compartment_sel = "comp2"
      )

    degs_path_sign <- list(
      degs_GSE5847,
      degs_GSE10797,
      degs_GSE14548,
      degs_GSE83591,
      degs_GSE68744,
      degs_GSE88715
    )

    save(degs_path_sign,
         file = "degs_hippo_KS_s.RData")


  incommon <- Reduce(intersect, list(
    names(degs_path_sign[[1]]$kExt1),
    names(degs_path_sign[[2]]$kExt1),
    names(degs_path_sign[[3]]$kExt1),
    names(degs_path_sign[[4]]$kExt1),
    names(degs_path_sign[[5]]$kExt1),
    names(degs_path_sign[[6]]$kExt1)
  ))


  require(biomaRt)
  ensembl.human <- useEnsembl(biomart = 'genes',
                              dataset = 'hsapiens_gene_ensembl')

  #extracellular region
  GO <- getBM(
    attributes = c('hgnc_symbol'),
    filters = 'go',
    values = "GO:0005576",
    mart = ensembl.human
  )

  extracellular <- intersect(incommon, paste(GO[, 1], "_tis1", sep = ""))

  all_sign <-
    (rank(-((degs_path_sign[[1]]$kExt1 / degs_path_sign[[1]]$kInt1)[extracellular]
    )) *
      rank(-((degs_path_sign[[2]]$kExt1 / degs_path_sign[[2]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[3]]$kExt1 / degs_path_sign[[3]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[4]]$kExt1 / degs_path_sign[[4]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[5]]$kExt1 / degs_path_sign[[5]]$kInt1)[extracellular]
      )) *
      rank(-((degs_path_sign[[6]]$kExt1 / degs_path_sign[[6]]$kInt1)[extracellular]
      )))

   names(all_sign) <- gsub("_tis1", "", rownames(all_sign))


  top_interactors<-names(sort(rank(all_sign))[1:5])

  write.csv(top_interactors, "top5hippo.csv")
  }
