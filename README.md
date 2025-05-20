# crossWGCNA_manuscript
Code to reproduce the results of the manuscript "Cross-tissue gene expression interactions from bulk, single cell and spatial transcriptomics with crossWGCNA".

## Description
To run the code, input data should be downloaded from Zenodo (10.5281/zenodo.8355624), unzipped and placed in a working directory containing the "scripts" folder from the current repository.
Note that the file "GSE83591_Correspondence_LCM.xlsx" is available upon request. CrossWGCNA networks are nevertheless made available in the Zenodo file, "results" folder.
The code is organized as follows:
- 00_main: wrapper to run all the analyses except for the spatial transcriptomics pipeline. Input parameters can be modified, but to reproduce the manuscript's results keep the parameters as set.
- 01_degs: runs crossWGCNA to compute inter-tissue degrees for six laser capture microdissection (LCM) datasets.
- 02_GSEA: Gene Set Enrichment Analysis (GSEA) of kRatio-ranked genes for the six LCM datasets.
- GSE161529_degs: computes crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset.
- GSE161529_degs_rand: computes crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset randomizing the matching between stroma and epithelium samples.
- 03_cocultures: tests whether genes with high kRatio in LCM are differentially expressed upon cocultures more often than by chance.
- 04_LCM_vs_singlecell: comparse LCM and single-cell kRatios.
- 04_LCM_vs_singlecell_rand: comparse LCM and single-cell kRatios with randomized single cell matching.
- 05_compare_versions: compares the two proposed algorithms.
- 06_datasets_coherence: compares the coherence of kratios between different LCM datasets.
- 07_mrnavsprot: runs crossWGCNA to infer inter-regulatory-layer interactions.
- ST_pipeline: runs crossWGCNA on a spatial transcriptomics dataset.
- sc_clusts_netdiff: calls GSE161529_degs_clust.R and makes plots showing that clustering parameters don't impact crossWGCNA degrees.
- GSE161529_degs_clust: computes crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset with different combinations of parameters for cell clustering.


