# crossWGCNA_manuscript
Code to reproduce the results of the manuscript "Cross-tissue gene expression interactions from bulk, single cell and spatial transcriptomics with crossWGCNA".

## Description
To run the code, input data should be downloaded from Zenodo (10.5281/zenodo.8355624), unzipped and placed in a working directory containing the "scripts" folder from the current repository.
The code is organized as follows:
- 00_main: wrapper to run all the analyses except for the spatial transcriptomics pipeline. Input parameters can be modified, but to reproduce the manuscript's results keep the parameters as set.
- 01_degs: runs crossWGCNA to compute inter-tissue degrees for six laser capture microdissection (LCM) datasets.
- 02_GSEA: Gene Set Enrichment Analysis (GSEA) of kRatio-ranked genes for the six LCM datasets.
- GSE161529_degs: computes crossWGCNA degrees on the GSE161529 single cell RNA-seq dataset.
- 03_cocultures: tests whether genes with high kRatio in LCM are differentially expressed upon cocultures more often than by chance.
- 04_LCM_vs_singlecell: comparse LCM and single-cell kRatios.
- 05_compare_versions: compares the two proposed algorithms.
- 06_datasets_coherence: compares the coherence of kratios between different LCM datasets.
- 07_mrnavsprot: runs crossWGCNA to infer inter-regulatory-layer interactions.
- ST_pipeline: runs crossWGCNA on a spatial transcriptomics dataset.


