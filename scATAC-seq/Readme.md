**Module5: scATACseq module**

The single cell ATAC-seq module accepts input in 10X CellRanger-atac format that comprises of :

1. Peak/Cell matrix (filtered_peak_bc_matrix.h5): This is analogous to the gene expression count matrix used to analyze single-cell RNA-seq. However, instead of genes, each row of the matrix represents a region of the genome (a peak), that is predicted to represent a region of open chromatin. Each value in the matrix represents the number of Tn5 integration sites for each single barcode (i.e. a cell) that map within each peak. 

2. Metadata

3. Fragments file: This represents a full list of all unique fragments across all single cells. This file contains all fragments associated with each single cell, as opposed to only fragments that map to peaks.

ezSingleCell employs the Signac package (https://stuartlab.org/signac/) and includes variety of analysis tasks such as Cell Clustering, Target Enrichment analysis (Link Peak to Genes), Gene Set Enrichment Analysis of scATAC-seq data using tools such as rGREAT and fgsea, Differential peak analysis. In addition, scATAC-seq supports interaction with other modules in ezSingleCell such as scRNA-seq or Data Integration module to annotate the cell types in scATAC-seq data. 

![image](https://github.com/JinmiaoChenLab/ezSingleCell2/assets/8286779/d292f2e4-e004-46c0-aa52-eeba52025f56)
