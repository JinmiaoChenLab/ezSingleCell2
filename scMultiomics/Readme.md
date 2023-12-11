**Module 4: scRNA-seq module**

The single cell Multiomics module can perform CITE-seq (RNA + Protein) and 10X Multiome (RNA + ATAC) analysis. It accepts input in 10X cellranger format. 

For CITE-seq data analysis, ezSingleCell supports data generated from TotalSeq A,B and C while for 10X Multiome analysis, the accepted formats is the output from cellranger-arc. 

ezSingleCell accepts input in the form of counts table of both RNA and Protein modalities (in .h5 format) in case of CITE-Seq data and counts table of both RNA and ATAC modalities as well as fragments file in case of 10X Multiome data.

The data for scMultiomics can be downloaded from:

a) CITE-seq : https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3 

b) 10X Multiome: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_unsorted_3k? 

Currently, we employ 2 methods to perform single cell Multiomics analysis:

a) Seurat (https://satijalab.org/seurat/articles/multimodal_vignette), and 

b) MOFA2 (https://biofam.github.io/MOFA2/)

**For CITE-Seq:**

![image](https://github.com/JinmiaoChenLab/ezSingleCell2/assets/8286779/1ce3fac0-b0b1-47b3-970d-f962d43c00e7)
