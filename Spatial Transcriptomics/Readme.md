**Module 3: Spatial Transcriptomics module**

The Spatial Transcriptomics module offers analysis for 10X Visium and Xenium data. 

For 10X Visium data, ezSingleCell offers data analysis such as spatial clustering and deconvolution using 2 methods namely:

1. Seurat (https://satijalab.org/seurat/articles/spatial_vignette.html) 

2. In-house developed algorithm named GraphST avaialable at https://github.com/JinmiaoChenLab/GraphST [Long Y et al, Nature Communications (2023)]

The acceptable input format is in the form of 10X SpaceRanger output. The test dataset used is 'Human Breast Cancer' data comprising of 3798 spots downloaded from https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0.  

ezSingleCell also offers analysis support for sub-cellular resolution data such as data acquired using the **Xenium** platform. Users can perform clustering analysis and interactively visualize the expression patterns at sub-cellular level. In addition, users can zoom into a certain section to have a closer look at the cellular composition and inter-cellular interactions. Users can also view the expression profile of each gene on the tissue slice.

The test dataset for Xenium can be downloaded using the link :

https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip

![image](https://github.com/JinmiaoChenLab/ezSingleCell2/assets/8286779/3878d12b-c808-4b0f-8c24-e46d36b39cba)

