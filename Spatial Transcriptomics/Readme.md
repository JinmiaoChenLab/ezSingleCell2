**Module 3: Spatial Transcriptomics module**

The Spatial Transcriptomics module offers analysis for 10X Visium and Xenium data.

For 10X Visium, ezSingleCell accepts data in the form of 10X SpaceRanger output. The test dataset used is Mouse Brain Sagittal Anterior data comprising of 2695 spots downloaded from https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0. 

ezSingleCell also offers analysis support for sub-cellular resolution data such as data acquired using the **Xenium** platform. Users can perform clustering analysis and interactively visualize the expression patterns at sub-cellular level. In addition, users can zoom into a certain section to have a closer look at the cellular composition and inter-cellular interactions. Users can also view the expression profile of each gene on the tissue slice.

The test dataset for Xenium can be downloaded using the following commands :

wget https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip

unzip Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip
