####### Install R packages automatically #########
source("install_R_packages.R")


######## Load R packages ###############
library(devtools)
require(usethis)
library(shiny)
library(servr)
library(ggplot2)
library(pheatmap)
library(M3C)
library(RUVSeq)
library(scales)
library(dtwclust)
library(dplyr)
library(DESeq2)
library(ggcorrplot)
library(tibble)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(clusterProfiler)
library(cowplot)
library(scater)
library(hdf5r)
library(MAST)
library(Seurat)
#########################

options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=10000*1024^2)

ui <- fluidPage(
  tags$img(src='https://avatars1.githubusercontent.com/u/8896007?s=400&u=b0029c2e64f405ea0a46d311239b674a430ec77c&v=4'
           ,height='60',width='60', align='left'),
  tags$head(includeHTML(("GoogleAnalytics.html"))),
  # App title ----
  titlePanel("ezsinglecell : An integrated computational toolbox"),


  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the number of Dataset ----
      selectInput("Module",
                  label = "Module",
                  choices = c("Single Cell RNASeq Analysis", "Single cell data integration Analysis", "Single cell multiomics Analysis", "Flow cytometry Analysis", "Imaging mass cytometry Analysis", "Single cell ATAC-seq Analysis", "Spatial Transcriptomics Analysis", "Nanostring DSP Analysis"),
                  selected = "Single Cell RNASeq Analysis"),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis'",
        selectInput("scInput",
                    label = "Select Data Input Type",
                    choices = c("Raw Counts Matrix", "H5", "R Object"),
                    selected = "Raw Counts Matrix")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix'",
        fileInput("scCounts",
                    label = "Upload Counts File (Accepted Format: tab delimited text)",
                    accept = ".txt")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'H5'",
        fileInput("scH5",
                    label = "Upload H5 output from Cellranger or other toolkits (Accepted Format: H5)",
                    accept = ".h5")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'R Object'",
        fileInput("scRobj",
                    label = "Upload Seurat R object (Set Object name to seurat.object before upload)",
                    accept = ".Robj")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5' || input.scInput == 'R Object'",
        selectInput("Species_singlecell",
                    label = "Select Species",
                    choices = c("human", "mouse"),
                    selected = "human")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        numericInput("scMinCells",
                    label = "Input Minimum number of cells to express all genes",
                    value = 3,
                    min = 0,
                    max = 200000),
        verbatimTextOutput("3")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        numericInput("scMinFeatures",
                    label = "Input Minimum number of features all cells should express",
                    value = 100,
                    min = 0,
                    max = 30000),
        verbatimTextOutput("100")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        selectInput("scNormalization",
                    label = "Select Normalization Method",
                    choices = c("LogNormalize", "SCTransform"),
                    selected = "LogNormalize")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        selectInput("scReduction",
                    label = "Select Dimension Reduction",
                    choices = c("umap", "tsne"),
                    selected = "umap")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        selectInput("scDETest",
                    label = "Select Differential Expression Test (Please see: some of these tests increase run time significantly)",
                    choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"),
                    selected = "wilcox")),

      #conditionalPanel(
       # condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        #selectInput("scCellCycle",
        #            label = "Regress Cell Cycle Effect",
        #            choices = c("Yes", "No"),
        #            selected = "Yes")
        #),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        numericInput("VarFeatures",
                    label = "Input Number of Variable Features to use",
                    value = 2000,
                    min = 100,
                    max = 10000),
        verbatimTextOutput("2000")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        numericInput("scDims",
                    label = "Input Number of Dimensions to use",
                    value = 10,
                    min = 1,
                    max = 100),
        verbatimTextOutput("10")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        numericInput("scRes",
                    label = "Input resolution for clustering",
                    value = 0.5,
                    min = 0,
                    max = 10),
        verbatimTextOutput("0.5")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5' || input.scInput == 'R Object'",
        selectInput("scVisualization",
                    label = "Select single cell Visualization",
                    choices = c("Gene Expression Plot", "Dimension Reduction Plot", "Top10 Markers Heatmap", "Violin Plot", "DotPlot", "QC Metrics Plot"),
                    selected = "Dimension Reduction Plot")),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5' || input.scInput == 'R Object' && input.scVisualization == 'Gene Expression Plot' || input.scVisualization == 'Violin Plot' || input.scVisualization == 'DotPlot'",
        selectInput("scGene",
                    label = "Select Genes",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE)),

      conditionalPanel(
        condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
        downloadButton('scRNAObjectDownload', 'Download scRNA-Seq Seurat Object'))

    ),
    mainPanel(
          tabsetPanel(type = "tabs",
              tabPanel("Overview", fluidRow(
                p(strong("ezSinglecell"), "is a RShiny application developed with an intention to empower researchers from wet and dry lab to perform downstream Bioinformatics analysis. ezsinglecell powered by RShiny is packed with 8 modules : ", strong("Single cell RNA-seq, Data Integration, Multiomics, Flow cytometry, Imaging mass cytometry, single cell ATAC-seq, Spatial Transcriptomics and DSP Nanostring."), " These modules are designed in order to help researchers design a hypothesis or answer research questions with little or no expertise in Bioinformatics. In future, ezSinglecell will also be available as a command line application and New modules and functionalities will be added periodically.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("ezSingle cell RShiny is available on", a("GitHub", href = "https://github.com/raman91/ezsinglecell", target = "_blank"), ", if interested in hosting on your own servers.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")),
                p("Please post issues, suggestions and improvements using", a("Issues/suggestions", href = "https://github.com/raman91/ezsinglecell", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p("To view other tools and contributions please visit", a("GitHub", href = "https://github.com/JinmiaoChenLab", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"))),
              tabPanel("Result Window", textOutput('DisplayText'), DT::dataTableOutput("result"), downloadButton('downloadResult', 'Download Results')),
              tabPanel("Visualization Window", textOutput('DisplayText1'), downloadButton('downloadPlot', 'Save Plot'), plotOutput("Plot")),
              tabPanel("Getting Started", fluidRow(
                p(strong("ezsinglecell"), "is easy to use and is packed with powerful modules to help you analyze your data. Results generated from the modules are loaded on the Result window whereas the Visualization plots are displyed on Visualization windows.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Description of various modules"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Single Cell RNASeq Analysis module"), "is designed to perform single cell RNA-Seq analysis. ezsinglecell uses ", a("Seurat", href = "https://www.cell.com/action/showPdf?pii=S0092-8674%2821%2900583-3", target = "_blank"), "for analyzing single cell RNA-Seq data. In this module, users can input raw counts matrix, H5 output from cellranger or a processed Seurat object. ezsinglecell allows users to download R Objects after processing.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Single Cell RNASeq Integration Analysis module"), "is designed to perform single cell RNA-Seq integration analysis. ezsinglecell uses ", a("Seurat", href = "https://www.cell.com/action/showPdf?pii=S0092-8674%2821%2900583-3", target = "_blank"), "for analyzing single cell RNA-Seq integration data. In this module, users can input raw counts matrix, H5 output from cellranger or a processed Seurat object. ezsinglecell allows users to download R Objects after processing.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Single Cell Multiomics RNASeq Analysis module"), "is designed to perform single cell multiomics RNA-Seq analysis. ezsinglecell uses ", a("Seurat", href = "https://www.cell.com/action/showPdf?pii=S0092-8674%2821%2900583-3", target = "_blank"), "for analyzing single cell multiomics RNA-Seq data. In this module, users can input raw counts matrix, H5 output from cellranger or a processed Seurat object. ezsinglecell allows users to download R Objects after processing.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Flow cytometry module"), "is designed to perform flow cytometry data analysis. ezsinglecell integrates state-of-the-art bioinformatics methods and in-house novel algorithms for data pre-processing, data visualization through linear or non-linear dimensionality reduction (UMAP/tSNE), cell clustering, automatic identification of cell subsets and inference of the relatedness between cell subsets. Our web server can analyse datasets with millions of cells and performs clustering of 20 million cells in 1.7 hours using a fast clustering algorithm called ", a("FastPG", href = "https://www.biorxiv.org/content/10.1101/2020.06.19.159749v1.full", target = "_blank"), ".", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                p(strong("Imaging mass cytometry module"), "is designed to perform imaging mass cytometry data analysis. ezsinglecell is able to analyse high-dimensional imaging mass cytometry data to unravel the cellular composition, spatial architecture from the datasets, interactive visualization of cell subpopulations and progression profiles of key markers. The module also allows for cell-cell interactions from the imaging mass cytometry data."),
                p(strong("Single cell ATAC-seq module"), "is designed to perform single cell ATAC-seq data analysis."),
                p(strong("Spatial Transcriptomics module"), "is designed to perform spatial transcriptomics data analysis."),
                p(strong("DSP Nanostring analysis module"), "is designed to perform DSP Nanostring data analysis."),
                )),
              tabPanel("What's New in ezsinglecell", fluidRow(
                p(strong("Changelog"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:25px"),
                p(strong("Version 1.0 Log:"), "ezsinglecell has 8 modules (single cell RNA-seq, Data Integration, Multiomics, Flow cytometry, Imaging mass cytometry, single cell ATAC-seq, Spatial Transcriptomics and DSP Nanostring) and is built with an idea to empower researchers in performing bioinformatics analysis with little or no bioinformatics expertise", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")))
          )
      )
  )
)


server <- function(input, output, session) {

 scData <- reactive({
    if(input$Module == "Single Cell RNASeq Analysis"){
      if(input$scInput == "Raw Counts Matrix"){
        inFile <- input$scCounts
        if (is.null(inFile))
          return(NULL)
        data <- c()
        data <- readSparseCounts(inFile$datapath, sep = "\t", row.names = TRUE, col.names = TRUE)
        return(data)
      }
      else if(input$scInput == "H5"){
        inFile <- input$scH5
        if (is.null(inFile))
          return(NULL)
        data <- c()
        data <- Read10X_h5(inFile$datapath, use.names = TRUE, unique.features = TRUE)
        return(data)
      }
    }
 })

 observe({
     if(input$Module == "Single Cell RNASeq Analysis"){
      if(input$scInput == "Raw Counts Matrix"){
        inFile <- input$scCounts
        if (is.null(inFile))
          return(NULL)
        data <- c()
        data <- scData()
        genes <- rownames(data)
        genes <- data.frame(genes)
        genes <- genes$genes
        updateSelectInput(session, "scGene",
                         label = "Select Genes",
                         choices = genes)
      }
      else if(input$scInput == "H5"){
        inFile <- input$scH5
        if (is.null(inFile))
          return(NULL)
        data <- c()
        data <- scData()
        genes <- rownames(data)
        genes <- data.frame(genes)
        genes <- genes$genes
        updateSelectInput(session, "scGene",
                         label = "Select Genes",
                         choices = genes)
      }
      else if(input$scInput == "R Object"){
        inFile <- input$scRobj
        if (is.null(inFile))
          return(NULL)
        load(inFile$datapath)
        genes <- rownames(seurat.object)
        genes <- data.frame(genes)
        genes <- genes$genes
        updateSelectInput(session, "scGene",
                         label = "Select Genes",
                         choices = genes)
      }
    }
 })

SingleCell <- reactive({
  if(input$Module == "Single Cell RNASeq Analysis"){
    if((input$scInput == "Raw Counts Matrix")||(input$scInput == "H5")){
      data <- c()
      data <- scData()
      if(is.null(data))
        return(NULL)
      seurat.object <- CreateSeuratObject(counts = data, project = "scAnalysis", min.cells = input$scMinCells, min.features = input$scMinFeatures)
      if(input$Species_singlecell == "human"){
        seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
        #cc.genes <- readLines(con = "data/CellCycle_Human.txt")
      }
      else if(input$Species_singlecell == "mouse"){
        seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^mt-")
        #cc.genes <- readLines(con = "data/CellCycle_Mouse.txt")
      }
      nFeature_RNA_cutoff <- quantile(seurat.object@meta.data$nFeature_RNA, .95)
      nFeature_RNA_cutoff <- data.frame(nFeature_RNA_cutoff)
      percentMT_cutoff <- quantile(seurat.object@meta.data$percent.mt, .95)
      percentMT_cutoff <- data.frame(percentMT_cutoff)
      if(percentMT_cutoff$percentMT_cutoff == 0){
        seurat.object <- subset(seurat.object, subset = nFeature_RNA > 0 & nFeature_RNA < nFeature_RNA_cutoff$nFeature_RNA_cutoff & percent.mt < 1)
      }
      else if(percentMT_cutoff$percentMT_cutoff > 0){
        seurat.object <- subset(seurat.object, subset = nFeature_RNA > 0 & nFeature_RNA < nFeature_RNA_cutoff$nFeature_RNA_cutoff & percent.mt < percentMT_cutoff$percentMT_cutoff)
      }
      if(input$scNormalization == 'LogNormalize'){
        seurat.object <- NormalizeData(seurat.object)
        seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = input$VarFeatures)
        #if(input$scCellCycle == "Yes"){
          #s.genes <- cc.genes[1:43]
          #g2m.genes <- cc.genes[44:97]
          #seurat.object <- CellCycleScoring(object = seurat.object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
          #seurat.object@meta.data$CC.Difference <- seurat.object@meta.data$S.Score - seurat.object@meta.data$G2M.Score
          #seurat.object <- ScaleData(seurat.object, vars.to.regress = c("percent.mt", "CC.Difference"))
        #}
        #else if(input$scCellCycle == "No"){
      #}
        all.genes <- rownames(seurat.object)
        seurat.object <- ScaleData(seurat.object, features = all.genes)

        seurat.object <- RunPCA(seurat.object)
        seurat.object <- FindNeighbors(seurat.object, dims = 1:input$scDims)
        seurat.object <- FindClusters(seurat.object, resolution = input$scRes)
        if(input$scReduction == "umap"){
          seurat.object <- RunUMAP(seurat.object, dims = 1:input$scDims)
        }
        else if(input$scReduction == "tsne"){
          seurat.object <- RunTSNE(seurat.object, dims = 1:input$scDims)
        }
        return(seurat.object)
      }
      else if(input$scNormalization == 'SCTransform'){
        #if(input$scCellCycle == "Yes"){
          #s.genes <- cc.genes[1:43]
          #g2m.genes <- cc.genes[44:97]
          seurat.object <- SCTransform(seurat.object, variable.features.n = input$VarFeatures)
          #seurat.object <- CellCycleScoring(object = seurat.object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
          #seurat.object@meta.data$CC.Difference <- seurat.object@meta.data$S.Score - seurat.object@meta.data$G2M.Score
          seurat.object <- SCTransform(seurat.object, vars.to.regress = "percent.mt", variable.features.n = input$VarFeatures)
        #}
        #else if(input$scCellCycle == "No"){
         # seurat.object <- SCTransform(seurat.object, vars.to.regress = "percent.mt", variable.features.n = input$VarFeatures)
        #}
        seurat.object <- RunPCA(seurat.object)
        seurat.object <- FindNeighbors(seurat.object, dims = 1:input$scDims)
        seurat.object <- FindClusters(seurat.object, resolution = input$scRes)
        if(input$scReduction == "umap"){
          seurat.object <- RunUMAP(seurat.object, dims = 1:input$scDims)
        }
        else if(input$scReduction == "tsne"){
          seurat.object <- RunTSNE(seurat.object, dims = 1:input$scDims)
        }
        return(seurat.object)
      }
    }
    else if(input$scInput == "R Object"){
      inFile <- input$scRobj
      if(is.null(inFile))
        return(NULL)
      load(inFile$datapath)
      Idents(seurat.object) <- seurat.object@meta.data$seurat_clusters
      return(seurat.object)
    }
  }
})

 SingleCellMarkers <- reactive({
  seurat.object <- SingleCell()
  if(is.null(seurat.object))
    return(NULL)
  markers_clusters <- FindAllMarkers(seurat.object, test.use = input$scDETest)
  return(markers_clusters)
 })


 PCAplot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
      return(NULL)
     data <- c()
     if(input$Module == "Normalization"){
      data <- as.matrix(NormalizationData())
     }
     else if(input$Module == "Visualization"){
      data <- as.matrix(VisualizationData())
     }
     data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
     data.t <- t(data)
     pca <- prcomp(data.t, center=T, scale. = T)
     pc1 <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
     pc2 <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
     PC1_use <- paste0("PC1", "(", pc1, "%)")
     PC2_use <- paste0("PC2", "(", pc2, "%)")
     Samples_temp <- rownames(data.t)
     Samples <- factor(Samples_temp)
     scores <- data.frame(Samples_temp, pca$x[,1:3])
     MIN_X <- min(scores$PC1)
     Max_X <- max(scores$PC1)
     header <- "Principal Component Analysis"
     qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=1) + geom_text(aes(label=Samples_temp), hjust=0, vjust=0) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
 })

 TSNEplot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
      return(NULL)
     inGroupFile <- input$GroupFile
     if (is.null(inGroupFile))
      return(NULL)
     data <- c()
     if(input$Module == "Normalization"){
      data <- as.matrix(NormalizationData())
     }
     else if(input$Module == "Visualization"){
      data <- as.matrix(VisualizationData())
     }
     Group <- as.matrix(read.table(inGroupFile$datapath, header=T, sep="\t", row.names=1, check.names=F))
     Group <- data.frame(Group)
     GroupUse <- as.factor(Group$group)
     tsne(data, labels=GroupUse, perplex=input$perplexity) + xlab("tSNE_1") + ylab("tSNE_2") + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
 })

 HEATMAPplot <- reactive({
     inFile <- input$File
     if (is.null(inFile))
      return(NULL)
     data <- c()
     if(input$Module == "Normalization"){
      data <- as.matrix(NormalizationData())
     }
     else if(input$Module == "Visualization"){
      data <- as.matrix(VisualizationData())
     }
     data <- data[apply(data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
     if(input$Cluster == "Rows"){
      pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
     }
     else if(input$Cluster == "Columns"){
        pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
     }
     else if(input$Cluster == "Rows and Columns"){
        pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
     }
     else if(input$Cluster == "None"){
        pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
     }
 })

 CorrelationPlot <- reactive({
     if(input$Correlation == "Samples"){
      inFile <- input$File
        if (is.null(inFile))
          return(NULL)
        titleuse <- paste0("Displaying Correlation Plot of ", input$Correlation)
        data <- CorrelationData()
        ggcorrplot(data, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
     }
     else if(input$Correlation == "Genes"){
        inFile <- input$File
        if (is.null(inFile))
          return(NULL)
        Genes_List <- input$GeneList
        if (is.null(Genes_List))
          return(NULL)
        titleuse <- paste0("Displaying Correlation Plot of ", input$Correlation)
        data <-   CorrelationData()
        ggcorrplot(data, hc.order = TRUE, outline.color = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726"), legend.title = "Correlation") + labs(title = titleuse)
     }
 })

 DifferentialExpressionPlot <- reactive({
     inFile <- input$Counts
     if (is.null(inFile))
      return(NULL)
     DEresult <- DEData()
     DEgenes <- rownames(DEresult)
     data <- as.matrix(read.table(inFile$datapath, sep="\t", header=T, row.names=1, check.names=F))
     datause <- subset(data, rownames(data) %in% DEgenes)
     datause <- log2(datause+1)
     if(input$PlotType == "pca"){
      data <- as.matrix(datause)
        data <- data[apply(data[,-1], 1, function(x) !all(x==0)),]
        data.t <- t(data)
        pca <- prcomp(data.t, center=T, scale. = T)
        pc1 <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
      pc2 <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
      PC1_use <- paste0("PC1", "(", pc1, "%)")
        PC2_use <- paste0("PC2", "(", pc2, "%)")
        Samples_temp <- rownames(data.t)
        Samples <- factor(Samples_temp)
        scores <- data.frame(Samples_temp, pca$x[,1:3])
        MIN_X <- min(scores$PC1)
        Max_X <- max(scores$PC1)
        header <- "Principal Component Analysis"
        qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=1) + geom_text(aes(label=Samples_temp), hjust=0, vjust=0) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header) + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
    }
    else if(input$PlotType == "tsne"){
      inGroupFile <- input$GroupFile
      if (is.null(inGroupFile))
        return(NULL)
      data <- as.matrix(datause)
      Group <- as.matrix(read.table(inGroupFile$datapath, header=T, sep="\t", row.names=1, check.names=F))
      Group <- data.frame(Group)
      GroupUse <- as.factor(Group$group)
      tsne(data, labels=GroupUse, perplex=input$perplexity) + xlab("tSNE_1") + ylab("tSNE_2") + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5)))
    }
    else if(input$PlotType == "heatmap"){
      data <- as.matrix(datause)
      data <- data[apply(data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
      if(input$Cluster == "Rows"){
        pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
      }
      else if(input$Cluster == "Columns"){
       pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
      }
      else if(input$Cluster == "Rows and Columns"){
       pheatmap(data, scale=input$Scaling, cluster_rows=TRUE, cluster_cols=TRUE, main="Heatmap", border_color = "NA")
      }
      else if(input$Cluster == "None"){
       pheatmap(data, scale=input$Scaling, cluster_rows=FALSE, cluster_cols=FALSE, main="Heatmap", border_color = "NA")
      }
    }
    else if(input$PlotType == "Functional and Pathway Enrichment"){
      if(input$SpeciesUse == "human"){
       data <- as.matrix(datause)
       genes_de <- rownames(data)
       genes_de_id <- mapIds(org.Hs.eg.db, genes_de, 'ENTREZID', 'SYMBOL')
       enrichemnt <- enrichPathway(gene = genes_de_id, pvalueCutoff = 0.05, readable=T, organism = "human", maxGSSize = 5000)
       emapplot(enrichemnt)
      }
      else if(input$SpeciesUse == "mouse"){
       data <- as.matrix(datause)
       genes_de <- rownames(data)
       genes_de_id <- mapIds(org.Mm.eg.db, genes_de, 'ENTREZID', 'SYMBOL')
       enrichemnt <- enrichPathway(gene = genes_de_id, pvalueCutoff = 0.05, readable=T, organism = "mouse", maxGSSize = 5000)
       emapplot(enrichemnt)
      }
    }
    #else if(input$PlotType == "Volcano Plots"){
    #   data <- DEresult
    #   EnhancedVolcano(data, lab = rownames(data), x = 'logFC', y = 'FDR')
    #}
 })

 SingleCellPlot <- reactive({
  seurat.object <- SingleCell()
  if(input$scVisualization == "Gene Expression Plot"){
    FeaturePlot(seurat.object, features = input$scGene, pt.size = 2)
  }
  else if(input$scVisualization == "Dimension Reduction Plot"){
    DimPlot(seurat.object, label = TRUE, pt.size = 2)
  }
  else if(input$scVisualization == "Top10 Markers Heatmap"){
    top10 <- SingleCellMarkers() %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)   ###### based on v4 Seurat
    DoHeatmap(seurat.object, features = top10$gene)
  }
  else if(input$scVisualization == "Violin Plot"){
    VlnPlot(seurat.object, features = input$scGene) + RotatedAxis()
  }
  else if(input$scVisualization == "DotPlot"){
    DotPlot(seurat.object, features = input$scGene) + RotatedAxis()
  }

  else if(input$scVisualization == "QC Metrics Plot"){
    Idents(seurat.object) <- seurat.object@meta.data$orig.ident
    VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  }
 })

 output$result <- DT::renderDataTable({
   if(input$Module == "Single Cell RNASeq Analysis"){
      DT::datatable(SingleCellMarkers())
    }
 })

 output$Plot <- renderPlot({
   if(input$Module == "Single Cell RNASeq Analysis"){
      SingleCellPlot()
    }
  }, width=1200, height=1000)

 output$downloadResult <- downloadHandler(
      filename = function() {
        paste0(input$Module, "_Analysis_Result", "-", Sys.Date(), ".txt")
      },
      content = function(file) {
       if(input$Module == "Single Cell RNASeq Analysis"){
          write.table(SingleCellMarkers(), file, sep="\t", quote=F)
        }
      }
 )

  output$DisplayText <- renderText({"Please Upload Input Data"})
  output$DisplayText1 <- renderText({"Please Upload Input Data"})

  observe({
    if(input$Module == "Single Cell RNASeq Analysis"){
      output$DisplayText <- renderText({"Please Upload Input Data"})
      inFile <- input$scCounts
      inFile1 <- input$scH5
      inFile2 <- input$scRobj
      if((!is.null(scData()))&&(is.null(SingleCell()))&&(is.null(SingleCellMarkers()))){
        output$DisplayText <- renderText({"ezsinglecell is processing your data"})
        output$DisplayText1 <- renderText({"ezsinglecell is generating selected visualization plot"})
      }
      else if((!is.null(scData()))&&(!is.null(SingleCell()))&&(!is.null(SingleCellMarkers()))){
        output$DisplayText <- renderText({"Processing Complete"})
        output$DisplayText1 <- renderText({"Plot Generated"})
      }
    }
  })

  output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0(input$Module, "_Plot", "-", Sys.Date(), ".png")
      },
      content = function(file) {
        if(input$Module == "Single Cell RNASeq Analysis"){
          SingleCellPlot()
          ggsave(file, width = 15, height = 15)
       }
    },
    contentType = 'image/png'
  )

  output$scRNAObjectDownload <- downloadHandler(
    filename = function() {
        paste0(input$Module, "_Analysis", "-", Sys.Date(), ".Robj")
      },
      content = function(file) {
        seurat.object <- SingleCell()
        save(seurat.object, file = file)
      }
  )
}
shinyApp(ui = ui, server = server)
