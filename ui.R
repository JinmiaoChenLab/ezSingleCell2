library(shiny)
library(plotly)
library(shinythemes)
library(shinydashboard)
library(shinyalert)
library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library(Seurat)
library(pbmcapply)
library(devtools)
library(harmony)
library(reticulate)
library(sceasy)
library(pheatmap)
library(heatmaply)
library(liana)
library(reticulate)
library(anndata)
library(EnhancedVolcano)
library(ggplot2)
library(lisi)
library(cowplot)
library(FastIntegration)
library(Signac)
library(cluster)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(purrr)
library(dplyr)
library(DT)
library(msigdbr)
library(fgsea)
library(MOFA2)
library(ktplots)
library(ggplot2)
library(SeuratWrappers)
#library(car)
#library(nortest)
library(shinyFiles)
library(bslib)

shinyUI(fluidPage(theme = shinytheme("cerulean"),
                  tagList(tags$head(tags$style(type = 'text/css','.navbar-brand{display:none;}')),             
    navbarPage("",
               tabPanel(icon("home"),
                        titlePanel(h1(p(strong("ezSinglecell : An integrated one-stop single-cell and spatial analysis toolbox for bench scientists")), align = "center")), 
                        fluidRow(
                                 column(
                                   br(),
                                   
                                   p(strong("ezSingleCell"), "is an integrated one-stop single-cell and spatial analysis toolbox developed by Chen Jinmiao's lab with an intention to empower bench scientists to perform downstream Bioinformatics analysis. In the current version, we incorporate 5 modules : Single cell RNA-seq, Single cell Data Integration, Single cell Multiomics, Single Cell ATAC-seq and Spatial Transcriptomics.", style="text-align:justify;color:black;font-size:15px"),
                                   column(6,
                                   imageOutput("demo_image"),
                                   ),
                                   column(6,
                                   imageOutput("demo_image1"),
                                   ),
                                   br(),
                                   br(),
                                   p("In this web server, we combine in-house novel algorithms such as CELLiD (for cell type identification), along with existing top performing methods for both basic and advanced downstream analyses such as batch effect removal, trajectory, cell-cell communication, differential abundance, and spatial deconvolution.", style="text-align:justify;color:black;font-size:15px"),
                                   p("Currently ezSingleCell supports inputs in different formats such as text, csv or 10X cell ranger output.  ezSingleCell is also available as a software package with a" , a("ShinyApp interface", href = "https://github.com/JinmiaoChenLab/ezSingleCell2", target = "_blank"), ",that can be run on a computer with basic memory requirements", style="text-align:justify;color:black;font-size:15px"),
                                   p("The toolkit uses example data from", a("DISCO", href = "https://www.immunesinglecell.org/", target="_blank"),"that contains data from 6141 samples, covering 107 tissues/cell lines/organoids, 158 diseases, and 20 platforms.", style="text-align:justify;color:black;font-size:15px"), width=12)),
                    
                            
                   # p(em("Developed by"),br("CJM Lab"),style="text-align:center; font-family: times")
                        ),
               tabPanel("Single cell RNA-Sequencing",
                        navlistPanel(widths=c(2,10),
                                     tabPanel("Overview",
                                              h2(p("Workflow for scRNA-Seq module")),
                                              br(),
                                              column(12,
                                              imageOutput("scrna_image1"),
                                              ),
                                              #column(6,
                                              #imageOutput("scrna_image2")       
                                              #),
                                              ),
                                     tabPanel("Upload your data",
                                              column(9,
                                                     column(5,
                                                            #h4('Load Data:'),
                                                            wellPanel(
                                                              titlePanel(h4(p("Load your input data"))),
                                                              br(),
                                                              selectInput("scInput",
                                                                          label = "Select Data Input Type",
                                                                          choices = c("Raw Counts Matrix", "10X cellranger", "rds object"),
                                                                          selected = "Raw Counts Matrix"),
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'Raw Counts Matrix'",
                                                                fileInput("tpmFiles",
                                                                          label = "Counts File (Accepted Format: text)",
                                                                          accept = ".txt"),
                                                                       actionButton("loadexample_tpm", "Load example and run", icon = icon("hand-o-right"))
                                                               ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == '10X cellranger'",
                                                                fileInput("scH5",
                                                                          label = "Cellranger output (Accepted Format: .h5)",
                                                                          accept = ".h5"),
                                                                actionButton("loadexample_scH5", "Load example and run", icon = icon("hand-o-right"))
                                                                ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'rds object'",
                                                                fileInput("rds",
                                                                          label = "Seurat Object (Accepted Format: .rds)",
                                                                          accept = ".rds"),
                                                                actionButton("loadexample_rds", "Load example and run", icon = icon("hand-o-right"))
                                                                ),
                                                              br(),
                                                                textInput(inputId = "projName",
                                                                          label = "Project Name",
                                                                          value = "scRNA"),
                                                              
                                                              fluidRow(
                                                                actionButton("loadButton", "Load data", icon = icon("hand-o-right")),
                                                                actionButton("reset_scRNA", "Reset", icon = icon("repeat"))
                                                              ),
                                                            )),
                                                     
                                                     column(3,
                                                            numericInput(inputId = "min.genes",
                                                                         label = "Min. genes",
                                                                         value = 200,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            numericInput(inputId = "min.cells",
                                                                         label = "Min. cells",
                                                                         value = 3,
                                                                         min = 1)
                                                     ),
                                                     
                                                     column(3,
                                                     actionButton("create_seurat", "Process", icon = icon("hand-o-right"))
                                                     ),
                                                     
                                                     
                                     ),
                                     
                                       column(12,
                                              withSpinner(dataTableOutput('countdataDT')),
                                       
                                     downloadButton('downloadCount', 'Download Table'))),
                                     
                                     tabPanel("Quality control Plot",
                                            tabsetPanel(id="qc_scRNA",
                                              tabPanel("Violin Plot", 
                                                column(3,
                                                       plotOutput("nFeature_RNAPlot")
                                                ),
                                                column(3,
                                                       plotOutput("mitoPlot")
                                                ),
                                                column(3,
                                                       plotOutput("nCount_RNAPlot")
                                                ),
                                                column(12,
                                                       column(3,
                                                       downloadButton('download_nFeature_RNA', 'Download nFeature (as png)'),
                                                       ),
                                                       column(3,
                                                       downloadButton('download_mito', 'Download mito (as png)'),
                                                       ),
                                                       column(3,
                                                       downloadButton('download_nCount_RNA', 'Download nCount (as png)'),   
                                                       ),
                                                      )
                                                  ),
                                                
                                     tabPanel("Feature Scatter Plot",
                                              column(5,
                                                     plotlyOutput("FeatureScatterPlot1")
                                              ),
                                              column(5,
                                                     plotlyOutput("FeatureScatterPlot2")
                                              ),
                                              br(),
                                              br(),
                                              column(12,
                                                     column(5,
                                                            downloadButton('download_FeatureScatterPlot1', 'Download (as png)'),
                                                     ),
                                                     column(5,
                                                            downloadButton('download_FeatureScatterPlot2', 'Download (as png)'),
                                                     ),
                                                  ),
                                                 )
                                                ),
                                              ),
                                     tabPanel("Normalization and Variable Feature Selection",
                                              
                                              selectInput("norm1",
                                                          label = "Normalization method",
                                                          choices = c("LogNormalize", "SCTransform")
                                              ),
                                              
                                              textOutput("nVarGenes"),
                                              fluidRow(
                                                column(3,
                                                       numericInput("var.genes",
                                                                    label = "Number of variable genes",
                                                                    value = 2000,
                                                                    min = 500,
                                                                    step = 500)
                                                ),
                                                column(3,
                                                       selectInput("selection.method",
                                                                   label = "Selection method",
                                                                   choices = c("vst", "dispersion"))
                                                ),
                                                column(4,
                                                       br(),
                                                       actionButton("findVarGenes", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                       #actionButton("doSCTransform", "Run SCTransform", icon = icon("hand-pointer-o"))
                                                       # actionButton("doVarplot", "Plot variable genes", icon = icon("hand-pointer-o"))
                                                )),
                                              plotOutput("VarGenes", width = "100%"),
                                            ),
                                     tabPanel("PCA",
                                              tabsetPanel(id="Pca",
                                                          tabPanel(title="PCA Plot", value="P_panel1",
                                                                   br(),
                                                                   column(3,
                                                                   selectInput("assays1",
                                                                               label = "Normalize by:",
                                                                               choices = c("LogNormalization", "SCTransform")
                                                                    ),
                                                                   ),
                                                                   br(),
                                                                   column(3,
                                                                          actionButton("doPCA", "Run PCA", icon = icon("hand-pointer-o"))
                                                                   ),
                                                                   br(),
                                                                   br(),
                                                                   column(6,
                                                                   plotlyOutput("PCA2DPlot", width = "100%")),
                                                                   column(12,
                                                                   column(3,
                                                                          downloadButton('download_PCA', 'Download PCA Plot (as png)'),
                                                                   ),
                                                                   column(3,
                                                                          downloadButton('download_PCA_embedding', 'Download PCA Embedding (as csv)'),
                                                                   ),
                                                                  ),
                                                          ),
                                                          tabPanel(title="PC Gene Visualisation", value="P_panel2",
                                                                   br(),
                                                                   selectInput("select.pc",
                                                                               label = "PC to plot",
                                                                               choices = c(1:50)
                                                                   ),
                                                                   fluidRow(
                                                                     column(4,
                                                                            plotOutput("vizPlot", width = "100%", height = "600px")
                                                                     ),
                                                                     column(8,
                                                                            plotOutput("PCHeatmap", width = "100%", height = "600px")
                                                                     ),
                                                                     column(12,
                                                                     column(3,
                                                                            downloadButton('download_vizPlot', 'Download vizPlot (as png)'),
                                                                     ),
                                                                     column(3,
                                                                            downloadButton('download_PCHeatmap', 'Download PCHeatmap (as png)'),
                                                                     ),
                                                                    ),
                                                                   ),
                                                                   br(),
                                                                   DT::dataTableOutput("PCtable"),
                                                                   column(3,
                                                                          downloadButton('download_PCTable', 'Download top genes (as csv)'),
                                                                   ),
                                                          ),
                                                          
                                                          tabPanel(title="Elbow", value="P_panel4",
                                                                   br(),
                                                                   #actionButton("doElbow", label = "Get Elbow Plot"),
                                                                   #br(),
                                                                   br(),
                                                                   plotOutput("Elbow", width = "100%"),
                                                                   column(12,
                                                                          column(3,
                                                                                 downloadButton('download_Elbow', 'Download ElbowPlot (as png)'),
                                                                          ),
                                                                   ),
                                                                   
                                                          )
                                              )),
                                     
                                     tabPanel("Clustering",
                                              fluidRow(
                                                column(3,
                                                       numericInput("clus.res",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       selectInput("dim.used",
                                                                   label = "PC to use",
                                                                   choices = c(10:50)),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                                                       textOutput("cluster.done"),
                                                       br()
                                                )),
                                              br(),
                                              plotlyOutput("Cluster2DPlot_1", width = "50%"),
                                              br(),
                                              column(12,
                                                     column(3,
                                                            downloadButton('download_Cluster', 'Download ClusterPlot (as png)'),
                                                     ),
                                                     column(3,
                                                            downloadButton('download_ClusterTable', 'Download Cluster table (as csv)'),
                                                     ),
                                              ),
                                     ),    
                                     
                                     tabPanel("UMAP",
                                              
                                              fluidRow(
                                                column(3,
                                                       numericInput("dim.used",
                                                                    label = "Dimensions used",
                                                                    value = 10)
                                                ),
                                                br(),
                                                column(3,
                                                       actionButton("doUmap", "Run UMAP", icon = icon("hand-pointer-o")),
                                                       textOutput("Umap.done"),
                                                       br()
                                                )),
                                              br(),
                                              plotlyOutput("Umap_2d_plot_1", width = "50%"),
                                              br(),
                                              column(12,
                                                     column(3,
                                                            downloadButton('download_UMAP', 'Download UMAP Plot (as png)'),
                                                     ),
                                                     column(3,
                                                            downloadButton('download_UMAP_embedding', 'Download UMAP Embeddings (as csv)'),
                                                     ),
                                              ),
                                     ),
                                        tabPanel("tSNE",
                                                 
                                                 fluidRow(
                                                   column(3,
                                                          numericInput("dim.used",
                                                                       label = "Dimensions used",
                                                                       value = 10)
                                                   ),
                                                   column(3,
                                                          uiOutput("perplex.option")
                                                   ),
                                                   column(3,
                                                          br(),
                                                          actionButton("doTsne", "Run TSNE", icon = icon("hand-pointer-o")),
                                                          textOutput("Tsne.done"),
                                                          br()
                                                   )),
                                                 br(),
                                                 plotlyOutput("Tsne_2d_plot_1", width = "50%"),
                                                 br(),
                                                 column(12,
                                                        column(3,
                                                               downloadButton('download_Tsne', 'Download tSNE Plot (as png)'),
                                                        ),
                                                        column(3,
                                                               downloadButton('download_Tsne_embedding', 'Download tSNE Embeddings (as csv)'),
                                                        ),
                                                 ),
                                        ),
                                     
                                     tabPanel("Cell type identification",
                                              column(3,
                                                     selectInput("cellatlas",
                                                                 label = "Reference Atlas",
                                                                 choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                 selected = "all")
                                              ),
                                              column(3,
                                                     br(),
                                                     actionButton("doCELLiD", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                     textOutput("CELLiD.done"),
                                                     br()
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              plotlyOutput("Umap_cellid", width = "50%"),
                                              plotlyOutput("Umap_cellid1", width = "50%"),
                                              br(),
                                              column(12,
                                                     column(4,
                                                            downloadButton('download_Umap_cellid', 'Download CELLiD predict1 (as png)'),
                                                     ),
                                                     column(4,
                                                            downloadButton('download_Umap_cellid1', 'Download CELLiD predict2 (as png)'),
                                                     ),
                                                     column(4,
                                                            downloadButton('download_cellid_prediction', 'Download CELLiD predictions (in csv)'),
                                                     ),
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              DT::dataTableOutput("ct.table")
                                     ),
                                     
                                     tabPanel("Cell-cell similarity",
                                              column(3,
                                                     br(),
                                                     selectInput("cell1",
                                                                 label = "Group by",
                                                                 choices = c("seurat_clusters", "primary.predict"),
                                                                 selected = "primary.predict"),
                                              ),
                                              column(3,
                                                     br(),
                                                     br(),
                                                     actionButton("cell_cell", "Run cell-cell similarity", icon = icon("hand-pointer-o")),
                                                     textOutput("CELL.done"),
                                                     br()
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              column(9,
                                              plotlyOutput("cell_cell_sim")
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              column(12,
                                                     column(4,
                                                            downloadButton('download_cell_cell_sim', 'Download Celltype similarity plot (as png)'),
                                                     ),
                                                     column(4,
                                                            downloadButton('download_cor.table', 'Download Celltype similarity table (in csv)'),
                                                     ),
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              DT::dataTableOutput("cor.table")
                                     ),
                                      tabPanel("DEGs",
                                               fluidRow(
                                                 column(3,
                                                 selectInput("deg1",
                                                             label = "Group by",
                                                             choices = c("seurat_clusters", "primary.predict"),
                                                             selected = "primary.predict"),
                                                 ),
                                                 column(3, numericInput("min_pct",
                                                                        label = "min.pct",
                                                                        value = 0.25,
                                                                        min = 0,
                                                                        step = 0.01)
                                                 ),
                                                 column(3, numericInput("logfc",
                                                                        label = "logfc.threshold",
                                                                        value = 0.25,
                                                                        min = 0,
                                                                        step = 0.01)
                                                 ),
                                                 column(3, selectInput("test.use",
                                                                       label = "Test use",
                                                                       choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                 ),
                                                 br(),
                                                 column(3,
                                                        actionButton("doDeg", "Run DEGs", icon = icon("hand-pointer-o"))
                                                 )),
                                               br(),
                                               column(12,
                                                      DT::dataTableOutput("Deg.table"),
                                                      br(),
                                                      br(),
                                                      br(),
                                                      plotOutput("Deg3.plot", width = "100%")
                                               )
                                               
                                      ),
                                     
                                     tabPanel("Data visualization",
                                              
                                              fluidRow(
                                                column(3,
                                                selectInput("deg2",
                                                            label = "Group by",
                                                            choices = c("seurat_clusters", "primary.predict"),
                                                            selected = "primary.predict"),
                                                ),
                                                column(6,
                                                       uiOutput("deg.gene.select"),
                                                       plotlyOutput("Deg.plot", width = "150%"),
                                                       br(),
                                                       plotlyOutput("Deg1.plot", width = "150%"),
                                                       br(),
                                                       plotOutput("Deg2.plot", width = "150%"),
                                                ),
                                              )
                                     ),
                                     
                                     tabPanel("Volcano Plot",
                                              br(),
                                              column(12,
                                              column(3,
                                                     selectInput("deg3",
                                                                 label = "Group by",
                                                                 choices = c("seurat_clusters", "primary.predict"),
                                                                 selected = "primary.predict"),
                                              ),
                                              column(3,
                                              uiOutput("gene1.select"),
                                              ),
                                              column(3,
                                              uiOutput("gene2.select"),
                                              ),
                                              br(),
                                              column(3,
                                                     actionButton("doVolcano", "Run Volcano plot", icon = icon("hand-pointer-o"))
                                              ),
                                            ),
                                            column(12,
                                                   column(3, numericInput("min_pct_a",
                                                                          label = "min.pct",
                                                                          value = 0.25,
                                                                          min = 0,
                                                                          step = 0.01)
                                                   ),
                                                   column(3, numericInput("logfc_a",
                                                                          label = "logfc.threshold",
                                                                          value = 0.25,
                                                                          min = 0,
                                                                          step = 0.01)
                                                   ),
                                                   column(3, selectInput("test.use_a",
                                                                         label = "Test use",
                                                                         choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                   )
                                                  ),
                                            column(12,
                                              plotOutput("volcano.plot", width = "75%"),
                                            ),
                                              br(),
                                            column(12,
                                              plotOutput("dega.plot", width = "100%")
                                            ),
                                              ),
                                     
                                     tabPanel("GSEA",
                                              column(12,
                                                     br(),
                                                     column(2, selectInput("species_gsea",
                                                                           label = "Species",
                                                                           choices = c("Homo sapiens", "Mus musculus"))
                                                     ),
                                                     column(2, selectInput("category_gsea",
                                                                           label = "Collection",
                                                                           choices = c("H", "C2", "C5", "C7", "C8"))
                                                     ),
                                                     column(2,
                                                            uiOutput("gsea.ct1.select"),
                                                     ),
                                                     column(2,
                                                            uiOutput("gsea.ct2.select"),
                                                     ),
                                              
                                                    
                                                     column(2, numericInput("min_pct1",
                                                                            label = "min.pct",
                                                                            value = 0.25,
                                                                            min = 0,
                                                                            step = 0.01)
                                                     ),
                                                     column(2, numericInput("logfc1",
                                                                            label = "logfc.threshold",
                                                                            value = 0.25,
                                                                            min = 0,
                                                                            step = 0.01)
                                                     ),
                                              ),
                                                     column(12,
                                                     column(2, selectInput("test.use1",
                                                                           label = "Test use",
                                                                           choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                     ),
                                                     column(4,
                                                            uiOutput("gsea.select"),
                                                     ),
                                                     br(),
                                                     column(3,
                                                     actionButton("gsea", "Run gene set enrichment analysis", icon = icon("hand-pointer-o")),
                                                     ),
                                                     textOutput("gsea.done"),
                                                     br()
                                              ),
                                              br(),
                                              column(7,
                                              plotOutput("gsea_plot", width = "100%"),
                                              ),
                                              br(),
                                              br(),
                                              column(12,
                                              plotOutput("gsea_plot1", width = "100%"),
                                              ),
                                              br(),
                                              br(),
                                              DT::dataTableOutput("gsea.table"),
                                              column(4,
                                                     downloadButton('download_gsea.table', 'Download Celltype similarity table (in csv)'),
                                              ),
                                     ),
                                     
                                     tabPanel("Cell-cell communication", 
                                              column(3,
                                                     br(),
                                                     actionButton("doCC", "Run analysis", icon = icon("hand-pointer-o")),
                                                     textOutput("cc.done"),
                                                     br()
                                              ),
                                              
                                              column(3, selectInput("cc_pval",
                                                                    label = "p-value",
                                                                    choices = c("0.01", "0.05"))
                                              ),
                                              column(3,
                                              uiOutput("CC.gene.select"),
                                              ),
                                              column(3,
                                              uiOutput("CC.gene1.select"),
                                              ),
                                              br(),
                                              br(),
                                              column(12,
                                              plotlyOutput("CC_plot1", width = "100%"),
                                              ),
                                              br(),
                                              column(12,
                                                     plotOutput("CC_plot2", width = "100%"),
                                              ),
                                              )
                                     ),
                       ),
               
               ##------------Data Integration module--------------------##
               
               tabPanel("Single cell data integration",
                        navlistPanel(widths=c(2,10),
                                     tabPanel("Overview",
                                              h2(p("Workflow for Data Integration module")),
                                              br(),
                                              imageOutput("integration_image"),         
                                     ),
                                     tabPanel("Upload your data",
                                              column(9,
                                                     column(5,
                                                            #h4('Load Data:'),
                                                            wellPanel(
                                                              titlePanel(h4(p("Load your input data"))),
                                                              br(),
                                                              selectInput("scInput1",
                                                                          label = "Select Data Input Type",
                                                                          choices = c("Raw Counts Matrix", "10X cellranger"),
                                                                          selected = "Raw Counts Matrix"),
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                fileInput("tpmFiles1",
                                                                          label = "Counts File (Accepted Format: text)",
                                                                          accept = ".txt", 
                                                                          multiple = T),
                                                                fileInput("cellAnnoFiles1",
                                                                          label = "Upload Metadata (Accepted Format: tab delimited text)",
                                                                          accept = ".txt",
                                                                          multiple = T),
                                                                actionButton("loadexample1", "Load example and run", icon = icon("hand-o-right")),
                                                                br(),
                                                                br(),
                                                                ),
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == '10X cellranger'",
                                                                fileInput("scH5_1",
                                                                          label = "Cellranger output (Accepted Format: .h5)",
                                                                          accept = ".h5",
                                                                          multiple = T),
                                                                actionButton("loadexample1a", "Load example and run", icon = icon("hand-o-right"))),
                                                                br(),
                                                                selectInput("scAnalysis_integ",
                                                                            label = "Analysis method",
                                                                            choices = c("Seurat", "Harmony", "fastMNN"),
                                                                            selected = "Seurat"),
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == 'Raw Counts Matrix' || input.scInput1 == '10X cellranger'",
                                                                textInput(inputId = "projName1",
                                                                          label = "Project Name",
                                                                          value = "Integration")),
                                                              
                                                              fluidRow(
                                                                
                                                                actionButton("loadButton1", "Load data", icon = icon("hand-o-right")),
                                                                actionButton("reset_intg", "Reset", icon = icon("repeat"))
                                                              ),
                                                            )),
                                                     
                                                     column(3,
                                                            numericInput(inputId = "min.genes1",
                                                                         label = "Min. genes",
                                                                         value = 200,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            numericInput(inputId = "min.cells1",
                                                                         label = "Min. cells",
                                                                         value = 3,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            actionButton("create_seurat1", "Process", icon = icon("hand-o-right"))
                                                     ),
                                                     
                                                     
                                              ),
                                              
                                              column(12,
                                                     withSpinner(dataTableOutput('countdataDT1'))
                                              )),
                                     
                                     tabPanel("Quality control",
                                          tabsetPanel(id="qc_intg",
                                            tabPanel("Violin Plot",            
                                              column(3,
                                                     plotOutput("nFeature_RNAPlot1")
                                              ),
                                              column(3,
                                                     plotOutput("mitoPlot1")
                                              ),
                                              column(3,
                                                     plotOutput("nCount_RNAPlot1")
                                              )),
                                            tabPanel("Feature Scatter", 
                                                     column(5,
                                                            plotlyOutput("FeatureScatterPlot1a")
                                                     ),
                                                     column(5,
                                                            plotlyOutput("FeatureScatterPlot2a")
                                                     ),        
                                            )),
                                          
                                          
                                     ),
                                     
                                     tabPanel("Normalization and Variable Gene Plot",
                                              textOutput("nVarGenes_bef_intg"),
                                              fluidRow(
                                                column(3,
                                                       numericInput("var.genes_bef_intg",
                                                                    label = "Number of variable genes",
                                                                    value = 2000,
                                                                    min = 500,
                                                                    step = 500)
                                                ),
                                                column(3,
                                                       selectInput("selection.method_bef_intg",
                                                                   label = "Selection method",
                                                                   choices = c("vst", "dispersion"))
                                                ),
                                                column(4,
                                                       br(),
                                                       actionButton("findVarGenes_bef_intg", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                       )),
                                              plotOutput("VarGenes_bef_intg", width = "100%")       
                                     ),
                                     
                                     tabPanel("Before Data Integration",
                                              
                                              tabsetPanel(id="bef_data_integration",
                                                          
                                                          tabPanel("PCA", 
                                                                   br(),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Seurat'",
                                                                     tabsetPanel(id="Pca_bef_intg_seurat",
                                                                                 tabPanel(title="PCA Plot",
                                                                                          br(),
                                                                                          fluidRow(
                                                                                            column(3,
                                                                                                   actionButton("runPCA_bef_intg_seurat", "Run PCA", icon = icon("hand-pointer-o"))
                                                                                            ),
                                                                                          ),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                                            plotlyOutput("PCAplot_bef_seurat_tpm1", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_seurat_tpm2", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_seurat_tpm3", width = "100%")),
                                                                                          br(),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == '10X cellranger'",
                                                                                            plotlyOutput("PCAplot_bef_seurat_h5_1", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_seurat_h5_2", width = "100%")),
                                                                                 ),
                                                                                 tabPanel(title="PC Gene Visualisation",
                                                                                          br(),
                                                                                          selectInput("select.pc_bef_intg_seurat",
                                                                                                      label = "PC to plot",
                                                                                                      choices = c(1:50)
                                                                                          ),
                                                                                          fluidRow(
                                                                                            column(4,
                                                                                                   plotOutput("vizPlot_bef_intg_seurat", width = "100%", height = "600px")
                                                                                            ),
                                                                                            column(8,
                                                                                                   plotOutput("PCHeatmap_bef_intg_seurat", width = "100%", height = "600px")
                                                                                            )
                                                                                          ),
                                                                                          DT::dataTableOutput("PCtable_bef_intg_seurat")
                                                                                 ),
                                                                                 
                                                                                 tabPanel(title="Elbow", 
                                                                                          br(),
                                                                                          #actionButton("doElbow", label = "Get Elbow Plot"),
                                                                                          #br(),
                                                                                          br(),
                                                                                          plotOutput("Elbow_bef_intg_seurat", width = "100%")
                                                                                          
                                                                                 ))),
                                                                   
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Harmony'",
                                                                     tabsetPanel(id="Pca_bef_intg_harmony",
                                                                                 tabPanel(title="PCA Plot",
                                                                                          br(),
                                                                                          fluidRow(
                                                                                            column(3,
                                                                                                   actionButton("runPCA_bef_intg_harmony", "Run PCA", icon = icon("hand-pointer-o"))
                                                                                            ),
                                                                                          ),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                                            plotlyOutput("PCAplot_bef_harmony_tpm1", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_harmony_tpm2", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_harmony_tpm3", width = "100%")),
                                                                                          br(),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == '10X cellranger'",
                                                                                            plotlyOutput("PCAplot_bef_harmony_h5_1", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_harmony_h5_2", width = "100%")),
                                                                                 ),
                                                                                 tabPanel(title="PC Gene Visualisation",
                                                                                          br(),
                                                                                          selectInput("select.pc_bef_intg_harmony",
                                                                                                      label = "PC to plot",
                                                                                                      choices = c(1:50)
                                                                                          ),
                                                                                          fluidRow(
                                                                                            column(4,
                                                                                                   plotOutput("vizPlot_bef_intg_harmony", width = "100%", height = "600px")
                                                                                            ),
                                                                                            column(8,
                                                                                                   plotOutput("PCHeatmap_bef_intg_harmony", width = "100%", height = "600px")
                                                                                            )
                                                                                          ),
                                                                                          DT::dataTableOutput("PCtable_bef_intg_harmony")
                                                                                 ),
                                                                                 
                                                                                 tabPanel(title="Elbow", 
                                                                                          br(),
                                                                                          br(),
                                                                                          plotOutput("Elbow_bef_intg_harmony", width = "100%")
                                                                                 ))),
                                                                      conditionalPanel(
                                                                            condition = "input.scAnalysis_integ == 'fastMNN'",
                                                                            tabsetPanel(id="Pca_bef_intg_fastmnn",
                                                                                 tabPanel(title="PCA Plot",
                                                                                          br(),
                                                                                          fluidRow(
                                                                                            column(3,
                                                                                                   actionButton("runPCA_bef_intg_fastmnn", "Run PCA", icon = icon("hand-pointer-o"))
                                                                                            ),
                                                                                          ),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                                            plotlyOutput("PCAplot_bef_fastmnn_tpm1", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_fastmnn_tpm2", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_fastmnn_tpm3", width = "100%")),
                                                                                          br(),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == '10X cellranger'",
                                                                                            plotlyOutput("PCAplot_bef_fastmnn_h5_1", width = "100%"),
                                                                                            plotlyOutput("PCAplot_bef_fastmnn_h5_2", width = "100%")),
                                                                                 ),
                                                                                 tabPanel(title="PC Gene Visualisation",
                                                                                          br(),
                                                                                          selectInput("select.pc_bef_intg_fastmnn",
                                                                                                      label = "PC to plot",
                                                                                                      choices = c(1:50)
                                                                                          ),
                                                                                          fluidRow(
                                                                                            column(4,
                                                                                                   plotOutput("vizPlot_bef_intg_fastmnn", width = "100%", height = "600px")
                                                                                            ),
                                                                                            column(8,
                                                                                                   plotOutput("PCHeatmap_bef_intg_fastmnn", width = "100%", height = "600px")
                                                                                            )
                                                                                          ),
                                                                                          DT::dataTableOutput("PCtable_bef_intg_fastmnn")
                                                                                 ),
                                                                                 
                                                                                 tabPanel(title="Elbow", 
                                                                                          br(),
                                                                                          br(),
                                                                                          plotOutput("Elbow_bef_intg_fastmnn", width = "100%")
                                                                                 )))),
                                                          
                                                          tabPanel("Clustering",
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Seurat'",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              numericInput("clus.res_bef_intg_seurat",
                                                                                           label = "Resolution used",
                                                                                           value = 0.6,
                                                                                           min = 0.1,
                                                                                           step = 0.1)
                                                                       ),
                                                                       column(3,
                                                                              selectInput("dim.used_bef_intg_seurat",
                                                                                          label = "PC to use",
                                                                                          choices = c(10:50)),
                                                                       ),
                                                                       column(3,
                                                                              br(),
                                                                              actionButton("findCluster_bef_intg_seurat", "Find Clusters", icon = icon("hand-pointer-o")),
                                                                              #textOutput("cluster1.done"),
                                                                       )),
                                                                     br(),
                                                                     plotlyOutput("Cluster2DPlot_bef_intg_seurat", width = "100%")
                                                                     
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Harmony'",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              numericInput("clus.res_bef_intg_harmony",
                                                                                           label = "Resolution used",
                                                                                           value = 0.6,
                                                                                           min = 0.1,
                                                                                           step = 0.1)
                                                                       ),
                                                                       column(3,
                                                                              selectInput("dim.used_bef_intg_harmony",
                                                                                          label = "PC to use",
                                                                                          choices = c(10:50)),
                                                                       ),
                                                                       column(3,
                                                                              br(),
                                                                              actionButton("findCluster_bef_intg_harmony", "Find Clusters", icon = icon("hand-pointer-o")),
                                                                              #textOutput("cluster2.done"),
                                                                       )),
                                                                     br(),
                                                                     plotlyOutput("Cluster2DPlot_bef_intg_harmony", width = "100%")
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'fastMNN'",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              numericInput("clus.res_bef_intg_fastmnn",
                                                                                           label = "Resolution used",
                                                                                           value = 0.6,
                                                                                           min = 0.1,
                                                                                           step = 0.1)
                                                                       ),
                                                                       column(3,
                                                                              selectInput("dim.used_bef_intg_fastmnn",
                                                                                          label = "PC to use",
                                                                                          choices = c(10:50)),
                                                                       ),
                                                                       column(3,
                                                                              br(),
                                                                              actionButton("findCluster_bef_intg_fastmnn", "Find Clusters", icon = icon("hand-pointer-o")),
                                                                              #textOutput("cluster1.done"),
                                                                       )),
                                                                     br(),
                                                                     plotlyOutput("Cluster2DPlot_bef_intg_fastmnn", width = "100%")
                                                                     
                                                                   )
                                                          ),
                                                          
                                                          tabPanel("UMAP",
                                                                   
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Seurat'",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg_seurat",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionButton("runUMAP_bef_intg_seurat", "Run UMAP", icon = icon("hand-pointer-o")),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("UMAPplot_bef_seurat_tpm1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_seurat_tpm2", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_seurat_tpm3", width = "100%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("UMAPplot_bef_seurat_h5_1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_seurat_h5_2", width = "100%")),
                                                                            
                                                                     )),
                                                                   
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Harmony'",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg_harmony",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionButton("runUMAP_bef_intg_harmony", "Run UMAP", icon = icon("hand-pointer-o")),
                                                                            
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("UMAPplot_bef_harmony_tpm1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_harmony_tpm2", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_harmony_tpm3", width = "100%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("UMAPplot_bef_harmony_h5_1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_harmony_h5_2", width = "100%")),
                                                                     )),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'fastMNN'",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg_fastmnn",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionButton("runUMAP_bef_intg_fastmnn", "Run UMAP", icon = icon("hand-pointer-o")),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("UMAPplot_bef_fastmnn_tpm1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_fastmnn_tpm2", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_fastmnn_tpm3", width = "100%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("UMAPplot_bef_fastmnn_h5_1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_fastmnn_h5_2", width = "100%")),
                                                                            
                                                                     ))
                                                                   
                                                          ),
                                                          
                                                          tabPanel("tSNE",
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Seurat'",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg_seurat",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionButton("runTSNE_bef_intg_seurat", "Run TSNE", icon = icon("hand-pointer-o")),
                                                                            #textOutput("Intg.done"),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("TSNEplot_bef_seurat_tpm1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_seurat_tpm2", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_seurat_tpm3", width = "100%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("TSNEplot_bef_seurat_h5_1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_seurat_h5_2", width = "100%")),
                                                                     )),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Harmony'",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg_harmony",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionButton("runTSNE_bef_intg_harmony", "Run TSNE", icon = icon("hand-pointer-o")),
                                                                            #textOutput("Intg.done"),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("TSNEplot_bef_harmony_tpm1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_harmony_tpm2", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_harmony_tpm3", width = "100%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("TSNEplot_bef_harmony_h5_1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_harmony_h5_2", width = "100%")),
                                                                     )),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'fastMNN'",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg_fastmnn",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionButton("runTSNE_bef_intg_fastmnn", "Run TSNE", icon = icon("hand-pointer-o")),
                                                                            #textOutput("Intg.done"),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("TSNEplot_bef_fastmnn_tpm1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_fastmnn_tpm2", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_fastmnn_tpm3", width = "100%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("TSNEplot_bef_fastmnn_h5_1", width = "100%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_fastmnn_h5_2", width = "100%")),
                                                                     ))
                                                                  ),
                                                          
                                                          tabPanel("Cell type identification",
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Seurat'",
                                                                     br(),
                                                                     
                                                                     column(3,
                                                                            selectInput("cellatlas1",
                                                                                        label = "Reference Atlas",
                                                                                        choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                                        selected = "all")
                                                                     ),
                                                                     column(3,
                                                                            br(),
                                                                            actionButton("doCELLiD_bef_intg_seurat", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                                            #textOutput("CELLiD.done"),
                                                                            br()
                                                                     ),
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     plotlyOutput("Umap_cellid_bef_intg_seurat", width = "50%"),
                                                                     plotlyOutput("Umap_cellid_bef_intg_seurat1", width = "50%"),
                                                                     br(),
                                                                     br(),
                                                                     DT::dataTableOutput("ct_bef_intg_seurat.table")
                                                                     ),
                                                                   
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'Harmony'",
                                                                     column(3,
                                                                            selectInput("cellatlas2",
                                                                                        label = "Reference Atlas",
                                                                                        choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                                        selected = "all")
                                                                     ),
                                                                     column(3,
                                                                            br(),
                                                                            actionButton("doCELLiD_bef_intg_harmony", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                                            #textOutput("CELLiD.done"),
                                                                            br()
                                                                     ),
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     plotlyOutput("Umap_cellid_bef_intg_harmony", width = "50%"),
                                                                     plotlyOutput("Umap_cellid_bef_intg_harmony1", width = "50%"),
                                                                     br(),
                                                                     br(),
                                                                     DT::dataTableOutput("ct_bef_intg_harmony.table")
                                                                     ),
                                                                   conditionalPanel(
                                                                     condition = "input.scAnalysis_integ == 'fastMNN'",
                                                                     column(3,
                                                                            selectInput("cellatlas3",
                                                                                        label = "Reference Atlas",
                                                                                        choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                                        selected = "all")
                                                                     ),
                                                                     column(3,
                                                                            br(),
                                                                            actionButton("doCELLiD_bef_intg_fastmnn", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                                            #textOutput("CELLiD.done"),
                                                                            br()
                                                                     ),
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     plotlyOutput("Umap_cellid_bef_intg_fastmnn", width = "100%"),
                                                                     plotlyOutput("Umap_cellid_bef_intg_fastmnn1", width = "100%"),
                                                                     br(),
                                                                     br(),
                                                                     DT::dataTableOutput("ct_bef_intg_fastmnn.table")
                                                                     )
                                                            ),
                                                          )
                                     ),
                                              
                                     tabPanel("Data Integration", 
                                       fluidRow(
                                         column(12,
                                                br(),
                                                conditionalPanel(
                                                  condition = "input.scAnalysis_integ == 'Seurat'",
                                                column(5,
                                                numericInput("nfeatures_intg_seurat",
                                                               label = "Integration features",
                                                               value = 2000,
                                                               min = 500,
                                                               step = 1)),
                                                br(),
                                                column(5,
                                                actionButton("doIntg_seurat", "Run Data Integration", icon = icon("hand-pointer-o")),
                                                ),
                                                br(),
                                         ),
                                         conditionalPanel(
                                           condition = "input.scAnalysis_integ == 'Harmony'",
                                           column(5,
                                           numericInput("nfeatures_intg_harmony",
                                                        label = "Integration features",
                                                        value = 2000,
                                                        min = 500,
                                                        step = 1)),
                                           br(),
                                           column(5,
                                           actionButton("doIntg_harmony", "Running Data Integration", icon = icon("hand-pointer-o")),
                                           ),
                                         ),
                                         conditionalPanel(
                                           condition = "input.scAnalysis_integ == 'fastMNN'",
                                           column(5,
                                           numericInput("nfeatures_intg_fastmnn",
                                                        label = "Integration features",
                                                        value = 2000,
                                                        min = 500,
                                                        step = 1)),
                                           br(),
                                           column(5,
                                           actionButton("doIntg_fastmnn", "Running Data Integration", icon = icon("hand-pointer-o")),
                                           )),
                                         
                                         ))),
                                     tabPanel("PCA (After Data Integration)", 
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'Seurat'",
                                              tabsetPanel(id="Pca1",
                                                          tabPanel(title="PCA Plot",
                                                                   br(),
                                                                   conditionalPanel(
                                                                     condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                   plotlyOutput("PCAplot_seurat_tpm1", width = "100%"),
                                                                   plotlyOutput("PCAplot_seurat_tpm2", width = "100%"),
                                                                   plotlyOutput("PCAplot_seurat_tpm3", width = "100%")),
                                                                   br(),
                                                                   conditionalPanel(
                                                                     condition = "input.scInput1 == '10X cellranger'",
                                                                     plotlyOutput("PCAplot_seurat_h5_1", width = "100%"),
                                                                     plotlyOutput("PCAplot_seurat_h5_2", width = "100%")),
                                                          ),
                                                          tabPanel(title="PC Gene Visualisation",
                                                                   br(),
                                                                   selectInput("select.pc_intg_seurat",
                                                                               label = "PC to plot",
                                                                               choices = c(1:50)
                                                                   ),
                                                                   fluidRow(
                                                                     column(4,
                                                                            plotOutput("vizPlot_intg_seurat", width = "100%", height = "600px")
                                                                     ),
                                                                     column(8,
                                                                            plotOutput("PCHeatmap_intg_seurat", width = "100%", height = "600px")
                                                                     )
                                                                   ),
                                                                   DT::dataTableOutput("PCtable_intg_seurat")
                                                          ),
                                                          
                                                          tabPanel(title="Elbow", value="P_panel4",
                                                                   br(),
                                                                   #actionButton("doElbow", label = "Get Elbow Plot"),
                                                                   #br(),
                                                                   br(),
                                                                   plotOutput("Elbow_intg_seurat", width = "100%")
                                                          ))),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'Harmony'",
                                                tabsetPanel(id="Pca2",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                     plotlyOutput("PCAplot_harmony_tpm1", width = "100%"),
                                                                     plotlyOutput("PCAplot_harmony_tpm2", width = "100%"),
                                                                     plotlyOutput("PCAplot_harmony_tpm3", width = "100%")),
                                                                     br(),
                                                                   conditionalPanel(
                                                                     condition = "input.scInput1 == '10X cellranger'",
                                                                     plotlyOutput("PCAplot_harmony_h5_1", width = "100%"),
                                                                     plotlyOutput("PCAplot_harmony_h5_2", width = "100%")),
                                                            ),
                                                            tabPanel(title="PC Gene Visualisation",
                                                                     br(),
                                                                     selectInput("select.pc_intg_harmony",
                                                                                 label = "PC to plot",
                                                                                 choices = c(1:50)
                                                                     ),
                                                                     fluidRow(
                                                                       column(4,
                                                                              plotOutput("vizPlot_intg_harmony", width = "100%", height = "600px")
                                                                       ),
                                                                       column(8,
                                                                              plotOutput("PCHeatmap_intg_harmony", width = "100%", height = "600px")
                                                                       )
                                                                     ),
                                                                     DT::dataTableOutput("PCtable_intg_harmony")
                                                            ),
                                                            
                                                            tabPanel(title="Elbow", value="P_panel4",
                                                                     br(),
                                                                     br(),
                                                                     plotOutput("Elbow_intg_harmony", width = "100%")
                                                            ))),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'fastMNN'",
                                                tabsetPanel(id="Pca2",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                       plotlyOutput("PCAplot_fastmnn_tpm1", width = "100%"),
                                                                       plotlyOutput("PCAplot_fastmnn_tpm2", width = "100%"),
                                                                       plotlyOutput("PCAplot_fastmnn_tpm3", width = "100%")),
                                                                     br(),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == '10X cellranger'",
                                                                       plotlyOutput("PCAplot_fastmnn_h5_1", width = "100%"),
                                                                       plotlyOutput("PCAplot_fastmnn_h5_2", width = "100%")),
                                                            ),
                                                            tabPanel(title="PC Gene Visualisation",
                                                                     br(),
                                                                     selectInput("select.pc_intg_fastmnn",
                                                                                 label = "PC to plot",
                                                                                 choices = c(1:50)
                                                                     ),
                                                                     fluidRow(
                                                                       column(4,
                                                                              plotOutput("vizPlot_intg_fastmnn", width = "100%", height = "600px")
                                                                       ),
                                                                     column(8,
                                                                     DT::dataTableOutput("PCtable_intg_fastmnn")
                                                                     )),
                                                                  ),
                                                            
                                                            ))
                                              ),
                                     
                                     tabPanel("Clustering (After Data Integration)",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'Seurat'",
                                              fluidRow(
                                                column(3,
                                                       numericInput("clus.res_intg_seurat",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       selectInput("dim.used_intg_seurat",
                                                                   label = "PC to use",
                                                                   choices = c(10:50)),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("findCluster_intg_seurat", "Find Clusters", icon = icon("hand-pointer-o")),
                                                       textOutput("cluster1.done"),
                                                )),
                                              br(),
                                              plotlyOutput("Cluster2DPlot_intg_seurat", width = "100%")
                                              
                                            ),
                                            conditionalPanel(
                                              condition = "input.scAnalysis_integ == 'Harmony'",
                                              fluidRow(
                                                column(3,
                                                       numericInput("clus.res_intg_harmony",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       selectInput("dim.used_intg_harmony",
                                                                   label = "PC to use",
                                                                   choices = c(10:50)),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("findCluster_intg_harmony", "Find Clusters", icon = icon("hand-pointer-o")),
                                                       textOutput("cluster2.done"),
                                                )),
                                              br(),
                                              plotlyOutput("Cluster2DPlot_intg_harmony", width = "100%")
                                              ),
                                            conditionalPanel(
                                              condition = "input.scAnalysis_integ == 'fastMNN'",
                                              fluidRow(
                                                column(3,
                                                       numericInput("clus.res_intg_fastmnn",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       selectInput("dim.used_intg_fastmnn",
                                                                   label = "PC to use",
                                                                   choices = c(10:50)),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("findCluster_intg_fastmnn", "Find Clusters", icon = icon("hand-pointer-o")),
                                                       #textOutput("cluster2.done"),
                                                )),
                                              br(),
                                              plotlyOutput("Cluster2DPlot_intg_fastmnn", width = "100%")
                                            )
                                     ),
                                     
                                     tabPanel("UMAP (After Data Integration)", 
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'Seurat'",
                                                column(3,
                                                       numericInput("dim.used_intg_seurat",
                                                                    label = "Dimensions used",
                                                                    value = 10)
                                                ),
                                       column(9,
                                              br(),
                                              actionButton("runUMAP_intg_seurat", "Run UMAP", icon = icon("hand-pointer-o")),
                                              br(),
                                              br(),
                                              conditionalPanel(
                                                condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                plotlyOutput("UMAPplot_seurat_tpm1", width = "100%"),
                                                br(),
                                                plotlyOutput("UMAPplot_seurat_tpm2", width = "100%"),
                                                br(),
                                                plotlyOutput("UMAPplot_seurat_tpm3", width = "100%"),
                                                br(),
                                                textOutput("UMAP_lisi_seurat")),
                                              br(),
                                              conditionalPanel(
                                                condition = "input.scInput1 == '10X cellranger'",
                                                plotOutput("UMAPplot_seurat_h5_1", width = "100%"),
                                                br(),
                                                plotOutput("UMAPplot_seurat_h5_2", width = "100%")),
                                      )),
                        
                        conditionalPanel(
                          condition = "input.scAnalysis_integ == 'Harmony'",
                          column(3,
                                 numericInput("dim.used_intg_harmony",
                                              label = "Dimensions used",
                                              value = 10)
                          ),
                          column(9,
                                 br(),
                                 actionButton("runUMAP_intg_harmony", "Run UMAP", icon = icon("hand-pointer-o")),
                          
                                 br(),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == 'Raw Counts Matrix'",
                                   plotlyOutput("UMAPplot_harmony_tpm1", width = "100%"),
                                   br(),
                                   plotlyOutput("UMAPplot_harmony_tpm2", width = "100%"),
                                   br(),
                                   plotlyOutput("UMAPplot_harmony_tpm3", width = "100%"),
                                   br(),
                                   textOutput("UMAP_lisi_harmony")),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == '10X cellranger'",
                                   plotlyOutput("UMAPplot_harmony_h5_1", width = "100%"),
                                   br(),
                                   plotlyOutput("UMAPplot_harmony_h5_2", width = "100%")),
                          )),
                        conditionalPanel(
                          condition = "input.scAnalysis_integ == 'fastMNN'",
                          column(3,
                                 numericInput("dim.used_intg_fastmnn",
                                              label = "Dimensions used",
                                              value = 10)
                          ),
                          column(9,
                                 br(),
                                 actionButton("runUMAP_intg_fastmnn", "Run UMAP", icon = icon("hand-pointer-o")),
                                 
                                 br(),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == 'Raw Counts Matrix'",
                                   plotlyOutput("UMAPplot_fastmnn_tpm1", width = "100%"),
                                   br(),
                                   plotlyOutput("UMAPplot_fastmnn_tpm2", width = "100%"),
                                   br(),
                                   plotlyOutput("UMAPplot_fastmnn_tpm3", width = "100%"),
                                   br(),
                                   textOutput("UMAP_lisi_fastmnn")),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == '10X cellranger'",
                                   plotlyOutput("UMAPplot_fastmnn_h5_1", width = "100%"),
                                   br(),
                                   plotlyOutput("UMAPplot_fastmnn_h5_2", width = "100%")),
                          ))
                        ),
                        
                        tabPanel("tSNE (After Data Integration)", 
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Seurat'",
                                   column(3,
                                          numericInput("dim.used_intg_seurat",
                                                       label = "Dimensions used",
                                                       value = 10)
                                   ),
                          column(9,
                                 br(),
                                 actionButton("runTSNE_intg_seurat", "Run TSNE", icon = icon("hand-pointer-o")),
                                 #textOutput("Intg.done"),
                                 br(),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == 'Raw Counts Matrix'",
                                   plotlyOutput("TSNEplot_seurat_tpm1", width = "100%"),
                                   br(),
                                   plotlyOutput("TSNEplot_seurat_tpm2", width = "100%"),
                                   br(),
                                   plotlyOutput("TSNEplot_seurat_tpm3", width = "100%")),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == '10X cellranger'",
                                   plotlyOutput("TSNEplot_seurat_h5_1", width = "100%"),
                                   br(),
                                   plotlyOutput("TSNEplot_seurat_h5_2", width = "100%")),
                          )),
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Harmony'",
                                   column(3,
                                          numericInput("dim.used_intg_harmony",
                                                       label = "Dimensions used",
                                                       value = 10)
                                   ),
                          column(9,
                                 br(),
                                 actionButton("runTSNE_intg_harmony", "Run TSNE", icon = icon("hand-pointer-o")),
                                 #textOutput("Intg.done"),
                                 br(),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == 'Raw Counts Matrix'",
                                   plotlyOutput("TSNEplot_harmony_tpm1", width = "100%"),
                                   br(),
                                   plotlyOutput("TSNEplot_harmony_tpm2", width = "100%"),
                                   br(),
                                   plotlyOutput("TSNEplot_harmony_tpm3", width = "100%")),
                                 br(),
                                 conditionalPanel(
                                   condition = "input.scInput1 == '10X cellranger'",
                                   plotlyOutput("TSNEplot_harmony_h5_1", width = "100%"),
                                   br(),
                                   plotlyOutput("TSNEplot_harmony_h5_2", width = "100%")),
                          )),
                          conditionalPanel(
                            condition = "input.scAnalysis_integ == 'fastMNN'",
                            column(3,
                                   numericInput("dim.used_intg_fastmnn",
                                                label = "Dimensions used",
                                                value = 10)
                            ),
                            column(9,
                                   br(),
                                   actionButton("runTSNE_intg_fastmnn", "Run TSNE", icon = icon("hand-pointer-o")),
                                   #textOutput("Intg.done"),
                                   br(),
                                   br(),
                                   conditionalPanel(
                                     condition = "input.scInput1 == 'Raw Counts Matrix'",
                                     plotlyOutput("TSNEplot_fastmnn_tpm1", width = "100%"),
                                     br(),
                                     plotlyOutput("TSNEplot_fastmnn_tpm2", width = "100%"),
                                     br(),
                                     plotlyOutput("TSNEplot_fastmnn_tpm3", width = "100%")),
                                   br(),
                                   conditionalPanel(
                                     condition = "input.scInput1 == '10X cellranger'",
                                     plotlyOutput("TSNEplot_fastmnn_h5_1", width = "100%"),
                                     br(),
                                     plotlyOutput("TSNEplot_fastmnn_h5_2", width = "100%")),
                            ))
                          ),
                        
                        tabPanel("Cell type identification (After Data Integration)",
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Seurat'",
                                   column(3,
                                          selectInput("cellatlas1a",
                                                      label = "Reference Atlas",
                                                      choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                      selected = "all")
                                   ),
                                 column(3,
                                        br(),
                                        actionButton("doCELLiD_intg_seurat", "Run CELLiD", icon = icon("hand-pointer-o")),
                                        #textOutput("CELLiD.done"),
                                        br()
                                 ),
                                 br(),
                                 br(),
                                 br(),
                                 plotlyOutput("Umap_cellid_intg_seurat", width = "50%"),
                                 plotlyOutput("Umap_cellid_intg_seurat1", width = "50%"),
                                 br(),
                                 br(),
                                 DT::dataTableOutput("ct_intg_seurat.table")
                                 ),
                                 
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Harmony'",
                                   column(3,
                                          selectInput("cellatlas1b",
                                                      label = "Reference Atlas",
                                                      choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                      selected = "all")
                                   ),
                                   column(3,
                                          br(),
                                          actionButton("doCELLiD_intg_harmony", "Run CELLiD", icon = icon("hand-pointer-o")),
                                          #textOutput("CELLiD.done"),
                                          br()
                                   ),
                                   br(),
                                   br(),
                                   br(),
                                   plotlyOutput("Umap_cellid_intg_harmony", width = "50%"),
                                   plotlyOutput("Umap_cellid_intg_harmony1", width = "50%"),
                                   br(),
                                   br(),
                                   DT::dataTableOutput("ct_intg_harmony.table")
                                   ),
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'fastMNN'",
                                   column(3,
                                          selectInput("cellatlas1c",
                                                      label = "Reference Atlas",
                                                      choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                      selected = "all")
                                   ),
                                   column(3,
                                          br(),
                                          actionButton("doCELLiD_intg_fastmnn", "Run CELLiD", icon = icon("hand-pointer-o")),
                                          #textOutput("CELLiD.done"),
                                          br()
                                   ),
                                   br(),
                                   br(),
                                   br(),
                                   plotlyOutput("Umap_cellid_intg_fastmnn", width = "50%"),
                                   plotlyOutput("Umap_cellid_intg_fastmnn1", width = "50%"),
                                   br(),
                                   br(),
                                   DT::dataTableOutput("ct_intg_fastmnn.table")
                                   )
                              ),
                        
                        tabPanel("DEGs",
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Seurat'",
                                 fluidRow(
                                   column(3, selectInput("min_pct_intg_seurat",
                                                         label = "min.pct",
                                                         choices = c("0.1", "0.25"))
                                   ),
                                   
                                   column(3, selectInput("logfc_intg_seurat",
                                                         label = "logfc.threshold",
                                                         choices = c("0.1", "0.25"))
                                   ),
                                   
                                   column(3, selectInput("test.use_intg_seurat",
                                                         label = "Test use",
                                                         choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                   ),
                                   
                                   column(4,
                                          actionButton("doDeg_intg_seurat", "Run DEGs", icon = icon("hand-pointer-o"))
                                   )),
                                 br(),
                                 column(12,
                                        DT::dataTableOutput("Deg.table_intg_seurat")
                                 )),
                                 
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Harmony'",
                                   fluidRow(
                                     column(3, selectInput("min_pct_intg_harmony",
                                                           label = "min.pct",
                                                           choices = c("0.1", "0.25"))
                                     ),
                                     
                                     column(3, selectInput("logfc_intg_harmony",
                                                           label = "logfc.threshold",
                                                           choices = c("0.1", "0.25"))
                                     ),
                                     
                                     column(3, selectInput("test.use_intg_harmony",
                                                           label = "Test use",
                                                           choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                     ),
                                     
                                     column(4,
                                            actionButton("doDeg_intg_harmony", "Run DEGs", icon = icon("hand-pointer-o"))
                                     )),
                                   br(),
                                   column(12,
                                          DT::dataTableOutput("Deg.table_intg_harmony")
                                   )),
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'fastMNN'",
                                   fluidRow(
                                     column(3, selectInput("min_pct_intg_fastmnn",
                                                           label = "min.pct",
                                                           choices = c("0.1", "0.25"))
                                     ),
                                     
                                     column(3, selectInput("logfc_intg_fastmnn",
                                                           label = "logfc.threshold",
                                                           choices = c("0.1", "0.25"))
                                     ),
                                     
                                     column(3, selectInput("test.use_intg_fastmnn",
                                                           label = "Test use",
                                                           choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                     ),
                                     
                                     column(4,
                                            actionButton("doDeg_intg_fastmnn", "Run DEGs", icon = icon("hand-pointer-o"))
                                     )),
                                   br(),
                                   column(12,
                                          DT::dataTableOutput("Deg.table_intg_fastmnn")
                                          # br(),
                                          # plotlyOutput("Deg3.plot_intg_seurat", width = "100%")
                                   ))
                        ),
                        
                        tabPanel("Data Visualization",
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Seurat'",
                                   
                                     column(6,
                                            uiOutput("deg.gene.select_intg_seurat"),
                                            plotlyOutput("Deg.plot_intg_seurat", width = "200%"),
                                            br(),
                                            plotlyOutput("Deg1.plot_intg_seurat", width = "100%"),
                                            br(),
                                            plotOutput("Deg2.plot_intg_seurat", width = "200%")
                                     )
                                   ),
                                 
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'Harmony'",
                                   
                                   column(6,
                                          uiOutput("deg.gene.select_intg_harmony"),
                                          plotlyOutput("Deg.plot_intg_harmony", width = "200%"),
                                          br(),
                                          plotlyOutput("Deg1.plot_intg_harmony", width = "100%"),
                                          br(),
                                          plotOutput("Deg2.plot_intg_harmony", width = "200%")
                                   )
                                 ),
                                 conditionalPanel(
                                   condition = "input.scAnalysis_integ == 'fastMNN'",
                                   
                                   column(6,
                                          uiOutput("deg.gene.select_intg_fastmnn"),
                                          plotlyOutput("Deg.plot_intg_fastmnn", width = "200%"),
                                          br(),
                                          plotlyOutput("Deg1.plot_intg_fastmnn", width = "100%"),
                                          br(),
                                          plotOutput("Deg2.plot_intg_fastmnn", width = "200%")
                                   )
                                 )
                        ),
                      )),
               
               ##------------Spatial Transcriptomics module--------------------##
               
               tabPanel("Spatial Transcriptomics",
                        
                        navlistPanel(widths=c(2,10),
                                     tabPanel("Overview",
                                              h2(p("Workflow for Spatial Transcriptomics module")),
                                              br(),
                                              imageOutput("spatial_image"),         
                                     ),
                                     tabPanel("Upload your data",
                                              column(9,
                                                     column(5,
                                                            #h4('Load Data:'),
                                                            wellPanel(
                                                              titlePanel(h4(p("Load your input data"))),
                                                              br(),
                                                              selectInput("scInput3",
                                                                          label = "Select Data Input Type",
                                                                          choices = c("SpaceRanger output"),
                                                                          selected = "SpaceRanger output"),
                                                              selectInput("scAnalysis_platform",
                                                                          label = "Platform",
                                                                          choices = c("Visium", "Xenium"),
                                                                          selected = "Visium"),
                                                              selectInput("scAnalysis_sp",
                                                                          label = "Analysis method",
                                                                          choices = c("Seurat", "GraphST"),
                                                                          selected = "Seurat"),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                                fileInput("tpmFile_spatial",
                                                                          label = "Upload H5 output from Spaceranger (Accepted Format: .h5)",
                                                                          accept = ".h5"),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                                shinyFiles::shinyDirButton(id = 'dir', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                                                                verbatimTextOutput("dir", placeholder = TRUE),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                                shinyFiles::shinyDirButton(id = 'dir_xenium', label = "Path to xenium files", title = "Sheets Folder Selector"),
                                                                verbatimTextOutput("dir_xenium", placeholder = TRUE),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                                actionButton("loadexample_seurat_spatial", "Load example and run", icon = icon("hand-o-right")),
                                                              
                                                              textInput(inputId = "projName3",
                                                                        label = "Project Name",
                                                                        value = "Spatial"),
                                                              fluidRow(
                                                                actionButton("loadButton3", "Load Data", icon = icon("hand-o-right")),
                                                                actionButton("reset3", "Reset Data", icon = icon("repeat")),
                                                              )),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                                actionButton("loadexample_xenium", "Load example and run", icon = icon("hand-o-right")),
                                                                
                                                                textInput(inputId = "projName_xenium",
                                                                          label = "Project Name",
                                                                          value = "Xenium"),
                                                                fluidRow(
                                                                actionButton("load_xenium", "Load Data", icon = icon("hand-pointer-o")),
                                                                actionButton("reset_xenium", "Reset Data", icon = icon("repeat")),
                                                                )),
                                                              )),
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_platform == 'Visium'",
                                                     column(3,
                                                            numericInput(inputId = "min.genes_sp",
                                                                         label = "Min. genes",
                                                                         value = 200,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            numericInput(inputId = "min.cells_sp",
                                                                         label = "Min. cells",
                                                                         value = 3,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            actionButton("create_seurat_sp", "Process", icon = icon("hand-o-right"))
                                                     )
                                                    ), 
                                                    conditionalPanel(
                                                      condition = "input.scAnalysis_platform == 'Xenium'",
                                                      column(3,
                                                             br(),
                                                             numericInput("count_xenium",
                                                                          label = "nCount_Xenium",
                                                                          value = 0,
                                                                          min = 0,
                                                                          step = 1)
                                                      ),
                                                    ), 
                                              ),
                                              
                                              column(12,
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_platform == 'Visium'",
                                                       h4(p("Spatial counts")),
                                                       withSpinner(dataTableOutput('countdataDT_spatial')),
                                                       h4(p("H&E Image")),
                                                       withSpinner(plotOutput("h_e", width = "100%")),
                                                     ),
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_platform == 'Xenium'",
                                                       h4(p("Spatial counts")),
                                                       withSpinner(dataTableOutput('countdataDT_xenium')),
                                                       h4(p("H&E Image")),
                                                       withSpinner(plotOutput("h_e_xenium", width = "100%")),
                                                     ),
                                                    ),
                                     ),
                                     tabPanel("Quality Control",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                column(5,
                                                       plotOutput("nCount_SpatialPlot")
                                                ),
                                                column(5,
                                                       plotOutput("SpatialFeaturePlot")
                                                ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                              column(12,
                                                     column(5,
                                                            plotOutput("nFeature_xenium")
                                                     ),
                                                     column(5,
                                                            plotOutput("nCount_xenium")
                                                     ),
                                              )),
                                            ),
                                     tabPanel("FeatureScatter Plot",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                column(6,
                                                       plotlyOutput("FeatureScatterPlot_sp")
                                                ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                column(6,
                                                       plotlyOutput("FeatureScatterPlot_xenium")
                                                ),
                                              ),
                                              
                                     ),
                                     tabPanel("Normalization and Variable Gene Plot",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                fluidRow(
                                                  column(3,
                                                         numericInput("var.genes_sp",
                                                                      label = "Number of variable genes",
                                                                      value = 2000,
                                                                      min = 500,
                                                                      step = 500)
                                                  ),
                                                  column(4,
                                                         br(),
                                                         actionButton("doSCTransform_sp", "Run scTransform", icon = icon("hand-pointer-o")),
                                                         #actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                         
                                                  )),
                                                plotOutput("VarGenes_sp", width = "100%"),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                fluidRow(
                                                  column(3,
                                                         numericInput("var.genes_xenium",
                                                                      label = "Number of variable genes",
                                                                      value = 2000,
                                                                      min = 500,
                                                                      step = 500)
                                                  ),
                                                  column(4,
                                                         br(),
                                                         actionButton("doSCTransform_xenium", "Run scTransform", icon = icon("hand-pointer-o")),
                                                         #actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                         
                                                  )),
                                                plotOutput("VarGenes_xenium", width = "100%"),
                                              ),
                                            ),
                                     tabPanel("Gene expression visualization",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                column(4,
                                                       actionButton("Vis_sp", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                
                                                fluidRow(
                                                  column(6,
                                                         uiOutput("sp.gene.select"),
                                                         plotOutput("sp.plot", width = "100%"),
                                                  ),
                                                )),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                              column(6,
                                                     uiOutput("xenium.gene.select"),
                                                     plotOutput("markerplot_xenium")
                                              ),
                                              column(6,
                                                     uiOutput("xenium.gene1.select"),
                                                     plotOutput("featureplot_xenium")
                                              ),
                                              
                                              br(),
                                              br(),
                                              column(2,
                                                     numericInput("x1_xenium",
                                                                  label = "x1",
                                                                  value = 1200,
                                                                  min = 0,
                                                                  step = 1),
                                              ),
                                              column(2,
                                                     numericInput("x2_xenium",
                                                                  label = "x2",
                                                                  value = 2800,
                                                                  min = 0,
                                                                  step = 1),
                                              ),
                                              column(2,
                                                     numericInput("y1_xenium",
                                                                  label = "y1",
                                                                  value = 3700,
                                                                  min = 0,
                                                                  step = 1),
                                              ),
                                              column(2,
                                                     numericInput("y2_xenium",
                                                                  label = "y2",
                                                                  value = 4500,
                                                                  min = 0,
                                                                  step = 1),
                                              ),
                                              column(2,
                                                     br(),
                                              actionButton("crop_xenium", "Crop", icon = icon("hand-o-right")),
                                                     ),
                                              fluidRow(
                                                column(6,
                                                       uiOutput("xenium.gene2.select"),
                                                       plotOutput("zoom_xenium")
                                                ),
                                              ),
                                          ),
                                     ),
                                     tabPanel("PCA",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                tabsetPanel(id="Pca_sp",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              actionButton("runPCA_spatial", "Run PCA", icon = icon("hand-pointer-o"))
                                                                       ),
                                                                     ),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_spatial", width = "50%"),
                                                                     
                                                            ),
                                                            tabPanel(title="PC Gene Visualisation",
                                                                     br(),
                                                                     selectInput("select.pc_spatial",
                                                                                 label = "PC to plot",
                                                                                 choices = c(1:50)
                                                                     ),
                                                                     fluidRow(
                                                                       column(4,
                                                                              plotOutput("vizPlot_spatial", width = "100%", height = "600px")
                                                                       ),
                                                                       column(8,
                                                                              plotOutput("PCHeatmap_spatial", width = "100%", height = "600px")
                                                                       )
                                                                     ),
                                                                     DT::dataTableOutput("PCtable_spatial")
                                                            ),
                                                            
                                                            tabPanel(title="Elbow", 
                                                                     br(),
                                                                     br(),
                                                                     plotOutput("Elbow_spatial", width = "100%")
                                                            ))),
                                                conditionalPanel(
                                                  condition = "input.scAnalysis_platform == 'Xenium'",
                                                  tabsetPanel(id="Pca_xenium",
                                                              tabPanel(title="PCA Plot",
                                                                       br(),
                                                                       fluidRow(
                                                                         column(3,
                                                                                actionButton("runPCA_xenium", "Run PCA", icon = icon("hand-pointer-o"))
                                                                         ),
                                                                       ),
                                                                       br(),
                                                                       plotlyOutput("PCAplot_xenium", width = "100%"),
                                                                       
                                                              ),
                                                             
                                                              
                                                              tabPanel(title="Elbow", 
                                                                       br(),
                                                                       br(),
                                                                       plotOutput("Elbow_xenium", width = "100%")
                                                              ))
                                                ),
                                              
                                              
                                     ),
                                     tabPanel("Clustering",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                              column(3,
                                                     selectInput("scAnalysis_sp1",
                                                     label = "Analysis method",
                                                     choices = c("Seurat", "GraphST"),
                                                     selected = "Seurat"),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'Seurat'",
                                                fluidRow(
                                                  column(3,
                                                         numericInput("clus.res_spatial",
                                                                      label = "Resolution used",
                                                                      value = 0.6,
                                                                      min = 0.1,
                                                                      step = 0.1)
                                                  ),
                                                  column(3,
                                                         numericInput("dim.used_spatial",
                                                                     label = "PC to use",
                                                                     value = 10,
                                                                     min = 1,
                                                                     step = 1),
                                                  ),
                                                  column(2,
                                                         br(),
                                                         actionButton("findCluster_spatial", "Find Clusters", icon = icon("hand-pointer-o")),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_spatial", width = "50%")
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'GraphST'",
                                                column(12,
                                                column(4,
                                                fileInput("tpmFile_graphst",
                                                          label = "Spaceranger output (Accepted Format: .h5)",
                                                          accept = ".h5"),
                                                ),
                                                column(3,
                                                shinyFiles::shinyDirButton(id = 'dir_graphst', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                                                verbatimTextOutput("dir_graphst", placeholder = TRUE),
                                                ),
                                                column(2,
                                                       numericInput("cluster_graphst",
                                                                    label = "Number of clusters",
                                                                    value = 15,
                                                                    min = 2,
                                                                    step = 1)
                                                  ),
                                                column(2,
                                                selectInput("cluster_method_graphst",
                                                            label = "Clustering method",
                                                            choices = c("mclust", "leiden", "louvain"),
                                                            selected = "mclust")
                                                ),
                                                column(12,
                                                       actionButton("findCluster_graphst", "Find Clusters", icon = icon("hand-pointer-o")),
                                                        ),
                                                ),
                                                br(),
                                                column(12,
                                                plotOutput("Cluster2DPlot_graphst", width = "100%")
                                                )
                                              ),
                                            ),
                                            
                                            conditionalPanel(
                                              condition = "input.scAnalysis_platform == 'Xenium'",
                                              fluidRow(
                                                column(3,
                                                       numericInput("clus.res_xenium",
                                                                    label = "Resolution used",
                                                                    value = 0.3,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       selectInput("dim.used_xenium",
                                                                   label = "Dimensions",
                                                                   choices = c(10:50), 
                                                                   selected = 30),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("findCluster_xenium", "Find Clusters", icon = icon("hand-pointer-o")),
                                                       br()
                                                )),
                                              br(),
                                              plotlyOutput("Cluster2DPlot_xenium", width = "100%")
                                            ),
                                     ),
                                     
                                     tabPanel("UMAP",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'Seurat' & input.scAnalysis_platform == 'Visium'",
                                                column(3,
                                                       numericInput("dim.used_spatial",
                                                                    label = "Dimensions used",
                                                                    value = 10)
                                                ),
                                                column(3,
                                                       numericInput("clus.res_spatial",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(9,
                                                       br(),
                                                       actionButton("runUMAP_spatial", "Run UMAP", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotlyOutput("DimPlot_spatial", width = "100%"),
                                                       br(),
                                                       plotOutput("SpatialDimPlot", width = "100%"),
                                                       br(),
                                                       textOutput("silhouette_seurat")
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'GraphST' & input.scAnalysis_platform == 'Visium'",
                                                plotOutput("SpatialDimPlot_GraphST", width = "100%"),
                                                br(),
                                                textOutput("silhouette_graphst")
                                                ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                
                                                column(3,
                                                       selectInput("dim.used_xenium",
                                                                   label = "Dimensions",
                                                                   choices = c(10:50), 
                                                                   selected = 30),
                                                      ),
                                               
                                                column(3,
                                                       br(),
                                                       actionButton("runUMAP_xenium", "Run UMAP", icon = icon("hand-pointer-o")),
                                                       ),
                                                br(),
                                                br(),
                                                column(12,
                                                       br(),
                                                       plotlyOutput("umap_xenium", width = "50%"),
                                                      ),
                                                br(),
                                                column(6,
                                                       uiOutput("xenium.gene3.select"),
                                                       plotlyOutput("xenium_feature.plot", width = "100%"),
                                                ),
                                                plotOutput("spatialumap_xenium", width = "50%"),
                                                ),
                                     ),
                                     tabPanel("Visualize Spatial Domains",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'Seurat' & input.scAnalysis_platform == 'Visium'",
                                                column(4,
                                                       actionButton("Vis_sp1", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                
                                                fluidRow(
                                                  column(6,
                                                         uiOutput("sp.cluster.select"),
                                                         plotOutput("sp1.plot", width = "100%"),
                                                  ),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'GraphST' & input.scAnalysis_platform == 'Visium'",
                                                column(4,
                                                       actionButton("Vis_sp2", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                
                                                fluidRow(
                                                  column(6,
                                                         uiOutput("sp.cluster1.select"),
                                                         plotOutput("sp2.plot", width = "100%"),
                                                  ),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                column(4,
                                                       actionButton("Vis_sp_xenium", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                
                                                fluidRow(
                                                  column(6,
                                                         uiOutput("xenium.cluster.select"),
                                                         plotOutput("xenium1.plot", width = "100%"),
                                                  ),
                                                )),
                                     ),
                                     tabPanel("Deconvolution",
                                              column(2,
                                                     selectInput("scAnalysis_sp2",
                                                                 label = "Analysis method",
                                                                 choices = c("Seurat", "GraphST"),
                                                                 selected = "Seurat"),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp2 == 'Seurat'",
                                                column(2,
                                                       fileInput(inputId = 'tpmFiles_scRNA',
                                                                 label = "scRNA Reference",
                                                                 multiple = FALSE,
                                                                 accept = ".rds")),
                                                
                                                column(2,
                                                       selectInput("dim.used_sc",
                                                                   label = "Dimensions",
                                                                   choices = c(1:50),
                                                                   selected = "30"),
                                                ),
                                                #column(2,
                                                #       numericInput("num.genes_sc",
                                                #                    label = "Number of genes",
                                                #                    value = 2000,
                                                #                    min = 500,
                                                #                    step = 500)
                                                #),
                                                #br(),
                                                column(3,
                                                       br(),
                                                       actionButton("loadexample_ref", "Load example reference dataset", icon = icon("hand-o-right"),
                                                       ),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("loadButton_process_sp", "Load scRNA-seq dataset", icon = icon("hand-o-right"),
                                                       ),
                                                ),
                                                
                                                column(12,
                                                       br(),
                                                       actionButton("process_scRNA", "Process reference dataset", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotOutput("scRNAPlot", width = "100%")),
                                                column(9,
                                                       br(),
                                                       actionButton("vis_spRNA", "Visualise spatial dataset", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotOutput("spRNAPlot", width = "100%")),
                                                column(9,
                                                       br(),
                                                       actionButton("doDeconv", "Deconvolute", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       uiOutput("ct.select"),
                                                       plotOutput("DeconvPlot", width = "100%")
                                                ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp2 == 'GraphST'",
                                                
                                                column(3,
                                                       fileInput(inputId = 'tpmFiles_scRNA_graphst1',
                                                                 label = "scRNA reference dataset",
                                                                 multiple = FALSE,
                                                                 accept = ".rds")),
                                                #column(2,
                                                #       selectInput("dim.used_sc_graphst1",
                                                #                   label = "Dimensions",
                                                #                   choices = c(1:50),
                                                #                   selected = "30"),
                                                #      ),
                                                #column(2,
                                                #       numericInput("num.genes_sc_graphst1",
                                                #                    label = "Number of genes",
                                                #                    value = 2000,
                                                #                    min = 500,
                                                #                    step = 500)
                                                #       ),
                                                #br(),
                                                column(3,
                                                       br(),
                                                       actionButton("loadexample_ref1", "Load example reference dataset", icon = icon("hand-o-right"),
                                                       ),
                                                ),
                                                column(3,
                                                       actionButton("loadButton_process_sp_graphst1", "Load scRNA-seq dataset", icon = icon("hand-o-right"),
                                                       ),
                                                       ),
                                                column(12,
                                                       br(),
                                                       actionButton("process_scRNA_graphst1", "Process reference dataset", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotOutput("scRNAPlot_graphst1", width = "100%")),
                                                column(4,
                                                       fileInput("tpmFile_graphst1",
                                                                 label = "Spaceranger output (Accepted Format: .h5)",
                                                                 accept = ".h5"),
                                                ),
                                                column(3,
                                                       shinyFiles::shinyDirButton(id = 'dir_graphst1', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                                                       verbatimTextOutput("dir_graphst1", placeholder = TRUE),
                                                ),
                                                column(9,
                                                       br(),
                                                       actionButton("vis_spRNA_graphst", "Visualise spatial dataset", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotOutput("spRNAPlot_graphst", width = "100%")),
                                                column(9,
                                                       br(),
                                                       actionButton("doDeconv_graphst", "Deconvolute", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       uiOutput("graphst_ct.select"),
                                                       plotOutput("DeconvPlot_graphst", width = "100%")
                                                       ),
                                                      ),
                                                ),
                                     ),
               ),
               
               tabPanel("Single cell multiomics",
                        navlistPanel(widths=c(2,10),
                                     tabPanel("Overview",
                                              h2(p("Workflow for scMultiomics module")),
                                              br(),
                                              imageOutput("multiomics_image"),         
                                     ),
                                     tabPanel("Upload your data",
                                              column(9,
                                                     column(5,
                                                            #h4('Load Data:'),
                                                            wellPanel(
                                                              titlePanel(h4(p("Load your input data"))),
                                                              br(),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'CITE-seq'",
                                                                fileInput("tpmFiles2",
                                                                          label = "Cellranger output (Accepted Format: .h5)",
                                                                          accept = ".h5",
                                                                          multiple = T)),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                                fileInput(inputId = 'tpmFiles3',
                                                                          label = "Cellranger output (Accepted Format: .h5)",
                                                                          accept = ".h5")),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                                shinyFiles::shinyFilesButton(id = 'dir_multi_atac', label = "Path to fragments file", title = "Sheets Folder Selector", multiple = T),
                                                                verbatimTextOutput("dir_multi_atac", placeholder = TRUE)
                                                              ),
                                                              
                                                              selectInput("scAnalysis_mult",
                                                                          label = "Analysis method",
                                                                          choices = c("Seurat", "MOFA2"),
                                                                          selected = "Seurat"),
                                                              selectInput("scAnalysis_type",
                                                                          label = "scAnalysis type",
                                                                          choices = c("CITE-seq", "Multiome"),
                                                                          selected = "H5"),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'CITE-seq'",
                                                                actionButton("loadexample_cite_seurat", "Load example and run", icon = icon("hand-o-right")),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                                actionButton("loadexample_multiome_seurat", "Load example and run", icon = icon("hand-o-right")),
                                                              ),
                                                              textInput(inputId = "projName2",
                                                                        label = "Project Name",
                                                                        value = "Multiomics"),
                                                              fluidRow(
                                                                actionButton("loadButton2", "Load data", icon = icon("hand-o-right")),
                                                                actionButton("reset_mult", "Reset", icon = icon("repeat"))
                                                              ),
                                                            )),
                                                     
                                                     column(3,
                                                            numericInput(inputId = "min.genes2",
                                                                         label = "Min. genes",
                                                                         value = 200,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            numericInput(inputId = "min.cells2",
                                                                         label = "Min. cells",
                                                                         value = 3,
                                                                         min = 1)
                                                     ),
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'CITE-seq'",
                                                       column(3,
                                                              actionButton("create_seurat2a", "Process", icon = icon("hand-o-right"))
                                                       )),
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'Multiome'",
                                                       column(3,
                                                              actionButton("create_seurat2b", "Process", icon = icon("hand-o-right"))
                                                       )),
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'CITE-seq'",   
                                                column(12,
                                                       h4(p("RNA Data")),
                                                       withSpinner(dataTableOutput('countdataDT2a')),
                                                       h4(p("ADT/ATAC Data")),
                                                       withSpinner(dataTableOutput('countdataDT2b'))
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                column(12,
                                                       h4(p("RNA Data")),
                                                       withSpinner(dataTableOutput('countdataDT2c')),
                                                       h4(p("ADT/ATAC Data")),
                                                       withSpinner(dataTableOutput('countdataDT2d'))
                                                )),
                                     ),
                                     
                                     
                                     tabPanel("Violin Plot",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'CITE-seq'",
                                                column(3,
                                                       plotOutput("nFeature_RNAPlot2a")
                                                ),
                                                column(3,
                                                       plotOutput("mitoPlot2a")
                                                ),
                                                column(3,
                                                       plotOutput("nCount_RNAPlot2a")
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                column(3,
                                                       plotOutput("nFeature_RNAPlot2b")
                                                ),
                                                column(3,
                                                       plotOutput("mitoPlot2b")
                                                ),
                                                column(3,
                                                       plotOutput("nCount_RNAPlot2b")
                                                )),
                                     ),
                                     
                                     tabPanel("Feature Scatter Plot",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'CITE-seq'",
                                                column(5,
                                                       plotlyOutput("FeatureScatterPlot_mult_a")
                                                ),
                                                column(5,
                                                       plotlyOutput("FeatureScatterPlot_mult_b")
                                                )),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                column(5,
                                                       plotlyOutput("FeatureScatterPlot_mult_c")
                                                ),
                                                column(5,
                                                       plotlyOutput("FeatureScatterPlot_mult_d")
                                                )
                                              ),
                                     ),
                                     
                                     tabPanel("Normalization and Variable Gene Plot",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                #textOutput("nVarGenes_mult"),
                                                fluidRow(
                                                  column(3,
                                                         numericInput("var.genes_mult",
                                                                      label = "Number of variable genes",
                                                                      value = 2000,
                                                                      min = 500,
                                                                      step = 500)
                                                  ),
                                                  column(3,
                                                         selectInput("selection.method_mult",
                                                                     label = "Selection method",
                                                                     choices = c("vst", "dispersion"))
                                                  ),
                                                  column(4,
                                                         br(),
                                                         actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                  )),
                                                plotOutput("VarGenes_mult", width = "100%")),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                                #textOutput("nVarGenes_mult"),
                                                fluidRow(
                                                  column(3,
                                                         numericInput("var.genes_mult1",
                                                                      label = "Number of variable genes",
                                                                      value = 2000,
                                                                      min = 500,
                                                                      step = 500)
                                                  ),
                                                  column(3,
                                                         selectInput("selection.method_mult1",
                                                                     label = "Selection method",
                                                                     choices = c("vst", "dispersion"))
                                                  ),
                                                  column(4,
                                                         br(),
                                                         actionButton("findVarGenes_mult1", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                  )),
                                                plotOutput("VarGenes_mult1aa", width = "100%")),
                                              
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                
                                                fluidRow(
                                                  
                                                  column(3,
                                                         numericInput("var.genes_mult",
                                                                      label = "Number of variable genes",
                                                                      value = 2000,
                                                                      min = 500,
                                                                      step = 500)
                                                  ),
                                                  column(4,
                                                         br(),
                                                         actionButton("doSCTransform_multi", "Run scTransform", icon = icon("hand-pointer-o")),
                                                         #actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                         
                                                  )),
                                                plotOutput("VarGenes_mult1", width = "100%")),
                                     ),
                                     
                                     tabPanel("PCA",
                                              br(),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                tabsetPanel(id="Pca_mult_seurat",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              actionButton("runPCA_mult_seurat", "Run PCA", icon = icon("hand-pointer-o"))
                                                                       ),
                                                                     ),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_mult_seurat_h5_1", width = "100%"),
                                                                     
                                                            ),
                                                            tabPanel(title="PC Gene Visualisation",
                                                                     br(),
                                                                     selectInput("select.pc_mult_seurat",
                                                                                 label = "PC to plot",
                                                                                 choices = c(1:50)
                                                                     ),
                                                                     fluidRow(
                                                                       column(4,
                                                                              plotOutput("vizPlot_mult_seurat", width = "100%", height = "600px")
                                                                       ),
                                                                       column(8,
                                                                              plotOutput("PCHeatmap_mult_seurat", width = "100%", height = "600px")
                                                                       )
                                                                     ),
                                                                     DT::dataTableOutput("PCtable_mult_seurat")
                                                            ),
                                                            
                                                            tabPanel(title="Elbow", 
                                                                     br(),
                                                                     br(),
                                                                     plotOutput("Elbow_mult_seurat", width = "100%")
                                                            ))),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'Multiome'",
                                                tabsetPanel(id="Pca_mult_seurat",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              actionButton("runPCA_mult_seurat1", "Run PCA", icon = icon("hand-pointer-o"))
                                                                       ),
                                                                     ),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_mult_seurat_h5_2", width = "100%"),
                                                                     
                                                            ),
                                                            tabPanel(title="PC Gene Visualisation",
                                                                     br(),
                                                                     selectInput("select.pc_mult_seurat1",
                                                                                 label = "PC to plot",
                                                                                 choices = c(1:20)
                                                                     ),
                                                                     fluidRow(
                                                                       column(4,
                                                                              plotOutput("vizPlot_mult_seurat1", width = "100%", height = "600px")
                                                                       ),
                                                                       column(8,
                                                                              plotOutput("PCHeatmap_mult_seurat1", width = "100%", height = "600px")
                                                                       )
                                                                     ),
                                                                     DT::dataTableOutput("PCtable_mult_seurat1")
                                                            ),
                                                            
                                                            tabPanel(title="Elbow", 
                                                                     br(),
                                                                     br(),
                                                                     plotOutput("Elbow_mult_seurat1", width = "100%")
                                                            ))),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                                tabsetPanel(id="Pca_mult_mofa2",
                                                            tabPanel(title="Data overview",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              selectInput("num_factor_mofa",
                                                                                          label = "Number of factors",
                                                                                          choices = c(1:15)
                                                                              )),
                                                                       column(3,
                                                                              br(),
                                                                              actionButton("runMOFA_mult", "Run MOFA2", icon = icon("hand-pointer-o"))
                                                                       ),
                                                                     ),
                                                                     br(),
                                                                     plotOutput("mofaplot_mult_1", width = "50%"),
                                                            ),
                                                            tabPanel(title="Correlation between factors",
                                                                     br(),
                                                                     plotOutput("mofaplot_mult_2", width = "50%"),
                                                            ),
                                                            tabPanel(title="Plot variance decomposition",
                                                                     br(),
                                                                     plotOutput("mofaplot_mult_3", width = "50%"),
                                                            ),
                                                            tabPanel(title="Inspection of combinations of Factors",
                                                                     br(),
                                                                     column(3,
                                                                            selectInput("factor1",
                                                                                        label = "1st Factor to plot",
                                                                                        choices = c(1:15),
                                                                                        selected = "1"
                                                                            )),
                                                                     column(3,
                                                                            selectInput("factor2",
                                                                                        label = "2nd Factor to plot",
                                                                                        choices = c(1:15),
                                                                                        selected = "2"
                                                                            )),
                                                                     plotOutput("mofaplot_mult_4", width = "50%"),
                                                            ),
                                                            tabPanel(title="Plot molecular signatures",
                                                                     br(),
                                                                     column(3,
                                                                            selectInput("selection.view",
                                                                                        label = "Select view",
                                                                                        choices =c("RNA", "ADT"),
                                                                                        selected = "RNA")
                                                                     ),
                                                                     column(3,
                                                                            selectInput("factor1a",
                                                                                        label = "Factor to plot",
                                                                                        choices = c(1:10),
                                                                                        selected = "RNA"
                                                                            )),
                                                                     column(3,
                                                                            selectInput("num.features.mofa",
                                                                                        label = "Number of features to plot",
                                                                                        choices = c(1:100),
                                                                                        selected = "25"
                                                                            )),
                                                                     br(),
                                                                     br(),
                                                                     column(12,
                                                                            plotOutput("mofaplot_mult_5", width = "50%"),
                                                                     ),
                                                            ),
                                                )),
                                     ),
                                     
                                     tabPanel("Clustering",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                br(),
                                                fluidRow(
                                                  column(3,
                                                         numericInput("clus.res_mult_seurat",
                                                                      label = "Resolution used",
                                                                      value = 0.6,
                                                                      min = 0.1,
                                                                      step = 0.1)
                                                  ),
                                                  column(3,
                                                         selectInput("dim.used_mult_seurat",
                                                                     label = "PC to use",
                                                                     choices = c(10:50)),
                                                  ),
                                                  column(3,
                                                         br(),
                                                         actionButton("findCluster_mult_seurat", "Find Clusters", icon = icon("hand-pointer-o")),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_mult_seurat", width = "100%")
                                              ),  
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'Multiome'",
                                                br(),
                                                fluidRow(
                                                  column(3,
                                                         numericInput("clus.res_mult_seurat1",
                                                                      label = "Resolution used",
                                                                      value = 0.6,
                                                                      min = 0.1,
                                                                      step = 0.1)
                                                  ),
                                                  column(3,
                                                         selectInput("dim.used_mult_seurat1",
                                                                     label = "PC to use",
                                                                     choices = c(10:50)),
                                                  ),
                                                  column(3,
                                                         br(),
                                                         actionButton("findCluster_mult_seurat1", "Find Clusters", icon = icon("hand-pointer-o")),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_mult_seurat1", width = "100%")
                                              ),  
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                                br(),
                                                fluidRow(
                                                  column(3,
                                                         numericInput("num.clus.mofa",
                                                                      label = "Number of clusters",
                                                                      value = 10,
                                                                      min = 1,
                                                                      step = 1)
                                                  ),
                                                  column(3,
                                                         br(),
                                                         actionButton("findCluster_mofa2", "Find Clusters", icon = icon("hand-pointer-o")),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotOutput("Cluster2DPlot_mofa", width = "50%")
                                              ),  
                                     ),
                                     
                                     tabPanel("UMAP",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                br(),
                                                column(3,
                                                       numericInput("dim.used_mult_seurat",
                                                                    label = "Dimensions used",
                                                                    value = 10)
                                                ),
                                                column(3,
                                                       numericInput("clus.res_mult_seurat2a",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("runUMAP_mult_seurat", "Run UMAP", icon = icon("hand-pointer-o")),
                                                ),
                                                column(9,
                                                       br(),
                                                       br(),
                                                       plotOutput("UMAPplot_mult_seurat_1", width = "100%"),
                                                       br(),
                                                       plotOutput("UMAPplot_mult_seurat_2", width = "100%"),
                                                       br(),
                                                       plotOutput("UMAPplot_mult_seurat_3", width = "100%"),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'Multiome'",
                                                br(),
                                                column(3,
                                                       numericInput("dim.used_mult_seurat1",
                                                                    label = "Dimensions used",
                                                                    value = 10)
                                                ),
                                                column(3,
                                                       numericInput("clus.res_mult_seurat2b",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(9,
                                                       br(),
                                                       actionButton("runUMAP_mult_seurat1", "Run UMAP", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotlyOutput("UMAPplot1_a", width = "100%"),
                                                       br(),
                                                       plotlyOutput("UMAPplot1_b", width = "100%"),
                                                       br(),
                                                       plotlyOutput("UMAPplot1_c", width = "100%"),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                                br(),
                                                column(3,
                                                       numericInput("num.neighbors_mofa2",
                                                                    label = "Number of neighbors",
                                                                    value = 15)
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("runUMAP_mofa2", "Run UMAP", icon = icon("hand-pointer-o")),
                                                ),
                                                column(9,
                                                       br(),
                                                       br(),
                                                       plotOutput("UMAPplot_mofa_1", width = "100%"),
                                                )),
                                     ),
                                     
                                     # tabPanel("tSNE",
                                     #          conditionalPanel(
                                     #            condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                     #            br(),
                                     #            column(3,
                                     #                   numericInput("dim.used_mult_seurat",
                                     #                                label = "Dimensions used",
                                     #                                value = 10)
                                     #            ),
                                     #            column(9,
                                     #                   br(),
                                     #                   actionButton("runTSNE_mult_seurat", "Run tSNE", icon = icon("hand-pointer-o")),
                                     #                   br(),
                                     #                   br(),
                                     #                   plotlyOutput("TSNEplot_mult_seurat_1", width = "100%"),
                                     #                   br(),
                                     #                   plotlyOutput("TSNEplot_mult_seurat_2", width = "100%"),
                                     #                   br(),
                                     #                   plotlyOutput("TSNEplot_mult_seurat_3", width = "100%"),
                                     #            )),
                                     #                    conditionalPanel(
                                     #                     condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'Multiome'",
                                     #                      br(),
                                     #                     column(3,
                                     #                            numericInput("dim.used_mult_seurat1",
                                     #                                         label = "Dimensions used",
                                     #                                         value = 10)
                                     #                     ),
                                     #                      column(9,
                                     #                            br(),
                                     #                            actionButton("runTSNE_mult_seurat1", "Run tSNE", icon = icon("hand-pointer-o")),
                                     #                            br(),
                                     #                            br(),
                                     #                            plotlyOutput("TSNEplot1_a", width = "100%"),
                                     #                            br(),
                                     #                            plotlyOutput("TSNEplot1_b", width = "100%"),
                                     #                            br(),
                                     #                            plotlyOutput("TSNEplot1_c", width = "100%"),
                                     #                           
                                     #                     ),
                                     #             ),
                                     #          conditionalPanel(
                                     #             condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                     #             br(),
                                     #              column(3,
                                     #                    br(),
                                     #                    actionButton("runTSNE_mofa2", "Run UMAP", icon = icon("hand-pointer-o")),
                                     #                     ),
                                     #                column(9,
                                     #                   br(),
                                     #                   br(),
                                     #                   plotlyOutput("TSNEplot_mofa_1", width = "100%"),
                                     #                      )),
                                     #                ),
                                     
                                     tabPanel("Cell type identification",
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                column(3,
                                                       selectInput("cellatlas_mult_cite_seurat",
                                                                   label = "Reference Atlas",
                                                                   choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                   selected = "all")
                                                ),
                                                column(3,
                                                       selectInput("assay_mult_cite_seurat",
                                                                   label = "Assay",
                                                                   choices = c("rna.umap", "adt.umap", "wnn.umap"),
                                                                   selected = "all")
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("doCELLiD_mult_cite_seurat", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                       #textOutput("CELLiD.done"),
                                                       br()
                                                ),
                                                br(),
                                                br(),
                                                br(),
                                                plotOutput("Umap_cellid_mult_cite_seurat", width = "50%"),
                                                plotOutput("Umap_cellid_mult_cite_seurat1", width = "50%"),
                                                br(),
                                                DT::dataTableOutput("ct_cite_seurat.table"),
                                                column(4,
                                                       downloadButton('download_cellid_cite_seurat_prediction', 'Download CELLiD predictions (in csv)'),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'Multiome'",
                                                column(3,
                                                       selectInput("cellatlas_mult_multiome_seurat",
                                                                   label = "Reference Atlas",
                                                                   choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                   selected = "all")
                                                ),
                                                column(3,
                                                       selectInput("assay_mult_multiome_seurat",
                                                                   label = "Assay",
                                                                   choices = c("umap.rna", "umap.atac", "wnn.umap"),
                                                                   selected = "all")
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("doCELLiD_multiome_seurat", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                       #textOutput("CELLiD.done"),
                                                       br()
                                                ),
                                                br(),
                                                br(),
                                                br(),
                                                plotlyOutput("Umap_cellid_multiome_seurat", width = "50%"),
                                                plotlyOutput("Umap_cellid_multiome_seurat1", width = "50%"),
                                                br(),
                                                DT::dataTableOutput("ct_multiome_seurat.table"),
                                                column(4,
                                                       downloadButton('download_cellid_multiome_seurat_prediction', 'Download CELLiD predictions (in csv)'),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                                column(3,
                                                       selectInput("cellatlas_mult_cite_mofa",
                                                                   label = "Reference Atlas",
                                                                   choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                   selected = "all")
                                                ),
                                                column(3,
                                                       br(),
                                                       actionButton("doCELLiD_mult_cite_mofa", "Run CELLiD", icon = icon("hand-pointer-o")),
                                                       #textOutput("CELLiD.done"),
                                                       br()
                                                ),
                                                br(),
                                                br(),
                                                br(),
                                                plotOutput("Umap_cellid_mult_cite_mofa", width = "100%"),
                                                br(),
                                                plotOutput("Umap_cellid_mult_cite_mofa1", width = "100%"),
                                                br(),
                                                DT::dataTableOutput("ct_cite_mofa.table"),
                                                column(4,
                                                       downloadButton('download_cellid_cite_mofa_prediction', 'Download CELLiD predictions (in csv)'),
                                                )),
                                     ),
                                     
                                     tabPanel("Data Visualization",
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                column(3,
                                                       br(),
                                                       actionButton("Vis3", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                column(3,
                                                       selectInput("assay_mult_cite_seurat1",
                                                                   label = "Reduction",
                                                                   choices = c("rna.umap", "adt.umap", "wnn.umap"),
                                                                   selected = "rna.umap")
                                                ),
                                                column(3,
                                                       selectInput("assay_mult_cite_seurat2",
                                                                   label = "Assay",
                                                                   choices = c("RNA", "ADT"),
                                                                   selected = "RNA")
                                                ),
                                                #fluidRow(
                                                column(3,
                                                       uiOutput("vis.gene.select"),
                                                ),
                                                fluidRow(
                                                  column(6,
                                                         plotlyOutput("vis.plot", width = "100%"),
                                                         #br(),
                                                         #uiOutput("vis.gene.select1"),
                                                         #plotlyOutput("vis1.plot", width = "100%"),
                                                  )),
                                              ),  
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA2' & input.scAnalysis_type == 'CITE-seq'",
                                                column(4,
                                                       actionButton("Vis3a", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                
                                                fluidRow(
                                                  column(6,
                                                         uiOutput("vis.gene.select_a"),
                                                         plotlyOutput("vis.plot_a", width = "100%"),
                                                         br(),
                                                         uiOutput("vis.gene.select1_a"),
                                                         plotlyOutput("vis1.plot_a", width = "100%"),
                                                  ),
                                                )
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'Multiome'",
                                                column(3,
                                                       br(),
                                                       actionButton("Vis3_multiome", "Visualize", icon = icon("hand-pointer-o"))
                                                ),
                                                column(3,
                                                       selectInput("assay_multiome_seurat1a",
                                                                   label = "Reduction",
                                                                   choices = c("umap.rna", "umap.atac", "wnn.umap"),
                                                                   selected = "rna.umap")
                                                ),
                                                column(3,
                                                       selectInput("assay_multiome_seurat2a",
                                                                   label = "Assay",
                                                                   choices = c("RNA", "ATAC"),
                                                                   selected = "RNA")
                                                ),
                                                
                                                
                                                  column(3,
                                                         uiOutput("vis.gene.select_multiome"),
                                                        ),
                                                fluidRow(
                                                    column(12,
                                                    plotlyOutput("vis.plot_multiome", width = "50%"),
                                                          ),
                                                         #br(),
                                                         #uiOutput("vis.gene.select_multiome1"),
                                                         #plotlyOutput("vis1.plot_multiome", width = "100%"),
                                                  ),
                                                )))
               ),
               tabPanel("Single cell ATAC-seq",
                        
                        navlistPanel(widths=c(2,10),
                                     tabPanel("Overview",
                                              h2(p("Workflow for scATAC-Seq module")),
                                              br(),
                                              imageOutput("atac_image"),         
                                     ),
                                     tabPanel("Upload your data",
                                              column(9,
                                                     column(5,
                                                            #h4('Load Data:'),
                                                            wellPanel(
                                                              titlePanel(h4(p("Load your input data"))),
                                                              br(),
                                                              
                                                              selectInput("scInput_atac",
                                                                          label = "Select Data Input Type",
                                                                          choices = c("Cellranger atac output"),
                                                                          selected = "Cellranger atac output"),
                                                              
                                                              fileInput("tpmFiles_atac",
                                                                        label = "Upload Peak/Cell matrix",
                                                                        accept = ".h5"),
                                                              
                                                              fileInput("meta_atac",
                                                                        label = "Upload Metadata",
                                                                        accept = ".csv"),
                                                              
                                                              shinyFiles::shinyFilesButton(id = 'dir_atac', label = "Path to fragments file", title = "Sheets Folder Selector", multiple = T),
                                                              verbatimTextOutput("dir_atac", placeholder = TRUE),
                                                              
                                                              selectInput("scAnalysis_atac",
                                                                          label = "Analysis method",
                                                                          choices = c("Signac"),
                                                                          selected = "Signac"),
                                                              actionButton("loadexample_atac", "Load example and run", icon = icon("hand-o-right")),
                                                              textInput(inputId = "projName4",
                                                                        label = "Project Name",
                                                                        value = "ATAC"),
                                                              fluidRow(
                                                                actionButton("loadButton_atac", "Load Data", icon = icon("hand-o-right")),
                                                                actionButton("reset_atac", "Reset Data", icon = icon("repeat")),
                                                              ))),
                                                     column(3,
                                                            numericInput(inputId = "min.genes_atac",
                                                                         label = "Min. genes",
                                                                         value = 200,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            numericInput(inputId = "min.cells_atac",
                                                                         label = "Min. cells",
                                                                         value = 3,
                                                                         min = 1)
                                                     ),
                                                     column(3,
                                                            actionButton("create_seurat_atac", "Process", icon = icon("hand-o-right"))
                                                     )
                                              ),
                                              
                                              column(12,
                                                     h4(p("ATAC counts")),
                                                     withSpinner(dataTableOutput('countdataDT_atac')),
                                              ),
                                     ),
                                     
                                     tabPanel(title="ATAC-seq QC Plot", value = "ATAC-seq_QC_panel1",
                                              br(),
                                              fluidRow(
                                                
                                                column(5,
                                                       plotlyOutput("TSSPlot", width = "100%"),
                                                ),
                                                column(5,
                                                       plotlyOutput("FragmentHistogram", width = "100%"),
                                                ),
                                                column(4,
                                                       plotOutput("Vlnplot_atac_1")
                                                ),
                                                column(4,
                                                       plotOutput("Vlnplot_atac_2")
                                                ),
                                                column(4,
                                                       plotOutput("Vlnplot_atac_3")
                                                ),
                                                column(4,
                                                       plotOutput("Vlnplot_atac_4")
                                                ),
                                                column(4,
                                                       plotOutput("Vlnplot_atac_5")
                                                ),
                                              )),
                                     
                                     tabPanel("Normalization", fluidPage(
                                       hr(),
                                       tabPanel(title="Normalization", value = "norm_ATAC",
                                                br(),
                                                fluidRow(
                                                  column(3,
                                                         actionButton("donorm_ATAC", "Run Normalization", icon = icon("hand-pointer-o")),
                                                         textOutput("normalize_atac.done"),
                                                         br()
                                                  )),
                                                br(),
                                                plotlyOutput("DepthCor_ATAC", width = "100%"),
                                       ))),
                                     
                                     ##------Clustering of ATAC-seq data ------------
                                     tabPanel("Clustering of ATAC-seq", fluidPage(
                                       hr(),
                                       fluidRow(
                                         column(3,
                                                numericInput("clus.res_atac",
                                                             label = "Resolution used",
                                                             value = 0.6,
                                                             min = 0.1,
                                                             step = 0.1)
                                                
                                         ),
                                         column(3,
                                                selectInput("dim.used_atac",
                                                            label = "PC to plot",
                                                            choices = c(2:50),
                                                            selected = 15 
                                                ),
                                         ),
                                         
                                         column(3,
                                                br(),
                                                actionButton("doCluster_ATAC", "Run Clustering", icon = icon("hand-pointer-o")),
                                                textOutput("cluster_atac.done"),
                                                
                                                
                                         )),
                                       br(),
                                       plotlyOutput("cluster_ATAC", width = "100%"),
                                     )),
                                     
                                     tabPanel("UMAP on ATAC-seq", fluidPage(
                                       hr(),
                                       fluidRow(
                                         column(3,
                                                numericInput("dim.used_atac",
                                                             label = "Dimensions used",
                                                             value = 30)
                                         ),
                                         br(),
                                         column(3,
                                                actionButton("doUMAP_ATAC", "Running UMAP", icon = icon("hand-pointer-o")),
                                                textOutput("umap_atac.done"),
                                                br()
                                         )),
                                       br(),
                                       plotlyOutput("Umap_ATAC", width = "100%"),
                                     )),
                                     
                                     tabPanel("TSNE on ATAC-seq", fluidPage(
                                       hr(),
                                       fluidRow(
                                         column(3,
                                                numericInput("dim.used_atac",
                                                             label = "Dimensions used",
                                                             value = 30)
                                         ),
                                         br(),
                                         column(3,
                                                actionButton("doTSNE_ATAC", "Running TSNE", icon = icon("hand-pointer-o")),
                                                textOutput("tsne_atac.done"),
                                                br()
                                         )),
                                       br(),
                                       plotlyOutput("Tsne_ATAC", width = "100%"),
                                     )),
                                     tabPanel("DE Peaks", fluidPage(
                                       fluidRow(
                                         
                                         column(3, selectInput("min_pct_atac",
                                                               label = "min.pct",
                                                               choices = c("0.1", "0.25"))
                                         ),
                                         
                                         column(3, selectInput("logfc_atac",
                                                               label = "logfc.threshold",
                                                               choices = c("0.1", "0.25"))
                                         ),
                                         
                                         column(3, selectInput("test.use_atac",
                                                               label = "Test use",
                                                               choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                         ),
                                         br(),
                                         column(3,
                                                actionButton("doDeg_atac", "Run DEGs", icon = icon("hand-pointer-o"))
                                         )),
                                       br(),
                                       column(9,
                                              DT::dataTableOutput("Deg_atac.table"),
                                              br(),
                                              plotlyOutput("Deg_atac.plot", width = "100%")
                                       )
                                     )),
                                     
                                     tabPanel("Data visualization",
                                              fluidRow(
                                                column(6,
                                                       uiOutput("deg.atac.select"),
                                                       plotlyOutput("Deg_atac1.plot", width = "100%"),
                                                       br(),
                                                       plotlyOutput("Deg_atac2.plot", width = "100%"),
                                                       br(),
                                                       plotOutput("Deg_atac3.plot", width = "100%")
                                                ),
                                              )
                                     ),
                                     tabPanel("Coverage Plot",
                                              uiOutput("coverage.atac.select"),
                                              plotOutput("coverage.plot", width = "100%")
                                     ),
                                     tabPanel("Motif Analysis",
                                              
                                              tabsetPanel(id="motif_analysis",
                                                          tabPanel(title="Motif analysis",
                                                                   br(),
                                                                   column(3, numericInput("min_pct_motif",
                                                                                          label = "min.pct",
                                                                                          value = 0.25,
                                                                                          min = 0,
                                                                                          step = 0.01)
                                                                   ),
                                                                   column(3, numericInput("logfc_motif",
                                                                                          label = "logfc.threshold",
                                                                                          value = 0.25,
                                                                                          min = 0,
                                                                                          step = 0.01)
                                                                   ),
                                                                   column(3, selectInput("test.use_motif",
                                                                                         label = "Test use",
                                                                                         choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                                   ),
                                                                   br(),
                                                                   column(3,
                                                                          actionButton("doDeg_motif", "Run DEGs", icon = icon("hand-pointer-o"))
                                                                   ),  
                                                                   br(),
                                                                   fluidRow(
                                                                     column(6,
                                                                            uiOutput("motif.select"),
                                                                            plotOutput("motif.plot", width = "100%")
                                                                     )
                                                                   ),
                                                          ),
                                                          tabPanel(title="Motif activities",
                                                                   br(),
                                                                   column(3,
                                                                          actionButton("calc_motif_activity", "Calculate Motif Activity", icon = icon("hand-pointer-o")),
                                                                   ),
                                                                   column(3,
                                                                          uiOutput("motif_feature.select"),
                                                                   ),
                                                                   br(),
                                                                   fluidRow(
                                                                     column(12,
                                                                            plotlyOutput("motif_feature.plot", width = "50%"),
                                                                     )),  
                                                                   br(),
                                                                   fluidRow(
                                                                     column(12,
                                                                            uiOutput("motif1.select"),
                                                                            plotOutput("motif1.plot", width = "50%")
                                                                     )),
                                                                   
                                                          ),
                                              )),
                                     
                                     
                        ),
               ),
               tabPanel("Help",
                        wellPanel(
                          HTML(
                            '
      <p align="center" width="4">Singapore Immunology Network, Agency for Science, Technology and Research (A*STAR)</p>
      <p align="center" width="4">Github: <a href="https://github.com/JinmiaoChenLab/">https://github.com/JinmiaoChenLab/</a></p>
      <p align="center" width="4">Created by: <a href="mailto:Chen_Jinmiao@immunol.a-star.edu.sg">Jinmiao Chen Lab</a> </p>'
                          ))
                        )
               )
)))
