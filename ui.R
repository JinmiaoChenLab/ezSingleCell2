library(shiny)
library(plotly)
library(shinythemes)
library(shinydashboard)
library(shinyalert)
library(shinyjs)
library(shinyBS)
library(CellChat)
library(stringr)
library(fontawesome)
#library(hdf5r)
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
library(rGREAT)
library(anndata)
library(clustree)
library(magrittr)
library(rhandsontable)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(lisi)
library(cowplot)
library(FastIntegration)
library(Signac)
library(cluster)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(purrr)
library(dplyr)
library(DT)
library(msigdbr)
library(fgsea)
library(MOFA2)
library(stats)
library(ktplots)
library(ggplot2)
library(SeuratWrappers)
library(shinyWidgets)
#library(car)
#library(nortest)
library(shinyFiles)
library(bslib)
#library(torch)

tags$style(type="text/css",
           ".shiny-output-error { visibility: hidden; }",
           ".shiny-output-error:before { visibility: hidden; }"
)

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
                                   p("The toolkit uses example data from", a("DISCO", href = "https://www.immunesinglecell.org/", target="_blank"),"that contains data from 13998 samples, covering 461 tissues/cell lines/organoids, 158 diseases, and 20 platforms.", style="text-align:justify;color:black;font-size:15px"), width=12)),
                    
                            
                   # p(em("Developed by"),br("CJM Lab"),style="text-align:center; font-family: times")
                        
               wellPanel(
                 HTML(
                   '<p align="center" width="4">Singapore Immunology Network, Agency for Science, Technology and Research (A*STAR)</p>
                           <p align="center" width="4">Github: <a href="https://github.com/JinmiaoChenLab/">https://github.com/JinmiaoChenLab/</a></p>
                           <p align="center" width="4">Created by <a href="https://www.a-star.edu.sg/sign/people/principal-investigators/jinmiao-chen">Jinmiao Chen Lab</a> and collaborated with <a href="https://www.vishuo.com/en/"> Vishuo Biomedical Pte. Ltd.</a></p>'
                 )),
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
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'Raw Counts Matrix'",
                                                                titlePanel(h4(p("Load your example data"))),
                                                                actionBttn("loadexample_tpm", "Load example and run", icon = icon("hand-o-right"), size = 'sm', onclick = "$(tab).removeClass('disabled')"),
                                                                #bsPopover("loadexample_tpm", "Load Example Data","Press to load example data (Raw counts matrix)", placement = "bottom", trigger = "hover", options = NULL)
                                                              ),       
                                                                
                                                              conditionalPanel(
                                                                condition = "input.scInput == '10X cellranger'",
                                                                titlePanel(h4(p("Load your example data"))),
                                                                actionBttn("loadexample_scH5", "Load example and run", icon = icon("hand-o-right"), size = 'sm', onclick = "$(tab).removeClass('disabled')"),
                                                                #bsPopover("loadexample_scH5", "Load Example Data","Press to load example data (cellRanger output)", placement = "bottom", trigger = "hover", options = NULL)
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'rds object'",
                                                                titlePanel(h4(p("Load your example data"))),
                                                                actionBttn("loadexample_rds", "Load example and run", icon = icon("hand-o-right"), size = 'sm', onclick = "$(tab).removeClass('disabled')"),
                                                                #bsPopover("loadexample_rds", "Load Example Data","Press to load example data (rds object)", placement = "bottom", trigger = "hover", options = NULL)
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'h5ad'",
                                                                titlePanel(h4(p("Load your example data"))),
                                                                actionBttn("loadexample_h5ad", "Load example and run", icon = icon("hand-o-right"), size = 'sm', onclick = "$(tab).removeClass('disabled')"),
                                                                #bsPopover("loadexample_h5ad", "Load Example Data","Press to load example data (h5ad object)", placement = "bottom", trigger = "hover", options = NULL)
                                                              ),  
                                                              titlePanel(h4(HTML("<b>Load your input data</b>"))),
                                                              br(),
                                                              selectInput("scInput",
                                                                          label = "Select Data Input Type",
                                                                          choices = c("Raw Counts Matrix", "10X cellranger", "rds object", "h5ad"),
                                                                          selected = "Raw Counts Matrix"),
                                                              #bsPopover("scInput", "Select Input Format","Users can select input format", placement = "bottom", trigger = "hover", options = NULL),
                                                              
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'Raw Counts Matrix'",
                                                                fileInput("tpmFiles",
                                                                          label = "Counts File (Accepted Format: text)",
                                                                          accept = ".txt"),
                                                               ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == '10X cellranger'",
                                                                fileInput("scH5",
                                                                          label = "Cellranger output (Accepted Format: .h5)",
                                                                          accept = ".h5"),
                                                                ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'rds object'",
                                                                fileInput("rds",
                                                                          label = "Seurat Object (Accepted Format: .rds)",
                                                                          accept = ".rds"),
                                                                ),
                                                              conditionalPanel(
                                                                condition = "input.scInput == 'h5ad'",
                                                                fileInput("h5ad",
                                                                          label = "AnnData Object (Accepted Format: .h5ad)",
                                                                          accept = ".h5ad"),
                                                              ),
                                                              #column(6,
                                                              #       numericInput(inputId = "min.genes",
                                                              #                   label = "Min. genes",
                                                              #                   value = 200,
                                                              #                   min = 1)
                                                              #       ),
                                                              #column(6,
                                                              #      numericInput(inputId = "min.cells",
                                                              #                   label = "Min. cells",
                                                              #                   value = 3,
                                                              #                   min = 1)
                                                              #       ),
                                                                textInput(inputId = "projName",
                                                                          label = "Project Name",
                                                                          value = "scRNA"),
                                                              
                                                              fluidRow(
                                                                actionBttn("loadButton", "Load data", icon = icon("hand-o-right"), size = 'sm', onclick = "$(tab).removeClass('disabled')"),
                                                                #bsPopover("loadButton", "Load Data","Press to load your data", placement = "bottom", trigger = "hover", options = NULL),
                                                                actionBttn("reset_scRNA", "Reset", icon = icon("repeat"), size = 'sm'),
                                                                #bsPopover("reset_scRNA", "Reload Data","Press to reanalyze your data", placement = "bottom", trigger = "hover", options = NULL)
                                                              ),
                                                            )),
                                                     chooseSliderSkin("Modern"),
                                                     titlePanel(h4(p("Quality control"))),
                                                                          
                                                                             column(6,
                                                                                    plotOutput("nFeature_RNAPlot", width = "200%")
                                                                                    ),
                                                                             column(4,
                                                                                    br(),
                                                                                    downloadBttn('download_nFeature_RNA', 'Download (as png)', size = 'sm')
                                                                             ),
                                                                            
                                                                
                                                                    column(12,
                                                                           column(3,
                                                                                  numericInput("ob1",
                                                                                               label = "Min nFeature:",
                                                                                               value = 200,
                                                                                               min = 0,
                                                                                               step = 1),
                                                                           ),
                                                                           column(3,
                                                                                  numericInput("ob2",
                                                                                               label = "Max nFeature:",
                                                                                               value = 2500,
                                                                                               min = 0,
                                                                                               step = 1),
                                                                           ),
                                                                           column(3,
                                                                                  numericInput("ob3",
                                                                                               label = "Mt%:",
                                                                                               value = 5,
                                                                                               min = 0,
                                                                                               step = 1),
                                                                           ),
                                                                           column(3,
                                                                                  br(),
                                                                                  actionBttn("filter_seurat", "Filter", icon = icon("hand-o-right"), size = 'sm'),
                                                                                  #bsPopover("filter_seurat", "Filter Data","Press to filter data based on parameters", placement = "bottom", trigger = "hover", options = NULL),
                                                                           ),
                                                                        ),
                                                                      ),
                                     
                                       column(12,
                                              withSpinner(dataTableOutput('countdataDT')),
                                       
                                     downloadBttn('downloadCount', 'Download Table'),
                                     )),
                                     
                                     #tabPanel("2. Quality control Plot",
                                     #        tabsetPanel(id="qc_scRNA",
                                     #          tabPanel("Violin Plot", 
                                     #            column(3,
                                     #                   plotOutput("nFeature_RNAPlot")
                                     #           ),
                                     #            column(3,
                                     #                  plotOutput("mitoPlot")
                                     #           ),
                                     #            column(3,
                                     #                  plotOutput("nCount_RNAPlot")
                                     #          ),
                                     #          column(12,
                                     #                  column(3,
                                     #                 downloadBttn('download_nFeature_RNA', 'Download nFeature (as png)', size = 'sm'),
                                     #                 ),
                                     #                 column(3,
                                     #                   downloadBttn('download_mito', 'Download mito (as png)', size = 'sm'),
                                     #                  ),
                                     #                 column(3,
                                     #                downloadBttn('download_nCount_RNA', 'Download nCount (as png)', size = 'sm'),   
                                     #                  ),
                                     #                )
                                    #            ),
                                                
                                    # tabPanel("Feature Scatter Plot",
                                    #         column(5,
                                    #                plotlyOutput("FeatureScatterPlot1")
                                    #         ),
                                    #         column(5,
                                    #                plotlyOutput("FeatureScatterPlot2")
                                    #         ),
                                    #         br(),
                                    #         br(),
                                    #         column(12,
                                    #                column(5,
                                    #                       downloadBttn('download_FeatureScatterPlot1', 'Download (as png)', size = 'sm'),
                                    #                ),
                                    #                column(5,
                                    #                       downloadBttn('download_FeatureScatterPlot2', 'Download (as png)', size = 'sm'),
                                    #                ),
                                    #             ),
                                    #            )
                                    #           ),
                                    #         ),
                                     tabPanel("Normalization and Variable Feature Selection", value = "test",
                                              
                                              tags$script(
                                               '
                                                  var tab = $(\'a[data-value="test"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                                  ),
                                              
                                              selectInput("norm1",
                                                          label = "Normalization method",
                                                          choices = c("LogNormalize", "SCTransform"), 
                                                          selected = "SCTransform"
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
                                                       actionBttn("findVarGenes", "Identify highly variable genes", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab1).removeClass('disabled')"),
                                                       #bsPopover("findVarGenes", "Feature Selection","Press to identify highly variable features", placement = "bottom", trigger = "hover", options = NULL),
                                                       #actionButton("doSCTransform", "Run SCTransform", icon = icon("hand-pointer-o"))
                                                       # actionButton("doVarplot", "Plot variable genes", icon = icon("hand-pointer-o"))
                                                )),
                                              plotOutput("VarGenes", width = "100%"),
                                            ),
                                     tabPanel( "PCA", value = "test1",
                                              
                                              tags$script(
                                                '
                                                  var tab1 = $(\'a[data-value="test1"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              tabsetPanel(id="Pca",
                                                          tabPanel(title="PCA Plot", value="P_panel1",
                                                                   br(),
                                                                   #column(3,
                                                                   #selectInput("assays1",
                                                                   #            label = "Normalization method:",
                                                                   #            choices = c("LogNormalization", "SCTransform"), 
                                                                   #            selected = "SCTransform"
                                                                   #  ),
                                                                   # ),
                                                                   br(),
                                                                   column(3,
                                                                          actionBttn("doPCA", "Run PCA", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab2).removeClass('disabled')"),
                                                                          #bsPopover("doPCA", "Run PCA","Press to run PCA Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                                   ),
                                                                   br(),
                                                                   br(),
                                                                   column(6,
                                                                          plotlyOutput("PCA2DPlot", width = "100%")
                                                                          ),
                                                                  
                                                                   column(12,
                                                                          column(3,
                                                                                 br(),
                                                                                 downloadBttn('download_PCA', 'Download PCA Plot (as png)', size = 'sm'),
                                                                          ),
                                                                          column(3,
                                                                                 br(),
                                                                                 downloadBttn('download_PCA_embedding', 'Download PCA Embedding (as csv)', size = 'sm'),
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
                                                                            br(),
                                                                            downloadBttn('download_vizPlot', 'Download vizPlot (as png)', size = 'sm'),
                                                                     ),
                                                                     column(3,
                                                                            br(),
                                                                            downloadBttn('download_PCHeatmap', 'Download PCHeatmap (as png)', size = 'sm'),
                                                                     ),
                                                                   #br(),
                                                                   DT::dataTableOutput("PCtable"),
                                                                   column(3,
                                                                          downloadBttn('download_PCTable', 'Download top genes (as csv)', size = 'sm'),
                                                                      ),
                                                                     ),
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
                                                                                 downloadBttn('download_Elbow', 'Download ElbowPlot (as png)', size = 'sm'),
                                                                          ),
                                                                   ),
                                                                   
                                                          )
                                              )),
                                    
                                    tabPanel( "UMAP", value = "test2",
                                             
                                             tags$script(
                                               '
                                                  var tab2 = $(\'a[data-value="test2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
                                             #titlePanel(h4(p("UMAP Analysis"))),
                                             br(),
                                             fluidRow(
                                               column(3,
                                                      numericInput("dim.used",
                                                                   label = "Dimensions used",
                                                                   value = 10)
                                               ),
                                               br(),
                                               column(3,
                                                      actionBttn("doUmap", "Run UMAP", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab3).removeClass('disabled')"),
                                                      #bsPopover("doUmap", "Run UMAP","Press to run UMAP Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                      textOutput("Umap.done"),
                                                      br()
                                               )),
                                             br(),
                                             plotlyOutput("Umap_2d_plot_1", width = "50%"),
                                             br(),
                                             column(12,
                                                    column(3,
                                                           downloadBttn('download_UMAP', 'Download UMAP Plot (as png)', size = 'sm'),
                                                    ),
                                                    column(3,
                                                           downloadBttn('download_UMAP_embedding', 'Download UMAP Embeddings (as csv)', size = 'sm'),
                                                    ),
                                             ),
                                          ),
                                    tabPanel( "tSNE", value = "test2",
                                             
                                             tags$script(
                                               '
                                                  var tab2 = $(\'a[data-value="test2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
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
                                                      actionBttn("doTsne", "Run TSNE", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab3).removeClass('disabled')"),
                                                      #bsPopover("doTsne", "Run tSNE Analysis","Press to run tSNE Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                      textOutput("Tsne.done"),
                                                      br()
                                               )),
                                             br(),
                                             plotlyOutput("Tsne_2d_plot_1", width = "50%"),
                                             br(),
                                             column(12,
                                                    column(3,
                                                           downloadBttn('download_Tsne', 'Download tSNE Plot (as png)', size = 'sm'),
                                                    ),
                                                    column(3,
                                                           downloadBttn('download_Tsne_embedding', 'Download tSNE Embeddings (as csv)', size = 'sm'),
                                                    ),
                                             ),
                                    ),
                                     
                                    tabPanel("Clustering", value = "test2",
                                             
                                             tags$script(
                                               '
                                                  var tab2 = $(\'a[data-value="test2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
                                             tabsetPanel(id="cluster",
                                                         tabPanel(title="Clustering", value="C_panel1",
                                                                  br(),
                                                                  fluidRow(
                                                                    column(12,
                                                                           column(2,
                                                                                  numericInput("clus.res",
                                                                                               label = "Resolution used",
                                                                                               value = 0.6,
                                                                                               min = 0.1,
                                                                                               step = 0.1)
                                                                           ),
                                                                           column(2,
                                                                                  selectInput("dim.used",
                                                                                              label = "PC to use",
                                                                                              choices = c(10:50)),
                                                                           ),
                                                                           column(2,
                                                                                  br(),
                                                                                  actionBttn("findCluster", "Find Clusters", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab3).removeClass('disabled')"),
                                                                                  #bsPopover("findCluster", "Run Clustering","Press to run Cluster Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                                                  textOutput("cluster.done"),
                                                                                  br()
                                                                           ),
                                                                    ),
                                                                           column(6,
                                                                                  plotlyOutput("Cluster2DPlot_1", width = "100%"),
                                                                                  br(),
                                                                           ),
                                                                    ),
                                                                  column(6,
                                                                         column(6,
                                                                                downloadBttn('download_Cluster', 'Download ClusterPlot (as png)', size = 'sm'),
                                                                         ),
                                                                         column(6,
                                                                                downloadBttn('download_ClusterTable', 'Download Cluster table (as csv)', size = 'sm'),
                                                                         ),
                                                                  ),
                                                                ),
                                                        
                                                         tabPanel(title="Determine cluster resolution", value="C_panel2",
                                                                  br(),
                                                                  column(9,
                                                                         h4(p("Determine cluster resolution:")),
                                                                       
                                                                         column(4,
                                                                                numericInput("clus.res_a",
                                                                                             label = "Resolution (from)",
                                                                                             value = 0.6,
                                                                                             min = 0.1,
                                                                                             step = 0.1)
                                                                         ),
                                                                         column(4,
                                                                                numericInput("clus.res_b",
                                                                                             label = "Resolution (to)",
                                                                                             value = 1,
                                                                                             min = 0.1,
                                                                                             step = 0.1)
                                                                         ),
                                                                         column(4,
                                                                                br(),
                                                                                actionBttn("findoptimumCluster", "Determine optimum resolution", icon = icon("hand-pointer-o"), size = 'sm'),
                                                                                #bsPopover("findoptimumCluster", "Determine optimum resolution","Press to determine optimum cluster resolution", placement = "bottom", trigger = "hover", options = NULL),
                                                                         ),
                                                                  br(),
                                                                  column(12,
                                                                         br(),
                                                                         br(),
                                                                         br(),
                                                                         br(),
                                                                         plotOutput("OptimumCluster2DPlot_1", width = "100%"),
                                                                         br(),
                                                                  ),
                                                                  column(12,
                                                                         column(6,
                                                                                downloadBttn('download_OptimumCluster', 'Download ClusterPlot (as png)', size = 'sm'),
                                                                         ),
                                                                         column(6,
                                                                                downloadBttn('download_OptimumClusterTable', 'Download Cluster table (as csv)', size = 'sm'),
                                                                         ),
                                                                    ),
                                                                  ),
                                                                  ),
                                                         
                                                         tabPanel(title="Subcluster Analysis", value="C_panel3",
                                                                  column(12,
                                                                         titlePanel(h4(p("Subcluster Analysis"))),
                                                                  ),
                                                                  column(3,
                                                                         br(),
                                                                         uiOutput("subcluster.gene.select"),
                                                                  ),
                                                                  column(3,
                                                                         br(),
                                                                         numericInput("subcluster.res",
                                                                                      label = "Resolution used",
                                                                                      value = 0.6,
                                                                                      min = 0.1,
                                                                                      step = 0.1),
                                                                         #bsTooltip("subcluster", "Perform Subcluster Analysis", placement = "bottom", trigger = "hover",
                                                                         #options = NULL)
                                                                  ),
                                                                  
                                                                  column(3,
                                                                         br(),
                                                                         br(),
                                                                         actionBttn("subcluster", "Run Subcluster analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                                         #bsPopover("subcluster", "Subcluster Analysis","Press to perform Subcluster Analysis for a specific cluster at any given resolution", placement = "bottom", trigger = "hover", options = NULL),
                                                                         textOutput("Subcluster.done"),
                                                                         br()
                                                                  ),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  column(9,
                                                                         plotlyOutput("subcluster_plot")
                                                                  ),
                                                                  column(12,
                                                                         column(4,
                                                                                downloadBttn('download_subcluster', 'Download Subcluster plot (as png)', size = 'sm'),
                                                                                br(),
                                                                         ),
                                                                  ),
                                                             ),
                                                          ),
                                                         ), 
                                     
                                     tabPanel("Cell type Identification", value = "test3",
                                              
                                              tags$script(
                                                '
                                                  var tab3 = $(\'a[data-value="test3"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              fluidRow(
                                                column(3,
                                                       selectInput("cellid_method",
                                                                   label = "Celltype annotation method",
                                                                   choices = c("CELLiD", "Celltypist"),
                                                                   selected = "CELLiD"),
                                                ),
                                              ),
                                              
                                              tabsetPanel(id="ct1",
                                                          tabPanel(title="Celltype Identification", value="Ct_panel1",
                                                                   conditionalPanel(
                                                                     condition = "input.cellid_method == 'CELLiD'",
                                                                     column(3,
                                                                            selectInput("cellatlas",
                                                                                        label = "Reference Atlas",
                                                                                        choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                                        selected = "all")
                                                                     ),
                                                                     column(3,
                                                                            br(),
                                                                            actionBttn("doCELLiD", "Run CELLiD", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab4).removeClass('disabled')"),
                                                                            #bsPopover("doCELLiD", "Perform Celltype Identification","Press to run celltype identification using CELLiD", placement = "bottom", trigger = "hover", options = NULL),
                                                                            textOutput("CELLiD.done"),
                                                                            br()
                                                                     ),
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     column(12,
                                                                            column(6,
                                                                                   plotlyOutput("Umap_cellid", width = "100%"),
                                                                            ),
                                                                            column(6,
                                                                                   plotlyOutput("Umap_cellid1", width = "100%"),
                                                                            ),
                                                                     ),
                                                                     DT::dataTableOutput("ct.table"),         
                                                    column(12,
                                                     column(4,
                                                            downloadBttn('download_Umap_cellid', 'Download CELLiD predict1 (as png)', size = 'sm'),
                                                            br(),
                                                     ),
                                                     column(4,
                                                            downloadBttn('download_Umap_cellid1', 'Download CELLiD predict2 (as png)', size = 'sm'),
                                                            br(),
                                                     ),
                                                     column(4,
                                                            downloadBttn('download_cellid_prediction', 'Download CELLiD predictions (in csv)', size = 'sm'),
                                                            br(),
                                                           ),
                                                    ), 
                                                     tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler3",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Spatial Transcriptomics)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                                     column(3,
                                                            br(),
                                                            actionBttn("dodeconv_spatial", "Go to Spatial Deconvolution", icon = icon("hand-pointer-o"), size = 'sm'),
                                                            #bsPopover("dodeconv_spatial", "Perform Celltype Deconvolution","Press to navigate to Spatial module and perform celltype deconvolution", placement = "bottom", trigger = "hover", options = NULL),
                                                            ),
                                                     
                                                     tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler5",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell ATAC-seq)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                                     column(3,
                                                            br(),
                                                            actionBttn("doct_atac", "Annotate cell types for ATAC-data", icon = icon("hand-pointer-o"), size = 'sm'),
                                                            #bsPopover("doct_atac", "Perform Celltype Label Transfer","Press to navigate to scATAC-seq module and perform celltype label transfer", placement = "bottom", trigger = "hover", options = NULL),
                                                           ),
                                                       ),
                                                     conditionalPanel(
                                                       condition = "input.cellid_method == 'Celltypist'",
                                                       column(3,
                                                              selectInput("celltypistatlas",
                                                                          label = "Reference Atlas",
                                                                          choices = c("Immune_All_Low.pkl", "Autopsy_COVID19_Lung.pkl", "Pan_Fetal_Human.pkl", "Nuclei_Lung_Airway.pkl", "Developing_Human_Thymus.pkl", "Human_Lung_Atlas.pkl", "Developing_Mouse_Brain.pkl", "Developing_Human_Brain.pkl", "Cells_Lung_Airway.pkl", "Healthy_COVID19_PBMC.pkl", "Human_IPF_Lung.pkl", "Adult_Mouse_Gut.pkl", "Immune_All_High.pkl", "COVID19_Immune_Landscape.pkl", "Human_PF_Lung.pkl", "COVID19_HumanChallenge_Blood.pkl", "Lethal_COVID19_Lung.pkl", "Cells_Fetal_Lung.pkl", "Cells_Intestinal_Tract.pkl"),
                                                                          selected = "Immune_All_Low.pkl")
                                                       ),
                                                       column(3,
                                                              br(),
                                                              actionBttn("doCelltypist", "Run Celltypist", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab4).removeClass('disabled')"),
                                                              #bsPopover("doCelltypist", "Perform Celltype Identification","Press to Press to run celltype identification using CellTypist", placement = "bottom", trigger = "hover", options = NULL),
                                                              textOutput("Celltypist.done"),
                                                              br()
                                                             ),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       column(12,
                                                              column(6,
                                                                     plotlyOutput("Umap_celltypist", width = "100%"),
                                                                     ),
                                                              column(6,
                                                                     plotlyOutput("Umap_celltypist1", width = "100%"),
                                                                     ),
                                                              ),
                                                       DT::dataTableOutput("celltypist.table"),
                                                       #DT::dataTableOutput(outputId = "recoding"),
                                                       #DT::dataTableOutput(outputId = "newVars")
                                                       tags$head(tags$script('
                                                          Shiny.addCustomMessageHandler("myCallbackHandler3a",
                                                          function(typeMessage) {console.log(typeMessage)
                                                          if(typeMessage == 2){
                                                          console.log("got here");
                                                          $("a:contains(Spatial Transcriptomics)").click();
                                                          }
                                                          });
                                                        ')), 
                                                       column(3,
                                                              br(),
                                                              actionBttn("dodeconv_spatial1", "Go to Spatial Deconvolution", icon = icon("hand-pointer-o"), size = 'sm'),
                                                              #bsPopover("dodeconv_spatial1", "Perform Celltype Deconvolution","Press to navigate to Spatial module and perform celltype deconvolution", placement = "bottom", trigger = "hover", options = NULL),
                                                       ),
                                                       tags$head(tags$script('
                                                          Shiny.addCustomMessageHandler("myCallbackHandler5a",
                                                          function(typeMessage) {console.log(typeMessage)
                                                          if(typeMessage == 2){
                                                          console.log("got here");
                                                          $("a:contains(Single cell ATAC-seq)").click();
                                                          }
                                                          });
                                                          ')), 
                                                       column(3,
                                                              br(),
                                                              actionBttn("doct_atac1", "Annotate cell types for ATAC-data", icon = icon("hand-pointer-o"), size = 'sm'),
                                                              #bsPopover("doct_atac1", "Perform Celltype Label Transfer","Press to navigate to scATAC-seq module and perform celltype label transfer", placement = "bottom", trigger = "hover", options = NULL),
                                                       ),
                                                    ),
                                                  ),
                                                          
                                                          tabPanel(title="Rename Clusters", value="Ct_panel2",
                                                                   conditionalPanel(
                                                                     condition = "input.cellid_method == 'CELLiD'",
                                                                     column(12,
                                                                            column(6,
                                                                                   br(),
                                                                                   h4(p("Please modify your celltypes below:")),
                                                                                   br(),
                                                                                   actionBttn("commitButton", "Rename clusters (if needed)", size = "sm"),
                                                                                   #bsPopover("commitButton", "Rename clusters","Press to rename clusters", placement = "bottom", trigger = "hover", options = NULL),
                                                                                   br(),
                                                                                   br(),
                                                                                   rHandsontableOutput("hot"),       
                                                                            ),
                                                                            column(6,
                                                                                   h4(p("UMAP with renamed celltypes")),
                                                                                   plotlyOutput("Umap_cellid2", width = "100%"),
                                                                                   br(),
                                                                             ),
                                                                           ),
                                                                       ),
                                                                   conditionalPanel(
                                                                    condition = "input.cellid_method == 'Celltypist'",
                                                                    br(),
                                                                    h3(p("Not supported yet"))
                                                                   #  column(12,
                                                                   # column(6,
                                                                   #         br(),
                                                                   #        h4(p("Please modify your celltypes below:")),
                                                                   #          br(),
                                                                   #         #actionBttn("commitButton1", "Rename clusters", size = "sm"),
                                                                   #          br(),
                                                                   #          br(),
                                                                            #rHandsontableOutput("cot"),       
                                                                   #    ),
                                                                   #   column(6,
                                                                   #          h4(p("UMAP with renamed celltypes")),
                                                                   #          plotOutput("Umap_celltypist2", width = "100%"),
                                                                   #          br(),
                                                                   #   ),
                                                                     #),
                                                                     
                                                                   # ),
                                                                   ),
                                                                ),
                                                          
                                                          tabPanel(title="Visualize", value="Ct_panel3",
                                                                   conditionalPanel(
                                                                     condition = "input.cellid_method == 'CELLiD'",
                                                                   column(12,
                                                                          column(6,
                                                                                 br(),
                                                                                 uiOutput("ct.gene.select"),
                                                                                 actionBttn("Vis_seurat1", "Visualize", icon = icon("hand-pointer-o"), size = 'sm'),
                                                                                 plotlyOutput("ct.gene1.plot", width = "100%"),
                                                                                 br(),
                                                                                 ),
                                                                          column(6,
                                                                                 br(),
                                                                                 br(),
                                                                                 plotlyOutput("ct.gene.plot", width = "150%"),
                                                                                 br(),
                                                                                 ),
                                                                              ),
                                                                   
                                                                  column(12,
                                                                   column(4,
                                                                          downloadBttn('download_violn1', 'Download Violin plot (as png)', size = 'sm'),
                                                                          br(),
                                                                    ),
                                                                   column(4,
                                                                          downloadBttn('download_feature1', 'Download Feature plot (as png)', size = 'sm'),
                                                                          br(),
                                                                     ),
                                                                   ),
                                                                            ),
                                                                   conditionalPanel(
                                                                     condition = "input.cellid_method == 'Celltypist'",
                                                                     column(6,
                                                                            br(),
                                                                            uiOutput("celltypist.gene.select"),
                                                                            plotlyOutput("celltypist.gene1.plot", width = "100%"),
                                                                            br(),
                                                                     ),
                                                                     column(6,
                                                                            plotlyOutput("celltypist.gene.plot", width = "150%"),
                                                                            br(),
                                                                     ),
                                                                   ),
                                                          ),
                                                       
                                              ),
                                           ),
                                     tabPanel("Cell type Similarity", value = "test4",
                                              
                                              tags$script(
                                                '
                                                  var tab4 = $(\'a[data-value="test4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              column(3,
                                                     br(),
                                                     selectInput("cell1",
                                                                 label = "Group by",
                                                                 choices = c("seurat_clusters", "primary.predict", "newID"),
                                                                 selected = "primary.predict"),
                                                          ),
                                    
                                                  column(3,
                                                         br(),
                                                          selectInput("corr_method",
                                                                  label = "Statistics",
                                                                  choices = c("pearson", "spearman", "kendall"),
                                                                  selected = "pearson"),
                                                        ),                          
                                              column(3,
                                                     br(),
                                                     br(),
                                                     actionBttn("cell_cell", "Run celltype similarity", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab5).removeClass('disabled')"),
                                                     #bsPopover("cell_cell", "Perform Celltype Similarity","Press to run Celltype Similarity Analysis", placement = "bottom", trigger = "hover", options = NULL),
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
                                                            downloadBttn('download_cell_cell_sim', 'Download Celltype similarity plot (as png)', size = 'sm'),
                                                     ),
                                                     column(4,
                                                            downloadBttn('download_cor.table', 'Download Celltype similarity table (in csv)', size = 'sm'),
                                                     ),
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              DT::dataTableOutput("cor.table")
                                     ),
                                      tabPanel("DEGs", value = "test5",
                                               
                                               tags$script(
                                                 '
                                                  var tab5 = $(\'a[data-value="test5"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                               ),
                                               
                                               fluidRow(
                                               column(3,
                                                      selectInput("deg_method",
                                                                  label = "Type of DEG analysis",
                                                                  choices = c("Celltype specific", "Pairwise DEGs"),
                                                                  selected = "Celltype specific"),
                                                      ),
                                                    ),
                                               conditionalPanel(
                                                 condition = "input.deg_method == 'Celltype specific'",
                                               fluidRow(
                                                 column(3,
                                                 selectInput("deg1",
                                                             label = "Group by",
                                                             choices = c("seurat_clusters", "primary.predict", "newID"),
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
                                                        actionBttn("doDeg", "Run DEGs", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab6).removeClass('disabled')"),
                                                        #bsPopover("doDeg", "Perform DEG Analysis","Press to run DEG Analysis", placement = "bottom", trigger = "hover", options = NULL),
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
                                            conditionalPanel(
                                              condition = "input.deg_method == 'Pairwise DEGs'",
                                              #br(),
                                              column(12,
                                                     column(3,
                                                            selectInput("deg3",
                                                                        label = "Group by",
                                                                        choices = c("seurat_clusters", "primary.predict", "newID"),
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
                                                            actionBttn("doVolcano", "Run Pairwise DEGs", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab6).removeClass('disabled')"),
                                                            #bsPopover("doVolcano", "Perform Pairwise DEG Analysis","Press to run Pairwise DEG Analysis", placement = "bottom", trigger = "hover", options = NULL),
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
                                      ),
                                     
                                     tabPanel("Data visualization", value = "test5",
                                              
                                              tags$script(
                                                '
                                                  var tab5 = $(\'a[data-value="test5"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              fluidRow(
                                                column(3,
                                                selectInput("deg2",
                                                            label = "Group by",
                                                            choices = c("seurat_clusters", "primary.predict", "newID"),
                                                            selected = "primary.predict"),
                                                ),
                                                column(6,
                                                       uiOutput("deg.gene.select"),
                                                       actionBttn("Vis_seurat", "Visualize", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab6).removeClass('disabled')"),
                                                       plotlyOutput("Deg.plot", width = "150%"),
                                                       br(),
                                                       plotlyOutput("Deg1.plot", width = "150%"),
                                                       br(),
                                                       plotOutput("Deg2.plot", width = "150%"),
                                                ),
                                                column(4,
                                                       downloadBttn('download_violn', 'Download Violin plot (as png)', size = 'sm'),
                                                       br(),
                                                      ),
                                                column(4,
                                                       downloadBttn('download_feature', 'Download Feature plot (as png)', size = 'sm'),
                                                       br(),
                                                      ),
                                                column(4,
                                                       downloadBttn('download_ridge', 'Download Feature plot (as png)', size = 'sm'),
                                                       br(),
                                                ),
                                              )
                                     ),
                                     tabPanel("GSEA", value = "test6",
                                              
                                              tags$script(
                                                '
                                                  var tab6 = $(\'a[data-value="test6"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
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
                                                     actionBttn("gsea", "Run gene set enrichment analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                     #bsPopover("gsea", "Perform Gene Set Enrichment Analysis","Press to run Gene Set Enrichment Analysis", placement = "bottom", trigger = "hover", options = NULL),
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
                                                     downloadBttn('download_gsea.table', 'Download GSEA Results (in csv)', size = 'sm'),
                                              ),
                                     ),
                                     
                                     tabPanel("Cell-cell communication", value = "test6",
                                              
                                              tags$script(
                                                '
                                                  var tab6 = $(\'a[data-value="test6"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              fluidRow(
                                              column(3,
                                                     br(),
                                                     selectInput("cc1a",
                                                     label = "Group by",
                                                     choices = c("seurat_clusters", "primary.predict"),
                                                     selected = "primary.predict"),
                                                    ),
                                              column(3, 
                                                     br(),
                                                     selectInput("cc_method",
                                                                    label = "Method",
                                                                    choices = c("natmi", "connectome", "logfc", "sca", "cellphonedb", "cytotalk"), selected = "cellphonedb")
                                                                 ),
                                                    
                                              column(3, 
                                                     br(),
                                                     selectInput("cc_resource",
                                                                 label = "Method",
                                                                 choices = c("Default", "Consensus",  "Baccin2019", "CellCall", "CellChatDB", "Cellinker", "CellPhoneDB", "CellTalkDB", "connectomeDB2020", "EMBRACE", "Guide2Pharma", "HPMR", "ICELLNET", "iTALK", "Kirouac2010", "LRdb", "Ramilowski2015", "OmniPath"), selected = "CellPhoneDB")
                                                             ),
                                              
                                              column(3,
                                                     br(),
                                                     br(),
                                                     actionBttn("doCC", "Run Analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                     #bsPopover("doCC", "Perform Cell-Cell Communication Analysis","Press to run Cell-Cell Communication Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                     textOutput("cc.done"),
                                                     br()
                                              ),
                                            ),
                                              br(),
                                              br(),
                                              column(12,
                                              plotOutput("CC_plot1", width = "100%"),
                                              ),
                                              column(3,
                                                     uiOutput("CC.gene.select"),
                                              ),
                                              column(3,
                                                     uiOutput("CC.gene1.select"),
                                              ),
                                              br(),
                                              column(12,
                                                     plotOutput("CC_plot2", width = "100%"),
                                              ),
                                              column(12,
                                              h4(p("Interacting partners")),
                                                      ),
                                              withSpinner(dataTableOutput('cc.table')),
                                              br(),
                                            #column(12,
                                             # h4(p("p-values for the all the interacting partners")),
                                           #       ),
                                           #   withSpinner(dataTableOutput('cc.table1')),
                                              ),
                                     #tabPanel("13. Annotate scATAC-seq data using scRNA-seq data",
                                     #         tags$head(tags$script('
                                     #                               Shiny.addCustomMessageHandler("myCallbackHandler1",
                                     #                                function(typeMessage) {console.log(typeMessage)
                                     #                               if(typeMessage == 2){
                                     #                              console.log("got here");
                                     #                             $("a:contains(Single cell ATAC-seq)").click();
                                     #                             }
                                     #                             });
                                     #                             ')), 
                                     #         
                                     #         column(5,
                                     #                br(),
                                     #                actionBttn("loadexample_atacseq", "Load example scATAC data and annotate", icon = icon("hand-pointer-o"), size = 'sm'),
                                     #                #bsPopover("loadexample_atacseq", "Load example scATAC data and annotate","Press to load example scATAC data and annotate", placement = "bottom", trigger = "hover", options = NULL),
                                     #                ),
                                     #         column(4,
                                     #                br(),
                                     #                actionBttn("process_atacseq", "Process scATAC-seq data", icon = icon("hand-pointer-o"), size = 'sm'),
                                                     #bsPopover("process_atacseq", "Process scATAC-seq data","Press to navigate to scATAC-seq module and process", placement = "bottom", trigger = "hover", options = NULL),
                                     #                ),
                                     #                br(),
                                     #         column(12,
                                     #                plotOutput("annotate_scRNA_ATAC_plot", width = "70%"),
                                     #                plotOutput("annotate_scRNA_ATAC_plot1", width = "70%")
                                     #          ),
                                     #      )
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
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                titlePanel(h4(p("Load your example data"))),
                                                                actionBttn("loadexample1", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_intg).removeClass('disabled')"),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == '10X cellranger'",
                                                                titlePanel(h4(p("Load your example data"))),
                                                              actionBttn("loadexample1a", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_intg).removeClass('disabled')"),
                                                              ),
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
                                                                ),
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == '10X cellranger'",
                                                                fileInput("scH5_1",
                                                                          label = "Cellranger output (Accepted Format: .h5)",
                                                                          accept = ".h5",
                                                                          multiple = T),
                                                                ),
                                                              column(6,
                                                                     numericInput(inputId = "min.genes1",
                                                                                  label = "Min. genes",
                                                                                  value = 200,
                                                                                  min = 1)
                                                              ),
                                                              column(6,
                                                                     numericInput(inputId = "min.cells1",
                                                                                  label = "Min. cells",
                                                                                  value = 3,
                                                                                  min = 1)
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scInput1 == 'Raw Counts Matrix' || input.scInput1 == '10X cellranger'",
                                                                textInput(inputId = "projName1",
                                                                          label = "Project Name",
                                                                          value = "Integration")),
                                                              
                                                              fluidRow(
                                                                
                                                                actionBttn("loadButton1", "Load data", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_intg).removeClass('disabled')"),
                                                                actionBttn("reset_intg", "Reset", icon = icon("repeat"), size = "sm")
                                                              ),
                                                            )),
                                                     chooseSliderSkin("Modern"),
                                                     titlePanel(h4(p("Quality control"))),
                                                     column(6,
                                                            plotOutput("nFeature_RNAPlot1", width = "200%")
                                                     ),
                                                     column(12,
                                                            column(3,
                                                                   numericInput("ob1a",
                                                                                label = "Min nFeature:",
                                                                                value = 200,
                                                                                min = 0,
                                                                                step = 1),
                                                            ),
                                                            column(3,
                                                                   numericInput("ob2a",
                                                                                label = "Max nFeature:",
                                                                                value = 2500,
                                                                                min = 0,
                                                                                step = 1),
                                                            ),
                                                            column(3,
                                                                   numericInput("ob3a",
                                                                                label = "Mt%:",
                                                                                value = 5,
                                                                                min = 0,
                                                                                step = 1),
                                                            ),
                                                     column(3,
                                                            br(),
                                                            actionBttn("filter_seurat1", "Filter", icon = icon("hand-o-right"), size = "sm")
                                                     ),
                                                    ), 
                                                     
                                              ),
                                              
                                              column(12,
                                                     withSpinner(dataTableOutput('countdataDT1'))
                                              )),
                                     
                                     #tabPanel("2. Quality control",
                                    #      tabsetPanel(id="qc_intg",
                                    #        tabPanel("Violin Plot",            
                                    #          column(3,
                                    #                 plotOutput("nFeature_RNAPlot1")
                                    #          ),
                                    #          column(3,
                                    #                 plotOutput("mitoPlot1")
                                    #          ),
                                    #          column(3,
                                    #                 plotOutput("nCount_RNAPlot1")
                                    #          )),
                                    #        tabPanel("Feature Scatter", 
                                    #                 column(5,
                                    #                        plotlyOutput("FeatureScatterPlot1a")
                                    #                 ),
                                    #                 column(5,
                                    #                        plotlyOutput("FeatureScatterPlot2a")
                                    #                 ),        
                                    #       )),
                                    # ),
                                     
                                     tabPanel("Normalization and Variable Feature Selection", value = "test_intg",
                                              
                                              tags$script(
                                                '
                                                  var tab_intg = $(\'a[data-value="test_intg"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
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
                                                       actionBttn("findVarGenes_bef_intg", "Identify highly variable genes", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg1).removeClass('disabled')"),
                                                       )),
                                              plotOutput("VarGenes_bef_intg", width = "100%")       
                                     ),
                                     
                                     tabPanel("Before Data Integration", value = "test_intg1",
                                              
                                              tags$script(
                                                '
                                                  var tab_intg1 = $(\'a[data-value="test_intg1"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              tabsetPanel(id="bef_data_integration",
                                                          
                                                          tabPanel("PCA", 
                                                                   br(),
                                                                  
                                                                     tabsetPanel(id="Pca_bef_intg_seurat",
                                                                                 tabPanel(title="PCA Plot",
                                                                                          br(),
                                                                                          fluidRow(
                                                                                            column(3,
                                                                                                   actionBttn("runPCA_bef_intg", "Run PCA", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg2).removeClass('disabled')")
                                                                                            ),
                                                                                          ),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                                            plotlyOutput("PCAplot_bef_tpm1", width = "75%"),
                                                                                            plotlyOutput("PCAplot_bef_tpm2", width = "75%"),
                                                                                            plotlyOutput("PCAplot_bef_tpm3", width = "75%")),
                                                                                          br(),
                                                                                          conditionalPanel(
                                                                                            condition = "input.scInput1 == '10X cellranger'",
                                                                                            plotlyOutput("PCAplot_bef_h5_1", width = "75%"),
                                                                                            plotlyOutput("PCAplot_bef_h5_2", width = "75%")),
                                                                                 ),
                                                                                 tabPanel(title="PC Gene Visualisation",
                                                                                          br(),
                                                                                          selectInput("select.pc_bef_intg",
                                                                                                      label = "PC to plot",
                                                                                                      choices = c(1:50)
                                                                                          ),
                                                                                          fluidRow(
                                                                                            column(4,
                                                                                                   plotOutput("vizPlot_bef_intg", width = "100%", height = "600px")
                                                                                            ),
                                                                                            column(8,
                                                                                                   plotOutput("PCHeatmap_bef_intg", width = "100%", height = "600px")
                                                                                            )
                                                                                          ),
                                                                                          DT::dataTableOutput("PCtable_bef_intg")
                                                                                 ),
                                                                                 
                                                                                 tabPanel(title="Elbow", 
                                                                                          br(),
                                                                                          br(),
                                                                                          plotOutput("Elbow_bef_intg", width = "100%")
                                                                                          
                                                                                 )),
                                                                            ),
                                                          
                                                          tabPanel("UMAP",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionBttn("runUMAP_bef_intg", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm"),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("UMAPplot_bef_tpm1", width = "75%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_tpm2", width = "75%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_tpm3", width = "75%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("UMAPplot_bef_h5_1", width = "75%"),
                                                                              br(),
                                                                              plotlyOutput("UMAPplot_bef_h5_2", width = "75%")),
                                                                            ),
                                                                          ),
                                                          
                                                          tabPanel("tSNE",
                                                                     br(),
                                                                     column(3,
                                                                            numericInput("dim.used_bef_intg",
                                                                                         label = "Dimensions used",
                                                                                         value = 10)
                                                                     ),
                                                                     column(9,
                                                                            br(),
                                                                            actionBttn("runTSNE_bef_intg", "Run TSNE", icon = icon("hand-pointer-o"), size = "sm"),
                                                                            #textOutput("Intg.done"),
                                                                            br(),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                              plotlyOutput("TSNEplot_bef_tpm1", width = "75%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_tpm2", width = "75%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_tpm3", width = "75%")),
                                                                            br(),
                                                                            conditionalPanel(
                                                                              condition = "input.scInput1 == '10X cellranger'",
                                                                              plotlyOutput("TSNEplot_bef_h5_1", width = "75%"),
                                                                              br(),
                                                                              plotlyOutput("TSNEplot_bef_h5_2", width = "75%")),
                                                                     ),
                                                                  ),
                                                          
                                                          #tabPanel("Clustering",
                                                          #        br(),
                                                          #        fluidRow(
                                                          #          column(3,
                                                          #                 numericInput("clus.res_bef_intg",
                                                          #                              label = "Resolution used",
                                                          #                              value = 0.6,
                                                          #                              min = 0.1,
                                                          #                              step = 0.1)
                                                          #          ),
                                                          #          column(3,
                                                          #                 selectInput("dim.used_bef_intg",
                                                          #                             label = "PC to use",
                                                          #                             choices = c(10:50)),
                                                          #          ),
                                                          #          column(3,
                                                          #                 br(),
                                                          #                 actionBttn("findCluster_bef_intg", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm"),
                                                          #                 #textOutput("cluster1.done"),
                                                          #          ),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          plotlyOutput("Cluster2DPlot_bef_intg", width = "75%")
                                                          #        ),
                                                                   
                                                          #        column(12,
                                                          #               h4(p("Determine cluster resolution:")),
                                                          #               column(3,
                                                          #                      numericInput("clus.res_a1",
                                                          #                                   label = "Resolution (from)",
                                                          #                                   value = 0.6,
                                                          #                                   min = 0.1,
                                                          #                                   step = 0.1)
                                                          #               ),
                                                          #               column(3,
                                                          #                      numericInput("clus.res_b1",
                                                          #                                   label = "Resolution (to)",
                                                          #                                   value = 1,
                                                          #                                   min = 0.1,
                                                          #                                   step = 0.1)
                                                          #               ),
                                                          #               column(4,
                                                          #                      br(),
                                                          #                      actionBttn("findoptimumCluster1", "Determine optimum resolution", icon = icon("hand-pointer-o"), size = 'sm'),
                                                          #                      #bsPopover("findoptimumCluster", "Determine optimum resolution","Press to determine optimum cluster resolution", placement = "bottom", trigger = "hover", options = NULL),
                                                          #               ),
                                                          #        ),
                                                          #        br(),
                                                          #        column(9,
                                                          #               br(),
                                                          #               plotOutput("OptimumCluster2DPlot_a1", width = "100%"),
                                                          #               br(),
                                                          #               ),
                                                          #        column(6,
                                                          #               column(6,
                                                          #                      downloadBttn('download_OptimumCluster1', 'Download ClusterPlot (as png)', size = 'sm'),
                                                          #               ),
                                                          #               column(6,
                                                          #                      downloadBttn('download_OptimumClusterTable1', 'Download Cluster table (as csv)', size = 'sm'),
                                                          #               ),
                                                          #           ),
                                                          #        column(12,
                                                          #               titlePanel(h4(p("Subcluster Analysis"))),
                                                          #        ),
                                                          #        column(3,
                                                          #               br(),
                                                          #               uiOutput("subcluster.gene.select1"),
                                                          #        ),
                                                          #        column(3,
                                                          #               br(),
                                                          #               numericInput("subcluster.res1",
                                                          #                            label = "Resolution used",
                                                          #                            value = 0.6,
                                                          #                            min = 0.1,
                                                          #                            step = 0.1),
                                                          #               #bsTooltip("subcluster", "Perform Subcluster Analysis", placement = "bottom", trigger = "hover",
                                                          #               #options = NULL)
                                                          #        ),
                                                                  
                                                          #         column(3,
                                                          #               br(),
                                                          #               br(),
                                                          #               actionBttn("subcluster1", "Run Subcluster analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                                         #bsPopover("subcluster", "Subcluster Analysis","Press to perform Subcluster Analysis for a specific cluster at any given resolution", placement = "bottom", trigger = "hover", options = NULL),
                                                          #                textOutput("Subcluster1.done"),
                                                          #               br()
                                                          #        ),
                                                          #        br(),
                                                          #        br(),
                                                          #        br(),
                                                          #        column(9,
                                                          #               plotlyOutput("subcluster_plot1")
                                                          #        ),
                                                          #        column(12,
                                                          #               column(4,
                                                          #                      downloadBttn('download_subcluster1', 'Download Subcluster plot (as png)', size = 'sm'),
                                                          #                      br(),
                                                          #               ),
                                                          #        ),
                                                          #),
                                                          
                                                          #tabPanel("Cell type identification",
                                                                  #conditionalPanel(
                                                                   #   condition = "input.scAnalysis_integ == 'Seurat'",
                                                          #           br(),
                                                          #          column(3,
                                                          #                 selectInput("cellid_method1",
                                                          #                             label = "Celltype annotation method",
                                                          #                             choices = c("CELLiD", "Celltypist"),
                                                          #                             selected = "CELLiD"),
                                                          #          ),
                                                          #          conditionalPanel(
                                                          #            condition = "input.cellid_method1 == 'CELLiD'",
                                                          #          column(3,
                                                          #                 selectInput("cellatlas1",
                                                          #                             label = "Reference Atlas",
                                                          #                             choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                          #                             selected = "all")
                                                          #          ),
                                                          #          column(3,
                                                          #                 br(),
                                                          #                 actionBttn("doCELLiD_bef_intg", "Run CELLiD", icon = icon("hand-pointer-o"), size = "sm"),
                                                          #                 textOutput("CELLiD1.done"),
                                                          #                 br()
                                                          #          ),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          plotlyOutput("Umap_cellid_bef_intg", width = "50%"),
                                                          #          plotlyOutput("Umap_cellid_bef_intg1", width = "50%"),
                                                          #          br(),
                                                          #          br(),
                                                          #          DT::dataTableOutput("ct_bef_intg.table")
                                                                    #),
                                                          #         ),
                                                                   
                                                          #         conditionalPanel(
                                                          #          condition = "input.cellid_method1 == 'Celltypist'",
                                                          #          column(3,
                                                          #                 selectInput("celltypistatlas1",
                                                          #                             label = "Reference Atlas",
                                                          #                             choices = c("Immune_All_Low.pkl", "Autopsy_COVID19_Lung.pkl", "Pan_Fetal_Human.pkl", "Nuclei_Lung_Airway.pkl", "Developing_Human_Thymus.pkl", "Human_Lung_Atlas.pkl", "Developing_Mouse_Brain.pkl", "Developing_Human_Brain.pkl", "Cells_Lung_Airway.pkl", "Healthy_COVID19_PBMC.pkl", "Human_IPF_Lung.pkl", "Adult_Mouse_Gut.pkl", "Immune_All_High.pkl", "COVID19_Immune_Landscape.pkl", "Human_PF_Lung.pkl", "COVID19_HumanChallenge_Blood.pkl", "Lethal_COVID19_Lung.pkl", "Cells_Fetal_Lung.pkl", "Cells_Intestinal_Tract.pkl"),
                                                          #                             selected = "Immune_All_Low.pkl")
                                                          #          ),
                                                          #          column(3,
                                                          #                 br(),
                                                          #                 actionBttn("doCelltypist_bef_intg", "Run Celltypist", icon = icon("hand-pointer-o"), size = 'sm'),
                                                          #                 textOutput("Celltypist1.done"),
                                                          #                 br()
                                                          #          ),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          br(),
                                                          #          plotlyOutput("Umap_celltypist_bef_intg", width = "50%"),
                                                          #          plotlyOutput("Umap_celltypist_bef_intg1", width = "50%"),
                                                          #          br(),
                                                          #          br(),
                                                          #          DT::dataTableOutput("ct_celltypist_bef_intg.table")
                                                          #        ),
                                                          #     )
                                                            )
                                                           ),
                                              
                                     tabPanel("Data Integration", value = "test_intg2",
                                              
                                              tags$script(
                                                '
                                                  var tab_intg2 = $(\'a[data-value="test_intg2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                       fluidRow(
                                         column(12,
                                                br(),
                                                
                                                column(3,
                                                selectInput("scAnalysis_integ",
                                                            label = "Analysis method",
                                                            choices = c("Seurat", "Harmony", "scVI", "fastMNN"),
                                                            selected = "Seurat"),
                                                ),
                                                #br(),
                                                
                                                conditionalPanel(
                                                  condition = "input.scAnalysis_integ == 'Seurat'",
                                                  column(4,
                                                         numericInput("nfeatures_intg_seurat",
                                                                      label = "Integration features",
                                                                      value = 2000,
                                                                      min = 500,
                                                                      step = 1)),
                                                column(4,
                                                       br(),
                                                       actionBttn("doIntg_seurat", "Run Data Integration", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg3).removeClass('disabled')"),
                                                      ),
                                                    ),
                                         conditionalPanel(
                                           condition = "input.scAnalysis_integ == 'Harmony'",
                                           column(4,
                                           numericInput("nfeatures_intg_harmony",
                                                        label = "Integration features",
                                                        value = 2000,
                                                        min = 500,
                                                        step = 1)),
                                           br(),
                                           column(4,
                                           actionBttn("doIntg_harmony", "Run Data Integration", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg3).removeClass('disabled')"),
                                           ),
                                         ),
                                         conditionalPanel(
                                           condition = "input.scAnalysis_integ == 'scVI'",
                                           column(3,
                                           numericInput("nfeatures_intg_scvi",
                                                        label = "Integration features",
                                                        value = 2000,
                                                        min = 500,
                                                        step = 1),
                                           ),
                                           column(3,
                                                  numericInput("scvi_latent",
                                                               label = "Number of dimensions",
                                                               value = 10,
                                                               min = 1,
                                                               step = 1),
                                           ),
                                           column(3,
                                           br(),
                                           actionBttn("doIntg_scvi", "Run Data Integration", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg3).removeClass('disabled')"),
                                                      ),
                                           column(12,
                                                  h4('This may take long time to run (around 30 mins to 1 hr)'),
                                                  ),
                                               br(),
                                             ),
                                         conditionalPanel(
                                           condition = "input.scAnalysis_integ == 'fastMNN'",
                                           column(4,
                                           numericInput("nfeatures_intg_fastmnn",
                                                        label = "Integration features",
                                                        value = 2000,
                                                        min = 500,
                                                        step = 1)),
                                           br(),
                                           column(4,
                                           actionBttn("doIntg_fastmnn", "Run Data Integration", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg3).removeClass('disabled')"),
                                           )),
                                         ))),
                                    
                                     tabPanel("PCA", value = "test_intg3",
                                              
                                              tags$script(
                                                '
                                                  var tab_intg3 = $(\'a[data-value="test_intg3"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'Seurat'",
                                              tabsetPanel(id="Pca1",
                                                          tabPanel(title="PCA Plot",
                                                                   br(),
                                                                   column(3,
                                                                          actionBttn("runPCA_intg_seurat", "Run PCA", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg4).removeClass('disabled')"),
                                                                   ),
                                                                   conditionalPanel(
                                                                     condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_seurat_tpm1", width = "75%"),
                                                                     plotlyOutput("PCAplot_seurat_tpm2", width = "75%"),
                                                                     plotlyOutput("PCAplot_seurat_tpm3", width = "75%")),
                                                                     br(),
                                                                   conditionalPanel(
                                                                     condition = "input.scInput1 == '10X cellranger'",
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_seurat_h5_1", width = "75%"),
                                                                     plotlyOutput("PCAplot_seurat_h5_2", width = "75%")),
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
                                                                   br(),
                                                                   plotOutput("Elbow_intg_seurat", width = "100%")
                                                          ))),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'Harmony'",
                                                tabsetPanel(id="Pca1",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     column(3,
                                                                            actionBttn("runPCA_intg_harmony", "Run PCA", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg4).removeClass('disabled')"),
                                                                     ),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                       br(),
                                                                       br(),
                                                                       br(),
                                                                       plotlyOutput("PCAplot_harmony_tpm1", width = "75%"),
                                                                       plotlyOutput("PCAplot_harmony_tpm2", width = "75%"),
                                                                       plotlyOutput("PCAplot_harmony_tpm3", width = "75%")),
                                                                     br(),
                                                                   conditionalPanel(
                                                                     condition = "input.scInput1 == '10X cellranger'",
                                                                     br(),
                                                                     br(),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_harmony_h5_1", width = "75%"),
                                                                     plotlyOutput("PCAplot_harmony_h5_2", width = "75%")),
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
                                                condition = "input.scAnalysis_integ == 'scVI'",
                                                tabsetPanel(id="Pca2",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     column(3,
                                                                            actionBttn("runPCA_intg_scvi", "Run PCA", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg4).removeClass('disabled')"),
                                                                     ),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                       br(),
                                                                       br(),
                                                                       br(),
                                                                       plotlyOutput("PCAplot_scvi_tpm1", width = "75%"),
                                                                       plotlyOutput("PCAplot_scvi_tpm2", width = "75%"),
                                                                       plotlyOutput("PCAplot_scvi_tpm3", width = "75%")),
                                                                      br(),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == '10X cellranger'",
                                                                       br(),
                                                                       br(),
                                                                       br(),
                                                                       plotlyOutput("PCAplot_scvi_h5_1", width = "75%"),
                                                                       plotlyOutput("PCAplot_scvi_h5_2", width = "75%")),
                                                            ),
                                                          )),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_integ == 'fastMNN'",
                                                tabsetPanel(id="Pca2",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     column(3,
                                                                            actionBttn("runPCA_intg_fastmnn", "Run PCA", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg4).removeClass('disabled')"),
                                                                     ),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                                       br(),
                                                                       br(),
                                                                       br(),
                                                                       plotlyOutput("PCAplot_fastmnn_tpm1", width = "100%"),
                                                                       plotlyOutput("PCAplot_fastmnn_tpm2", width = "100%"),
                                                                       plotlyOutput("PCAplot_fastmnn_tpm3", width = "100%")),
                                                                      br(),
                                                                     conditionalPanel(
                                                                       condition = "input.scInput1 == '10X cellranger'",
                                                                       br(),
                                                                       br(),
                                                                       br(),
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
                                    
                                    tabPanel("UMAP", value = "test_intg4",
                                             
                                             tags$script(
                                               '
                                                  var tab_intg4 = $(\'a[data-value="test_intg4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
                                               column(3,
                                                      numericInput("dim.used_intg",
                                                                   label = "Dimensions used",
                                                                   value = 10, step = 1)
                                               ),
                                               column(9,
                                                      br(),
                                                      actionBttn("runUMAP_intg", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg5).removeClass('disabled')"),
                                                      br(),
                                                      br(),
                                                      conditionalPanel(
                                                        condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                        plotlyOutput("UMAPplot_intg_tpm1", width = "100%"),
                                                        br(),
                                                        plotlyOutput("UMAPplot_intg_tpm2", width = "100%"),
                                                        br(),
                                                        plotlyOutput("UMAPplot_intg_tpm3", width = "100%"),
                                                        br(),
                                                        textOutput("UMAP_lisi_intg")),
                                                      br(),
                                                      conditionalPanel(
                                                        condition = "input.scInput1 == '10X cellranger'",
                                                        plotOutput("UMAPplot_intg_h5_1", width = "100%"),
                                                        br(),
                                                        plotOutput("UMAPplot_intg_h5_2", width = "100%")),
                                               ),
                                            ),
                                    
                                    tabPanel("tSNE", value = "test_intg4",
                                             
                                             tags$script(
                                               '
                                                  var tab_intg4 = $(\'a[data-value="test_intg4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
                                               column(3,
                                                      numericInput("dim.used_intg1",
                                                                   label = "Dimensions used",
                                                                   value = 10)
                                               ),
                                               column(9,
                                                      br(),
                                                      actionBttn("runTSNE_intg", "Run TSNE", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg5).removeClass('disabled')"),
                                                      #textOutput("Intg.done"),
                                                      br(),
                                                      br(),
                                                      conditionalPanel(
                                                        condition = "input.scInput1 == 'Raw Counts Matrix'",
                                                        plotlyOutput("TSNEplot_intg_tpm1", width = "100%"),
                                                        br(),
                                                        plotlyOutput("TSNEplot_intg_tpm2", width = "100%"),
                                                        br(),
                                                        plotlyOutput("TSNEplot_intg_tpm3", width = "100%")),
                                                      br(),
                                                      conditionalPanel(
                                                        condition = "input.scInput1 == '10X cellranger'",
                                                        plotlyOutput("TSNEplot_intg_h5_1", width = "100%"),
                                                        br(),
                                                        plotlyOutput("TSNEplot_intg_h5_2", width = "100%")),
                                               ),
                                             ),
                                     
                                     tabPanel("Clustering", value = "test_intg4",
                                              
                                              tags$script(
                                                '
                                                  var tab_intg4 = $(\'a[data-value="test_intg4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              fluidRow(
                                                column(3,
                                                       numericInput("clus.res_intg",
                                                                    label = "Resolution used",
                                                                    value = 0.6,
                                                                    min = 0.1,
                                                                    step = 0.1)
                                                ),
                                                column(3,
                                                       selectInput("dim.used_intg2",
                                                                   label = "PC to use", selected = "30",
                                                                   choices = c(1:50)),
                                                ),
                                                column(3,
                                                       br(),
                                                       actionBttn("findCluster_intg", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg5).removeClass('disabled')"),
                                                       textOutput("cluster1.done"),
                                                ),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              plotlyOutput("Cluster2DPlot_intg", width = "50%")
                                            ),
                                            column(12,
                                                   h4(p("Determine cluster resolution:")),
                                                   column(3,
                                                          numericInput("clus.res_a2",
                                                                       label = "Resolution (from)",
                                                                       value = 0.6,
                                                                       min = 0.1,
                                                                       step = 0.1)
                                                   ),
                                                   column(3,
                                                          numericInput("clus.res_b2",
                                                                       label = "Resolution (to)",
                                                                       value = 1,
                                                                       min = 0.1,
                                                                       step = 0.1)
                                                   ),
                                                   column(4,
                                                          br(),
                                                          actionBttn("findoptimumCluster2", "Determine optimum resolution", icon = icon("hand-pointer-o"), size = 'sm'),
                                                          #bsPopover("findoptimumCluster", "Determine optimum resolution","Press to determine optimum cluster resolution", placement = "bottom", trigger = "hover", options = NULL),
                                                   ),
                                                   column(12,
                                                   column(9,
                                                          plotOutput("OptimumCluster2DPlot_b1", width = "100%"),
                                                          br(),
                                                          ),
                                                  ),
                                                   column(12,
                                                          column(6,
                                                                 downloadBttn('download_OptimumCluster2', 'Download ClusterPlot (as png)', size = 'sm'),
                                                          ),
                                                          column(6,
                                                                 downloadBttn('download_OptimumClusterTable2', 'Download Cluster table (as csv)', size = 'sm'),
                                                          ),
                                                       ),
                                                  column(12,
                                                         titlePanel(h4(p("Subcluster Analysis"))),
                                                  ),
                                                  column(3,
                                                         br(),
                                                         uiOutput("subcluster.gene.select2"),
                                                  ),
                                                  column(3,
                                                         br(),
                                                         numericInput("subcluster.res2",
                                                                      label = "Resolution used",
                                                                      value = 0.6,
                                                                      min = 0.1,
                                                                      step = 0.1),
                                                         #bsTooltip("subcluster", "Perform Subcluster Analysis", placement = "bottom", trigger = "hover",
                                                         #options = NULL)
                                                  ),
                                                  
                                                  column(3,
                                                         br(),
                                                         br(),
                                                         actionBttn("subcluster2", "Run Subcluster analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                         #bsPopover("subcluster", "Subcluster Analysis","Press to perform Subcluster Analysis for a specific cluster at any given resolution", placement = "bottom", trigger = "hover", options = NULL),
                                                         textOutput("Subcluster2.done"),
                                                         br()
                                                  ),
                                                  br(),
                                                  br(),
                                                  br(),
                                                  column(9,
                                                         plotlyOutput("subcluster_plot2")
                                                  ),
                                                  column(12,
                                                         column(4,
                                                                downloadBttn('download_subcluster2', 'Download Subcluster plot (as png)', size = 'sm'),
                                                                br(),
                                                         ),
                                                  ),
                                                 ),
                                              ),
                        
                        tabPanel("Cell type identification",value = "test_intg5",
                                 
                                 tags$script(
                                   '
                                                  var tab_intg5 = $(\'a[data-value="test_intg5"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                 ),
                                 
                                 fluidRow(
                                   column(3,
                                          selectInput("cellid_method5",
                                                      label = "Celltype annotation method",
                                                      choices = c("CELLiD", "Celltypist"),
                                                      selected = "CELLiD"),
                                   )
                                 ),
                                      conditionalPanel(
                                          condition = "input.cellid_method5 == 'CELLiD'",
                                          column(3,
                                                 #br(),
                                            selectInput("cellatlas1a",
                                                      label = "Reference Atlas",
                                                      choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                      selected = "all")
                                            ),
                                 column(3,
                                        br(),
                                        actionBttn("doCELLiD_intg", "Run CELLiD", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg6).removeClass('disabled')"),
                                        textOutput("CELLiD2.done"),
                                        br()
                                 ),
                                 br(),
                                 br(),
                                 br(),
                                 plotlyOutput("Umap_cellid_intg", width = "50%"),
                                 plotlyOutput("Umap_cellid_intg1", width = "50%"),
                                 br(),
                                 br(),
                                 DT::dataTableOutput("ct_intg.table")
                                 ),
                                
                                
                                conditionalPanel(
                                  condition = "input.cellid_method5 == 'Celltypist'",
                                  column(3,
                                         #br(),
                                         selectInput("celltypistatlas5",
                                                     label = "Reference Atlas",
                                                     choices = c("Immune_All_Low.pkl", "Autopsy_COVID19_Lung.pkl", "Pan_Fetal_Human.pkl", "Nuclei_Lung_Airway.pkl", "Developing_Human_Thymus.pkl", "Human_Lung_Atlas.pkl", "Developing_Mouse_Brain.pkl", "Developing_Human_Brain.pkl", "Cells_Lung_Airway.pkl", "Healthy_COVID19_PBMC.pkl", "Human_IPF_Lung.pkl", "Adult_Mouse_Gut.pkl", "Immune_All_High.pkl", "COVID19_Immune_Landscape.pkl", "Human_PF_Lung.pkl", "COVID19_HumanChallenge_Blood.pkl", "Lethal_COVID19_Lung.pkl", "Cells_Fetal_Lung.pkl", "Cells_Intestinal_Tract.pkl"),
                                                     selected = "Immune_All_Low.pkl")
                                  ),
                                  column(3,
                                         br(),
                                         actionBttn("doCelltypist_intg", "Run Celltypist", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg6).removeClass('disabled')"),
                                         textOutput("Celltypist5.done"),
                                         br()
                                  ),
                                  br(),
                                  br(),
                                  br(),
                                  #br(),
                                  #br(),
                                  plotlyOutput("Umap_celltypist_intg", width = "50%"),
                                  plotlyOutput("Umap_celltypist_intg1", width = "50%"),
                                  br(),
                                  br(),
                                  DT::dataTableOutput("ct_celltypist_intg.table")
                                ),
                                
                                tags$head(tags$script('
                                                      Shiny.addCustomMessageHandler("myCallbackHandler6",
                                                      function(typeMessage) {console.log(typeMessage)
                                                      if(typeMessage == 2){
                                                      console.log("got here");
                                                      $("a:contains(Spatial Transcriptomics)").click();
                                                      }
                                                      });
                                                    ')),
                                column(3,
                                       br(),
                                       actionBttn("dodeconv_spatial_intg", "Go to Spatial Deconvolution", icon = icon("hand-pointer-o"), size = 'sm'),
                                       #bsPopover("dodeconv_spatial", "Perform Celltype Deconvolution","Press to navigate to Spatial module and perform celltype deconvolution", placement = "bottom", trigger = "hover", options = NULL),
                                ),
                                #br(),
                                
                                
                                tags$head(tags$script('
                                                       Shiny.addCustomMessageHandler("myCallbackHandler7",
                                                       function(typeMessage) {console.log(typeMessage)
                                                       if(typeMessage == 2){
                                                       console.log("got here");
                                                       $("a:contains(Single cell ATAC-seq)").click();
                                                       }
                                                      });
                                                    ')), 
                                column(3,
                                       br(),
                                       actionBttn("doct_atac_intg", "Annotate cell types for ATAC-data", icon = icon("hand-pointer-o"), size = 'sm'),
                                       #bsPopover("doct_atac", "Perform Celltype Label Transfer","Press to navigate to scATAC-seq module and perform celltype label transfer", placement = "bottom", trigger = "hover", options = NULL),
                                ),
                              ),
                        
                        tabPanel("Cell type similarity", value = "test_intg6",
                                 
                                 tags$script(
                                   '
                                                  var tab_intg6 = $(\'a[data-value="test_intg6"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                 ),
                                 
                                 column(3,
                                        br(),
                                        selectInput("cell1a",
                                                    label = "Group by",
                                                    choices = c("seurat_clusters", "primary.predict", "celltype"),
                                                    selected = "primary.predict"),
                                 ),
                                 
                                 column(3,
                                        br(),
                                        selectInput("corr_method1",
                                                    label = "Statistics",
                                                    choices = c("pearson", "spearman", "kendall"),
                                                    selected = "pearson"),
                                 ),   
                                 column(3,
                                        br(),
                                        br(),
                                        actionBttn("cell_cell1", "Run celltype similarity", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg7).removeClass('disabled')"),
                                        #bsPopover("cell_cell", "Perform Celltype Similarity","Press to run Celltype Similarity Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                        textOutput("CELL1.done"),
                                        br()
                                 ),
                                 br(),
                                 br(),
                                 br(),
                                 column(9,
                                        plotlyOutput("cell_cell_sim1")
                                 ),
                                 br(),
                                 br(),
                                 br(),
                                 column(12,
                                        column(4,
                                               downloadBttn('download_cell_cell_sim1', 'Download Celltype similarity plot (as png)', size = 'sm'),
                                        ),
                                        column(4,
                                               downloadBttn('download_cor.table1', 'Download Celltype similarity table (in csv)', size = 'sm'),
                                        ),
                                 ),
                                 br(),
                                 br(),
                                 br(),
                                 br(),
                                 br(),
                                 DT::dataTableOutput("cor.table1")
                        ),
                        
                        tabPanel("DEGs", value = "test_intg7",
                                 
                                 tags$script(
                                   '
                                                  var tab_intg7 = $(\'a[data-value="test_intg7"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                 ),
                                 
                                   fluidRow(
                                   column(3,
                                          selectInput("deg_method1",
                                                      label = "Type of DEG analysis",
                                                      choices = c("Celltype specific", "Pairwise DEGs"),
                                                      selected = "Celltype specific"),
                                          ),
                                   ),
                                   
                                   conditionalPanel(
                                     condition = "input.deg_method1 == 'Celltype specific'",
                                 fluidRow(
                                   column(3,
                                          selectInput("deg1a",
                                                      label = "Group by",
                                                      choices = c("seurat_clusters", "Predicted celltype", "Metadata celltype"),
                                                      selected = "Predicted celltype"),
                                          ),
                                 
                                   column(3, numericInput("min_pct_intg",
                                                          label = "min.pct",
                                                          value = 0.25,
                                                          min = 0,
                                                          step = 0.01)
                                          ),
                                   column(3, numericInput("logfc_intg",
                                                          label = "logfc.threshold",
                                                          value = 0.25,
                                                          min = 0,
                                                          step = 0.01)
                                         ),
                                       
                                   column(3, selectInput("test.use_intg",
                                                         label = "Test use",
                                                         choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                         ),
                                   column(4,
                                          actionBttn("doDeg_intg", "Run DEGs", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_intg8).removeClass('disabled')")
                                        )),
                                 br(),
                                 column(12,
                                        DT::dataTableOutput("Deg.table_intg")
                                 )),
                                    
                                 conditionalPanel(
                                   condition = "input.deg_method1 == 'Pairwise DEGs'",
                                   #br(),
                                   column(12,
                                          column(3,
                                                 selectInput("deg3a",
                                                             label = "Group by",
                                                             choices = c("seurat_clusters", "Predicted celltype", "Metadata celltype"),
                                                             selected = "primary.predict"),
                                          ),
                                          column(3,
                                                 uiOutput("gene1a.select"),
                                          ),
                                          column(3,
                                                 uiOutput("gene2a.select"),
                                          ),
                                          br(),
                                          column(3,
                                                 actionBttn("doVolcano1", "Run Pairwise DEGs", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg8).removeClass('disabled')"),
                                                 #bsPopover("doVolcano", "Perform Pairwise DEG Analysis","Press to run Pairwise DEG Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                          ),
                                   ),
                                   column(12,
                                          column(3, numericInput("min_pct_a1",
                                                                 label = "min.pct",
                                                                 value = 0.25,
                                                                 min = 0,
                                                                 step = 0.01)
                                          ),
                                          column(3, numericInput("logfc_a1",
                                                                 label = "logfc.threshold",
                                                                 value = 0.25,
                                                                 min = 0,
                                                                 step = 0.01)
                                          ),
                                          column(3, selectInput("test.use_a1",
                                                                label = "Test use",
                                                                choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                          )
                                   ),
                                   column(12,
                                          plotOutput("volcano.plot1", width = "75%"),
                                   ),
                                   br(),
                                   column(12,
                                          plotOutput("dega1.plot", width = "100%")
                                   ),
                                 ),
                               ),
                        
                        tabPanel("Data Visualization", value = "test_intg7",
                                 
                                 tags$script(
                                   '
                                                  var tab_intg7 = $(\'a[data-value="test_intg7"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                 ),
                                 
                                 column(3,
                                     selectInput("deg1b",
                                               label = "Group by",
                                               choices = c("seurat_clusters", "Predicted celltype", "Metadata celltype"),
                                               selected = "Predicted celltype"),
                                        ),
                                     column(6,
                                            uiOutput("deg.gene.select_intg"),
                                            actionBttn("Vis_intg", "Visualize", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_intg8).removeClass('disabled')"),
                                            plotlyOutput("Deg.plot_intg", width = "200%"),
                                            br(),
                                            plotlyOutput("Deg1.plot_intg", width = "100%"),
                                            br(),
                                            plotOutput("Deg2.plot_intg", width = "200%")
                                     ),
                                 column(4,
                                        downloadBttn('download_violn_intg', 'Download Violin plot (as png)', size = 'sm'),
                                        br(),
                                 ),
                                 column(4,
                                        downloadBttn('download_feature_intg', 'Download Feature plot (as png)', size = 'sm'),
                                        br(),
                                 ),
                                 column(4,
                                        downloadBttn('download_ridge_intg', 'Download Feature plot (as png)', size = 'sm'),
                                        br(),
                                 ),
                                  ),
                        
                        tabPanel("GSEA", value = "test_intg8",
                                 
                                 tags$script(
                                   '
                                                  var tab_intg8 = $(\'a[data-value="test_intg8"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                 ),
                                 
                                 column(12,
                                        br(),
                                        column(2, selectInput("species_gsea1",
                                                              label = "Species",
                                                              choices = c("Homo sapiens", "Mus musculus"))
                                        ),
                                        column(2, selectInput("category_gsea1",
                                                              label = "Collection",
                                                              choices = c("H", "C2", "C5", "C7", "C8"))
                                        ),
                                        column(2,
                                               uiOutput("gsea.ct1a.select"),
                                        ),
                                        column(2,
                                               uiOutput("gsea.ct2a.select"),
                                        ),
                                        
                                        
                                        column(2, numericInput("min_pct1a",
                                                               label = "min.pct",
                                                               value = 0.25,
                                                               min = 0,
                                                               step = 0.01)
                                        ),
                                        column(2, numericInput("logfc1a",
                                                               label = "logfc.threshold",
                                                               value = 0.25,
                                                               min = 0,
                                                               step = 0.01)
                                        ),
                                 ),
                                 column(12,
                                        column(2, selectInput("test.use1a",
                                                              label = "Test use",
                                                              choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                        ),
                                        column(4,
                                               uiOutput("gsea1.select"),
                                        ),
                                        br(),
                                        column(3,
                                               actionBttn("gsea1", "Run gene set enrichment analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                               #bsPopover("gsea", "Perform Gene Set Enrichment Analysis","Press to run Gene Set Enrichment Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                        ),
                                        textOutput("gsea1.done"),
                                        br()
                                 ),
                                 br(),
                                 column(7,
                                        plotOutput("gsea_plot1a", width = "100%"),
                                 ),
                                 br(),
                                 br(),
                                 column(12,
                                        plotOutput("gsea_plot1b", width = "100%"),
                                 ),
                                 br(),
                                 br(),
                                 DT::dataTableOutput("gsea.table1"),
                                 column(4,
                                        downloadBttn('download_gsea.table1', 'Download GSEA Results (in csv)', size = 'sm'),
                                 ),
                           ),
                        tabPanel("Cell-cell communication", value = "test_intg8",
                                 
                                 tags$script(
                                   '
                                                  var tab_intg8 = $(\'a[data-value="test_intg8"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_intg.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                 ),
                                 
                                 fluidRow(
                                   column(3,
                                          br(),
                                          selectInput("cc1b",
                                                      label = "Group by",
                                                      choices = c("seurat_clusters", "primary.predict"),
                                                      selected = "primary.predict"),
                                           ),
                                 
                                 
                                 column(3, 
                                        br(),
                                        selectInput("cc_method1",
                                                    label = "Method",
                                                    choices = c("natmi", "connectome", "logfc", "sca", "cellphonedb", "cytotalk"), selected = "cellphonedb")
                                 ),
                                 
                                 column(3, 
                                        br(),
                                        selectInput("cc_resource1",
                                                    label = "Method",
                                                    choices = c("Default", "Consensus",  "Baccin2019", "CellCall", "CellChatDB", "Cellinker", "CellPhoneDB", "CellTalkDB", "connectomeDB2020", "EMBRACE", "Guide2Pharma", "HPMR", "ICELLNET", "iTALK", "Kirouac2010", "LRdb", "Ramilowski2015", "OmniPath"), selected = "CellPhoneDB")
                                 ),
                                 column(3,
                                        br(),
                                        br(),
                                        actionBttn("doCC1", "Run analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                        #bsPopover("doCC", "Perform Cell-Cell Communication Analysis","Press to run Cell-Cell Communication Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                        textOutput("cc1.done"),
                                        br()
                                  ),
                                 ),
                                 
                                 column(12,
                                        plotOutput("CC_plot1a", width = "100%"),
                                 ),
                                 
                                 column(3,
                                        uiOutput("CC.gene.select1"),
                                 ),
                                 column(3,
                                        uiOutput("CC.gene1.select1"),
                                 ),
                                 br(),
                                 br(),
                                 column(12,
                                        plotOutput("CC_plot2a", width = "100%"),
                                 ),
                                 column(12,
                                 h4(p("Interacting partners")),
                                 ),
                                 withSpinner(dataTableOutput('cc.table1a')),
                                 
                          ),
                        #tabPanel("13. Annotate scATAC-seq data using scRNA-seq data",
                        #         tags$head(tags$script('
                        #                                Shiny.addCustomMessageHandler("myCallbackHandler1a",
                        #                                function(typeMessage) {console.log(typeMessage)
                        #                                if(typeMessage == 2){
                        #                                  console.log("got here");
                        #                                  $("a:contains(Single cell ATAC-seq)").click();
                        #                                }
                        #                              });
                        #                            ')), 
                                 
                        #         column(5,
                        #                br(),
                        #                actionBttn("loadexample_atacseq1", "Load example scATAC data and annotate", icon = icon("hand-pointer-o"), size = 'sm'),
                        #                #bsPopover("loadexample_atacseq", "Load example scATAC data and annotate","Press to load example scATAC data and annotate", placement = "bottom", trigger = "hover", options = NULL),
                        #         ),
                        #         column(4,
                        #                br(),
                        #                actionBttn("process_atacseq1", "Process scATAC-seq data", icon = icon("hand-pointer-o"), size = 'sm'),
                        #                #bsPopover("process_atacseq", "Process scATAC-seq data","Press to navigate to scATAC-seq module and process", placement = "bottom", trigger = "hover", options = NULL),
                        #         ),
                        #        br(),
                        #     column(12,
                        #             plotOutput("annotate_scRNA_ATAC_plot_a", width = "70%"),
                        #              plotOutput("annotate_scRNA_ATAC_plot_1a", width = "70%")
                        #        ),
                        #)
                        ),
                        
                       
                    ),
               
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
                                                              conditionalPanel(
                                                                titlePanel(h4(p("Load your example data"))),
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'Seurat'",
                                                                actionBttn("loadexample_seurat_spatial", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_sp).removeClass('disabled')"),
                                                              ),
                                                              conditionalPanel(
                                                                titlePanel(h4(p("Load your example data"))),
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'GraphST'",
                                                              actionBttn("loadexample_graphst_spatial", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_sp).removeClass('disabled')"),
                                                              ),
                                                              conditionalPanel(
                                                                titlePanel(h4(p("Load your example data"))),
                                                                condition = "input.scAnalysis_platform == 'Xenium' & input.scAnalysis_sp == 'Seurat'",
                                                                actionBttn("loadexample_xenium", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_sp).removeClass('disabled')"),
                                                              ),
                                                              titlePanel(h4(HTML("<b>Load your input data</b>"))),
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
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'Seurat'",
                                                                fileInput("tpmFile_spatial",
                                                                          label = "Upload H5 output from Spaceranger (Accepted Format: .h5)",
                                                                          accept = ".h5"),
                                                                ),
                                                              conditionalPanel(
                                                              condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'GraphST'",
                                                              fileInput("tpmFile_graphst",
                                                                        label = "Upload H5 output from Spaceranger (Accepted Format: .h5)",
                                                                        accept = ".h5"),
                                                                ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'Seurat'",
                                                                shinyFiles::shinyDirButton(id = 'dir', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                                                                verbatimTextOutput("dir", placeholder = TRUE),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'GraphST'",
                                                                shinyFiles::shinyDirButton(id = 'dir_graphst', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                                                                verbatimTextOutput("dir_graphst", placeholder = TRUE),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                                shinyFiles::shinyDirButton(id = 'dir_xenium', label = "Path to xenium files", title = "Sheets Folder Selector"),
                                                                verbatimTextOutput("dir_xenium", placeholder = TRUE),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'Seurat'",
                                                              textInput(inputId = "projName3",
                                                                        label = "Project Name",
                                                                        value = "Spatial_seurat"),
                                                              fluidRow(
                                                                actionBttn("loadButton3", "Load Data", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_sp).removeClass('disabled')"),
                                                                actionBttn("reset3", "Reset Data", icon = icon("repeat"), size = "sm"),
                                                              )),
                                                              
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'GraphST'",
                                                                textInput(inputId = "projName3a",
                                                                          label = "Project Name",
                                                                          value = "Spatial_graphST"),
                                                                fluidRow(
                                                                  actionBttn("loadButton3a", "Load Data", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_sp).removeClass('disabled')"),
                                                                  actionBttn("reset3a", "Reset Data", icon = icon("repeat"), size = "sm"),
                                                                )),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_platform == 'Xenium'",
                                                                #actionButton("loadexample_xenium", "Load example and run", icon = icon("hand-o-right")),
                                                                
                                                                textInput(inputId = "projName_xenium",
                                                                          label = "Project Name",
                                                                          value = "Xenium"),
                                                                fluidRow(
                                                                actionBttn("load_xenium", "Load Data", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp).removeClass('disabled')"),
                                                                actionBttn("reset_xenium", "Reset Data", icon = icon("repeat"), size = "sm"),
                                                                )),
                                                              )),
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_platform == 'Visium'",
                                                       chooseSliderSkin("Modern"),
                                                       titlePanel(h4(p("Quality control"))),
                                                     column(7,
                                                     plotOutput("nCount_SpatialPlot", width = "200%"),
                                                     ),
                                                     column(7,
                                                     plotOutput("SpatialFeaturePlot", width = "200%"),
                                                     column(4,
                                                            numericInput("obs",
                                                                         label = "nMolecules:",
                                                                         value = 200,
                                                                         min = 0,
                                                                         step = 1),
                                                                        ),
                                                     column(4,
                                                            numericInput("obs1",
                                                                         label = "nGenes:",
                                                                         value = 200,
                                                                         min = 0,
                                                                         step = 1),
                                                                         ),
                                                     column(3,
                                                            numericInput("obs2",
                                                                         label = "Mt%:",
                                                                         value = 30,
                                                                         min = 0,
                                                                         step = 1),
                                                            ),
                                                     column(6,
                                                            actionBttn("filter_spatial", "Filter", icon = icon("hand-pointer-o"), size = "sm"),
                                                            #actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                            
                                                     )
                                                     #column(12,
                                                     #       sliderInput("obs", "Number of molecules:", min = 0, max = 30000, value = 200),
                                                     #),
                                                     #column(12,
                                                     #       sliderInput("obs1", "Number of genes:", min = 0, max = 30000, value = 200),
                                                     #),
                                                     #column(12,
                                                     #       sliderInput("obs2", "Mitochondrial percentage:", min = 0, max = 100, value = 30),
                                                     #),
                                                     ),
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
                                                       condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'Seurat'",
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
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_platform == 'Visium' & input.scAnalysis_sp == 'GraphST'",
                                                       h4(p("Spatial counts")),
                                                       withSpinner(dataTableOutput('countdataDT_spatial_graphst')),
                                                       h4(p("H&E Image")),
                                                       withSpinner(plotOutput("h_e_graphst", width = "100%")),
                                                     ),
                                                    ),
                                     ),
                                     #tabPanel("2. Quality Control",
                                              #conditionalPanel(
                                                   # condition = "input.scAnalysis_platform == 'Xenium'",
                                                #column(12,
                                                     #column(5,
                                                     #       plotOutput("nFeature_xenium")
                                                     #),
                                                     # column(5,
                                                     #        plotOutput("nCount_xenium")
                                                     #),
                                                #)
                                                #),
                                            #),
                                     #tabPanel("3. FeatureScatter Plot",
                                     #         conditionalPanel(
                                     #           condition = "input.scAnalysis_platform == 'Visium'",
                                     #            column(6,
                                     #                  plotlyOutput("FeatureScatterPlot_sp")
                                     #           ),
                                     #         ),
                                     #         conditionalPanel(
                                     #           condition = "input.scAnalysis_platform == 'Xenium'",
                                     #           column(6,
                                     #                  plotlyOutput("FeatureScatterPlot_xenium")
                                     #           ),
                                     #         ),
                                     #),
                                     tabPanel("Normalization and Variable Feature Selection", value = "test_sp",
                                              
                                              tags$script(
                                                '
                                                  var tab_sp = $(\'a[data-value="test_sp"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
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
                                                         actionBttn("doSCTransform_sp", "Run scTransform", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp1).removeClass('disabled')"),
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
                                                         actionBttn("doSCTransform_xenium", "Run scTransform", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp1).removeClass('disabled')"),
                                                         #actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                         
                                                  )),
                                                plotOutput("VarGenes_xenium", width = "100%"),
                                              ),
                                            ),
                                     #tabPanel("Gene expression visualization", value = "test_sp1",
                                              
                                     #         tags$script(
                                     #          '
                                     #            var tab_sp1 = $(\'a[data-value="test_sp1"]\').parent().addClass("disabled");
                                     #            $(function(){
                                     #              $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                     #                e.preventDefault();
                                     #                return false;
                                     #              });
                                     #            });
                                     #            '
                                     #        ),
                                              
                                     #        conditionalPanel(
                                     #          condition = "input.scAnalysis_platform == 'Visium'",
                                     #         fluidRow(
                                     #            column(6,
                                     #                   uiOutput("sp.gene.select"),
                                     #                   plotOutput("sp.plot", width = "100%"),
                                     #            ),
                                     #          )),
                                     #        conditionalPanel(
                                     #          condition = "input.scAnalysis_platform == 'Xenium'",
                                     #        column(6,
                                     #               uiOutput("xenium.gene.select"),
                                     #               plotOutput("markerplot_xenium")
                                     #        ),
                                     #        column(6,
                                     #               uiOutput("xenium.gene1.select"),
                                     #               plotOutput("featureplot_xenium")
                                     #        ),
                                              
                                     #        br(),
                                     #        br(),
                                     #        column(2,
                                     #               numericInput("x1_xenium",
                                     #                            label = "x1",
                                     #                            value = 1200,
                                     #                            min = 0,
                                     #                            step = 1),
                                     #        ),
                                     #        column(2,
                                     #               numericInput("x2_xenium",
                                     #                            label = "x2",
                                     #                            value = 2800,
                                     #                            min = 0,
                                     #                            step = 1),
                                     #        ),
                                     #        column(2,
                                     #               numericInput("y1_xenium",
                                     #                            label = "y1",
                                     #                            value = 3700,
                                     #                            min = 0,
                                     #                            step = 1),
                                     #        ),
                                     #        column(2,
                                     #               numericInput("y2_xenium",
                                     #                            label = "y2",
                                     #                            value = 4500,
                                     #                            min = 0,
                                     #                            step = 1),
                                     #        ),
                                     #        column(2,
                                     #               br(),
                                     #        actionBttn("crop_xenium", "Crop", icon = icon("hand-o-right")),
                                     #               ),
                                     #        fluidRow(
                                     #          column(6,
                                     #                 uiOutput("xenium.gene2.select"),
                                     #                 plotOutput("zoom_xenium")
                                     #          ),
                                     #        ),
                                     #    ),
                                     #),
                                     tabPanel("PCA", value = "test_sp1",
                                              
                                              tags$script(
                                                '
                                                  var tab_sp1 = $(\'a[data-value="test_sp1"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                tabsetPanel(id="Pca_sp",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              actionBttn("runPCA_spatial", "Run PCA", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp2).removeClass('disabled')")
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
                                                                                actionBttn("runPCA_xenium", "Run PCA", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp2).removeClass('disabled')")
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
                                     
                                     
                                     tabPanel("UMAP",  value = "test_sp2",
                                              
                                              tags$script(
                                                '
                                                  var tab_sp2 = $(\'a[data-value="test_sp2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_platform == 'Visium'",
                                                column(3,
                                                       numericInput("dim.used_spatial",
                                                                    label = "Dimensions used",
                                                                    value = 10)
                                                ),
                                                column(3,
                                                       br(),
                                                       actionBttn("runUMAP_spatial", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp3).removeClass('disabled')")
                                                ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'Seurat' & input.scAnalysis_platform == 'Visium'",
                                                column(12,
                                                       br(),
                                                       br(),
                                                       plotlyOutput("DimPlot_spatial", width = "75%")
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp1 == 'GraphST' & input.scAnalysis_platform == 'Visium'",
                                                column(12,
                                                       br(),
                                                       br(),
                                                       plotlyOutput("DimPlot_GraphST", width = "75%")
                                                ),
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
                                                       actionBttn("runUMAP_xenium", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp3).removeClass('disabled')"),
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
                                     
                                     tabPanel("Clustering",  value = "test_sp2",
                                              
                                              tags$script(
                                                '
                                                  var tab_sp2 = $(\'a[data-value="test_sp2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
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
                                                           actionBttn("findCluster_spatial", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp3).removeClass('disabled')"),
                                                           #textOutput("cluster1.done"),
                                                    )),
                                                  br(),
                                                  plotlyOutput("Cluster2DPlot_spatial", width = "50%"),
                                                  br(),
                                                  #br(),
                                                  plotOutput("SpatialDimPlot", width = "100%"),
                                                  br(),
                                                  textOutput("silhouette_seurat")
                                                ),
                                                conditionalPanel(
                                                  condition = "input.scAnalysis_sp1 == 'GraphST'",
                                                  
                                                  #column(4,
                                                  #      fileInput("tpmFile_graphst",
                                                  #                label = "Spaceranger output (Accepted Format: .h5)",
                                                  #                accept = ".h5"),
                                                  #),
                                                  #column(3,
                                                  #       shinyFiles::shinyDirButton(id = 'dir_graphst', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                                                  #       verbatimTextOutput("dir_graphst", placeholder = TRUE),
                                                  #),
                                                  fluidRow(
                                                    column(3,
                                                           numericInput("cluster_graphst",
                                                                        label = "Number of clusters",
                                                                        value = 15,
                                                                        min = 2,
                                                                        step = 1)
                                                    ),
                                                    column(3,
                                                           selectInput("cluster_method_graphst",
                                                                       label = "Clustering method",
                                                                       choices = c("mclust", "leiden", "louvain"),
                                                                       selected = "mclust")
                                                    ),
                                                    column(2,
                                                           br(),
                                                           actionBttn("findCluster_graphst", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp3).removeClass('disabled')"),
                                                    ),
                                                  ),
                                                  br(),
                                                  column(12,
                                                         plotOutput("Cluster2DPlot_graphst", width = "100%")
                                                  ),
                                                  br(),
                                                  plotOutput("SpatialDimPlot_GraphST", width = "100%"),
                                                  br(),
                                                  textOutput("silhouette_graphst")
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
                                                         actionBttn("findCluster_xenium", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp3).removeClass('disabled')"),
                                                         br()
                                                  )),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_xenium", width = "100%")
                                              ),
                                     ),
                                     
                                     tabPanel("DEGs", value = "test_sp3",
                                              
                                              tags$script(
                                                '
                                                  var tab_sp3 = $(\'a[data-value="test_sp3"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp2 == 'Seurat'",
                                                fluidRow(
                                                  column(3,
                                                         selectInput("deg_method_sp",
                                                                     label = "Type of DEG analysis",
                                                                     choices = c("Celltype specific", "Pairwise DEGs"),
                                                                     selected = "Celltype specific"),
                                                  ),
                                                ),
                                                conditionalPanel(
                                                  condition = "input.deg_method_sp == 'Celltype specific'",
                                                  fluidRow(
                                                    column(3,
                                                           selectInput("deg_spatial",
                                                                       label = "Group by",
                                                                       choices = c("seurat_clusters"),
                                                                       selected = "seurat_clusters"),
                                                    ),
                                                    column(3, numericInput("min_pct_spatial",
                                                                           label = "min.pct",
                                                                           value = 0.75,
                                                                           min = 0,
                                                                           step = 0.01)
                                                    ),
                                                    column(3, numericInput("logfc_spatial",
                                                                           label = "logfc.threshold",
                                                                           value = 0.75,
                                                                           min = 0,
                                                                           step = 0.01)
                                                    ),
                                                    column(3, selectInput("test.use_spatial",
                                                                          label = "Test use",
                                                                          choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                    ),
                                                    column(4,
                                                           actionBttn("doDeg_spatial", "Run DEGs", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp4).removeClass('disabled')")
                                                    )),
                                                  br(),
                                                  column(12,
                                                         DT::dataTableOutput("Deg.table_spatial")
                                                  ),
                                                  br(),
                                                  br(),
                                                  plotOutput("Deg.heatmap_spatial", width = "100%"),
                                                ),
                                                
                                                conditionalPanel(
                                                  condition = "input.deg_method_sp == 'Pairwise DEGs'",
                                                  #br(),
                                                  column(12,
                                                         column(3,
                                                                selectInput("deg_sp1",
                                                                            label = "Group by",
                                                                            choices = c("seurat_clusters"),
                                                                            selected = "seurat_clusters"),
                                                         ),
                                                         column(3,
                                                                uiOutput("gene_sp.select"),
                                                         ),
                                                         column(3,
                                                                uiOutput("gene_sp1.select"),
                                                         ),
                                                         br(),
                                                         column(3,
                                                                actionBttn("doVolcano_sp", "Run Pairwise DEGs", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_sp4).removeClass('disabled')"),
                                                                #bsPopover("doVolcano", "Perform Pairwise DEG Analysis","Press to run Pairwise DEG Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                         ),
                                                  ),
                                                  column(12,
                                                         column(3, numericInput("min_pct_sp1",
                                                                                label = "min.pct",
                                                                                value = 0.25,
                                                                                min = 0,
                                                                                step = 0.01)
                                                         ),
                                                         column(3, numericInput("logfc_sp1",
                                                                                label = "logfc.threshold",
                                                                                value = 0.25,
                                                                                min = 0,
                                                                                step = 0.01)
                                                         ),
                                                         column(3, selectInput("test.use_sp1",
                                                                               label = "Test use",
                                                                               choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                         )
                                                  ),
                                                  column(12,
                                                         plotOutput("volcano.plot_sp1", width = "75%"),
                                                  ),
                                                  br(),
                                                  column(12,
                                                         plotOutput("dega_sp.plot", width = "100%")
                                                  ),
                                                ),
                                              ),
                                           ),
                                     
                                     tabPanel("Deconvolution", value = "test_sp4",
                                              
                                              tags$script(
                                                '
                                                  var tab_sp4 = $(\'a[data-value="test_sp4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                             
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp2 == 'Seurat' | input.scAnalysis_sp2 == 'GraphST'",
                                                
                                                column(4,
                                                       br(),
                                                       actionBttn("loadexample_ref", "Load example reference dataset", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_sp5).removeClass('disabled')"),
                                                ),
                                                tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler2",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 1){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell RNA-Sequencing)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                                column(4,
                                                       br(),
                                                       actionBttn("loaduser_ref", "Load and process user reference dataset", icon = icon("hand-o-right"), size = "sm"),
                                                ),
                                                column(12,
                                                       #br(),
                                                       #actionButton("process_scRNA", "Process reference dataset", icon = icon("hand-pointer-o")),
                                                       br(),
                                                       br(),
                                                       plotOutput("scRNAPlot", width = "75%")),
                                                column(9,
                                                       br(),
                                                       actionBttn("vis_spRNA", "Visualise spatial dataset", icon = icon("hand-pointer-o"), size = "sm"),
                                                       br(),
                                                       br(),
                                                       plotOutput("spRNAPlot", width = "100%")),
                                                ),
                                              column(4,
                                                     br(),
                                                     selectInput("scAnalysis_sp2",
                                                                 label = "Analysis method",
                                                                 choices = c("Seurat", "GraphST"),
                                                                 selected = "Seurat"),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp2 == 'Seurat'",
                                                column(4,
                                                       br(),
                                                       br(),
                                                       actionBttn("doDeconv", "Deconvolute", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp5).removeClass('disabled')"),
                                                ),
                                                       br(),
                                                column(9,
                                                       uiOutput("ct.select"),
                                                       plotOutput("DeconvPlot", width = "100%")
                                                ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_sp2 == 'GraphST'",
                                               
                                                column(4,
                                                       br(),
                                                       br(),
                                                       actionBttn("doDeconv_graphst", "Deconvolute", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_sp5).removeClass('disabled')"),
                                                      ),
                                                       br(),
                                                column(9,
                                                       uiOutput("graphst_ct.select"),
                                                       plotOutput("DeconvPlot_graphst", width = "100%")
                                                       ),
                                                      )),
                                            
                                                 tabPanel("Data Visualization", value = "test_sp5",
                                                          
                                                  tags$script(
                                                            '
                                                  var tab_sp5 = $(\'a[data-value="test_sp5"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                                  ),
                                                          
                                                          tabsetPanel(id="vis_sp",
                                                            tabPanel(title="Visualize genes", value="Vis_sp_panel1",
                                                              fluidRow(
                                                                column(3,
                                                                       br(),
                                                                       selectInput("deg_sp",
                                                                               label = "Group by",
                                                                               choices = c("seurat_clusters"),
                                                                               selected = "seurat_clusters"),
                                                                      ),
                                                            column(6,
                                                                   br(),
                                                                   uiOutput("deg.gene.select_spatial"),
                                                                   actionBttn("Vis_spatial", "Visualize", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_sp5).removeClass('disabled')"),
                                                                   plotlyOutput("Deg.plot_spatial", width = "150%"),
                                                                   br(),
                                                                   plotOutput("Deg1.plot_spatial", width = "150%"),
                                                                   br(),
                                                                   plotOutput("Deg2.plot_spatial", width = "150%"),
                                                            ),
                                                            column(4,
                                                                   downloadBttn('download_violn_spatial', 'Download Violin plot (as png)', size = 'sm'),
                                                                   br(),
                                                            ),
                                                            column(4,
                                                                   downloadBttn('download_feature_spatial', 'Download Feature plot (as png)', size = 'sm'),
                                                                   br(),
                                                            ),
                                                            column(4,
                                                                   downloadBttn('download_ridge_spatial', 'Download Feature plot (as png)', size = 'sm'),
                                                                   br(),
                                                            ),
                                                          )
                                                        ),
                                                        tabPanel(title="Visualize Spatial Domains", value="Vis_sp_panel2",
                                                                 conditionalPanel(
                                                                   condition = "input.scAnalysis_sp1 == 'Seurat' & input.scAnalysis_platform == 'Visium'",
                                                                   
                                                                   fluidRow(
                                                                     column(6,
                                                                            br(),
                                                                            uiOutput("sp.cluster.select"),
                                                                            plotOutput("sp1.plot", width = "100%"),
                                                                           ),
                                                                         )),
                                                                 
                                                                 conditionalPanel(
                                                                   condition = "input.scAnalysis_sp1 == 'GraphST' & input.scAnalysis_platform == 'Visium'",
                                                                   
                                                                   fluidRow(
                                                                     column(6,
                                                                            br(),
                                                                            uiOutput("sp.cluster1.select"),
                                                                            plotOutput("sp2.plot", width = "100%"),
                                                                           ),
                                                                         )),
                                                                 
                                                                 conditionalPanel(
                                                                   condition = "input.scAnalysis_platform == 'Xenium'",
                                                                   column(4,
                                                                          actionBttn("Vis_sp_xenium", "Visualize", icon = icon("hand-pointer-o"), size = "sm")
                                                                          ),
                                                                   
                                                                   fluidRow(
                                                                     column(6,
                                                                            uiOutput("xenium.cluster.select"),
                                                                            plotOutput("xenium1.plot", width = "100%"),
                                                                             ),
                                                                           )),
                                                                        ),
                                                                      ) 
                                                                    ),
                                                                 tabPanel("GSEA", value = "test_sp5",
                                                                          
                                                                      tags$script(
                                                                                                      '
                                                                            var tab_sp5 = $(\'a[data-value="test_sp5"]\').parent().addClass("disabled");
                                                                            $(function(){
                                                                              $(tab_sp.parent()).on("click", "li.disabled", function(e) {
                                                                                e.preventDefault();
                                                                                return false;
                                                                              });
                                                                            });
                                                                            '
                                                                        ),
                                                                          
                                                                          column(12,
                                                                                 br(),
                                                                                 column(2, selectInput("species_gsea_sp",
                                                                                                       label = "Species",
                                                                                                       choices = c("Homo sapiens", "Mus musculus"))
                                                                                 ),
                                                                                 column(2, selectInput("category_gsea_sp",
                                                                                                       label = "Collection",
                                                                                                       choices = c("H", "C2", "C5", "C7", "C8"))
                                                                                 ),
                                                                                 column(2,
                                                                                        uiOutput("gsea.ct_sp1.select"),
                                                                                 ),
                                                                                 column(2,
                                                                                        uiOutput("gsea.ct_sp2.select"),
                                                                                 ),
                                                                                 
                                                                                 
                                                                                 column(2, numericInput("min_pct_sp1",
                                                                                                        label = "min.pct",
                                                                                                        value = 0.25,
                                                                                                        min = 0,
                                                                                                        step = 0.01)
                                                                                 ),
                                                                                 column(2, numericInput("logfc_sp1",
                                                                                                        label = "logfc.threshold",
                                                                                                        value = 0.25,
                                                                                                        min = 0,
                                                                                                        step = 0.01)
                                                                                 ),
                                                                          ),
                                                                          column(12,
                                                                                 column(2, selectInput("test.use_sp1",
                                                                                                       label = "Test use",
                                                                                                       choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                                                 ),
                                                                                 column(4,
                                                                                        uiOutput("gsea_sp.select"),
                                                                                 ),
                                                                                 br(),
                                                                                 column(3,
                                                                                        actionBttn("gsea_sp", "Run gene set enrichment analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                                                        #bsPopover("gsea", "Perform Gene Set Enrichment Analysis","Press to run Gene Set Enrichment Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                                                 ),
                                                                                 textOutput("gsea_sp.done"),
                                                                                 br()
                                                                          ),
                                                                          br(),
                                                                          column(7,
                                                                                 plotOutput("gsea_sp_plot", width = "100%"),
                                                                          ),
                                                                          br(),
                                                                          br(),
                                                                          column(12,
                                                                                 plotOutput("gsea_sp_plot1", width = "100%"),
                                                                          ),
                                                                          br(),
                                                                          br(),
                                                                          DT::dataTableOutput("gsea_sp.table"),
                                                                          column(4,
                                                                                 downloadBttn('download_gsea_sp.table', 'Download GSEA Results (in csv)', size = 'sm'),
                                                                          ),
                                                                        ),
                                                                      ),
                       
               ),
               
               ##---------------Single cell multiomics module-----------##
               
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
                                                              
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'CITE-seq'",
                                                                titlePanel(h4(p("Load example data"))),
                                                                actionBttn("loadexample_cite", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_mult).removeClass('disabled')"),
                                                              ),
                                                              conditionalPanel(
                                                                condition = "input.scAnalysis_type == 'Multiome'",
                                                                titlePanel(h4(p("Load example data"))),
                                                                actionBttn("loadexample_multiome", "Load example and run", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_mult).removeClass('disabled')"),
                                                              ),
                                                              
                                                              titlePanel(h4(HTML("<b>Load your input data</b>"))),
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
                                                              
                                                              selectInput("scAnalysis_type",
                                                                          label = "scAnalysis type",
                                                                          choices = c("CITE-seq", "Multiome"),
                                                                          selected = "H5"),
                                                              selectInput("scAnalysis_mult",
                                                                          label = "Analysis method",
                                                                          choices = c("Seurat", "MOFA+"),
                                                                          selected = "Seurat"),
                                                              column(6,
                                                                     numericInput(inputId = "min.genes2",
                                                                                  label = "Min. genes",
                                                                                  value = 200,
                                                                                  min = 1)
                                                              ),
                                                              column(6,
                                                                     numericInput(inputId = "min.cells2",
                                                                                  label = "Min. cells",
                                                                                  value = 3,
                                                                                  min = 1)
                                                              ),
                                                              textInput(inputId = "projName2",
                                                                        label = "Project Name",
                                                                        value = "Multiomics"),
                                                              fluidRow(
                                                                actionBttn("loadButton2", "Load data", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_mult).removeClass('disabled')"),
                                                                actionBttn("reset_mult", "Reset", icon = icon("repeat"), size = "sm")
                                                              ),
                                                            )),
                                                     
                                                     chooseSliderSkin("Modern"),
                                                     titlePanel(h4(p("Quality control"))),
                                                     conditionalPanel(
                                                         condition = "input.scAnalysis_type == 'CITE-seq'",
                                                     column(6,
                                                            h4(p("RNA QC")),
                                                            plotOutput("nFeature_RNAPlot2a", width = "200%")
                                                            ),
                                                     column(6,
                                                            h4(p("ADT QC")),
                                                            plotOutput("nFeature_RNAPlot2c", width = "150%")
                                                      ),
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'Multiome'",
                                                       column(6,
                                                              h4(p("RNA QC")),
                                                              plotOutput("nFeature_RNAPlot2b", width = "200%")
                                                       ),
                                                       column(6,
                                                              h4(p("ATAC QC")),
                                                              plotOutput("nFeature_RNAPlot2d", width = "150%")
                                                       ),
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'CITE-seq'",
                                                       column(4,
                                                              numericInput("obsa",
                                                                           label = "Min nFeature - RNA:",
                                                                           value = 200,
                                                                           min = 0,
                                                                           step = 1),
                                                       ),
                                                       column(4,
                                                              numericInput("obsa1",
                                                                           label = "Max nFeature - RNA:",
                                                                           value = 2500,
                                                                           min = 0,
                                                                           step = 1),
                                                       ),
                                                       column(3,
                                                              numericInput("obsa2",
                                                                           label = "Mt%:",
                                                                           value = 5,
                                                                           min = 0,
                                                                           step = 1),
                                                       ),
                                                       column(4,
                                                              numericInput("obsa3",
                                                                           label = "Min nCount - ADT:",
                                                                           value = 200,
                                                                           min = 0,
                                                                           step = 1),
                                                       ),
                                                       column(4,
                                                              numericInput("obsa4",
                                                                           label = "Max nCount - ADT:",
                                                                           value = 50000,
                                                                           min = 0,
                                                                           step = 1),
                                                       ),
                                                     ), 
                                                     
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'Multiome'",
                                                     column(4,
                                                            numericInput("obss",
                                                                         label = "Min nFeature - RNA:",
                                                                         value = 200,
                                                                         min = 0,
                                                                         step = 1),
                                                     ),
                                                     column(4,
                                                            numericInput("obss1",
                                                                         label = "Max nFeature - RNA:",
                                                                         value = 2500,
                                                                         min = 0,
                                                                         step = 1),
                                                     ),
                                                     column(3,
                                                            numericInput("obss2",
                                                                         label = "Mt%:",
                                                                         value = 5,
                                                                         min = 0,
                                                                         step = 1),
                                                     ),
                                                     column(4,
                                                            numericInput("obss3",
                                                                         label = "Min nCount - ATAC:",
                                                                         value = 200,
                                                                         min = 0,
                                                                         step = 1),
                                                     ),
                                                     column(4,
                                                            numericInput("obss4",
                                                                         label = "Max nCount - ATAC:",
                                                                         value = 50000,
                                                                         min = 0,
                                                                         step = 1),
                                                     ),
                                                    ),
                                                    
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'CITE-seq'",
                                                       column(3,
                                                              br(),
                                                              actionBttn("filter_seurat2a", "Filter", icon = icon("hand-o-right"), size = "sm")
                                                       )),
                                                     conditionalPanel(
                                                       condition = "input.scAnalysis_type == 'Multiome'",
                                                       column(3,
                                                              br(),
                                                              actionBttn("filter_seurat2b", "Filter", icon = icon("hand-o-right"), size = "sm")
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
                                     
                                     
                                     #tabPanel("2. Violin Plot",
                                     #        conditionalPanel(
                                     #          condition = "input.scAnalysis_type == 'CITE-seq'",
                                     #          column(3,
                                     #                 plotOutput("nFeature_RNAPlot2a")
                                     #          ),
                                     #          column(3,
                                     #                 plotOutput("mitoPlot2a")
                                     #          ),
                                     #          column(3,
                                     #                 plotOutput("nCount_RNAPlot2a")
                                     #          )),
                                     #        
                                     #        conditionalPanel(
                                     #          condition = "input.scAnalysis_type == 'Multiome'",
                                     #          column(3,
                                     #                 plotOutput("nFeature_RNAPlot2b")
                                     #          ),
                                     #          column(3,
                                     #                 plotOutput("mitoPlot2b")
                                     #          ),
                                     #          column(3,
                                     #                 plotOutput("nCount_RNAPlot2b")
                                     #          )),
                                     #),
                                     
                                     #tabPanel("Feature Scatter Plot",
                                     #          conditionalPanel(
                                     #          condition = "input.scAnalysis_type == 'CITE-seq'",
                                     #          column(5,
                                     #                 plotlyOutput("FeatureScatterPlot_mult_a")
                                     #          ),
                                     #          column(5,
                                     #                 plotlyOutput("FeatureScatterPlot_mult_b")
                                     #          )),
                                     #        conditionalPanel(
                                     #          condition = "input.scAnalysis_type == 'Multiome'",
                                     #          column(5,
                                     #                 plotlyOutput("FeatureScatterPlot_mult_c")
                                     #          ),
                                     #          column(5,
                                     #                 plotlyOutput("FeatureScatterPlot_mult_d")
                                     #          )
                                     #        ),
                                     #),
                                     
                                     tabPanel("Normalization and Variable Feature Selection", value = "test_mult",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult = $(\'a[data-value="test_mult"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_type == 'CITE-seq'",
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
                                                         actionBttn("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult1).removeClass('disabled')"),
                                                  )),
                                                plotOutput("VarGenes_mult", width = "100%")),
                                              
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
                                                         actionBttn("doSCTransform_multi", "Run scTransform", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult1).removeClass('disabled')"),
                                                         #actionButton("findVarGenes_mult", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                                         
                                                  )),
                                                plotOutput("VarGenes_mult1", width = "100%")),
                                     ),
                                     
                                     tabPanel("PCA", value = "test_mult1",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult1 = $(\'a[data-value="test_mult1"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              br(),
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                tabsetPanel(id="Pca_mult_seurat",
                                                            tabPanel(title="PCA Plot",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              actionBttn("runPCA_mult_seurat", "Run PCA", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult2).removeClass('disabled')")
                                                                       ),
                                                                     ),
                                                                     br(),
                                                                     plotlyOutput("PCAplot_mult_seurat_h5_1", width = "50%"),
                                                                     
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
                                                                              actionBttn("runPCA_mult_seurat1", "Run PCA", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult2).removeClass('disabled')")
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
                                                condition = "input.scAnalysis_mult == 'MOFA+' & input.scAnalysis_type == 'CITE-seq'",
                                                tabsetPanel(id="Pca_mult_mofa2",
                                                            tabPanel(title="Data overview",
                                                                     br(),
                                                                     fluidRow(
                                                                       column(3,
                                                                              selectInput("num_factor_mofa",
                                                                                          label = "Number of factors",
                                                                                          choices = c(1:15),
                                                                                          selected = "10"
                                                                              )),
                                                                       column(3,
                                                                              br(),
                                                                              actionBttn("runMOFA_mult", "Run MOFA+", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult2).removeClass('disabled')")
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
                                     
                                     tabPanel("UMAP", value = "test_mult2",
                                     
                                              tags$script(
                                                '
                                                  var tab_mult2 = $(\'a[data-value="test_mult2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
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
                                                       actionBttn("runUMAP_mult_seurat", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult3).removeClass('disabled')"),
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
                                                column(3,
                                                       br(),
                                                       actionBttn("runUMAP_mult_seurat1", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult3).removeClass('disabled')"),
                                                ),
                                                column(9,
                                                       br(),
                                                       br(),
                                                       plotlyOutput("UMAPplot1_a", width = "100%"),
                                                       br(),
                                                       plotlyOutput("UMAPplot1_b", width = "100%"),
                                                       br(),
                                                       plotlyOutput("UMAPplot1_c", width = "100%"),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA+' & input.scAnalysis_type == 'CITE-seq'",
                                                br(),
                                                column(3,
                                                       numericInput("num.neighbors_mofa2",
                                                                    label = "Number of neighbors",
                                                                    value = 10)
                                                ),
                                                column(3,
                                                       br(),
                                                       actionBttn("runUMAP_mofa2", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult3).removeClass('disabled')"),
                                                ),
                                                column(9,
                                                       br(),
                                                       br(),
                                                       plotOutput("UMAPplot_mofa_1", width = "100%"),
                                                )),
                                     ),
                                     
                                     tabPanel("Clustering", value = "test_mult3",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult3 = $(\'a[data-value="test_mult3"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
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
                                                         actionBttn("findCluster_mult_seurat", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult4).removeClass('disabled')"),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_mult_seurat_a", width = "75%"),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_mult_seurat_b", width = "75%"),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_mult_seurat_c", width = "75%")
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
                                                         actionBttn("findCluster_mult_seurat1", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult4).removeClass('disabled')"),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotlyOutput("Cluster2DPlot_mult_seurat1", width = "100%")
                                              ),  
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA+' & input.scAnalysis_type == 'CITE-seq'",
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
                                                         actionBttn("findCluster_mofa2", "Find Clusters", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult4).removeClass('disabled')"),
                                                         #textOutput("cluster1.done"),
                                                  )),
                                                br(),
                                                plotOutput("Cluster2DPlot_mofa", width = "50%")
                                              ),  
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
                                     
                                     tabPanel("Cell type identification", value = "test_mult4",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult4 = $(\'a[data-value="test_mult4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                column(3,
                                                       selectInput("cellid_method9",
                                                                   label = "Celltype annotation method",
                                                                   choices = c("CELLiD", "Celltypist"),
                                                                   selected = "CELLiD"),
                                                ),
                                                conditionalPanel(
                                                  condition = "input.cellid_method9 == 'CELLiD'",
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
                                                       actionBttn("doCELLiD_mult_cite_seurat", "Run CELLiD", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult5).removeClass('disabled')"),
                                                       #textOutput("CELLiD.done"),
                                                       br()
                                                ),
                                                br(),
                                                br(),
                                                br(),
                                                plotOutput("Umap_cellid_mult_cite_seurat", width = "50%"),
                                                br(),
                                                plotOutput("Umap_cellid_mult_cite_seurat1", width = "50%"),
                                                br(),
                                                DT::dataTableOutput("ct_cite_seurat.table"),
                                                column(4,
                                                       downloadButton('download_cellid_cite_seurat_prediction', 'Download CELLiD predictions (in csv)'),
                                                )),
                                                conditionalPanel(
                                                  condition = "input.cellid_method9 == 'Celltypist'",
                                                  column(3,
                                                         #br(),
                                                         selectInput("celltypistatlas9",
                                                                     label = "Reference Atlas",
                                                                     choices = c("Immune_All_Low.pkl", "Autopsy_COVID19_Lung.pkl", "Pan_Fetal_Human.pkl", "Nuclei_Lung_Airway.pkl", "Developing_Human_Thymus.pkl", "Human_Lung_Atlas.pkl", "Developing_Mouse_Brain.pkl", "Developing_Human_Brain.pkl", "Cells_Lung_Airway.pkl", "Healthy_COVID19_PBMC.pkl", "Human_IPF_Lung.pkl", "Adult_Mouse_Gut.pkl", "Immune_All_High.pkl", "COVID19_Immune_Landscape.pkl", "Human_PF_Lung.pkl", "COVID19_HumanChallenge_Blood.pkl", "Lethal_COVID19_Lung.pkl", "Cells_Fetal_Lung.pkl", "Cells_Intestinal_Tract.pkl"),
                                                                     selected = "Immune_All_Low.pkl")
                                                  ),
                                                  column(3,
                                                         selectInput("assay_mult1_cite_seurat",
                                                                     label = "Assay",
                                                                     choices = c("rna.umap", "adt.umap", "wnn.umap"),
                                                                     selected = "all")
                                                  ),
                                                  column(3,
                                                         br(),
                                                         actionBttn("doCelltypist_mult_cite_seurat", "Run Celltypist", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult5).removeClass('disabled')"),
                                                         textOutput("Celltypist9.done"),
                                                         br()
                                                  ),
                                                  br(),
                                                  br(),
                                                  br(),
                                                  #br(),
                                                  #br(),
                                                  plotlyOutput("Umap_celltypist_mult_cite_seurat", width = "50%"),
                                                  br(),
                                                  plotlyOutput("Umap_celltypist_mult_cite_seurat1", width = "50%"),
                                                  br(),
                                                  br(),
                                                  DT::dataTableOutput("ct_celltypist_mult_cite_seurat.table"),
                                                  column(4,
                                                         downloadButton('download_celltypist_cite_seurat_prediction', 'Download CELLiD predictions (in csv)'),
                                                  ),
                                                ),
                                                
                                                ),
                                              
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
                                                       actionBttn("doCELLiD_multiome_seurat", "Run CELLiD", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult5).removeClass('disabled')"),
                                                       #textOutput("CELLiD.done"),
                                                       br()
                                                ),
                                                br(),
                                                br(),
                                                br(),
                                                plotlyOutput("Umap_cellid_multiome_seurat", width = "50%"),
                                                br(),
                                                plotlyOutput("Umap_cellid_multiome_seurat1", width = "50%"),
                                                br(),
                                                DT::dataTableOutput("ct_multiome_seurat.table"),
                                                column(4,
                                                       downloadButton('download_cellid_multiome_seurat_prediction', 'Download CELLiD predictions (in csv)'),
                                                )),
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA+' & input.scAnalysis_type == 'CITE-seq'",
                                                column(3,
                                                       selectInput("cellatlas_mult_cite_mofa",
                                                                   label = "Reference Atlas",
                                                                   choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                   selected = "all")
                                                ),
                                                column(3,
                                                       br(),
                                                       actionBttn("doCELLiD_mult_cite_mofa", "Run CELLiD", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_mult5).removeClass('disabled')"),
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
                                     
                                     tabPanel("Cell type Similarity", value = "test_mult5",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult5 = $(\'a[data-value="test_mult5"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              column(3,
                                                     br(),
                                                     selectInput("cell_m1",
                                                                 label = "Group by",
                                                                 choices = c("seurat_clusters", "primary.predict"),
                                                                 selected = "primary.predict"),
                                              ),
                                              
                                              column(3,
                                                     br(),
                                                     selectInput("corr_method2",
                                                                 label = "Statistics",
                                                                 choices = c("pearson", "spearman", "kendall"),
                                                                 selected = "pearson"),
                                              ),  
                                              column(3,
                                                     br(),
                                                     br(),
                                                     actionBttn("cell_cell_mult", "Run celltype similarity", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_mult6).removeClass('disabled')"),
                                                     #bsPopover("cell_cell", "Perform Celltype Similarity","Press to run Celltype Similarity Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                     textOutput("CELL_mult.done"),
                                                     br()
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              column(9,
                                                     plotlyOutput("cell_cell_mult_sim")
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              column(12,
                                                     column(4,
                                                            downloadBttn('download_cell_cell_mult_sim', 'Download Celltype similarity plot (as png)', size = 'sm'),
                                                     ),
                                                     column(4,
                                                            downloadBttn('download_cor_mult.table', 'Download Celltype similarity table (in csv)', size = 'sm'),
                                                     ),
                                              ),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              br(),
                                              DT::dataTableOutput("cor_mult.table")
                                     ),
                                     
                                     tabPanel("DEGs", value = "test_mult6",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult6 = $(\'a[data-value="test_mult6"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              fluidRow(
                                                column(3,
                                                       selectInput("deg_method_mult",
                                                                   label = "Type of DEG analysis",
                                                                   choices = c("Celltype specific", "Pairwise DEGs"),
                                                                   selected = "Celltype specific"),
                                                ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.deg_method_mult == 'Celltype specific'",
                                                fluidRow(
                                                  column(3,
                                                         selectInput("deg_mult1",
                                                                     label = "Group by",
                                                                     choices = c("seurat_clusters", "primary.predict", "newID"),
                                                                     selected = "primary.predict"),
                                                  ),
                                                  column(3, numericInput("min_pct_mult",
                                                                         label = "min.pct",
                                                                         value = 0.25,
                                                                         min = 0,
                                                                         step = 0.01)
                                                  ),
                                                  column(3, numericInput("logfc_mult",
                                                                         label = "logfc.threshold",
                                                                         value = 0.25,
                                                                         min = 0,
                                                                         step = 0.01)
                                                  ),
                                                  column(3, selectInput("test.use_mult",
                                                                        label = "Test use",
                                                                        choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                  ),
                                                  br(),
                                                  column(3,
                                                         actionBttn("doDeg_mult", "Run DEGs", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_mult7).removeClass('disabled')"),
                                                         #bsPopover("doDeg", "Perform DEG Analysis","Press to run DEG Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                  )),
                                                br(),
                                                column(12,
                                                       DT::dataTableOutput("Deg_mult.table"),
                                                       br(),
                                                       br(),
                                                       br(),
                                                       plotOutput("Deg3_mult.plot", width = "100%")
                                                )
                                              ),
                                              conditionalPanel(
                                                condition = "input.deg_method_mult == 'Pairwise DEGs'",
                                                #br(),
                                                column(12,
                                                       column(3,
                                                              selectInput("deg3_mult",
                                                                          label = "Group by",
                                                                          choices = c("seurat_clusters", "primary.predict"),
                                                                          selected = "primary.predict"),
                                                       ),
                                                       column(3,
                                                              uiOutput("gene1_mult.select"),
                                                       ),
                                                       column(3,
                                                              uiOutput("gene2_mult.select"),
                                                       ),
                                                       br(),
                                                       column(3,
                                                              actionBttn("doVolcano_mult", "Run Pairwise DEGs", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_mult7).removeClass('disabled')"),
                                                              #bsPopover("doVolcano", "Perform Pairwise DEG Analysis","Press to run Pairwise DEG Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                       ),
                                                ),
                                                column(12,
                                                       column(3, numericInput("min_pct_a_mult",
                                                                              label = "min.pct",
                                                                              value = 0.25,
                                                                              min = 0,
                                                                              step = 0.01)
                                                       ),
                                                       column(3, numericInput("logfc_a_mult",
                                                                              label = "logfc.threshold",
                                                                              value = 0.25,
                                                                              min = 0,
                                                                              step = 0.01)
                                                       ),
                                                       column(3, selectInput("test.use_a_mult",
                                                                             label = "Test use",
                                                                             choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                       )
                                                ),
                                                column(12,
                                                       plotOutput("volcano_mult.plot", width = "75%"),
                                                ),
                                                br(),
                                                column(12,
                                                       plotOutput("dega_mult.plot", width = "100%")
                                                ),
                                              ),
                                           ),
                                      
                                     tabPanel("Data Visualization", value = "test_mult6",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult6 = $(\'a[data-value="test_mult6"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              fluidRow(
                                                column(3,
                                                       selectInput("deg_mult2",
                                                                   label = "Group by",
                                                                   choices = c("seurat_clusters", "primary.predict"),
                                                                   selected = "primary.predict"),
                                                ),
                                              ),  
                                              
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
                                                #column(3,
                                                #       br(),
                                                #       actionBttn("Vis3", "Visualize", icon = icon("hand-pointer-o"))
                                                #),
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
                                                         actionBttn("Vis_mult", "Visualize", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_sp5).removeClass('disabled')"),
                                                         #br(),
                                                         #uiOutput("vis.gene.select1"),
                                                         #plotlyOutput("vis1.plot", width = "100%"),
                                                  )),
                                                column(6,
                                                       #uiOutput("deg.gene_mult.select"),
                                                       plotlyOutput("Deg_mult.plot", width = "150%"),
                                                       br(),
                                                       plotlyOutput("vis.plot", width = "100%"),
                                                       #br(),
                                                       #plotlyOutput("Deg_mult1.plot", width = "150%"),
                                                       br(),
                                                       plotOutput("Deg_mult2.plot", width = "150%"),
                                                ),
                                                column(12,
                                                column(4,
                                                       downloadBttn('download_violn_mult', 'Download Violin plot (as png)', size = 'sm'),
                                                       br(),
                                                  ),
                                                column(4,
                                                       downloadBttn('download_feature_mult', 'Download Feature plot (as png)', size = 'sm'),
                                                       br(),
                                                  ),
                                                column(4,
                                                       downloadBttn('download_ridge_mult', 'Download Feature plot (as png)', size = 'sm'),
                                                       br(),
                                                  ),
                                                ),
                                              ),  
                                              conditionalPanel(
                                                condition = "input.scAnalysis_mult == 'MOFA+' & input.scAnalysis_type == 'CITE-seq'",
                                                #column(4,
                                                #       actionBttn("Vis3a", "Visualize", icon = icon("hand-pointer-o"), size = "sm")
                                                #),
                                                
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
                                                #column(3,
                                                #       br(),
                                                #       actionBttn("Vis3_multiome", "Visualize", icon = icon("hand-pointer-o"), size = "sm")
                                                #),
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
                                                )
                                              ),
                                     
                                           
                                     tabPanel("GSEA", value = "test_mult7",
                                              
                                              tags$script(
                                                '
                                                  var tab_mult7 = $(\'a[data-value="test_mult7"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_mult.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              column(12,
                                                     br(),
                                                     column(2, selectInput("species_gsea_mult",
                                                                           label = "Species",
                                                                           choices = c("Homo sapiens", "Mus musculus"))
                                                     ),
                                                     column(2, selectInput("category_gsea_mult",
                                                                           label = "Collection",
                                                                           choices = c("H", "C2", "C5", "C7", "C8"))
                                                     ),
                                                     column(2,
                                                            uiOutput("gsea.ct_mult1.select"),
                                                     ),
                                                     column(2,
                                                            uiOutput("gsea.ct_mult2.select"),
                                                     ),
                                                     
                                                     
                                                     column(2, numericInput("min_pct_mult1",
                                                                            label = "min.pct",
                                                                            value = 0.25,
                                                                            min = 0,
                                                                            step = 0.01)
                                                     ),
                                                     column(2, numericInput("logfc_mult1",
                                                                            label = "logfc.threshold",
                                                                            value = 0.25,
                                                                            min = 0,
                                                                            step = 0.01)
                                                     ),
                                              ),
                                              column(12,
                                                     column(2, selectInput("test.use_mult1",
                                                                           label = "Test use",
                                                                           choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                                            ),
                                                     column(4,
                                                            uiOutput("gsea_mult.select"),
                                                            ),
                                                     br(),
                                                     column(3,
                                                            actionBttn("gsea_mult", "Run gene set enrichment analysis", icon = icon("hand-pointer-o"), size = 'sm'),
                                                            #bsPopover("gsea", "Perform Gene Set Enrichment Analysis","Press to run Gene Set Enrichment Analysis", placement = "bottom", trigger = "hover", options = NULL),
                                                            ),
                                                     textOutput("gsea_mult.done"),
                                                     br()
                                                    ),
                                              br(),
                                              column(7,
                                                     plotOutput("gsea_mult_plot", width = "100%"),
                                                     ),
                                              br(),
                                              br(),
                                              column(12,
                                                     plotOutput("gsea_mult_plot1", width = "100%"),
                                                    ),
                                              br(),
                                              br(),
                                              DT::dataTableOutput("gsea_mult.table"),
                                              column(4,
                                                     downloadBttn('download_gsea_mult.table', 'Download GSEA Results (in csv)', size = 'sm'),
                                                    ),
                                            #),
                                          ),
                                        ),
                        
               ),
               
               ##------------Single cell ATAC-seq----------------##
               
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
                                                              titlePanel(h4(p("Load your example data"))),
                                                              actionBttn("loadexample_atac", "Load example and run", icon = icon("hand-o-right"), size = 'sm', onclick = "$(tab_atac).removeClass('disabled')"),
                                                              titlePanel(h4(HTML("<b>Load your input data</b>"))),
                                                              
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
                                                              
                                                              textInput(inputId = "projName4",
                                                                        label = "Project Name",
                                                                        value = "scATAC"),
                                                              fluidRow(
                                                                actionBttn("loadButton_atac", "Load Data", icon = icon("hand-o-right"), size = "sm", onclick = "$(tab_atac).removeClass('disabled')"),
                                                                actionBttn("reset_atac", "Reset Data", icon = icon("repeat"), size = "sm"),
                                                              ))),
                                                     
                                                               fluidRow(
                                                                 column(6,
                                                                        plotlyOutput("TSSPlot", width = "150%"),
                                                                 ),
                                                                 column(6,
                                                                      plotlyOutput("FragmentHistogram", width = "150%"),
                                                                 ),
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
                                                                        br(),
                                                                        actionBttn("create_seurat_atac", "Process", icon = icon("hand-o-right"), size = "sm")
                                                                 ),
                                                                 column(12,
                                                                 column(4,
                                                                      plotOutput("Vlnplot_atac_1", width = "150%")
                                                                 ),
                                                                 column(4,
                                                                      plotOutput("Vlnplot_atac_2", width = "150%")
                                                                 ),
                                                                 column(4,
                                                                      plotOutput("Vlnplot_atac_3", width = "150%")
                                                                 ),
                                                                 column(4,
                                                                      plotOutput("Vlnplot_atac_4", width = "150%")
                                                                 ),
                                                                 column(4,
                                                                      plotOutput("Vlnplot_atac_5", width = "150%")
                                                                 )),
                                                               #),
                                                     
                                              ),
                                              
                                              column(12,
                                                     h4(p("ATAC counts")),
                                                     withSpinner(dataTableOutput('countdataDT_atac')),
                                              ),
                                     ),
                                     
                                     #tabPanel(title="ATAC-seq QC Plot", value = "ATAC-seq_QC_panel1",
                                     #         br(),
                                    ),
                                     
                                     tabPanel("Normalization", value = "test_atac",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac = $(\'a[data-value="test_atac"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                       tabPanel(title="Normalization", value = "norm_ATAC",
                                                fluidRow(
                                                  column(3,
                                                         selectInput("top_features_atac",
                                                                      label = "Top features to use",
                                                                     choices = c("q0", "q5", "q10", "q25", "q50", "q75"),
                                                                     selected = "q5")),
                                                  column(3,
                                                         numericInput("npc_atac",
                                                                      label = "Singular values to compute",
                                                                      value = 50,
                                                                      min = 1,
                                                                      max = 50,
                                                                      step = 1)
                                                         ),
                                                 br(),
                                                  column(3,
                                                         actionBttn("donorm_ATAC", "Run Normalization", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac1).removeClass('disabled')"),
                                                         textOutput("normalize_atac.done"),
                                                         br()
                                                  )),
                                                br(),
                                                plotlyOutput("DepthCor_ATAC", width = "75%"),
                                       )),
                                     
                                     
                                     
                                     tabPanel("UMAP", value = "test_atac1",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac1 = $(\'a[data-value="test_atac1"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                       fluidRow(
                                         column(3,
                                                numericInput("dim.used_atac",
                                                             label = "Dimensions used",
                                                             value = 30)
                                         ),
                                         br(),
                                         column(3,
                                                actionBttn("doUMAP_ATAC", "Run UMAP", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac2).removeClass('disabled')"),
                                                textOutput("umap_atac.done"),
                                                br()
                                         )),
                                       br(),
                                       plotlyOutput("Umap_ATAC", width = "50%"),
                                     ),
                                     
                                     tabPanel("tSNE", value = "test_atac1",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac1 = $(\'a[data-value="test_atac1"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                       fluidRow(
                                         column(3,
                                                numericInput("dim.used_atac",
                                                             label = "Dimensions used",
                                                             value = 30)
                                         ),
                                         br(),
                                         column(3,
                                                actionBttn("doTSNE_ATAC", "Run TSNE", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac4).removeClass('disabled')"),
                                                textOutput("tsne_atac.done"),
                                                br()
                                         )),
                                       br(),
                                       plotlyOutput("Tsne_ATAC", width = "50%"),
                                     ),
                                    
                                    tabPanel("Clustering", value = "test_atac2",
                                             
                                             tags$script(
                                               '
                                                  var tab_atac2 = $(\'a[data-value="test_atac2"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
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
                                                                  selected = 30 
                                                      ),
                                               ),
                                               
                                               column(3,
                                                      br(),
                                                      actionBttn("doCluster_ATAC", "Run Clustering", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac4).removeClass('disabled')"),
                                                      textOutput("cluster_atac.done"),
                                                      
                                                      
                                               )),
                                             br(),
                                             plotlyOutput("cluster_ATAC", width = "50%"),
                                    ),
                        
                                     tabPanel("Calculate gene activity", value = "test_atac4",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac4 = $(\'a[data-value="test_atac4"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              br(),
                                              column(3,
                                                     actionBttn("calc_gene_activity_atac", "Calculate gene activity", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac5).removeClass('disabled')")
                                              ),
                                              br(),
                                              column(6,
                                                     uiOutput("gene_activity.select"),
                                                     plotlyOutput("gene_activity.plot", width = "100%")
                                              ),
                                              column(4,
                                                     downloadBttn('download_gene_activity', 'Download Gene Activity plot (as png)', size = 'sm'),
                                                     br(),
                                              ),
                                     ),
                                     tabPanel("Cell type identification", value = "test_atac5",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac5 = $(\'a[data-value="test_atac5"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              br(),
                                              fluidRow(
                                              column(3,
                                                     selectInput("module_1",
                                                                 label = "From module",
                                                                 choices = c("Single cell RNA-seq", "Data Integration"),
                                                                 selected = "Single cell RNA-seq"),
                                                ),
                                              ),  
                                              
                                              conditionalPanel(
                                                condition = "input.module_1 == 'Single cell RNA-seq'",
                                              tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler4",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell RNA-Sequencing)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                              tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler4a",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell RNA-Sequencing)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                              column(3,
                                                     selectInput("cellid_method10",
                                                                 label = "Celltype annotation method",
                                                                 choices = c("CELLiD", "Celltypist"),
                                                                 selected = "CELLiD"),
                                              ),
                                              conditionalPanel(
                                                condition = "input.cellid_method10 == 'CELLiD'",
                                              column(3,
                                                     selectInput("cellatlas_atac",
                                                                 label = "Reference Atlas",
                                                                 choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                 selected = "all")
                                              ),
                                              column(3,
                                                     br(),
                                                     actionBttn("load_eg_scRNA", "Load example scRNA-seq data and get celltype", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac6).removeClass('disabled')"),
                                              ),
                                              tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler4",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell RNA-Sequencing)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                              column(3,
                                                     br(),
                                                     actionBttn("load_user_scRNA", "Load user scRNA-seq data and get celltype", icon = icon("hand-pointer-o"), size = "sm")
                                              ),
                                              br(),
                                              column(12,
                                                     plotOutput("cellanno_atac.plot", width = "50%"),
                                              ),
                                              column(12,
                                                     actionBttn("annotate_atacseq", "Annotate scATAC-seq", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac6).removeClass('disabled')")
                                                     ),
                                              column(12,
                                                     plotOutput("cellanno_atac1.plot", width = "50%")
                                              ),
                                              column(3,
                                                     downloadBttn('download_cellanno_plot', 'Download cell annotation Plot (as png)', size = 'sm'),
                                              ),
                                              ),
                                              conditionalPanel(
                                                condition = "input.cellid_method10 == 'Celltypist'",
                                                column(3,
                                                       #br(),
                                                       selectInput("celltypistatlas_atac",
                                                                   label = "Reference Atlas",
                                                                   choices = c("Immune_All_Low.pkl", "Autopsy_COVID19_Lung.pkl", "Pan_Fetal_Human.pkl", "Nuclei_Lung_Airway.pkl", "Developing_Human_Thymus.pkl", "Human_Lung_Atlas.pkl", "Developing_Mouse_Brain.pkl", "Developing_Human_Brain.pkl", "Cells_Lung_Airway.pkl", "Healthy_COVID19_PBMC.pkl", "Human_IPF_Lung.pkl", "Adult_Mouse_Gut.pkl", "Immune_All_High.pkl", "COVID19_Immune_Landscape.pkl", "Human_PF_Lung.pkl", "COVID19_HumanChallenge_Blood.pkl", "Lethal_COVID19_Lung.pkl", "Cells_Fetal_Lung.pkl", "Cells_Intestinal_Tract.pkl"),
                                                                   selected = "Immune_All_Low.pkl")
                                                ),
                                                column(3,
                                                       br(),
                                                       actionBttn("load_eg1_scRNA", "Load example scRNA-seq data and get celltype", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac6).removeClass('disabled')"),
                                                ),
                                                tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler4a",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell RNA-Sequencing)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                                column(3,
                                                       br(),
                                                       actionBttn("load_user1_scRNA", "Load user scRNA-seq data and get celltype", icon = icon("hand-pointer-o"), size = "sm")
                                                ),
                                                br(),
                                                column(12,
                                                       plotOutput("cellanno1_atac.plot", width = "50%"),
                                                ),
                                                column(12,
                                                       actionBttn("annotate1_atacseq", "Annotate scATAC-seq", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac6).removeClass('disabled')")
                                                ),
                                                column(12,
                                                       plotOutput("cellanno1_atac1.plot", width = "50%")
                                                )
                                              )),
                                              conditionalPanel(
                                                condition = "input.module_1 == 'Data Integration'",
                                                
                                                column(3,
                                                       selectInput("cellid_method10a",
                                                                   label = "Celltype annotation method",
                                                                   choices = c("CELLiD", "Celltypist"),
                                                                   selected = "CELLiD"),
                                                ),
                                                conditionalPanel(
                                                  condition = "input.cellid_method10a == 'CELLiD'",
                                                  column(3,
                                                         selectInput("cellatlas_atac_a",
                                                                     label = "Reference Atlas",
                                                                     choices = c("all", "adipose", "adrenal_gland", "blood", "bone_marrow", "brain", "breast", "breast_milk", "eye", "gut", "heart", "kidney", "liver", "lung", "pancreas", "PDAC", "skin", "testis", "thymus", "tonsil"),
                                                                     selected = "all")
                                                  ),
                                                  column(3,
                                                         br(),
                                                         actionBttn("load_eg_scRNA_a", "Load example scRNA-seq data and get celltype", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac6).removeClass('disabled')"),
                                                  ),
                                                  tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler_4",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell data integration)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                                  column(3,
                                                         br(),
                                                         actionBttn("load_user_scRNA_a", "Load user multiple scRNA-seq data, perform batch correction and get celltype", icon = icon("hand-pointer-o"), size = "sm")
                                                  ),
                                                  br(),
                                                  column(12,
                                                         plotOutput("cellanno_atac.plot_a", width = "50%"),
                                                  ),
                                                  column(12,
                                                         actionBttn("annotate_atacseq_a", "Annotate scATAC-seq", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac6).removeClass('disabled')")
                                                  ),
                                                  column(12,
                                                         plotOutput("cellanno_atac1.plot_a", width = "50%")
                                                  )),
                                                conditionalPanel(
                                                  condition = "input.cellid_method10a == 'Celltypist'",
                                                  column(3,
                                                         #br(),
                                                         selectInput("celltypistatlas_atac_a",
                                                                     label = "Reference Atlas",
                                                                     choices = c("Immune_All_Low.pkl", "Autopsy_COVID19_Lung.pkl", "Pan_Fetal_Human.pkl", "Nuclei_Lung_Airway.pkl", "Developing_Human_Thymus.pkl", "Human_Lung_Atlas.pkl", "Developing_Mouse_Brain.pkl", "Developing_Human_Brain.pkl", "Cells_Lung_Airway.pkl", "Healthy_COVID19_PBMC.pkl", "Human_IPF_Lung.pkl", "Adult_Mouse_Gut.pkl", "Immune_All_High.pkl", "COVID19_Immune_Landscape.pkl", "Human_PF_Lung.pkl", "COVID19_HumanChallenge_Blood.pkl", "Lethal_COVID19_Lung.pkl", "Cells_Fetal_Lung.pkl", "Cells_Intestinal_Tract.pkl"),
                                                                     selected = "Immune_All_Low.pkl")
                                                  ),
                                                  column(3,
                                                         br(),
                                                         actionBttn("load_eg1_scRNA_a", "Load example scRNA-seq data and get celltype", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac6).removeClass('disabled')"),
                                                  ),
                                                  tags$head(tags$script('
                                                                    Shiny.addCustomMessageHandler("myCallbackHandler_4a",
                                                                    function(typeMessage) {console.log(typeMessage)
                                                                    if(typeMessage == 2){
                                                                    console.log("got here");
                                                                    $("a:contains(Single cell data integration)").click();
                                                                    }
                                                                    });
                                                                    ')), 
                                                  column(3,
                                                         br(),
                                                         actionBttn("load_user1_scRNA_a", "Load user multiple scRNA-seq data, perform batch correction and get celltype", icon = icon("hand-pointer-o"), size = "sm")
                                                  ),
                                                  br(),
                                                  column(12,
                                                         plotOutput("cellanno1_atac.plot_a", width = "50%"),
                                                  ),
                                                  column(12,
                                                         actionBttn("annotate1_atacseq_a", "Annotate scATAC-seq", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac6).removeClass('disabled')")
                                                  ),
                                                  column(12,
                                                         plotOutput("cellanno1_atac1.plot_a", width = "50%")
                                                  )
                                                )),
                                              
                                              ),
                                       
                                     tabPanel("DE Peaks", value = "test_atac6",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac6 = $(\'a[data-value="test_atac6"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                       fluidRow(
                                         column(3, numericInput("min_pct_atac",
                                                                label = "min.pct",
                                                                value = 0.75,
                                                                min = 0,
                                                                step = 0.01)
                                         ),
                                         
                                         column(3,numericInput("logfc_atac",
                                                               label = "logfc.threshold",
                                                               value = 0.75,
                                                               min = 0,
                                                               step = 0.01)
                                         ),
                                         
                                         column(3, selectInput("test.use_atac",
                                                               label = "Test use",
                                                               choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                         ),
                                         br(),
                                         column(3,
                                                actionBttn("doDeg_atac", "Run DE Peaks", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac7).removeClass('disabled')")
                                         )),
                                       br(),
                                       column(9,
                                              DT::dataTableOutput("Deg_atac.table"),
                                              br(),
                                              plotlyOutput("Deg_atac.plot", width = "100%")
                                       )
                                     ),
                                     
                                    tabPanel("Target Predict (Link Peaks to Genes)", value = "test_atac7",
                                             
                                             tags$script(
                                               '
                                                  var tab_atac7 = $(\'a[data-value="test_atac7"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                      
                                      column(6,
                                             uiOutput("peak_gene.atac.select"),
                                            ),
                                      plotOutput("coverage1.plot", width = "100%"),
                                      br(),
                                      br(),
                                      br(),
                                      br(),
                                      column(12,
                                      actionBttn("link_peak_genes", "Link peaks to genes", icon = icon("hand-pointer-o"), size = "sm", onclick = "$(tab_atac8).removeClass('disabled')")
                                      ),
                                      br(),
                                      #column(12,
                                      #plotOutput("coverage2.plot", width = "100%"),
                                      #      ),
                                      br(),
                                      DT::dataTableOutput("peaks.table"),
                                      column(4,
                                             downloadBttn('download_coverage_atac', 'Download Coverage plot (as png)', size = 'sm'),
                                             br(),
                                      ),
                                    ),
                                    
                                    
                                    tabPanel("Gene Set Enrichment Analysis for scATAC-seq", value = "test_atac8",
                                             
                                             tags$script(
                                               '
                                                  var tab_atac8 = $(\'a[data-value="test_atac8"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
                                      fluidRow(
                                        column(2,
                                               selectInput("gsea_atac_method",
                                                           label = "Enrichment method",
                                                           choices = c("fgsea", "GREAT"),
                                                           selected = "fgsea"),
                                        ),
                                      ),
                                      conditionalPanel(
                                        condition = "input.gsea_atac_method == 'fgsea'",
                                      column(12,
                                             br(),
                                             column(2, selectInput("species_gsea_atac",
                                                                   label = "Species",
                                                                   choices = c("Homo sapiens", "Mus musculus"))
                                             ),
                                             column(2, selectInput("category_gsea_atac",
                                                                   label = "Collection",
                                                                   choices = c("H", "C2", "C5", "C7", "C8"))
                                             ),
                                             column(2,
                                                    uiOutput("gsea_atac.ct1.select"),
                                             ),
                                             column(2,
                                                    uiOutput("gsea_atac.ct2.select"),
                                             ),
                                             
                                             
                                             column(2, numericInput("min_pct_atac1",
                                                                    label = "min.pct",
                                                                    value = 0.25,
                                                                    min = 0,
                                                                    step = 0.01)
                                             ),
                                             column(2, numericInput("logfc_atac1",
                                                                    label = "logfc.threshold",
                                                                    value = 0.25,
                                                                    min = 0,
                                                                    step = 0.01)
                                             ),
                                      ),
                                      column(12,
                                             column(2, selectInput("test.use_atac1",
                                                                   label = "Test use",
                                                                   choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                             ),
                                             column(4,
                                                    uiOutput("gsea_atac.select"),
                                             ),
                                             br(),
                                             column(3,
                                                    actionBttn("gsea_atac", "Run gene set enrichment analysis", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac9).removeClass('disabled')"),
                                             ),
                                             textOutput("gsea_atac.done"),
                                             br()
                                      ),
                                      br(),
                                      column(7,
                                             plotOutput("gsea_atac_plot", width = "100%"),
                                      ),
                                      br(),
                                      br(),
                                      column(12,
                                             plotOutput("gsea_atac_plot1", width = "100%"),
                                      ),
                                      br(),
                                      br(),
                                      DT::dataTableOutput("gsea_atac.table"),
                                      column(4,
                                             downloadBttn('download_gsea_atac.table', 'Download Celltype similarity table (in csv)', size = 'sm'),
                                      )),
                                      conditionalPanel(
                                        condition = "input.gsea_atac_method == 'GREAT'",
                                        column(12,
                                               br(),
                                               column(2, selectInput("species_gsea_atac",
                                                                     label = "Species",
                                                                     choices = c("Homo sapiens", "Mus musculus"))
                                               ),
                                               
                                               column(2,
                                                      uiOutput("great_atac.ct1.select"),
                                               ),
                                               column(2,
                                                      uiOutput("great_atac.ct2.select"),
                                               ),
                                               
                                               column(2, selectInput("gene_set_atac",
                                                                     label = "Gene set",
                                                                     choices = c("BP", "CC", "MP", "H", "C2", "C5", "C7", "C8"))
                                               ),
                                               
                                               column(3, selectInput("tss_atac",
                                                                     label = "TSS source",
                                                                     choices = c("txdb:hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene", "RefSeq:hg19", "GREAT:hg19", "Gencode_v19"))
                                               ),
                                               column(2, numericInput("min_pct_great",
                                                                      label = "min.pct",
                                                                      value = 0.25,
                                                                      min = 0,
                                                                      step = 0.01)
                                                      ),
                                               column(2, numericInput("logfc_great",
                                                                      label = "logfc.threshold",
                                                                      value = 0.25,
                                                                      min = 0,
                                                                      step = 0.01)
                                                     ),
                                               column(2, selectInput("test.use_atac2",
                                                                     label = "Test use",
                                                                     choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                               ),
                                              column(3,
                                                      actionBttn("gsea_atac1", "Run gene set enrichment analysis", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac9).removeClass('disabled')"),
                                                     ),
                                               textOutput("gsea_atac1.done"),
                                               br(),
                                               plotOutput("great.plot", width = "100%"),
                                               br(),
                                               plotlyOutput("great1.plot", width = "100%", height = "800px"),
                                        )),
                                    ),
                                    
                                    tabPanel("Data visualization", value = "test_atac9",
                                             
                                             tags$script(
                                               '
                                                  var tab_atac9 = $(\'a[data-value="test_atac9"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                             ),
                                             
                                              fluidRow(
                                                column(6,
                                                       uiOutput("deg.atac.select"),
                                                       actionBttn("Vis_atac", "Visualize", icon = icon("hand-pointer-o"), size = 'sm', onclick = "$(tab_atac10).removeClass('disabled')"),
                                                       plotlyOutput("Deg_atac1.plot", width = "100%"),
                                                       br(),
                                                       plotlyOutput("Deg_atac2.plot", width = "100%"),
                                                       br(),
                                                       plotOutput("Deg_atac3.plot", width = "100%")
                                                ),
                                              ),
                                             column(12,
                                                column(4,
                                                       downloadBttn('download_violn_atac', 'Download Violin plot (as png)', size = 'sm'),
                                                       br(),
                                                ),
                                                column(4,
                                                       downloadBttn('download_feature_atac', 'Download Feature plot (as png)', size = 'sm'),
                                                       br(),
                                                ),
                                                column(4,
                                                       downloadBttn('download_ridge_atac', 'Download Ridge plot (as png)', size = 'sm'),
                                                       br(),
                                                ),
                                              ),
                                     ),
                                     tabPanel("Coverage Plot", value = "test_atac10",
                                              
                                              tags$script(
                                                '
                                                  var tab_atac10 = $(\'a[data-value="test_atac10"]\').parent().addClass("disabled");
                                                  $(function(){
                                                    $(tab_atac.parent()).on("click", "li.disabled", function(e) {
                                                      e.preventDefault();
                                                      return false;
                                                    });
                                                  });
                                                  '
                                              ),
                                              
                                              column(3,
                                              uiOutput("coverage.atac.select"),
                                              ),
                                              column(3,
                                              uiOutput("coverage.atac_feature.select"),
                                              ),
                                              plotOutput("coverage.plot", width = "100%")
                                     ),
                                    
                                     #tabPanel("Integrate with scRNA-seq",
                                     #         tags$head(tags$script('
                                     #                                Shiny.addCustomMessageHandler("myCallbackHandler",
                                     #                                function(typeMessage) {console.log(typeMessage)
                                     #                                if(typeMessage == 1){
                                     #                                console.log("got here");
                                     #                                $("a:contains(Single cell RNA-Sequencing)").click();
                                     #                                }
                                     #                                });
                                     #                                ')),                                   
                                     #          column(4,
                                     #                 br(),
                                     #                 actionBttn("integrate_scATAC_scRNA", "Integrate scATAC and scRNA-seq data", icon = icon("hand-pointer-o")),
                                     #          ),
                                     #          column(4,
                                     #                 br(),
                                     #                 actionBttn("process_scrnaseq", "Process scRNA-seq data", icon = icon("hand-pointer-o")),
                                     #          ),
                                     #),
                                     #tabPanel("Motif Analysis",
                                              
                                     #         tabsetPanel(id="motif_analysis",
                                     #                      tabPanel(title="Motif analysis",
                                     #                               br(),
                                     #                               column(3, numericInput("min_pct_motif",
                                     #                                                      label = "min.pct",
                                     #                                                      value = 0.25,
                                     #                                                      min = 0,
                                     #                                                      step = 0.01)
                                     #                               ),
                                     #                               column(3, numericInput("logfc_motif",
                                     #                                                      label = "logfc.threshold",
                                     #                                                      value = 0.25,
                                     #                                                      min = 0,
                                     #                                                      step = 0.01)
                                     #                               ),
                                     #                               column(3, selectInput("test.use_motif",
                                     #                                                     label = "Test use",
                                     #                                                     choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                     #                               ),
                                     #                               br(),
                                     #                               column(3,
                                     #                                      actionButton("doDeg_motif", "Run DEGs", icon = icon("hand-pointer-o"))
                                     #                               ),  
                                     #                               br(),
                                     #                               fluidRow(
                                     #                                 column(6,
                                     #                                        uiOutput("motif.select"),
                                     #                                        plotOutput("motif.plot", width = "100%")
                                     #                                 )
                                     #                               ),
                                     #                      ),
                                     #                      tabPanel(title="Motif activities",
                                     #                               br(),
                                     #                               column(3,
                                     #                                      actionButton("calc_motif_activity", "Calculate Motif Activity", icon = icon("hand-pointer-o")),
                                     #                               ),
                                     #                               column(3,
                                     #                                      uiOutput("motif_feature.select"),
                                     #                               ),
                                     #                               br(),
                                     #                               fluidRow(
                                     #                                 column(12,
                                     #                                        plotlyOutput("motif_feature.plot", width = "50%"),
                                     #                                 )),  
                                     #                               br(),
                                     #                               fluidRow(
                                     #                                 column(12,
                                     #                                        uiOutput("motif1.select"),
                                     #                                        plotOutput("motif1.plot", width = "50%")
                                     #                                 )),
                                                                   
                                     #                     ),
                                     #          )),
                                     
                                     
                        ),
                        
               ),
               
              
               
               ##--------------Help page--------------##
               
               tabPanel("Help",
                        column(6,
                        h1("User manual"),
                        h4("Manual for scRNA-seq analysis in ezSingleCell: ", tags$a(href="https://drive.google.com/file/d/1CWjcJaNjMh6f2OcvLtakcD2wJqu-wJZa/view?usp=share_link", "Manual")),
                        h4("Manual for Data Integration analysis in ezSingleCell: ", tags$a(href="https://drive.google.com/file/d/1ghpi5ZR9195LHa4J77XFE6XN_nG253sc/view?usp=share_link", "Manual")),
                        h4("Manual for Spatial Transcriptomics analysis in ezSingleCell: ", tags$a(href="https://drive.google.com/file/d/1az1bM_7i4dZopYJSYpKpiqFkSo_Uo_CN/view?usp=share_link", "Manual")),
                        h4("Manual for scMultiomics analysis in ezSingleCell: ", tags$a(href="https://drive.google.com/file/d/1EL7334xnfZyobJ1MhXm8brOG1zp-dmMo/view?usp=drive_link", "Manual")),
                        h4("Manual for scATAC-seq analysis in ezSingleCell:", tags$a(href="https://drive.google.com/file/d/1nncXRZX3YfySgk46oj3sMzBTVyjd860E/view?usp=drive_link", "Manual")),
                        ),
                        column(6,
                               imageOutput("demo_image2"),     
                        ),
                        br(),
                        br(),
                        br(),
                        br(),
                        column(12,
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                               br(),
                        
                        ))
               )
)))
