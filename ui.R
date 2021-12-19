#source("install_R_packages.R")
library(shiny)
library(plotly)
library(ggplot2)
library(shinyjs)
library(DT)
library(devtools)
library(assertthat)
library(Seurat)
#library(SeuratData)
library(cowplot)
library(sctransform)
library(SeuratWrappers)
library(dplyr)
library(cytomapper)
library(data.table)
library(selectr)
library(readr)
library(hdf5r)
library(future)
library(CellChat)
library(tidyverse)
library(pbmcapply)
library(monocle3)
library(leidenbase)
library(Biobase)
library(biovizBase)
library(SingleCellExperiment)
library(harmony)
library(MOFA2)
library(flowCore)
library(cytofkit2)
#library(rliger)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Signac)
library(BayesSpace)
library(SPOTlight)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)
library(reshape2)
library(shinydashboard)
library(shinyalert)
library(shinyFiles)
library(shinyWidgets)
library(shinythemes)

#import('h5py')

shiny_one_panel = fluidPage(
                   theme = shinytheme("united"),
    titlePanel("ezSinglecell : An integrated one-stop single-cell analysis toolbox for bench scientists"),
    hr(),

    fluidRow(
        ##------Sidebar---------
        column(3,
               #h4('Load Data:'),
               wellPanel(
                   selectInput("Module",
                               label = "Module",
                               choices = c("Homepage", "Single Cell RNASeq Analysis", "Single cell data integration Analysis", "Single cell multiomics Analysis", "Flow cytometry Analysis", "Imaging mass cytometry Analysis", "Single cell ATAC-seq Analysis", "Spatial Transcriptomics Analysis"),
                               selected = "Single Cell RNASeq Analysis"),

                   fluidRow(
                       conditionalPanel(
                           condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'Raw Counts Matrix' & input.Species == 'human'",
                           column(12,
                                  actionButton("loadexample_tpm_human", "Load example and run", icon = icon("hand-o-right"))
                           )),

                       conditionalPanel(
                           condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'H5' & input.Species == 'human'",
                           column(12,
                                  actionButton("loadexample_scH5_human", "Load example and run", icon = icon("hand-o-right"))
                           )),

                       conditionalPanel(
                           condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'Raw Counts Matrix' & input.Species == 'mouse'",
                           column(12,
                                  actionButton("loadexample_tpm_mouse", "Load example and run", icon = icon("hand-o-right"))
                           )),

                       conditionalPanel(
                           condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'H5' & input.Species == 'mouse'",
                           column(12,
                                  actionButton("loadexample_scH5_mouse", "Load example and run", icon = icon("hand-o-right"))
                           )),

                       #br(),

                   ),

                   conditionalPanel(
                       condition = "input.Module == 'Single Cell RNASeq Analysis'",
                       selectInput("scInput",
                                   label = "Select Data Input Type",
                                   choices = c("Raw Counts Matrix", "H5"),
                                   selected = "Raw Counts Matrix")),

                   conditionalPanel(
                       condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'Raw Counts Matrix'",
                       fileInput("tpmFiles",
                                 label = "Upload Counts File (Accepted Format: tab delimited text)",
                                 accept = ".txt")),

                   conditionalPanel(
                       condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'Raw Counts Matrix'",
                       fileInput("cellAnnoFiles",
                                 label = "Upload Metadata (Accepted Format: tab delimited text)",
                                 accept = ".txt")),

                   conditionalPanel(
                       condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'H5'",
                       fileInput("scH5",
                                 label = "Upload H5 output from Cellranger (Accepted Format: .h5)",
                                 accept = ".h5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix' || input.scInput1 == 'H5' & input.Species == 'human'",
                       fluidRow(
                           actionButton("loadexample1", "Load example and run", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'CITE-seq' & input.scInput2 == 'H5' & input.Species == 'human'",
                       fluidRow(
                           actionButton("loadexample2", "Load example data", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'mRNA+scATAC-seq' & input.scInput2 == 'H5' & input.scAnalysis_mult == 'Seurat' & input.Species == 'human'",
                       fluidRow(
                           actionButton("loadexample2_b", "Load example data", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis'",
                       selectInput("scInput1",
                                   label = "Select Data Input Type",
                                   choices = c("Raw Counts Matrix", "H5"),
                                   selected = "Raw Counts Matrix")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis'",
                       selectInput("scAnalysis_integ",
                                   label = "Analysis method",
                                   choices = c("Seurat", "Harmony"),
                                   selected = "Seurat")),

                   conditionalPanel(
                       condition = "input.Module == 'Single Cell RNASeq Analysis' || input.Module == 'Single cell data integration Analysis'",
                       selectInput("Species",
                                   label = "Select Species",
                                   choices = c("human", "mouse"),
                                   selected = "human")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                       fileInput("tpmFiles1",
                                 label = "Upload Counts File (Accepted Format: tab delimited text)",
                                 accept = ".txt",
                                 multiple = T)),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                       fileInput("cellAnnoFiles1",
                                 label = "Upload Metadata (Accepted Format: tab delimited text)",
                                 accept = ".txt",
                                 multiple = T)),

                   #conditionalPanel(
                   #    condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'Raw Counts Matrix' & input.demo_analyse_sc == 'Demo'",
                   #    selectInput("tpmFiles_demo",
                   #                label = "Upload Counts File (Accepted Format: tab delimited text)",
                   #                choices = c(" ", "TPM.txt"),
                   #                selected = " ")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                       fileInput("scH5_1",
                                 label = "Upload H5 output from Cellranger (Accepted Format: .h5)",
                                 accept = ".h5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis'",
                       selectInput("scInput2",
                                   label = "Select Data Input Type",
                                   choices = "H5",
                                   selected = "H5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis'",
                       selectInput("scAnalysis_mult",
                                   label = "Analysis method",
                                   choices = c("Seurat", "MOFA"),
                                   selected = "Seurat")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis'",
                       selectInput("scAnalysis_type",
                                   label = "scAnalysis type",
                                   choices = c("CITE-seq", "mRNA+scATAC-seq"),
                                   selected = "H5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'CITE-seq'",
                       fileInput("tpmFiles2",
                                 label = "Upload H5 output from Cellranger (Accepted Format: .h5)",
                                 multiple = FALSE,
                                 accept = ".h5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_mult == 'MOFA' & input.scAnalysis_type == 'CITE-seq'",
                       fileInput("meta_mofa",
                                 label = "Upload Metadata (Accepted Format: .csv)",
                                 multiple = FALSE,
                                 accept = c(".csv", ".txt"))),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'mRNA+scATAC-seq'",
                       fileInput(inputId = 'tpmFiles3b',
                                 label = "Gene expression file",
                                 accept = ".h5"),

                       conditionalPanel(
                           condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'mRNA+scATAC-seq'",
                           shinyFiles::shinyFilesButton(id = 'dir_multi_atac', label = "Path to fragments file", title = "Sheets Folder Selector", multiple = T),
                           verbatimTextOutput("dir_multi_atac", placeholder = TRUE)
                       )),

                   #fileInput(inputId = 'fragFiles',
                   #             label = "Fragments file",
                   #             multiple = FALSE,
                   #             accept = ".gz"),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell ATAC-seq Analysis'",
                       fluidRow(
                           actionButton("loadexample6", "Load example and run", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell ATAC-seq Analysis'",
                       selectInput("scInput_atac",
                                   label = "Select Data Input Type",
                                   choices = "H5",
                                   selected = "H5")),

                   conditionalPanel(
                       condition = "input.Module == 'Spatial Transcriptomics Analysis'",
                       fluidRow(
                           actionButton("loadexample5", "Load example and run", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Spatial Transcriptomics Analysis'",
                       selectInput("scInput3",
                                   label = "Select Data Input Type",
                                   choices = "H5",
                                   selected = "H5")),

                   conditionalPanel(
                       condition = "input.Module == 'Spatial Transcriptomics Analysis'",
                       shinyFiles::shinyDirButton(id = 'dir', label = "Path to spaceRanger output file", title = "Sheets Folder Selector"),
                       verbatimTextOutput("dir", placeholder = TRUE)
                   ),

                   conditionalPanel(
                       condition = "input.Module == 'Flow cytometry Analysis'",
                       selectInput("demo_analyse_sc",
                                   label = "Perform demo or analyse",
                                   choices = c("Demo", "Analysis"),
                                   selected = "Demo")),

                   conditionalPanel(
                       condition = "input.Module == 'Flow cytometry Analysis' & input.demo_analyse_sc == 'Demo'",
                       fluidRow(
                           actionButton("loadexample3", "Load example data", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Imaging mass cytometry Analysis'",
                       fluidRow(
                           actionButton("loadexample4", "Load example data", icon = icon("hand-o-right"))
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell ATAC-seq Analysis' & input.scInput_atac == 'H5'",
                       fileInput("tpmFiles_atac",
                                 label = "Upload Peak/Cell matrix",
                                 accept = ".h5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell ATAC-seq Analysis'",
                       fileInput("meta_atac",
                                 label = "Upload Metadata",
                                 accept = ".csv")),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell ATAC-seq Analysis'",
                       shinyFiles::shinyFilesButton(id = 'dir_atac', label = "Path to fragments file", title = "Sheets Folder Selector", multiple = T),
                       verbatimTextOutput("dir_atac", placeholder = TRUE)
                   ),

                   #conditionalPanel(
                   #    condition = "input.Module == 'Single cell ATAC-seq Analysis'",
                   #    fileInput("frag_atac",
                   #              label = "Upload Fragments",
                   #              accept = ".gz")),

                   conditionalPanel(
                       condition = "input.Module == 'Spatial Transcriptomics Analysis'",
                       fileInput("tpmFiles3",
                                 label = "Upload H5 output from Spaceranger (Accepted Format: .h5)",
                                 accept = ".h5")),

                   conditionalPanel(
                       condition = "input.Module == 'Single Cell RNASeq Analysis' & input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'" ,

                       column(6,
                              numericInput(inputId = "min.genes",
                                           label = "Min. genes",
                                           value = 200,
                                           min = 1)
                       ),
                       column(6,
                              numericInput(inputId = "min.cells",
                                           label = "Min. cells",
                                           value = 3,
                                           min = 1)
                       ),

                       textInput(inputId = "projName",
                                 label = "Project Name",
                                 value = "Seurat_analysis"),

                       fluidRow(
                           column(6,
                                  actionButton("loadButton", "Create Seurat Object", icon = icon("hand-o-right"))
                           ),
                           column(6,
                                  actionButton("reset", "Reset", icon = icon("repeat"))
                           ),
                       )),

                   conditionalPanel(
                       condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix' || input.scInput1 == 'H5'" ,

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

                       textInput(inputId = "projName1",
                                 label = "Project Name",
                                 value = "Seurat_analysis"),
                       fluidRow(
                           column(6,
                                  actionButton("loadButton1", "Create Seurat object", icon = icon("hand-o-right"))
                           ),
                           column(6,
                                  actionButton("reset1", "Reset Data", icon = icon("repeat"))
                           ),


                       )
                   )),

               conditionalPanel(
                   condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'CITE-seq'" ,

                   textInput(inputId = "projName2",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton2_a", "Create  Seurat Object", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset2", "Reset Data", icon = icon("repeat"))
                       ),


                   )
               ),

               conditionalPanel(
                   condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_type == 'mRNA+scATAC-seq'" ,

                   textInput(inputId = "projName2",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton2_b", "Create  Seurat Object", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset2", "Reset Data", icon = icon("repeat"))
                       ),


                   )
               ),

               conditionalPanel(
                   condition = "input.Module == 'Imaging mass cytometry Analysis'",

                   fluidRow(
                       ##------Sidebar---------
                       h4('Load Data:'),
                       wellPanel(
                           fileInput(inputId = 'sce',
                                     label = "Input single cell data",
                                     multiple = FALSE,
                                     accept = c(".rds",
                                                ".RData",
                                                ".Robj")),
                           fileInput(inputId = 'mask',
                                     label = "Input cell mask",
                                     multiple = FALSE,
                                     accept = c(".rds",
                                                ".RData",
                                                ".Robj")),
                           fluidRow(
                               column(6,
                                      actionButton("loadButton_img", "Load Data", icon = icon("hand-o-right"))
                               ),
                               column(6,
                                      actionButton("reset_img", "Reset Data", icon = icon("repeat"))
                               )
                           )
                       )
                   )),

               conditionalPanel(
                   condition = "input.Module == 'Single cell ATAC-seq Analysis' & input.scInput_atac == 'H5'" ,


                   column(6,
                          numericInput(inputId = "min.cells_atac",
                                       label = "Min. cells",
                                       value = 3,
                                       min = 1)
                   ),
                   column(6,
                          numericInput(inputId = "min.features_atac",
                                       label = "Min. Features",
                                       value = 200,
                                       min = 1)
                   ),

                   textInput(inputId = "projName",
                             label = "Project Name",
                             value = "Seurat_ATACseq_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton_atac", "Create  Seurat Object", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset_atac", "Reset Data", icon = icon("repeat"))
                       ),
                   )),

               conditionalPanel(
                   condition = "input.Module == 'Spatial Transcriptomics Analysis'" ,

                   textInput(inputId = "projName3",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                       column(6,
                              actionButton("loadButton3", "Create  Seurat Object", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset3", "Reset Data", icon = icon("repeat"))
                       ),


                   )
               ),

               ##------Plot download---------

               conditionalPanel(
                   condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",


                   h4("Export to PDF:"),
                   wellPanel(
                       ## Conditional panel for different plots
                       conditionalPanel(" input.QC == 'QC_panel1' && input.tabs == 'QC plots' ",
                                        actionButton("PDFa", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.QC == 'QC_panel2' && input.tabs == 'QC plots' ",
                                        actionButton("PDFb", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.tabs == 'Normalization and Variable Gene Plot' ",
                                        actionButton("PDFc", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.Pca == 'P_panel1' && input.tabs == 'PCA' ",
                                        actionButton("PDFd", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.Pca == 'P_panel2' && input.tabs == 'PCA' ",
                                        actionButton("PDFe", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.Pca == 'P_panel4' && input.tabs == 'PCA' ",
                                        actionButton("PDFh", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.tabs == 'Clustering' ",
                                        actionButton("PDFf", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.tabs == 'UMAP' ",
                                        actionButton("PDFi", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.tabs == 'TSNE' ",
                                        actionButton("PDFj", "Download", icon = icon("download"))
                       ),
                       conditionalPanel(" input.tabs == 'DEGs' ",
                                        actionButton("PDFk", "Download", icon = icon("download"))
                       ),

                       conditionalPanel(
                           condition = "input.Module == 'Spatial Transcriptomics Analysis' && input.scInput3 == 'H5'",


                           h4("Export to PDF:"),
                           wellPanel(
                               ## Conditional panel for different plots
                               conditionalPanel(" input.Spatial_QC == 'Spatial_QC_panel1' && input.tabs == 'Spatial QC plots' ",
                                                actionButton("PDFa", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.Spatial_QC == 'Spatial_QC_panel2' && input.tabs == 'Spatial QC plots' ",
                                                actionButton("PDFb", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.tabs == 'Normalization and Variable Gene Plot' ",
                                                actionButton("PDFc", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.Pca == 'P_panel1' && input.tabs == 'PCA' ",
                                                actionButton("PDFd", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.Pca == 'P_panel2' && input.tabs == 'PCA' ",
                                                actionButton("PDFe", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.Pca == 'P_panel4' && input.tabs == 'PCA' ",
                                                actionButton("PDFh", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.tabs == 'Clustering' ",
                                                actionButton("PDFf", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.tabs == 'UMAP' ",
                                                actionButton("PDFi", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.tabs == 'TSNE' ",
                                                actionButton("PDFj", "Download", icon = icon("download"))
                               ),
                               conditionalPanel(" input.tabs == 'DEGs' ",
                                                actionButton("PDFk", "Download", icon = icon("download"))
                               ),
                           )),

                       ## ensure no spill over in button text
                       tags$head(
                           tags$style(HTML('
                                   .btn {
                                   white-space: normal;
                                   }'
                           )
                           )
                       ),
                       ## Conditional is separate from pdf options
                       hr(),
                       fluidRow(
                           column(6,
                                  sliderInput(inputId="pdf_w", label = "PDF width(in):",
                                              min=3, max=20, value=8, width=100, ticks=F)
                           ),
                           column(6,
                                  sliderInput(inputId="pdf_h", label = "PDF height(in):",
                                              min=3, max=20, value=8, width=100, ticks=F)
                           )),

                       #actionButton("OpenDir", "Open download folder", icon = icon("folder"))
                   )),

               ##------Save Data---------
               conditionalPanel(
                   condition = "input.Module == 'Single Cell RNASeq Analysis' && input.scInput == 'Raw Counts Matrix' || input.scInput == 'H5'",
                   hr(),
                   actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),

                   hr(),
                   h4(tags$a(href="mailto:Chen_Jinmiao@immunol.a-star.edu.sg?subject=[cytof-question]",
                             "Contact Us")),
                   imageOutput("logo", height = "60px")
               )),
        ##------Main area---------

        conditionalPanel(
            condition = "input.Module == 'Homepage'",
            column(9,
                   tabsetPanel(type = "pills",
                               tabPanel("Overview", fluidRow(
                                   br(),
                                   #img(src='a.jpg',height='300',width='250', align='center'),
                                   p(strong("ezSinglecell"), "is a web server developed by Jinmiao chen's lab with an intention to empower bench scientists to perform downstream Bioinformatics analysis. ezsinglecell contains 7 modules : ", strong("Single cell RNA-seq, Data Integration, Multiomics, Flow cytometry, Imaging mass cytometry, single cell ATAC-seq and Spatial Transcriptomics."), " These modules are designed in order to help researchers with little or no expertise in Bioinformatics design a hypothesis or answer biological questions.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p(strong("ezSinglecell is available on", a("GitHub", href = "https://github.com/raman91/ezsinglecell", target = "_blank"), ", if interested in hosting on your own servers.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")),
                                   p("Please post issues, suggestions and improvements using", a("Issues/suggestions", href = "https://github.com/raman91/ezsinglecell", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p("To view other tools and contributions please visit", a("GitHub", href = "https://github.com/JinmiaoChenLab", target = "_blank"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"))),

                               tabPanel("Getting Started", fluidRow(
                                   p(strong("ezsinglecell"), "is easy to use and is packed with powerful modules to help you analyze your data. Results generated from the modules are loaded on the Result window whereas the Visualization plots are displyed on Visualization windows.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p(strong("Description of various modules"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p(strong("Single Cell RNASeq Analysis module"), "is designed to perform single cell RNA-Seq analysis. ezsinglecell uses ", a("Seurat", href = "https://www.cell.com/action/showPdf?pii=S0092-8674%2821%2900583-3", target = "_blank"), "for analyzing single cell RNA-Seq data. In this module, users can input raw counts matrix, H5 output from cellranger or a processed Seurat object. ezsinglecell allows users to download R Objects after processing.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p(strong("Single Cell RNASeq Integration Analysis module"), "is designed to perform single cell RNA-Seq integration analysis. ezsinglecell uses ", a("Seurat", href = "https://www.cell.com/action/showPdf?pii=S0092-8674%2821%2900583-3", target = "_blank"), "for analyzing single cell RNA-Seq integration data. In this module, users can input raw counts matrix, H5 output from cellranger or a processed Seurat object. ezsinglecell allows users to download R Objects after processing.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p(strong("Single Cell Multiomics RNASeq Analysis module"), "is designed to perform single cell multiomics RNA-Seq analysis. ezsinglecell uses ", a("Seurat", href = "https://www.cell.com/action/showPdf?pii=S0092-8674%2821%2900583-3", target = "_blank"), "for analyzing single cell multiomics RNA-Seq data. In this module, users can input raw counts matrix, H5 output from cellranger or a processed Seurat object. ezsinglecell allows users to download R Objects after processing.", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px"),
                                   p(strong("Flow cytometry module"), "is designed to perform flow cytometry data analysis. ezsinglecell integrates state-of-the-art bioinformatics methods and in-house novel algorithms for data pre-processing, data visualization through linear or non-linear dimensionality reduction (UMAP/tSNE), cell clustering, automatic identification of cell subsets and inference of the relatedness between cell subsets."),
                                   p(strong("Imaging mass cytometry module"), "is designed to perform imaging mass cytometry data analysis. ezsinglecell is able to analyse high-dimensional imaging mass cytometry data to unravel the cellular composition, spatial architecture from the datasets, interactive visualization of cell subpopulations and progression profiles of key markers."),
                                   p(strong("Single cell ATAC-seq module"), "is designed to perform single cell ATAC-seq data analysis."),
                                   p(strong("Spatial Transcriptomics module"), "is designed to perform spatial transcriptomics data analysis."),
                               )),
                               tabPanel("What's New in ezsinglecell", fluidRow(
                                   p(strong("Changelog"), style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:25px"),
                                   p(strong("Version 1.0 Log:"), "ezsinglecell has 7 modules (single cell RNA-seq, Data Integration, Multiomics, Flow cytometry, Imaging mass cytometry, single cell ATAC-seq and Spatial Transcriptomics) and is built with an idea to empower researchers in performing bioinformatics analysis with little or no bioinformatics expertise", style="text-align:justify;color:black;background-color:white;padding:20px;border-radius:10px;font-size:15px")))
                   )
            )
        ),

        conditionalPanel(
            condition = "input.Module == 'Single Cell RNASeq Analysis'",
            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------QC plots---------

                               tabPanel("QC plots", fluidPage(
                                   hr(),
                                   tabsetPanel(id="QC",
                                               tabPanel(title="Violin Plots", value = "QC_panel1",
                                                        br(),
                                                        fluidRow(
                                                            column(6,
                                                                   plotlyOutput("nFeature_RNAPlot", width = "100%"),
                                                                   br(),
                                                                   plotlyOutput("mitoPlot", width = "100%"),
                                                                   br(),
                                                                   plotlyOutput("nCount_RNAPlot", width = "100%")
                                                            ),
                                                            #column(6,
                                                            #       verbatimTextOutput("name")
                                                            #)
                                                        )

                                               ),
                                               tabPanel(title="Feature Scatter Plots", value="QC_panel2",
                                                        br(),
                                                        fluidRow(
                                                            column(6,
                                                                   plotlyOutput("FeatureScatterPlot1", width = "100%"),
                                                                   br(),
                                                                   plotlyOutput("FeatureScatterPlot2", width = "100%")
                                                            ))
                                               )))
                               ),
                               #------Normalization and Variable Genes---------
                               tabPanel("Normalization and Variable Gene Plot", fluidPage(
                                   hr(),

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
                                              actionButton("findVarGenes", "Identify highly variable genes", icon = icon("hand-pointer-o")),
                                              #actionButton("doSCTransform", "Run SCTransform", icon = icon("hand-pointer-o"))
                                              # actionButton("doVarplot", "Plot variable genes", icon = icon("hand-pointer-o"))
                                       )),
                                   plotOutput("VarGenes", width = "100%")
                               )),
                               ##------PCA---------
                               tabPanel("PCA", fluidPage(
                                   hr(),
                                   tabsetPanel(id="Pca",
                                               tabPanel(title="PCA Plot", value="P_panel1",
                                                        br(),
                                                        fluidRow(
                                                            column(3,
                                                                   actionButton("doPCA", "Run PCA", icon = icon("hand-pointer-o"))
                                                            ),
                                                        ),
                                                        selectInput("assays1",
                                                                    label = "Normalize by:",
                                                                    choices = c("RNA", "SCT")
                                                        ),
                                                        plotlyOutput("PCA2DPlot", width = "100%"),
                                               ),
                                               tabPanel(title="PC Gene Visualisation", value="P_panel2",
                                                        br(),
                                                        selectInput("select.pc",
                                                                    label = "PC to plot",
                                                                    choices = c(1:20)
                                                        ),
                                                        fluidRow(
                                                            column(4,
                                                                   plotOutput("vizPlot", width = "100%", height = "600px")
                                                            ),
                                                            column(8,
                                                                   plotOutput("PCHeatmap", width = "100%", height = "600px")
                                                            )
                                                        ),
                                                        DT::dataTableOutput("PCtable")
                                               ),

                                               tabPanel(title="Elbow", value="P_panel4",
                                                        br(),
                                                        #actionButton("doElbow", label = "Get Elbow Plot"),
                                                        #br(),
                                                        br(),
                                                        plotOutput("Elbow", width = "100%")

                                               )
                                   )
                               )),
                               tabPanel("Clustering", fluidPage(
                                   hr(),
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
                                   plotlyOutput("Cluster2DPlot_1", width = "100%")
                               )),

                               tabPanel("UMAP", fluidPage(
                                   hr(),
                                   fluidRow(
                                       column(3,
                                              numericInput("dim.used",
                                                           label = "Dimensions used",
                                                           value = 10)
                                       ),
                                       br(),
                                       column(3,
                                              br(),
                                              actionButton("doUmap", "Run UMAP", icon = icon("hand-pointer-o")),
                                              textOutput("Umap.done"),
                                              br()
                                       )),
                                   br(),
                                   plotlyOutput("Umap_2d_plot_1", width = "100%")

                               )),
                               tabPanel("TSNE", fluidPage(
                                   hr(),
                                   fluidRow(
                                       column(3,
                                              numericInput("dim.used",
                                                           label = "Dimensions used",
                                                           value = 10)
                                       ),
                                       br(),
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
                                   plotlyOutput("Tsne_2d_plot_1", width = "100%")

                               )),
                               tabPanel("Cell type identification", fluidPage(
                                   hr(),
                                   column(3,
                                          br(),
                                          actionButton("doCELLiD", "Run CELLiD", icon = icon("hand-pointer-o")),
                                          textOutput("CELLiD.done"),
                                          br()
                                   )),
                                   br(),
                                   plotlyOutput("Umap_cellid", width = "100%")

                               ),
                               tabPanel("DEGs", fluidPage(
                                   hr(),
                                   fluidRow(

                                       column(3, selectInput("min_pct",
                                                             label = "min.pct",
                                                             choices = c("0.1", "0.25"))
                                       ),

                                       column(3, selectInput("logfc",
                                                             label = "logfc.threshold",
                                                             choices = c("0.1", "0.25"))
                                       ),

                                       column(3, selectInput("test.use",
                                                             label = "Test use",
                                                             choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
                                       ),

                                       column(4,
                                              actionButton("doDeg", "Run DEGs", icon = icon("hand-pointer-o"))
                                       )),
                                   br(),
                                   column(9,
                                          DT::dataTableOutput("Deg.table"),
                                          br(),
                                          plotlyOutput("Deg3.plot", width = "100%")
                                   ))
                                   #)
                               ),

                               tabPanel("Data visualization", fluidPage(
                                   fluidRow(
                                       column(6,
                                              uiOutput("deg.gene.select"),
                                              plotlyOutput("Deg.plot", width = "100%"),
                                              br(),
                                              plotlyOutput("Deg1.plot", width = "100%")
                                       ),
                                   ))),

                               tabPanel("Trajectory", fluidPage(
                                   fluidRow(
                                       column(3,
                                              br(),
                                              actionButton("doMonocle3", "Run Monocle3", icon = icon("hand-pointer-o")),
                                              textOutput("Trajectory.done"),
                                              br()
                                       )),
                                   br(),
                                   plotlyOutput("Monocle3_plot", width = "100%")
                               )),

                               tabPanel("Cell-cell communication", fluidPage(
                                   br(),
                                   tabsetPanel(id="Cell-cell communication",
                                               tabPanel(title="Aggregated cell-cell communication network",
                                                        fluidRow(
                                                            column(3,
                                                                   br(),
                                                                   actionButton("doCC", "Run Cell-cell communication", icon = icon("hand-pointer-o")),
                                                                   textOutput("CC.done"),
                                                                   br()
                                                            )),
                                                        br(),
                                                        plotOutput("CC_plot1", width = "50%"),
                                                        br(),
                                                        plotOutput("CC_plot2", width = "50%")
                                               ),

                                               tabPanel(title="Visualize each signaling pathway",
                                                        fluidRow(
                                                            column(3,
                                                                   br(),
                                                                   actionButton("doCC1", "Visualize each signaling pathway", icon = icon("hand-pointer-o")),
                                                                   textOutput("CC1.done"),
                                                                   br()
                                                            )),
                                                        br(),
                                                        uiOutput("CC.gene.select"),
                                                        column(3, selectInput("layout_cc",
                                                                              label = "Layout",
                                                                              choices = c("hierarchy", "circle", "chord"))
                                                        ),
                                                        br(),
                                                        plotOutput("CC_plot3", width = "50%"),
                                                        br(),
                                                        plotOutput("CC_plot4", width = "50%")

                                               ),

                                               tabPanel(title="Contribution of each ligand-receptor pair",
                                                        fluidRow(
                                                            plotOutput("CC_plot5", width = "50%"),
                                                            br(),
                                                            column(3,
                                                                   br(),
                                                                   actionButton("doCC2", "Contribution of each ligand-receptor pair", icon = icon("hand-pointer-o")),
                                                                   textOutput("CC2.done"),
                                                                   br()
                                                            )),
                                                        br(),
                                                        uiOutput("LR.gene.select"),
                                                        column(3, selectInput("layout_cc1",
                                                                              label = "Layout",
                                                                              choices = c("hierarchy", "circle", "chord"))
                                                        ),
                                                        plotOutput("CC_plot6", width = "50%"),
                                                        br(),
                                                        uiOutput("cell_group.select"),
                                                        br(),
                                                        plotOutput("CC_plot7", width = "50%"),
                                                        br(),
                                                        plotOutput("CC_plot8", width = "50%")
                                               ),
                                               tabPanel(title="Plot the gene expression",
                                                        br(),
                                                        plotOutput("CC_plot9", width = "50%"),
                                               ),
                                               tabPanel(title="Other plots",
                                                        br(),
                                                        fluidRow(
                                                            column(3,
                                                                   br(),
                                                                   actionButton("doCC3", "Visualize each signaling pathway", icon = icon("hand-pointer-o")),
                                                                   textOutput("CC3.done"),
                                                            )),
                                                        br(),
                                                        plotOutput("CC_plot10", width = "50%"),
                                                        br(),
                                                        plotOutput("CC_plot11", width = "50%"),
                                                        br(),
                                                        plotOutput("CC_plot12", width = "50%"),
                                                        br(),
                                                        plotOutput("CC_plot13", width = "50%")
                                               ),
                                   ))),

                               tabPanel("Differential abundance", fluidPage(
                                   br(),
                                   tabsetPanel(id="Milo analysis",
                                               tabPanel(title="Visualize the data",
                                                        fluidRow(
                                                            column(3,
                                                                   br(),
                                                                   selectInput("group.by_milo",
                                                                               label = "Group by",
                                                                               choices = c("Genotype", "Day")),
                                                            ),
                                                            column(3,
                                                                   br(),
                                                                   actionButton("do_milo_preprocess", "Visualize", icon = icon("hand-pointer-o")),
                                                                   br()
                                                            )),
                                                        br(),
                                                        plotOutput("vis_milo_plot1", width = "50%")
                                                        #br(),
                                                        #plotOutput("vis_milo_plot2", width = "50%")
                                               ),

                                               tabPanel(title="Milo analysis",
                                                        fluidRow(

                                                            column(3,
                                                                   br(),
                                                                   selectInput("dim.used_milo",
                                                                               label = "Dimensions to use",
                                                                               choices = c(10:50)),
                                                            ),
                                                            br(),
                                                            br(),
                                                            actionButton("do_milo", "Analyse", icon = icon("hand-pointer-o")),
                                                            br()
                                                        ),
                                                        br(),
                                                        plotOutput("vis_milo_plot3", width = "50%"),
                                                        br(),
                                                        plotOutput("vis_milo_plot4", width = "50%")
                                                        #br(),
                                                        #plotOutput("vis_milo_plot5", width = "50%")
                                                        #br(),
                                                        #plotOutput("vis_milo_plot6", width = "50%"),
                                                        #br(),
                                                        #plotOutput("vis_milo_plot7", width = "50%")
                                               ),
                                   ),

                               )),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        "TPM file (Human): ", tags$a(href="https://drive.google.com/file/d/1i78xdIhsUOj6s7czNOMTuFJRDowb5M51/view?usp=sharing", "Human TPM file"),
                                        br(),
                                        br(),
                                        "Metadata file (Human): ", tags$a(href="https://drive.google.com/file/d/1Jb2noiAKeXLc96K0f7uZ3Rr3Y_PK3gry/view?usp=sharing", "Human Metadata file"),
                                        br(),
                                        br(),
                                        "Cellranger H5 output file for human: ", tags$a(href="https://drive.google.com/file/d/1A1bp8pWAR4Y2qQpFrMM3XRFckNkeY1mg/view?usp=sharing", "Human H5 file"),
                                        br(),
                                        br(),
                                        "Cellranger H5 output file for mouse: ", tags$a(href="https://drive.google.com/file/d/1EsNoJg7_E86KQrRwXPASGNaS4r1Jtzqw/view?usp=sharing", "Mouse H5 file")),







                               ##------END---------
                   ))
        ),

        conditionalPanel(
            condition = "input.Module == 'Single cell data integration Analysis' & input.scAnalysis_integ == 'Seurat'",
            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------Data Integration using Seurat---------

                               tabPanel("Data Integration and PCA", fluidPage(
                                   hr(),
                                   fluidRow(
                                       column(3,
                                              br(),
                                              actionButton("doIntg_seurat", "Running Data Integration", icon = icon("hand-pointer-o")),
                                              br(),
                                       )),

                                   column(9,
                                          br(),
                                          actionButton("runPCA_intg_seurat", "Run PCA", icon = icon("hand-pointer-o")),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("PCAplot_seurat_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_seurat_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_seurat_tpm3", width = "100%")),

                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("PCAplot_seurat_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_seurat_h5_2", width = "100%"))),
                               )),
                               tabPanel("UMAP", fluidPage(
                                   column(9,
                                          br(),
                                          actionButton("runUMAP_intg_seurat", "Run UMAP", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("UMAPplot_seurat_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_seurat_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_seurat_tpm3", width = "100%")),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("UMAPplot_seurat_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_seurat_h5_2", width = "100%")),
                                   ))),

                               tabPanel("TSNE", fluidPage(
                                   column(9,
                                          br(),
                                          actionButton("runTSNE_intg_seurat", "Run TSNE", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("TSNEplot_seurat_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_seurat_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_seurat_tpm3", width = "100%")),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("TSNEplot_seurat_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_seurat_h5_2", width = "100%")),
                                   ))),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(3, "Human TPM data:",
                                               br(),
                                               br(),
                                               "TPM Matrix File 1 : ", tags$a(href="https://drive.google.com/file/d/11-PaWmvOTytKNF7NIUfrkN0I-M-ANiYl/view?usp=sharing", "TPM Matrix File 1"),
                                               br(),
                                               br(),
                                               "TPM Matrix File 2 : ", tags$a(href="https://drive.google.com/file/d/1l7elllspoZ95xuGfw4NKIBPRzVgy4bp1/view?usp=sharing", "TPM Matrix File 2"),
                                               br(),
                                               br(),
                                               "TPM Matrix File 3: ", tags$a(href="https://drive.google.com/file/d/1fN4hb6OVg_U8FCVg-kqE_GR2g9Ooqnu1/view?usp=sharing", "TPM Matrix File 3"),
                                               br(),
                                               br(),
                                               "Metadata File 1 : ", tags$a(href="https://drive.google.com/file/d/1Qv4Mhheog-fgvR_xXjSltrOyXna5Xrqt/view?usp=sharing", "Metadata File 1"),
                                               br(),
                                               br(),
                                               "Metadata File 2 : ", tags$a(href="https://drive.google.com/file/d/1zwhq8n8MASJjsoRInuIf51pZA82Yv829/view?usp=sharing", "Metadata File 2"),
                                               br(),
                                               br(),
                                               "Metadata File 3: ", tags$a(href="https://drive.google.com/file/d/19tlKX7X1m9rfYgSSrOLOxhmPdcvV414M/view?usp=sharing", "Metadata File 3")),
                                        br(),
                                        br(),
                                        column(4, "Human cellranger output:",
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 1): ", tags$a(href="https://drive.google.com/file/d/1KkvEd4x_Q-9ITboj7ZC3EV7PSllxh8Rn/view?usp=sharing", "Cellranger output 1"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 2): ", tags$a(href="https://drive.google.com/file/d/1mEJS8XqieZO_T4MQRPdtdZ1FV_hgg1_e/view?usp=sharing", "Cellranger output 2"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 3): ", tags$a(href="https://drive.google.com/file/d/14wYhoNFJP09I5aOaopjqX_xEP3TNgujI/view?usp=sharing", "Cellranger output 3")),
                                        br(),
                                        br(),
                                        column(4, "Mouse cellranger output:",
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 1): ", tags$a(href="https://drive.google.com/file/d/1GkvtSzjIurMsCP5nOZIMv6iMDpajJErs/view?usp=sharing", "Cellranger output 1"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 2): ", tags$a(href="https://drive.google.com/file/d/1Ylj3SZblbWplk9D7eLY238JgdY8GQm7w/view?usp=sharing", "Cellranger output 2"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 3): ", tags$a(href="https://drive.google.com/file/d/1qgQfrBjC-X_wRoPvb8QRfxqAzavqtVRS/view?usp=sharing", "Cellranger output 3"))),
                               #))
                               #))
                   ))),

        conditionalPanel(
            condition = "input.Module == 'Single cell data integration Analysis' & input.scAnalysis_integ == 'Harmony'",
            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------Data Integration using Harmony---------

                               tabPanel("Data Integration and PCA", fluidPage(
                                   hr(),
                                   fluidRow(
                                       column(3,
                                              br(),
                                              actionButton("doIntg_harmony", "Running Data Integration", icon = icon("hand-pointer-o")),
                                              br(),
                                       )),

                                   column(9,
                                          br(),
                                          actionButton("runPCA_intg_harmony", "Run PCA", icon = icon("hand-pointer-o")),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("PCAplot_harmony_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_harmony_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_harmony_tpm3", width = "100%")),

                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("PCAplot_harmony_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_harmony_h5_2", width = "100%"))),
                               )),
                               tabPanel("UMAP", fluidPage(
                                   column(9,
                                          br(),
                                          actionButton("runUMAP_intg_harmony", "Run UMAP", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("UMAPplot_harmony_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_harmony_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_harmony_tpm3", width = "100%")),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("UMAPplot_harmony_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_harmony_h5_2", width = "100%")),
                                   ))),

                               tabPanel("TSNE", fluidPage(
                                   column(9,
                                          br(),
                                          actionButton("runTSNE_intg_harmony", "Run TSNE", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("TSNEplot_harmony_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_harmony_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_harmony_tpm3", width = "100%")),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("TSNEplot_harmony_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_harmony_h5_2", width = "100%")),
                                   ))),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(3, "Human TPM data:",
                                               br(),
                                               br(),
                                               "TPM Matrix File 1 : ", tags$a(href="https://drive.google.com/file/d/11-PaWmvOTytKNF7NIUfrkN0I-M-ANiYl/view?usp=sharing", "TPM Matrix File 1"),
                                               br(),
                                               br(),
                                               "TPM Matrix File 2 : ", tags$a(href="https://drive.google.com/file/d/1l7elllspoZ95xuGfw4NKIBPRzVgy4bp1/view?usp=sharing", "TPM Matrix File 2"),
                                               br(),
                                               br(),
                                               "TPM Matrix File 3: ", tags$a(href="https://drive.google.com/file/d/1fN4hb6OVg_U8FCVg-kqE_GR2g9Ooqnu1/view?usp=sharing", "TPM Matrix File 3"),
                                               br(),
                                               br(),
                                               "Metadata File 1 : ", tags$a(href="https://drive.google.com/file/d/1Qv4Mhheog-fgvR_xXjSltrOyXna5Xrqt/view?usp=sharing", "Metadata File 1"),
                                               br(),
                                               br(),
                                               "Metadata File 2 : ", tags$a(href="https://drive.google.com/file/d/1zwhq8n8MASJjsoRInuIf51pZA82Yv829/view?usp=sharing", "Metadata File 2"),
                                               br(),
                                               br(),
                                               "Metadata File 3: ", tags$a(href="https://drive.google.com/file/d/19tlKX7X1m9rfYgSSrOLOxhmPdcvV414M/view?usp=sharing", "Metadata File 3")),
                                        br(),
                                        br(),
                                        column(4, "Human cellranger output:",
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 1): ", tags$a(href="https://drive.google.com/file/d/1KkvEd4x_Q-9ITboj7ZC3EV7PSllxh8Rn/view?usp=sharing", "Cellranger output 1"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 2): ", tags$a(href="https://drive.google.com/file/d/1mEJS8XqieZO_T4MQRPdtdZ1FV_hgg1_e/view?usp=sharing", "Cellranger output 2"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 3): ", tags$a(href="https://drive.google.com/file/d/14wYhoNFJP09I5aOaopjqX_xEP3TNgujI/view?usp=sharing", "Cellranger output 3")),
                                        br(),
                                        br(),
                                        column(4, "Mouse cellranger output:",
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 1): ", tags$a(href="https://drive.google.com/file/d/1GkvtSzjIurMsCP5nOZIMv6iMDpajJErs/view?usp=sharing", "Cellranger output 1"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 2): ", tags$a(href="https://drive.google.com/file/d/1Ylj3SZblbWplk9D7eLY238JgdY8GQm7w/view?usp=sharing", "Cellranger output 2"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 3): ", tags$a(href="https://drive.google.com/file/d/1qgQfrBjC-X_wRoPvb8QRfxqAzavqtVRS/view?usp=sharing", "Cellranger output 3"))),
                               #))
                               #))
                   ))),

        conditionalPanel(
            condition = "input.Module == 'Single cell data integration Analysis' & input.scAnalysis_integ == 'RLiger'",
            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(3, "Human TPM data:",
                                               br(),
                                               br(),
                                               "TPM Matrix File 1 : ", tags$a(href="https://drive.google.com/file/d/11-PaWmvOTytKNF7NIUfrkN0I-M-ANiYl/view?usp=sharing", "TPM Matrix File 1"),
                                               br(),
                                               br(),
                                               "TPM Matrix File 2 : ", tags$a(href="https://drive.google.com/file/d/1l7elllspoZ95xuGfw4NKIBPRzVgy4bp1/view?usp=sharing", "TPM Matrix File 2"),
                                               br(),
                                               br(),
                                               "TPM Matrix File 3: ", tags$a(href="https://drive.google.com/file/d/1fN4hb6OVg_U8FCVg-kqE_GR2g9Ooqnu1/view?usp=sharing", "TPM Matrix File 3"),
                                               br(),
                                               br(),
                                               "Metadata File 1 : ", tags$a(href="https://drive.google.com/file/d/1Qv4Mhheog-fgvR_xXjSltrOyXna5Xrqt/view?usp=sharing", "Metadata File 1"),
                                               br(),
                                               br(),
                                               "Metadata File 2 : ", tags$a(href="https://drive.google.com/file/d/1zwhq8n8MASJjsoRInuIf51pZA82Yv829/view?usp=sharing", "Metadata File 2"),
                                               br(),
                                               br(),
                                               "Metadata File 3: ", tags$a(href="https://drive.google.com/file/d/19tlKX7X1m9rfYgSSrOLOxhmPdcvV414M/view?usp=sharing", "Metadata File 3")),
                                        br(),
                                        br(),
                                        column(4, "Human cellranger output:",
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 1): ", tags$a(href="https://drive.google.com/file/d/1KkvEd4x_Q-9ITboj7ZC3EV7PSllxh8Rn/view?usp=sharing", "Cellranger output 1"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 2): ", tags$a(href="https://drive.google.com/file/d/1mEJS8XqieZO_T4MQRPdtdZ1FV_hgg1_e/view?usp=sharing", "Cellranger output 2"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 3): ", tags$a(href="https://drive.google.com/file/d/14wYhoNFJP09I5aOaopjqX_xEP3TNgujI/view?usp=sharing", "Cellranger output 3")),
                                        br(),
                                        br(),
                                        column(4, "Mouse cellranger output:",
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 1): ", tags$a(href="https://drive.google.com/file/d/1GkvtSzjIurMsCP5nOZIMv6iMDpajJErs/view?usp=sharing", "Cellranger output 1"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 2): ", tags$a(href="https://drive.google.com/file/d/1Ylj3SZblbWplk9D7eLY238JgdY8GQm7w/view?usp=sharing", "Cellranger output 2"),
                                               br(),
                                               br(),
                                               "Cellranger output (Sample 3): ", tags$a(href="https://drive.google.com/file/d/1qgQfrBjC-X_wRoPvb8QRfxqAzavqtVRS/view?usp=sharing", "Cellranger output 3"))),

                               ##------Data Integration using RLiger---------

                               tabPanel("Data Integration and PCA", fluidPage(
                                   hr(),
                                   fluidRow(
                                       column(3,
                                              br(),
                                              actionButton("doIntg_rliger", "Running Data Integration", icon = icon("hand-pointer-o")),
                                              br(),
                                       )),

                                   column(9,
                                          br(),
                                          actionButton("runPCA_intg_rliger", "Run PCA", icon = icon("hand-pointer-o")),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("PCAplot_rliger_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_rliger_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_rliger_tpm3", width = "100%")),

                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("PCAplot_rliger_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("PCAplot_rliger_h5_2", width = "100%"))),
                               )),
                               tabPanel("UMAP", fluidPage(
                                   column(9,
                                          br(),
                                          actionButton("runUMAP_intg_rliger", "Run UMAP", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("UMAPplot_rliger_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_rliger_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_rliger_tpm3", width = "100%")),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("UMAPplot_rliger_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot_rliger_h5_2", width = "100%")),
                                   ))),

                               tabPanel("TSNE", fluidPage(
                                   column(9,
                                          br(),
                                          actionButton("runTSNE_intg_rliger", "Run TSNE", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'Raw Counts Matrix'",
                                              plotlyOutput("TSNEplot_rliger_tpm1", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_rliger_tpm2", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_rliger_tpm3", width = "100%")),
                                          br(),
                                          conditionalPanel(
                                              condition = "input.Module == 'Single cell data integration Analysis' & input.scInput1 == 'H5'",
                                              plotlyOutput("TSNEplot_rliger_h5_1", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot_rliger_h5_2", width = "100%")),
                                   ))),

                   ))),

        conditionalPanel(
            condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'CITE-seq'",
            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------CITE-seq (Seurat)---------

                               tabPanel("Dimension Reduction (UMAP)", fluidPage(
                                   hr(),

                                   column(3,
                                          numericInput("clus.res1",
                                                       label = "Resolution used",
                                                       value = 0.6,
                                                       min = 0.1,
                                                       step = 0.1)
                                   ),

                                   column(3,
                                          selectInput("dim.used1",
                                                      label = "PC to plot",
                                                      choices = c(5:50)),
                                   ),
                                   fluidRow(
                                       column(9,
                                              br(),
                                              actionButton("runUMAP3", "Run UMAP", icon = icon("hand-pointer-o")),
                                              #textOutput("Intg.done"),
                                              br(),
                                              #),
                                              #column(6,
                                              br(),
                                              plotlyOutput("UMAPplot4_a", width = "100%")),
                                   ),
                               )),


                               tabPanel("Dimension Reduction (TSNE)", fluidPage(
                                   hr(),
                                   column(9,
                                          br(),
                                          actionButton("runTSNE3", "Run TSNE", icon = icon("hand-pointer-o")),
                                          #textOutput("Intg.done"),
                                          br(),
                                          br(),
                                          plotlyOutput("TSNEplot4_a", width = "100%"),
                                          br(),
                                   ),
                               )),
                               tabPanel("Data visualization", fluidPage(
                                   hr(),
                                   column(4,
                                          actionButton("Vis3", "Visualize", icon = icon("hand-pointer-o"))
                                   ),

                                   fluidRow(
                                       column(6,
                                              uiOutput("vis.gene.select"),
                                              plotlyOutput("vis.plot", width = "100%"),
                                              br(),
                                              uiOutput("vis.gene.select1"),
                                              plotlyOutput("vis1.plot", width = "100%"),
                                       ),

                                   ))),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(4, "CITE-Seq:",
                                               br(),
                                               br(),
                                               "Cell ranger output: ", tags$a(href="https://drive.google.com/file/d/16MjQIIFqvgZoL3KhKre4kC6hoIWhsKIE/view?usp=sharing", "CITE-seq cell ranger output")),
                               ),
                   ))),

        conditionalPanel(
            condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_mult == 'Seurat' & input.scAnalysis_type == 'mRNA+scATAC-seq'",

            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               tabPanel("Quality control", fluidPage(
                                   hr(),
                                   fluidRow(

                                       column(9,
                                              plotlyOutput("nCount_ATAC.plot", width = "100%"),
                                              br(),
                                              plotlyOutput("nCount_RNA.plot", width = "100%"),
                                              br(),
                                              plotlyOutput("percent.mt.plot", width = "100%")),
                                   ),
                               )),

                               tabPanel("UMAP", fluidPage(

                                   fluidRow(
                                       column(5,
                                              br(),
                                              numericInput("clus.res2",
                                                           label = "Resolution used",
                                                           value = 0.6,
                                                           min = 0.1,
                                                           step = 0.1)
                                       ),


                                       column(4,
                                              br(),
                                              selectInput("dim.used2",
                                                          label = "PC to use",
                                                          choices = c(5:50))
                                       )),

                                   fluidRow(
                                       column(9,
                                              actionButton("doSCTransform_multi", "Run scTransform", icon = icon("hand-pointer-o")),
                                              br(),
                                              br(),
                                              actionButton("runUMAP4", "Run UMAP", icon = icon("hand-pointer-o")),
                                              #textOutput("Intg.done"),
                                              br(),
                                              #),
                                              #column(6,
                                              br(),
                                              plotlyOutput("UMAPplot5_a", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot5_b", width = "100%"),
                                              br(),
                                              plotlyOutput("UMAPplot5_c", width = "100%")),
                                   ),
                               )),

                               tabPanel("TSNE", fluidPage(
                                   fluidRow(
                                       column(9,
                                              br(),
                                              actionButton("runTSNE4", "Run TSNE", icon = icon("hand-pointer-o")),
                                              #textOutput("Intg.done"),
                                              br(),
                                              br(),
                                              plotlyOutput("TSNEplot5_a", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot5_b", width = "100%"),
                                              br(),
                                              plotlyOutput("TSNEplot5_c", width = "100%")),

                                       br(),
                                   ))),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(4, "mRNA + scATAC-seq:",
                                               br(),
                                               br(),
                                               "Cell ranger output: ", tags$a(href="https://drive.google.com/file/d/1JtZNl_euNHVJgzwBlRlKrW6vyVKPlLoP/view?usp=sharing", "Cell ranger output"),
                                               br(),
                                               br(),
                                               "Fragments file: ",  tags$a(href="https://drive.google.com/file/d/1aPK9EsROxO3RFnTC39K1UfXqBvn8DfDI/view?usp=sharing", "Fragments file"),
                                               br(),
                                               br(),
                                               "Fragments index file: ",  tags$a(href="https://drive.google.com/file/d/1DxMsbyMcIi8K2nmxL3GvGtgQLWDEyaUI/view?usp=sharing", "Fragments index file")),
                               ),


                   ))),

        conditionalPanel(
            condition = "input.Module == 'Single cell multiomics Analysis' & input.scAnalysis_mult == 'MOFA' & input.scAnalysis_type == 'CITE-seq'",

            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               tabPanel("Data overview", fluidPage(

                                   fluidRow(
                                       br(),
                                       column(3,
                                              selectInput("factor",
                                                          label = "Number of factors",
                                                          choices = c(2:50))),
                                       #column(3,
                                       #       selectInput("nfeatures",
                                       #                   label = "Number of features",
                                       #                   choices = c(2:15))),
                                       #column(3,
                                       #       selectInput("view",
                                       #                   label = "View",
                                       #                   choices = c("RNA", "ADT"))),

                                       #  column(3,
                                       #       selectInput("factors",
                                       #                   label = "Color by",
                                       #                   choices = c("Phenotype", "Combined"))),

                                       column(3,
                                              br(),
                                              actionButton("data_overview", "Plot Data Overview", icon = icon("hand-pointer-o")),
                                       )),
                                   br(),
                                   plotOutput("mofa_a", width = "50%"),
                                   br(),
                                   plotOutput("mofa_b", width = "100%")
                                   #br(),
                                   #plotOutput("factor_a", width = "100%"),
                                   #br(),
                                   #plotOutput("factor_b", width = "100%")

                                   #br(),
                                   #plotOutput("mofa_c", width = "100%")
                               )),

                               tabPanel("UMAP", fluidPage(

                                   fluidRow(
                                       br(),
                                       column(3,
                                              selectInput("n_neighbors",
                                                          label = "Number of neighbours",
                                                          choices = c(5:50))),
                                       column(3,
                                              br(),
                                              actionButton("umap_mofa", "Run UMAP", icon = icon("hand-pointer-o")),
                                       )),
                                   plotOutput("umap_mofa_mult", width = "100%"),

                               )),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(4, "CITE-Seq:",
                                               br(),
                                               br(),
                                               "Cell ranger output: ", tags$a(href="https://drive.google.com/file/d/11Nxq1eE7f_QFQxAkxydTqg-d6USTDte3/view?usp=sharing", "CITE-seq cell ranger output")),
                               ),
                   ))),

        conditionalPanel(
            condition = "input.Module == 'Flow cytometry Analysis'",

            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               tabPanel("Data Analysis",

                                        fluidPage(
                                            br(),
                                            box(
                                                title = "Load data",
                                                status = "primary",
                                                width = 3,
                                                solidHeader = TRUE,
                                                collapsible = FALSE,
                                                collapsed = FALSE,
                                                fluidPage(
                                                    tags$div(title='The fcs files to be analyzed. One or multiple fcs files are allowed. When multiple fcs files are selected, cells from each fcs file are combined for analysis.'
                                                             , fileInput("rawfcs", "Raw FCS files", multiple = TRUE
                                                                         , accept = c("FCSfile/fcs", '.fcs'))
                                                             , fileInput("sample_anno", "Meta data", multiple = FALSE, accept = c("Txtfile/txt", '.txt'))
                                                    )
                                                    , tags$div(title="Select the list of makers to be used for analysis."
                                                               , selectInput('markers', 'Markers', choices = NULL, selected = NULL
                                                                             , selectize = TRUE
                                                                             , multiple = TRUE)
                                                    )
                                                    , tags$div(title="A prefix that will be added to the names of result files."
                                                               , textInput('project_name', 'Project name', value = 'cytofkit')
                                                    )
                                                    , tags$div(title="When multiple fcs files are selected, cell expression data can be merged using one of the four different methods including \"ceil\",\"all\", \"min\",\"fixed\". \n\n\"ceil\" (the default option): up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each fcs file and combined for analysis. \n\n\"all\": all cells from each fcs file are combined for analysis. \n\n\"min\": The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. \n\n\"fixed\": a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each fcs file and combined for analysis."
                                                               , selectInput('merge_method', 'Merge method', choices = c('all', 'min', 'ceil', 'fixed')
                                                                             , selected = "ceil")
                                                    )
                                                    , tags$div(title="Up to fixedNum of cells from each fcs file are used for analysis."
                                                               , numericInput('fix_number', 'Fixed Number', value = 5000)
                                                    )
                                                    , tags$div(title="Data Transformation method, including \"cytofAsinh\"(Customized Asinh transformation for CyTOF data), \"autoLgcl\"(automatic logicle transformation for CyTOF data), \"logicle\"(customize your own parameters for logicle transformation) and \"none\"(if your data is already transformed)."
                                                               , selectInput('transform_method', 'Tranformation Method'
                                                                             , choices = c('autoLgcl', 'cytofAsinh', 'logicle', 'arcsinh', 'none'))
                                                    )
                                                    , fluidRow(
                                                        tags$div(title=''
                                                                 , downloadButton('download_analysis_res', 'Save result'))
                                                    )
                                                    , br()
                                                    , fluidRow(actionButton('reset', 'Reset')
                                                               , actionButton('submit', 'Submit')
                                                               , actionButton('quit', 'Quit'))
                                                )),
                                            box(
                                                title = "Dimensionality Reduction",
                                                status = "primary",
                                                width = 3,
                                                height = 630,
                                                # background = "orange",
                                                solidHeader = TRUE,
                                                collapsible = FALSE,
                                                collapsed = FALSE,
                                                fluidPage(
                                                    tags$div(title="The method(s) used for visualizing the clustering results, multiple selections are allowed. Including \"pca\", \"isomap\", \"tsne\". \n\nWARNING: \"tsne\" is the default selection, \"isomap\" may take long time."
                                                             , selectInput('dr_method', 'Dimensionality reduction methods'
                                                                           , choices = c('PCA', 'isoMAP', 'tSNE', 'UMAP')
                                                                           , multiple = TRUE
                                                                           , selected = 'tSNE')
                                                    )
                                                    , tags$div(title=''
                                                               , numericInput('tsne_perplexity', 'tSNE perplexity', value = 30)
                                                    )
                                                    , tags$div(title=''
                                                               , numericInput('tsne_interation', 'tSNE Max Iterations', value = 1000)
                                                    )
                                                    , tags$div(title=''
                                                               , numericInput('seed', 'Seed', value = 42))
                                                )),
                                            box(
                                                title = "Clustering",
                                                status = "primary",
                                                width = 3,
                                                height = 630,
                                                # background = "orange",
                                                solidHeader = TRUE,
                                                collapsible = FALSE,
                                                collapsed = FALSE,
                                                fluidPage(
                                                    tags$div(title="The method(s) for clustering, including \"DensVM\", \"ClusterX\", \"Rphenograph\", and \"FlowSOM\". "
                                                             , selectInput('cluster_method', 'Cluster Method(s)'
                                                                           , choices = c('Rphenograph', 'ClusterX', 'DensVM', 'FlowSOM')
                                                                           , multiple = TRUE
                                                                           , selected = c('Rphenograph'))
                                                    )
                                                    , tags$div(title="Number of nearest neighbours to pass to Rphenograph."
                                                               , numericInput('rphenograph_k', 'Rphenograph_k', value = 30)
                                                    )
                                                    , tags$div(title="Number of clusters for meta clustering in FlowSOM."
                                                               , numericInput('flowsom_k', 'FlowSOM_k', value = 40))
                                                )),
                                            box(
                                                title = "Pseudo-time",
                                                status = "primary",
                                                width = 3,
                                                height = 630,
                                                # background = "orange",
                                                solidHeader = TRUE,
                                                collapsible = FALSE,
                                                collapsed = FALSE,
                                                fluidPage(
                                                    tags$div(title="The method used for cellular progression analysis including \"diffusion map\" and \"isomap\"\n\nIf \"NULL\" was selected, no progression estimation will be performed."
                                                             , selectInput('progressionMethods', 'Cellular Progression'
                                                                           , choices = c('NULL', 'diffusionmap', 'isomap')))
                                                )),

                                        )),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(4,
                                               "Example 1 : ", tags$a(href="https://github.com/JinmiaoChenLab/cytofkit2/tree/master/inst/extdata", "Example data 1")),
                                        br(),
                                        br(),
                                        column(4,
                                               "Example 2 : ", tags$a(href="https://drive.google.com/drive/folders/1EM79nsoRey8ZU6WklyihUmJOFZoLxhIC?usp=sharing", "Example data 2")),
                               ),
                   ))
        ),

        conditionalPanel(
            condition = "input.Module == 'Imaging mass cytometry Analysis'",

            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------Plot cells---------
                               tabPanel("Cell Segmentation", fluidPage(
                                   br(),
                                   fluidRow(
                                       column(10,
                                              plotOutput("plot_cell1", width = "100%")
                                       )
                                       #column(6,
                                       #       verbatimTextOutput("name")
                                       #)
                                   ),

                               )),

                               tabPanel("Cell type classification", fluidPage(
                                   br(),
                                   fluidRow(
                                       column(10,
                                              plotOutput("plot_cell3", width = "100%")
                                       )
                                   ),

                               )),

                               tabPanel("Marker expression", fluidPage(
                                   br(),
                                   fluidRow(
                                       column(10,
                                              uiOutput("sce.select"),
                                              plotOutput("plot_cell5", width = "100%"),
                                              br(),
                                              plotOutput("plot_cell4", width = "100%")
                                       )
                                   ),

                               )),

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(4,
                                               "Single cell data : ", tags$a(href="https://drive.google.com/file/d/1XgY3WC5ykJkNt5RdGtAq5LMDweugSTEq/view?usp=sharing", "Single cell data")),
                                        br(),
                                        br(),
                                        column(4,
                                               "Cell mask : ", tags$a(href="https://drive.google.com/file/d/1sWDsaD1i7hRWEkD5NB4814AjL2tfBOj4/view?usp=sharing", "Cell mask")),
                               ),
                   ))),

        conditionalPanel(
            condition = "input.Module == 'Single cell ATAC-seq Analysis'",

            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------ATAC QC plots---------
                               tabPanel("Computing QC Metrics", fluidPage(

                                   tabPanel(title="ATAC-seq QC Plot", value = "ATAC-seq_QC_panel1",
                                            br(),
                                            fluidRow(

                                                column(10,
                                                       plotlyOutput("TSSPlot", width = "100%"),
                                                       br(),
                                                       plotlyOutput("FragmentHistogram", width = "100%"),
                                                       #br(),
                                                       #plotlyOutput("VlnPlot_atac", width = "100%"),
                                                )))
                               )),

                               ##------Normalization---------
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
                                       br(),
                                       selectInput("dim.used_atac",
                                                   label = "PC to plot",
                                                   choices = c(2:50)
                                       ),
                                       column(3,
                                              actionButton("doCluster_ATAC", "Run Clustering", icon = icon("hand-pointer-o")),
                                              textOutput("cluster_atac.done"),
                                              br()
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
                                                           value = 10)
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
                                                           value = 10)
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

                               tabPanel("Test Data",
                                        br(),
                                        br(),
                                        column(4,
                                               "Peak/Cell matrix : ", tags$a(href="https://drive.google.com/file/d/1W-_7GkVG6h_VDcFo7erSYb7OmCmRsUAh/view?usp=sharing", "Peak/Cell matrix")),
                                        br(),
                                        br(),
                                        column(4,
                                               "Fragment file : ", tags$a(href="https://drive.google.com/file/d/1M6rMCM5usMIWMplIB1_K8ng0Us3kfIYS/view?usp=sharing", "Fragment file")),
                                        br(),
                                        br(),
                                        column(4,
                                               "Fragment index file : ", tags$a(href="https://drive.google.com/file/d/1btpXCtekZGWnIOTQA68JIRS53iO5o87D/view?usp=sharing", "Fragments index file")),
                                        br(),
                                        br(),
                                        column(4,
                                               "Metadata : ", tags$a(href="https://drive.google.com/file/d/1xDfv2smekXGG8_UqsJakCJpS1pQf6IDk/view?usp=sharing", "Metadata"))),
                   ))),

        conditionalPanel(
            condition = "input.Module == 'Spatial Transcriptomics Analysis'",
            column(9,
                   tabsetPanel(type = "pills", id = "tabs",
                               ## add file preview tab

                               ##------Spatial QC plots---------
                               tabPanel("Spatial QC plots", fluidPage(
                                   hr(),
                                   tabsetPanel(id="Spatial_QC",
                                               tabPanel(title="Spatial QC Plot", value = "Spatial_QC_panel1",
                                                        br(),
                                                        fluidRow(

                                                            column(10,
                                                                   plotOutput("nCount_SpatialPlot", width = "50%"),
                                                                   br(),
                                                                   plotOutput("SpatialFeaturePlot", width = "100%"),
                                                                   #plotlyOutput("nCount_RNAPlot", width = "100%")
                                                            ),
                                                        )

                                               )

                                   )
                               )),

                               ##------Dimension Reduction---------
                               tabPanel("Dimension Reduction", fluidPage(
                                   hr(),
                                   tabsetPanel(id="Pca_spatial",
                                               tabPanel(title="Dimension Reduction Plot",
                                                        br(),
                                                        fluidRow(
                                                            column(3,
                                                                   actionButton("doPCA_spatial", "Run", icon = icon("hand-pointer-o"))
                                                            ),

                                                            column(3,
                                                                   numericInput("clus.res2",
                                                                                label = "Resolution used",
                                                                                value = 0.6,
                                                                                min = 0.1,
                                                                                step = 0.1)
                                                            ),
                                                        ),
                                                        br(),
                                                        plotlyOutput("DimPlot_spatial", width = "100%"),
                                                        br(),
                                                        plotOutput("SpatialDimPlot", width = "100%"),
                                                        br(),
                                                        fluidRow(
                                                            column(6,
                                                                   uiOutput("deg1.gene.select"),
                                                                   plotOutput("Deg1_spatial.plot", width = "100%"),
                                                                   #plotlyOutput("Deg1.plot", width = "100%")
                                                            ),
                                                        ),
                                               ),
                                               tabPanel(title="PC Gene Visualisation",
                                                        br(),
                                                        selectInput("select.pc1",
                                                                    label = "PC to plot",
                                                                    choices = c(1:20)
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
                                               )
                                   )
                               )),

                               ##------DEGs---------
                               tabPanel("DEGs for spatial", fluidPage(
                                   hr(),
                                   fluidRow(

                                       column(3, selectInput("min_pct_spatial",
                                                             label = "min.pct",
                                                             choices = c("0.1", "0.25"))
                                       ),

                                       column(3, selectInput("logfc_spatial",
                                                             label = "logfc.threshold",
                                                             choices = c("0.1", "0.25"))
                                       ),

                                       column(4,
                                              actionButton("doDeg_spatial", "Run DEGs", icon = icon("hand-pointer-o"))
                                       )),
                                   fluidRow(
                                       column(6,
                                              uiOutput("deg2.gene.select"),
                                              plotOutput("Deg2_spatial.plot", width = "100%"),
                                       ),
                                       column(6,
                                              DT::dataTableOutput("Deg_spatial.table"),
                                       )),
                                   br(),
                                   fluidRow(column(4,
                                                   actionButton("doDegn", "Find Spatially Variable Features", icon = icon("hand-pointer-o"))
                                   )),
                                   br(),
                                   fluidRow(
                                       column(6,
                                              uiOutput("degn.gene.select"),
                                              plotOutput("Degn.plot", width = "100%"),
                                       )),
                               )),

                               ##------Deconvolution---------

                               tabPanel("Deconvolution using Spotlight", fluidPage(
                                   hr(),
                                   fluidRow(
                                       column(3,
                                              fileInput(inputId = 'tpmFiles_scRNA',
                                                        label = "scRNA-seq dataset",
                                                        multiple = FALSE,
                                                        accept = ".rds")),
                                       br(),
                                       column(3,
                                              actionButton("loadButton4", "Load scRNA-seq dataset", icon = icon("hand-o-right"),
                                              ),
                                       ),
                                       column(3,
                                              actionButton("loadexample_deconv", "Load example scRNA-seq dataset", icon = icon("hand-o-right"),
                                              ),
                                       ),
                                       #br(),
                                       fluidRow(
                                           #selectInput("dim.used",
                                           #            label = "Dimensions",
                                           #            choices = c(1:50)
                                           #),
                                           #br(),
                                           #fluidRow(
                                           column(9,
                                                  br(),
                                                  actionButton("process_scRNA", "Process scRNA-seq dataset", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  #column(6,
                                                  br(),
                                                  plotOutput("scRNAPlot", width = "100%")),
                                           br(),
                                           column(9,
                                                  br(),
                                                  actionButton("vis_spRNA", "Visualise spatial dataset", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  plotOutput("spRNAPlot", width = "100%")),
                                           br(),
                                           #),
                                           column(9,
                                                  br(),
                                                  actionButton("get_markers", "Get Marker Genes", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  actionButton("doDeconv", "Deconvolute", icon = icon("hand-pointer-o")),
                                                  br(),
                                                  br(),
                                                  plotOutput("DeconvPlot", width = "100%")),
                                           br(),
                                       )))),

                                   tabPanel("Test Data",
                                            br(),
                                            br(),
                                            column(6, "Test Data:",
                                                   br(),
                                                   br(),
                                                   "Spatial data download: ", tags$a(href="https://drive.google.com/drive/folders/1sw6Jgn-voHmomhkm8Xm2E37hiqn4ChAT?usp=sharing", "Spatial Data download"))),




                              ))),

    )
)



dbHeader <- dashboardHeader(title = "ezSinglecell")
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/JinmiaoChenLab',
                                           tags$img(src='https://avatars1.githubusercontent.com/u/8896007?s=400&u=b0029c2e64f405ea0a46d311239b674a430ec77c&v=4'
                                                    ,height='60',width='60', align='left')
                                           , tags$div('ezsinglecell', style='color:white;font-family:arial rounded MT bold'))

dashboardPage(title = "ezSinglecell : An integrated one-stop single-cell analysis toolbox for bench scientists",
	      skin = "yellow",
              dbHeader,
              dashboardSidebar(collapsed = TRUE,
                               sidebarMenu(id = "sbm",
                                           menuItem(tags$p(style = "display:inline;font-size: 20px;", ""), tabName = "seurat", icon = icon('cog'))


                               )# end of sidebarMenu
              ),#end of dashboardSidebar
              dashboardBody(
                  #includeCSS("www/custom.css")
                  useShinyalert()
                  , shinyjs::useShinyjs()
                  , tabItem(
                      tabName = "seurat"
                      , shiny_one_panel
                  ) # End of tabItem

              )# end of dashboard body
)# end of dashboard page
