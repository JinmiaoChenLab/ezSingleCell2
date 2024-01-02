source('ui.R')
#source('cellphonedb.R')
#source('plot_heatmaps.R')

## max data size
#options(shiny.maxRequestSize = 1024^10)
options(shiny.maxRequestSize = 4096*4096*100*100)
options(shiny.launch.browser = T)
options(bitmapType = 'cairo')
options(timeout=1000000)
options(show.error.messages = F)

shinyServer(function(input, output, session) {
  
  #use_python('/home/ezsinglecell/miniconda3/bin/python', required = NULL)
  
  v <- reactiveValues(scData = NULL,
                      idents = NULL,
                      isPCAdone = NULL,
                      isUMAPdone = NULL,
                      isTSNEdone = NULL,
                      isTrajectorydone = NULL,
                      isCELLiDdone = NULL,
                      isCELLdone = NULL,
                      isGSEAdone = NULL,
                      pcGenes = NULL,
                      plotlySelection = NULL,
                      ips.markers = NULL,
                      markers = NULL,
                      selected_markers = NULL)
  
  #celltypes <- NULL
  prePlot <- function(){
    while(names(dev.cur()) != "null device"){
      dev.off()
    }
  }
  observe({
    #s <- event_data("plotly_selected")
    #cells <- s[["key"]]
    v$plotlySelection <- event_data("plotly_selected")[["key"]]
  })
  
  output$demo_image <- renderImage({
    list(src = "www/Picture1.png",
         width = 500,
         height = 400)
  }, deleteFile = FALSE)
  
  output$demo_image1 <- renderImage({
    list(src = "www/Picture2.png",
         width = 330,
         height = 400)
  }, deleteFile = FALSE)
  
  output$demo_image2 <- renderImage({
    list(src = "www/overview1.png",
         width = 500,
         height = 600)
  }, deleteFile = FALSE)
  
  output$scrna_image1 <- renderImage({
    list(src = "www/scRNA-seq.png",
         height = 500)
  }, deleteFile = FALSE)
  
  #output$scrna_image2 <- renderImage({
  #  list(src = "www/scRNA-seq_2.png",
  #       height = 500)
  #}, deleteFile = FALSE)
  
  ##-------------------scRNA-seq module-------------------
  
  observeEvent(input$loadexample_tpm, {
      withProgress(message="Loading example data...", value=0.5, {
        tpmFiles_demo <- read.table("scRNA/pbmc.txt", header = T, row.names = 1, check.names = F)
        sObj <- CreateSeuratObject(tpmFiles_demo,
                                   project = input$projName)
                                   #min.genes = input$min.genes,
                                   #min.cells = input$min.cells)
        
        sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
      })
      if (is.null(tpmFiles_demo)){
        v$scData <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0.5, {
          print(tpmFiles_demo)
          #print(tpmFiles$name)
          print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample_tpm", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
          v$scData <- tpmFiles_demo
          v$scData1 <- sObj
        })
      }
  })
  
  observeEvent(input$loadexample_scH5, {
    withProgress(message="Loading example data...", value=0.5, {
      tpmFiles_demo <- Read10X_h5("scRNA/filtered_feature_bc_matrix_human.h5")
      sObj <- CreateSeuratObject(tpmFiles_demo,
                                 project = input$projName)
                                 #min.genes = input$min.genes,
                                 #min.cells = input$min.cells)
      
      sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
    })
    if (is.null(tpmFiles_demo)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0.5, {
        #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample_scH5", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        v$scData <- tpmFiles_demo
        v$scData1 <- sObj
      })
    }
  })
  
  observeEvent(input$loadexample_rds, {
    withProgress(message="Loading example data...", value=0.5, {
      tpmFiles_demo <- readRDS("scRNA/pbmc.rds")
      sObj <- CreateSeuratObject(tpmFiles_demo@assays$RNA@counts,
                                 project = input$projName)
      print(sObj)
      sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
    })
    if (is.null(tpmFiles_demo)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0.5, {
        #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample_rds", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        v$scData <- tpmFiles_demo@assays$RNA@counts
        v$scData1 <- sObj
      })
    }
  })
  
  observeEvent(input$loadexample_h5ad, {
    withProgress(message="Loading example data...", value=0.5, {
      tpmFiles_demo <- read_h5ad("scRNA/pbmc.h5ad")
      tpmFiles_demo = t(tpmFiles_demo$to_df())
      sObj <- CreateSeuratObject(tpmFiles_demo, project = input$projName)
      #min.genes = input$min.genes,
      #min.cells = input$min.cells)
      print(sObj)
      sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
      print(sObj@meta.data)
    })
    if (is.null(tpmFiles_demo)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0.5, {
        #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample_h5ad", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        v$scData <- tpmFiles_demo
        v$scData1 <- sObj
      })
    }
  })
  
  observeEvent(input$loadButton, {
    if(input$scInput == "Raw Counts Matrix"){
      tpmFiles <- input$tpmFiles
      if (is.null(tpmFiles)){
        v$scData <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else
      {
        withProgress(message="Loading and Processing Data...", value=0, {
          print(tpmFiles$datapath)
          print(tpmFiles$name)
          print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
          exp.data <- read.table(tpmFiles$datapath,
                                 sep="\t", header=TRUE, row.names=1, stringsAsFactors = F, check.names = F)
          sObj <- CreateSeuratObject(exp.data,
                                     project = input$projName)
                                     #min.genes = input$min.genes,
                                     #min.cells = input$min.cells)
          
          sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
          v$scData <- exp.data
          v$scData1 <- sObj
          label1 <- "Data loaded"
          updateActionButton(inputId = "loadButton", label = label1)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
    else if(input$scInput == "10X cellranger"){
      scH5 <- input$scH5
      if (is.null(scH5)){
        v$scData <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(scH5$datapath)
          print(scH5$name)
          print(file.exists(paste(scH5$datapath[1], "/", scH5$name[1], sep="")))
          exp.data <- Read10X_h5(scH5$datapath)
          sObj <- CreateSeuratObject(exp.data,
                                     project = input$projName)
                                     #min.genes = input$min.genes,
                                     #min.cells = input$min.cells)
          sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
          v$scData <- exp.data
          v$scData1 <- sObj
          label1 <- "Data loaded"
          updateActionButton(inputId = "loadButton", label = label1)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
    else if(input$scInput == "rds object"){
      rds <- input$rds
      if (is.null(rds)){
        v$scData <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(rds$datapath)
          print(rds$name)
          print(file.exists(paste(rds$datapath[1], "/", rds$name[1], sep="")))
          exp.data <- readRDS(rds$datapath)
          sObj <- CreateSeuratObject(exp.data@assays$RNA@counts,
                                     project = input$projName)
                                     #min.genes = input$min.genes,
                                     #min.cells = input$min.cells)
          sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
          v$scData <- exp.data@assays$RNA@counts
          v$scData1 <- sObj
          label1 <- "Data loaded"
          updateActionButton(inputId = "loadButton", label = label1)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
    else if(input$scInput == "h5ad"){
      h5ad <- input$h5ad
      if (is.null(h5ad)){
        v$scData <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(h5ad$datapath)
          print(h5ad$name)
          print(file.exists(paste(h5ad$datapath[1], "/", h5ad$name[1], sep="")))
          exp.data <- read_h5ad(h5ad$datapath)
          exp.data = t(exp.data$to_df())
          sObj <- CreateSeuratObject(exp.data, project = input$projName)
                                     
          #min.genes = input$min.genes,
          #min.cells = input$min.cells)
          sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
          v$scData <- exp.data
          v$scData1 <- sObj
          label1 <- "Data loaded"
          updateActionButton(inputId = "loadButton", label = label1)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
  })
  
  output$downloadCount <- downloadHandler(
    
    filename = function(){"counts_table.csv"}, 
    content = function(fname){
      withProgress(message="Downloading counts data...", value=0, {
      write.csv(v$scData, fname)
      })
    }
  )
    
  observeEvent(input$filter_seurat, {
    if(input$scInput == "Raw Counts Matrix"){
    tpmFiles <- input$tpmFiles
    tpmFiles <- v$scData
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message="Loading and Processing Data...", value=0, {
      print(v$scData)
      v$scData1 <- subset(v$scData1, subset = nFeature_RNA > input$ob1 & nFeature_RNA < input$ob2 & percent.mt < input$ob3)
      print(v$scData1)
      shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  }
    else if(input$scInput == "10X cellranger"){
      scH5 <- input$scH5
      scH5 <- v$scData
      if (is.null(scH5)){
        v$scData <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(v$scData)
          v$scData1 <- subset(v$scData1, subset = nFeature_RNA > input$ob1 & nFeature_RNA < input$ob2 & percent.mt < input$ob3)
          print(v$scData1)
          shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }  
    else if(input$scInput == "rds object"){
      rds <- input$rds
      rds <- v$scData
      if (is.null(rds)){
        v$scData <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(v$scData)
          v$scData1 <- subset(v$scData1, subset = nFeature_RNA > input$ob1 & nFeature_RNA < input$ob2 & percent.mt < input$ob3)
          print(v$scData1)
          shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
   else if(input$scInput == "h5ad"){
      h5ad <- input$h5ad
      h5ad <- v$scData
      if (is.null(h5ad)){
        v$scData <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(v$scData)
          v$scData1 <- subset(v$scData1, subset = nFeature_RNA > input$ob1 & nFeature_RNA < input$ob2 & percent.mt < input$ob3)
          print(v$scData1)
          shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }  
})
  
  observeEvent(input$reset_scRNA, {
    session$reload()
    print("Reset done")
  })
  #})
  
  output$logo <- renderImage({
    return(list(
      src = "inst/extdata/logo.png",
      contentType = "image/png",
      alt = "Singapore Immunology Network"
    ))
  }, deleteFile = FALSE)
  
  opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
      shell.exec(dir)
    } else {
      system(paste(Sys.getenv("R_BROWSER"), dir))
    }
  }
  
  ##---------------QC of scRNA-seq-------------------
  
  output$countdataDT <- renderDataTable({
    if(!is.null(v$scData))
    {
      if(ncol(v$scData) > 20 )
        as.matrix(v$scData[,1:20])
    }
  }, server = FALSE)
  
  plotInput <- reactive({
    p <- VlnPlot(v$scData1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend()
  })
  
  output$nFeature_RNAPlot <- renderPlot({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      print(plotInput())
    }
  })
  
  output$download_nFeature_RNA <- downloadHandler(
    filename = function(){"QC violin plots.png"}, 
    content = function(fname){
      ggsave(fname,plotInput(), width = 15, height = 10)
    }
  )
  
  ## FeatureScatter plot
  
  plotInput3 <- reactive({
    p <- FeatureScatter(v$scData1, "nCount_RNA", "nFeature_RNA")
  })
  
  output$FeatureScatterPlot1 <- renderPlotly({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      print(plotInput3())
    }
  })
  
  output$download_FeatureScatterPlot1 <- downloadHandler(
    filename = function(){"FeatureScatterPlot1.png"}, 
    content = function(fname){
      ggsave(fname,plotInput3(), height = 7, width = 7)
    }
  )
  
  plotInput4 <- reactive({
    p <- FeatureScatter(v$scData1, "nCount_RNA", "percent.mt")
  })
  
  output$FeatureScatterPlot2 <- renderPlotly({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      print(plotInput4())
    }
  })
  
  output$download_FeatureScatterPlot2 <- downloadHandler(
    filename = function(){"FeatureScatterPlot2.png"}, 
    content = function(fname){
      ggsave(fname,plotInput4(), height = 7, width = 7)
    }
  )
  
  observeEvent(input$findVarGenes, {
    tpmFiles <- v$scData1
      if (is.null(tpmFiles)){
        v$scData <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
    withProgress(message = "Finding variable genes...", value = 0, {
      if(input$norm1 == "LogNormalize"){
        v$scData1 <- NormalizeData(v$scData1)
        v$scData1 <- FindVariableFeatures(v$scData1,
                                          mean.function = ExpMean,
                                          dispersion.function = LogVMR,
                                          nfeatures = input$var.genes,
                                          selection.method = input$selection.method)
        all.genes <- rownames(v$scData1)
        v$scData1 <- ScaleData(v$scData1, features = all.genes)
        incProgress(0.5)
        #VarGeneText <- paste0("Number of variable genes: ", length(v$scData1@assays$RNA@var.features))
        #output$nVarGenes <- renderText(VarGeneText)
        varGenePlotInput <- function(){
          if(is.null(v$scData1)){
            return(NULL)
          }else{
            withProgress(message="Plotting variable genes...", value=0, {
              top10 <- head(VariableFeatures(v$scData1), 10)
              variable_feature1 <- VariableFeaturePlot(v$scData1)
              variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
              print (variable_feature1)
              print (variable_feature2)
              shinyalert("Highly variable features identified", "Highly variable features identified, please perform PCA", type = "success", imageWidth = 10, imageHeight = 10)
              #dev.off()
            })
          }
        }
        output$VarGenes <- renderPlot({
          varGenePlotInput()
        }, height = 500, width = 600)
        observeEvent(input$PDFc, {
          if(!is.null(v$scData1)){
            withProgress(message="Downloading plot PDF files...", value=0, {
              print(getwd())
              pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
              if(!dir.exists(pdfDir)){
                dir.create(pdfDir)
              }
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
              i = 0
              while(file.exists(filename2)){
                filename2 <- paste0(pdfDir, .Platform$file.sep,
                                    "Var_genes_plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
              }
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              plot1 <- VariableFeaturePlot(v$scData1)
              print(plot1)
              dev.off()
              txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
              txtfile <- sub(".pdf", ".txt", txtfile)
              write(v$scData1@assays$RNA@var.features, file = txtfile)
              
            })
          }
        })
      }
      
      else if(input$norm1 == "SCTransform"){
        v$scData1 <- SCTransform(v$scData1, variable.features.n = input$var.genes, vars.to.regress = "percent.mt", verbose = FALSE, conserve.memory = T)
        incProgress(0.5)
        #VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
        #output$nVarGenes <- renderText(VarGeneText)
        varGenePlotInput <- function(){
          if(is.null(v$scData1)){
            return(NULL)
          }else{
            withProgress(message="Plotting variable genes...", value=0, {
              top10 <- head(VariableFeatures(v$scData1), 10)
              variable_feature1 <- VariableFeaturePlot(v$scData1)
              variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
              print (variable_feature1)
              print (variable_feature2)
              shinyalert("SCTransform done", "SCTransform done, please perform PCA", type = "success", imageWidth = 10, imageHeight = 10)
            })
          }
        }
        
        output$VarGenes <- renderPlot({
          varGenePlotInput()
        }, height = 500, width = 600)
        observeEvent(input$PDFc, {
          if(!is.null(v$scData1)){
            withProgress(message="Downloading plot PDF files...", value=0, {
              print(getwd())
              pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
              if(!dir.exists(pdfDir)){
                dir.create(pdfDir)
              }
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
              i = 0
              while(file.exists(filename2)){
                filename2 <- paste0(pdfDir, .Platform$file.sep,
                                    "Var_genes_plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
              }
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              plot1 <- VariableFeaturePlot(v$scData1)
              print(plot1)
              dev.off()
              txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
              txtfile <- sub(".pdf", ".txt", txtfile)
              write(v$scData1@assays$RNA@var.features, file = txtfile)
              
            })
          }
        })
      }
    })
      }    
  })
  
  output$name <- renderPrint({
    s <- event_data("plotly_selected")
    c(s[["key"]], class(s[["key"]]))
  })
  
  
  
  ##---------------PCA of scRNA-seq-------------------
  # PCA plot
  observeEvent(input$doPCA, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Scaling Data...", value = 0,{
      incProgress(0.5, message = "Running PCA...")
      #if(input$assays1 == "LogNormalization"){
      if(length(colnames(v$scData1@assays$RNA))<5000){
      v$scData1 <- RunPCA(v$scData1, features = rownames(v$scData1))
      }
      else if(length(colnames(v$scData1@assays$RNA))>=5000){
        v$cells <- sketchData(v$scData1, percent = 0.25, idents = NULL, do.PCA = TRUE, dimPC = 30)
        v$scData1 <- v$scData1[,v$cells]
        v$scData1 <- RunPCA(v$scData1, features = rownames(v$scData1))
      }
      #print(v$scData1[["pca"]], dims = 1:5, nfeatures = 5)
      #else if(input$assays1 == "SCTransform"){
      #  v$scData1 <- RunPCA(v$scData1, features = rownames(v$scData1), assay = "SCT")
      #  print(v$scData1[["pca"]], dims = 1:5, nfeatures = 5)
      #}
      v$isPCAdone <- TRUE
      PCA_plot <- DimPlot(v$scData1, reduction = "pca", label = T)
      print(PCA_plot)
      incProgress(0.4, message = "Getting list of PC genes...")
      pc.table <- list()
      for(i in 1:20){
        pcg <- TopFeatures(v$scData1)
        pc.table[[i]] <- pcg
      }
      pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
      v$pcGenes <- pc.table
      shinyalert("PCA performed", "PCA performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
      label1 <- "PCA Done"
      updateActionButton(inputId = "doPCA", label = label1)
      })
    }
  })
  
  plotPCA <- reactive({
    if(is.null(v$scData1) || is.null(v$isPCAdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    p <- DimPlot(v$scData1, reduction = "pca", label = T)
      })
    }
  })
  
  output$PCA2DPlot <- renderPlotly({
    plotPCA()
  })
  
  output$download_PCA <- downloadHandler(
    filename = function(){"PCAPlot.png"}, 
    content = function(fname){
      ggsave(fname,PCAPlot(v$scData1), height = 7, width = 7)
    }
  )
  
  output$download_PCA_embedding <- downloadHandler(
    
    filename = function(){"pca_embeddings.csv"}, 
    content = function(fname){
      withProgress(message="Downloading counts data...", value=0, {
        write.csv(v$scData1@reductions$pca@cell.embeddings, fname)
      })
    }
  )
  
  # Viz plot
  
  plotViz <- reactive({
    if(is.null(v$scData1) || is.null(v$isPCAdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    p <- VizDimLoadings(v$scData1, dims = as.numeric(input$select.pc))
      })
    }
  })
  
  output$vizPlot <- renderPlot({
    print(plotViz())
  })
  
  output$download_vizPlot <- downloadHandler(
    filename = function(){"VizLoading_Plot.png"}, 
    content = function(fname){
      ggsave(fname,plotViz(), height = 7, width = 7)
    }
  )
  
  plotPCHeatmap <- reactive({
    if(is.null(v$scData1) || is.null(v$isPCAdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    p <- DimHeatmap(v$scData1, dims = as.numeric(input$select.pc))
      })
    }
  })
  
  output$PCHeatmap <- renderPlot({
    print(plotPCHeatmap())
  })
  
  output$download_PCHeatmap <- downloadHandler(
    filename = function(){"PC Heatmap.png"}, 
    content = function(fname){
      png(fname)
      DimHeatmap(v$scData1, dims = as.numeric(input$select.pc))
      dev.off()
    }
  )
  
  output$PCtable <- DT::renderDataTable({
    if(is.null(v$scData1) ){
      return(NULL)
    }else{
      v$pcGenes
    }
  }, server = FALSE, options = list(scrollX = TRUE))
  
  output$download_PCTable <- downloadHandler(
    filename = function(){"PC_table.csv"}, 
    content = function(fname){
      withProgress(message="Downloading PC Table...", value=0, {
        write.csv(v$pcGenes, fname)
      })
    }
  )
  
  plotElbow <- reactive({
    if(is.null(v$scData1) || is.null(v$isPCAdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    p <- ElbowPlot(v$scData1, ndims = 50)
      })
    }
  })
  
  output$Elbow <- renderPlot({
    print(plotElbow())
  })
  
  output$download_Elbow <- downloadHandler(
    filename = function(){"Elbow plot.png"}, 
    content = function(fname){
      ggsave(fname,plotElbow(), height = 7, width = 7)
    }
  )
  
  ##---------------Clustering of scRNA-seq-------------------
  
  observeEvent(input$findCluster, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Finding clusters...", value = 0.3, {
      #if(input$assays1 == "LogNormalization"){
      #DefaultAssay(v$scData1) <- "RNA"
      v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, nn.method = "rann")
      v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
      print(v$scData1@meta.data)
      output$cluster.done <- renderText(paste0("Clustering done!"))
      #}
      #else if(input$assays1 == "SCTransform"){
      #  DefaultAssay(v$scData1) <- "SCT"
      #  v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, assay = "SCT", nn.method = "rann")
      #  v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
      #  output$cluster.done <- renderText(paste0("Clustering done!"))
      #}
      v$isClusterdone <- TRUE
      shinyalert("Clustering performed", "Clustering performed, please identify optimum number of clusters", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }  
  })
  
  plotCluster <- reactive({
    p <- DimPlot(v$scData1, reduction = "umap", label = T) + NoLegend()
  })
  
  output$Cluster2DPlot_1 <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isClusterdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    print(plotCluster())
      })
    }
  })
  
  output$download_Cluster <- downloadHandler(
    filename = function(){"Cluster plot.png"}, 
    content = function(fname){
      ggsave(fname,plotCluster(), height = 7, width = 7)
    }
  )
  
  output$download_ClusterTable <- downloadHandler(
    
    filename = function(){"Cluster_table.csv"}, 
    content = function(fname){
      withProgress(message="Downloading Cluster Table...", value=0, {
        write.csv(v$scData1$seurat_clusters, fname)
      })
    }
  )
  
  ##---------Find Optimum Cluster Resolution---------##
  
  observeEvent(input$findoptimumCluster, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Finding optimum resolution...", value = 0.3, {
      #if(input$assays1 == "LogNormalization"){
      #DefaultAssay(v$scData1) <- "RNA"
      reses<-seq(from = input$clus.res_a, to = input$clus.res_b, by = 0.1)
      for (res in reses){
        v$scData1<-FindClusters(v$scData1, resolution = res)
        nores<-gsub(pattern = ".", replacement = "", res, fixed = T)
      }
      print(v$scData1@meta.data)
      output$optimumcluster.done <- renderText(paste0("Clustering done!"))
      #}
      #else if(input$assays1 == "SCTransform"){
      #  DefaultAssay(v$scData1) <- "SCT"
      #  v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, assay = "SCT", nn.method = "rann")
      #  v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
      #  output$cluster.done <- renderText(paste0("Clustering done!"))
      #}
      v$isOptimumClusterdone <- TRUE
      shinyalert("Number of clusters identified at different resolution", "Number of clusters identified at different resolution, please perform UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  plotOptimumCluster <- reactive({
    if(is.null(v$scData1) || is.null(v$isOptimumClusterdone)){
      plotly_empty()
    }else{
      withProgress(message="Plotting Clustering tree...", value=0, {
    p <- clustree(v$scData1)
    })
  }
})
  
  output$OptimumCluster2DPlot_1 <- renderPlot({
    print(plotOptimumCluster())
  })
  
  output$download_OptimumCluster <- downloadHandler(
    filename = function(){"Clustree plot.png"}, 
    content = function(fname){
      ggsave(fname,plotOptimumCluster(), height = 7, width = 7)
    }
  )
  
  output$download_OptimumClusterTable <- downloadHandler(
    
    filename = function(){"Clustering_different_resolution.csv"}, 
    content = function(fname){
      withProgress(message="Downloading Cluster Table at different resolution...", value=0, {
        write.csv(v$scData1@meta.data, fname)
      })
    }
  )
  
  ##---------------UMAP of scRNA-seq-------------------
  
  observeEvent(input$doUmap, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Running UMAP...", value = 0.3, {
      #if(input$assays1 == "LogNormalization"){
      v$scData1 <- RunUMAP(v$scData1, dims = 1:input$dim.used, spread = 1)
      output$Umap.done <- renderText(paste0("UMAP done!"))
      #}
      #else if(input$assays1 == "SCTransform"){
      #  v$scData1 <- RunUMAP(v$scData1, dims = 1:input$dim.used, assay = "SCT", spread = 1)
      # output$Umap.done <- renderText(paste0("UMAP done!"))
      #}
      v$isUMAPdone <- TRUE
      shinyalert("UMAP performed", "UMAP performed, please perform tSNE", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  plotUMAP <- reactive({
    if(is.null(v$scData1) || is.null(v$isUMAPdone)){
      return(NULL)
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    p <- DimPlot(v$scData1, reduction = "umap", label = T, label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Umap_2d_plot_1 <- renderPlotly({
    plotUMAP()
  })
  
  output$download_UMAP <- downloadHandler(
    filename = function(){"UMAP plot.png"}, 
    content = function(fname){
      ggsave(fname,plotUMAP(), height = 7, width = 7)
    }
  )
  
  output$download_UMAP_embedding <- downloadHandler(
    
    filename = function(){"UMAP Embeddings.csv"}, 
    content = function(fname){
      withProgress(message="Downloading UMAP Embeddings...", value=0, {
        write.csv(v$scData1@reductions$umap@cell.embeddings, fname)
      })
    }
  )
  
  ##---------------TSNE of scRNA-seq-------------------
  output$perplex.option <- renderUI({
    if(is.null(v$isPCAdone)){
      return(NULL)
    }else{
      ##perplexity test
      n.cells <- isolate(nrow(v$scData1@reductions$pca@cell.embeddings))
      max.perplex <- as.integer((n.cells - 1)/3)
      numericInput("perplexity",
                   label = "Perplexity",
                   value = if(max.perplex <30) max.perplex else 30,
                   min = 0,
                   max = max.perplex)
    }
  })
  
  observeEvent(input$doTsne, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Running tSNE...", value = 0.3, {
      #if(input$assays1 == "LogNormalization"){
      v$scData1 <- RunTSNE(v$scData1, dims = 1:input$dim.used, perplexity = input$perplexity)
      output$Tsne.done <- renderText(paste0("TSNE done!"))
      #}
      #else if(input$assays1 == "SCTransform"){
      #  v$scData1 <- RunTSNE(v$scData1, dims = 1:input$dim.used, perplexity = input$perplexity, assay = "SCT")
      #  output$Tsne.done <- renderText(paste0("TSNE done!"))
      #}
      v$isTSNEdone <- TRUE
      shinyalert("tSNE performed", "tSNE performed, please perform celltype annotation", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  plotTsne <- reactive({
    if(is.null(v$scData1) || is.null(v$isTSNEdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "tsne", label = T, label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Tsne_2d_plot_1 <- renderPlotly({
    plotTsne()
  })
  
  output$download_Tsne <- downloadHandler(
    filename = function(){"tSNE plot.png"}, 
    content = function(fname){
      ggsave(fname,plotTsne(), height = 7, width = 7)
    }
  )
  
  output$download_Tsne_embedding <- downloadHandler(
    filename = function(){"tSNE Embeddings.csv"}, 
    content = function(fname){
      withProgress(message="Downloading tSNE Embeddings...", value=0, {
        write.csv(v$scData1@reductions$tsne@cell.embeddings, fname)
      })
    }
  )
  
  ##---------------Run CELLiD of scRNA-seq-------------------
  
  observeEvent(input$doCELLiD, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Running CELLiD...", value = 0.3, {
      ref = readRDS('ref.rds')
      v$scData1.rna.data.average = AverageExpression(v$scData1)
      v$scData1.rna.data.average = round(v$scData1.rna.data.average$RNA, 2)
      #v$scData1.rna.data.average = data.frame(v$scData1.rna.data.average$RNA)
      if(input$cellatlas == "all"){
      v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, ref)
      print(v$res)
      v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
      v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
      newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
      colnames(v$res) <- newheaders
      print(v$scData1@meta.data)
      output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
      write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "adipose"){
        adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
        adipose1 <- ref[,adipose]
        colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, adipose1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "adrenal_gland"){
        adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
        adrenal_gland1 <- ref[,adrenal_gland]
        colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, adrenal_gland1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "blood"){
        blood <- colnames(ref)[grepl("blood",colnames(ref))] 
        blood1 <- ref[,blood]
        colnames(blood1) <- gsub("--blood","",colnames(blood1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, blood1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "bone_marrow"){
        bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
        bone_marrow1 <- ref[,bone_marrow]
        colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, bone_marrow1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "brain"){
        brain <- colnames(ref)[grepl("brain",colnames(ref))] 
        brain1 <- ref[,brain]
        colnames(brain1) <- gsub("--brain","",colnames(brain1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, brain1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "breast"){
        breast <- colnames(ref)[grepl("breast",colnames(ref))] 
        breast1 <- ref[,breast]
        colnames(breast1) <- gsub("--breast","",colnames(breast1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, breast1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "breast_milk"){
        breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
        breast_milk1 <- ref[,breast_milk]
        colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, breast_milk1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "eye"){
        eye <- colnames(ref)[grepl("eye",colnames(ref))] 
        eye1 <- ref[,eye]
        colnames(eye1) <- gsub("--eye","",colnames(eye1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, eye1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "gut"){
        gut <- colnames(ref)[grepl("gut",colnames(ref))] 
        gut1 <- ref[,gut]
        colnames(gut1) <- gsub("--gut","",colnames(gut1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, gut1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "heart"){
        heart <- colnames(ref)[grepl("heart",colnames(ref))] 
        heart1 <- ref[,heart]
        colnames(heart1) <- gsub("--heart","",colnames(heart1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, heart1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "kidney"){
        kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
        kidney1 <- ref[,kidney]
        colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, kidney1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "liver"){
        liver <- colnames(ref)[grepl("liver",colnames(ref))] 
        liver1 <- ref[,liver]
        colnames(liver1) <- gsub("--liver","",colnames(liver1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, liver1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "lung"){
        lung <- colnames(ref)[grepl("lung",colnames(ref))] 
        lung1 <- ref[,lung]
        colnames(lung1) <- gsub("--lung","",colnames(lung1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, lung1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "pancreas"){
        pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
        pancreas1 <- ref[,pancreas]
        colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, pancreas1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "PDAC"){
        PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
        PDAC1 <- ref[,PDAC]
        colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, PDAC1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "skin"){
        skin <- colnames(ref)[grepl("skin",colnames(ref))] 
        skin1 <- ref[,skin]
        colnames(skin1) <- gsub("--skin","",colnames(skin1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, skin1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "testis"){
        testis <- colnames(ref)[grepl("testis",colnames(ref))] 
        testis1 <- ref[,testis]
        colnames(testis1) <- gsub("--testis","",colnames(testis1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, testis1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "thymus"){
        thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
        thymus1 <- ref[,thymus]
        colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, thymus1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      if(input$cellatlas == "tonsil"){
        tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
        tonsil1 <- ref[,tonsil]
        colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
        v$res = FastIntegration::CELLiD(v$scData1.rna.data.average, tonsil1)
        print(v$res)
        v$scData1$primary.predict = v$res[as.numeric(v$scData1$seurat_clusters),1]
        v$scData1$secondary.predict = v$res[as.numeric(v$scData1$seurat_clusters),2]
        newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
        colnames(v$res) <- newheaders
        print(v$scData1@meta.data)
        output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
        write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
      }
      v$isCELLiDdone <- TRUE
      shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform cell-cell similarity", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  output$Umap_cellid <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP from CELLiD...", value=0, {
        DimPlot(v$scData1, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Umap_cellid1 <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP from CELLiD...", value=0, {
        DimPlot(v$scData1, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
      })
    }
  })
  
  output$ct.table <- DT::renderDataTable(
    v$res, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
  
  #onclick("commitButton",{
  #  proxy=dataTableProxy("ct.table")
  #  replaceData(ct.table,ds)
  #})
  
  plotCELLiD <- reactive({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Umap_cellid <- renderPlotly({
    plotCELLiD() 
  })
  
  output$download_Umap_cellid <- downloadHandler(
    filename = function(){"Celltype identification plot (Primary prediction).png"}, 
    content = function(fname){
      ggsave(fname,plotCELLiD(), height = 7, width = 7)
    }
  )
  
  plotCELLiD1 <- reactive({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Umap_cellid1 <- renderPlotly({
    plotCELLiD1()
  })
  
  output$download_Umap_cellid1 <- downloadHandler(
    filename = function(){"Celltype identification plot (Secondary prediction).png"}, 
    content = function(fname){
      ggsave(fname,plotCELLiD1(), height = 7, width = 7)
    }
  )
  
  output$download_cellid_prediction <- downloadHandler(
    
    filename = function(){"CELLiD_predictions.csv"}, 
    content = function(fname){
      withProgress(message="Downloading CELLiD predictions...", value=0, {
        write.csv(v$res, fname)
      })
    }
  )
  
  recodeValues <- reactive({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)) {return()}
    uniqueVars <- v$res[,1]
    print(uniqueVars)
    convert <- cbind(uniqueVars, uniqueVars)
    print(uniqueVars)
    print(convert)
    colnames(convert) <- c("Original ID", "newID")
    convert
  })
  
  output$recoding <- DT::renderDataTable({
    recodeValues()
  }, server = FALSE)
  
  plotCELLiD2 <- reactive({
    if(is.null(v$scData1) || is.null(v$isRenameCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "newID", label = T,  label.size = 3) + NoLegend()
        p
      })
    }
  })
  
  output$Umap_cellid2 <- renderPlotly({
    plotCELLiD2()
  })
  
  fname = 'CELLiD_predictions.csv'
  values = list()
  setHot = function(x) values[["hot"]] <<- x
  
  observeEvent(input$commitButton, {
    if (!is.null(values[["hot"]])) {
      write.csv(values[["hot"]], fname)
      print(fname)
      
      ##### use the recoded values to rename values throughout the entire data frame
      recodes <- read.csv(fname, stringsAsFactors = FALSE)
      
      for(i in 1:nrow(recodes)){
        if(recodes[i,2] != recodes[i,3]){
          v$res[v$res == recodes[i,2]] <- recodes[i,3]
          }
        }
      test <- read.csv('CELLiD_predictions.csv', header = T, row.names = 1)
      as.numeric(rownames(test)) -> test$seurat_clusters
      test[3] <- test[3]-1
      test <- test[,-1]
      print(test)
      test$seurat_clusters <- as.factor(test$seurat_clusters)
      meta <- v$scData1@meta.data
      rownames(meta) -> meta$cellID
      new <- merge(x=meta,y=test,by="seurat_clusters")
      colnames(new) <- gsub('.x','',names(new))
      colnames(new) <- gsub('.y','',names(new))
      new$cellID -> rownames(new)
      print(new)
      v$scData1 <- AddMetaData(v$scData1, metadata = new)
      print(v$scData1@meta.data)
      v$isRenameCELLiDdone <- TRUE
      #output$newVars <- DT::renderDataTable({
      #  v$res
      #})
    }  
  })
  
  output$hot = renderRHandsontable({
    if (!is.null(input$hot)) {
      DF = hot_to_r(input$hot)
    } else {
      if(is.null(v$res)) {return()}
      DF = recodeValues()
    }
    
    setHot(DF)
    
    rhandsontable(DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })
  
  output$ct.gene.select <- renderUI({
    if(is.null(v$scData1)|| is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      selectInput("ct.gene", label = "Gene to visualise",
                  choices = rownames(v$scData1@assays$RNA@counts))
    }
  })
  
  observeEvent(input$Vis_seurat1, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData1 <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
      withProgress(message="Visualizing...", value=0, {
        v$isVisCELLiDdone <- TRUE
      })
    }
  })
  
  plotViolin1 <- reactive({
    if(is.null(v$scData1)|| is.null(v$isVisCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating Violin Plot...", value=0, {
          print(v$scData1)
          v$scData1$primary.predict -> Idents(v$scData1)
          VlnPlot(v$scData1, input$ct.gene)
      })
    }
  })
  
  output$ct.gene.plot <- renderPlotly({
    plotViolin1()
  })
  
  output$download_violn1 <- downloadHandler(
    filename = function(){"Violin plot (CELLiD).png"}, 
    content = function(fname){
      ggsave(fname,plotViolin(), height = 7, width = 7)
    }
  )
  
  plotFeature1 <- reactive({
    if(is.null(v$scData1)|| is.null(v$isVisCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating Violin Plot...", value=0, {
        print(v$scData1)
        v$scData1$primary.predict -> Idents(v$scData1)
        FeaturePlot(v$scData1, input$ct.gene)
      })
    }
  })
  
  output$ct.gene1.plot <- renderPlotly({
    plotFeature1()
  })
  
  output$download_feature1 <- downloadHandler(
    filename = function(){"Feature plot (CELLiD).png"}, 
    content = function(fname){
      ggsave(fname,plotFeature1(), height = 7, width = 7)
    }
  )
  
  #onclick("BRefresh",{
  #  proxy=dataTableProxy("OPreview")
  #  replaceData(proxy,ds)
  #})
  
  ##---------------Run celltypist-----------------##
  
  observeEvent(input$doCelltypist, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Running Celltypist...", value = 0.3, {
      sc <- reticulate::import("scanpy", convert = FALSE)
      ct <- reticulate::import("celltypist", convert = FALSE)
      sceasy::convertFormat(v$scData1, from = "seurat", to = "anndata", outFile = 'ct_scrna.h5ad')
      v$adata = sc$read_h5ad('ct_scrna.h5ad')
      v$res = ct$annotate(filename = 'ct_scrna.h5ad', model = input$celltypistatlas, majority_voting=T)
      print(v$res)
      v$adata = v$res$to_adata()
      print("fff")
      v$adata$obs$to_csv('celltypist_predict.csv')
      v$meta <- read.csv('celltypist_predict.csv', header = T, row.names = 1)
      v$scData1 <- AddMetaData(v$scData1, metadata = v$meta)
      v$scData1$primary.predict <- v$scData1$majority_voting
      v$scData1$secondary.predict <- v$scData1$predicted_labels
      v$res1 <- v$scData1@meta.data
      print(v$scData1@meta.data)
      v$isCelltypistdone <- TRUE
      shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform cell-cell similarity", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  output$Umap_celltypist <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isCelltypistdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP from Celltypist...", value=0, {
        DimPlot(v$scData1, reduction = "umap", group.by = "majority_voting", label = T,  label.size = 3) 
      })
    }
  })
  
  output$Umap_celltypist1 <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isCelltypistdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP from Celltypist...", value=0, {
        DimPlot(v$scData1, reduction = "umap", group.by = "predicted_labels", label = T,  label.size = 3)
      })
    }
  })
  
  output$celltypist.table <- DT::renderDataTable(
    v$meta, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
  
  plotCELLTypist <- reactive({
    if(is.null(v$scData1) || is.null(v$isCelltypistdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "majority_voting", label = T,  label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Umap_celltypist <- renderPlotly({
    plotCELLTypist() 
  })
  
  output$download_Umap_celltypist <- downloadHandler(
    filename = function(){"Celltypist plot (Majority voting).png"}, 
    content = function(fname){
      ggsave(fname,plotCELLTypist(), height = 7, width = 7)
    }
  )
  
  plotCELLTypist1 <- reactive({
    if(is.null(v$scData1) || is.null(v$isCelltypistdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "predicted_labels", label = T,  label.size = 3) + NoLegend()
      })
    }
  })
  
  output$Umap_celltypist1 <- renderPlotly({
    plotCELLTypist1()
  })
  
  output$download_Umap_celltypist1 <- downloadHandler(
    filename = function(){"Celltypist plot (Predicted labels).png"}, 
    content = function(fname){
      ggsave(fname,plotCELLTypist1(), height = 7, width = 7)
    }
  )
  
  output$download_celltypist_prediction <- downloadHandler(
    filename = function(){"Celltypist predictions.csv"}, 
    content = function(fname){
      withProgress(message="Downloading Celltypist predictions...", value=0, {
        write.csv(v$meta, fname)
      })
    }
  )
  
  plotCELLTypist2 <- reactive({
    if(is.null(v$scData1) || is.null(v$isRenameCELLTypistdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        test <- read.csv('CELLTypist predictions.csv', header = T, row.names = 1)
        test <- test[c(10,8,5)]
        #test[4] <- test[4]
        #test <- test[,-1]
        test$seurat_clusters <- as.factor(test$seurat_clusters)
        meta <- v$scData1@meta.data
        rownames(meta) -> meta$cellID
        new <- merge(x=meta,y=test,by="seurat_clusters")
        new$cellID -> rownames(new)
        v$scData1 <- AddMetaData(v$scData1, metadata = new)
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "newID", label = T,  label.size = 3) + NoLegend()
        p
      })
    }
  })
  
  output$Umap_celltypist2 <- renderPlotly({
    plotCELLTypist2()
  })
  
  
  
  output$celltypist.gene.select <- renderUI({
    if(is.null(v$scData1)|| is.null(v$isCelltypistdone)){
      return(NULL)
    }else{
      selectInput("celltypist.gene", label = "Gene to visualise",
                  choices = rownames(v$scData1@assays$RNA@counts))
    }
  })
  
  output$celltypist.gene.plot <- renderPlotly({
    if(is.null(v$scData1)|| is.null(v$isCelltypistdone)){
      return(NULL)
    }else{
      withProgress(message="Generating Plot...", value=0, {
        print(v$scData1)
        v$scData1$majority_voting -> Idents(v$scData1)
        VlnPlot(v$scData1, input$celltypist.gene)
      })
    }
  })
  
  output$celltypist.gene1.plot <- renderPlotly({
    if(is.null(v$scData1)|| is.null(v$isCelltypistdone)){
      return(NULL)
    }else{
      withProgress(message="Generating Plot...", value=0, {
        FeaturePlot(v$scData1, input$celltypist.gene)
      })
    }
  })
  
  observe({
    if(input$dodeconv_spatial > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler3", "2")
    }
  })
  
  observe({
    if(input$dodeconv_spatial1 > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler3a", "2")
    }
  })
  
  observe({
    if(input$dodeconv_spatial_intg > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler6", "2")
    }
  })
  
  observe({
    if(input$doct_atac_intg > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler7", "2")
    }
  })
  
  observe({
    if(input$doct_atac > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler5", "2")
    }
  })
  
  observe({
    if(input$doct_atac1 > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler5a", "2")
    }
  })
  
  ##---------------Cell-cell similarity of scRNA-seq-------------------
  
  observeEvent(input$cell_cell, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message="Generating heatmap...", value=0, {
      if(input$cell1 == "primary.predict"){
      v$scData1$primary.predict -> Idents(v$scData1)
      v$scData1.rna.data.average1 = AverageExpression(v$scData1)
      v$scData1.rna.data.average1 = data.frame(v$scData1.rna.data.average1$RNA)
      print(v$scData1.rna.data.average1)
      v$cor <- cor(v$scData1.rna.data.average1, method = input$corr_method)
      print(v$cor)
      output$CELL.done <- renderText(paste0("Celltype similarity done!"))
      v$isCELLdone <- TRUE
      }
      if(input$cell1 == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        v$scData1.rna.data.average1 = AverageExpression(v$scData1)
        v$scData1.rna.data.average1 = data.frame(v$scData1.rna.data.average1$RNA)
        print(v$scData1.rna.data.average1)
        v$cor <- cor(v$scData1.rna.data.average1, method = input$corr_method)
        rownames(v$cor) <- substr(rownames(v$cor),2,nchar(rownames(v$cor)))
        colnames(v$cor) <- substr(colnames(v$cor),2,nchar(colnames(v$cor)))
        print(v$cor)
        output$CELL.done <- renderText(paste0("Celltype similarity done!"))
        v$isCELLdone <- TRUE
      }
      if(input$cell1 == "newID"){
        v$scData1$newID -> Idents(v$scData1)
        v$scData1.rna.data.average1 = AverageExpression(v$scData1)
        v$scData1.rna.data.average1 = data.frame(v$scData1.rna.data.average1$RNA)
        print(v$scData1.rna.data.average1)
        v$cor <- cor(v$scData1.rna.data.average1, method = input$corr_method)
        print(v$cor)
        output$CELL.done <- renderText(paste0("Celltype similarity done!"))
        v$isCELLdone <- TRUE
        }
      })
    }
  })
  
  plotCELL <- reactive({
    if(is.null(v$scData1.rna.data.average1) || is.null(v$isCELLdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating Celltype similarity plot...", value=0, {
        p <- heatmaply(as.matrix(v$cor), cexRow = 0.8, cexCol = 0.8, margins = c(10,10), k_col =2, k_row = 2, colors = rev(RColorBrewer::brewer.pal(9, "RdBu")))
      })
    }
  })
  
  output$cell_cell_sim <- renderPlotly({
    plotCELL()
  })
  
  output$download_cell_cell_sim <- downloadHandler(
    filename = "Celltype similarity.png",
    content = function(file) {
      png(file)
      heatmap(as.matrix(v$cor), col = RColorBrewer::brewer.pal(9, "RdBu"), cexRow = 0.8, cexCol = 0.8, margins = c(10,10))
      dev.off()
    }
  )
  
  output$cor.table <- DT::renderDataTable(
    v$cor, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
  
  output$download_cor.table <- downloadHandler(
    filename = function(){"Celltype_similarity.csv"}, 
    content = function(fname){
      withProgress(message="Downloading celltype_similarity...", value=0, {
        write.csv(v$cor, fname)
      })
    }
  )
  
  ##---------------Sub-cluster analysis--------------##
  
  output$subcluster.gene.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      selectInput("cluster_subcluster", label = "Cluster to subcluster",
                  choices = unique(v$scData1$seurat_clusters), selected = unique(v$scData1$seurat_clusters)[1])
    }
  })
  
  observeEvent(input$subcluster, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message="Performing subcluster analysis...", value=0, {
      v$scData1$seurat_clusters -> Idents(v$scData1)
      v$scData1.subset <- FindSubCluster(v$scData1, input$cluster_subcluster, "RNA_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res, algorithm = 1)
      print(v$scData1.subset)
      output$Subcluster.done <- renderText(paste0("Subclustering done!"))
      v$isSubclusterdone <- TRUE
      })
    }
  })
  
  plotSubclster <- reactive({
    if(is.null(v$scData1.subset) || is.null(v$isSubclusterdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating Celltype similarity plot...", value=0, {
        p <- DimPlot(v$scData1.subset, reduction = "umap", group.by = "sub.cluster", label = T, label.size = 2.5)
      })
    }
  })
  
  output$subcluster_plot <- renderPlotly({
    plotSubclster()
  })
  
  output$download_subcluster <- downloadHandler(
    filename = function(){"Subcluster plot.png"}, 
    content = function(fname){
      ggsave(fname,plotSubclster(), height = 7, width = 7)
    }
  )
  
  ##---------------DEGs of scRNA-seq-------------------
  
  observeEvent(input$doDeg, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
      withProgress(message="Finding DEGs...", value=0, {
        if(input$deg1 == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        if (input$norm1 == "LogNormalize"){ 
        ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "RNA", test.use = input$test.use)
        v$ips.markers <- ips.markers
        v$isDEGdone <- TRUE
        }
        if (input$norm1 == "SCTransform"){ 
        ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "SCT", test.use = input$test.use)
        v$ips.markers <- ips.markers
        v$isDEGdone <- TRUE
        }
        shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        }
        if(input$deg1 == "primary.predict"){
          v$scData1$primary.predict -> Idents(v$scData1)
          if (input$norm1 == "LogNormalize"){ 
          ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "RNA", test.use = input$test.use)
          v$ips.markers <- ips.markers
          v$isDEGdone <- TRUE
          }
          if (input$norm1 == "SCTransform"){ 
          ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "SCT", test.use = input$test.use)
          v$ips.markers <- ips.markers
          v$isDEGdone <- TRUE
          }
          shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        }
        if(input$deg1 == "newID"){
          v$scData1$newID -> Idents(v$scData1)
          if (input$norm1 == "LogNormalize"){ 
            ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "RNA", test.use = input$test.use)
            v$ips.markers <- ips.markers
            v$isDEGdone <- TRUE
          }
          if (input$norm1 == "SCTransform"){ 
            ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "SCT", test.use = input$test.use)
            v$ips.markers <- ips.markers
            v$isDEGdone <- TRUE
          }
          shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    }
  })
  
  observeEvent(input$Vis_seurat, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData1 <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
      withProgress(message="Visualizing...", value=0, {
        v$isVisdone <- TRUE
      })
    }
  })
  
  output$deg.gene.select <- renderUI({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      selectInput("deg.gene", label = "Gene to visualise",
                  choices = rownames(v$scData1@assays$RNA@counts), selected = rownames(v$scData1@assays$RNA@counts)[50])
    }
  })
  
  plotViolin <- reactive({
    if(is.null(v$scData1) || is.null(v$isVisdone)){
      plotly_empty()
    }else{
      withProgress(message="Visualizing...", value=0, {
        if(input$deg2 == "seurat_clusters"){
          print(v$scData1)
          v$scData1$seurat_clusters -> Idents(v$scData1)
          print(v$scData1)
          VlnPlot(v$scData1, input$deg.gene)
        }
        else if(input$deg2 == "primary.predict"){
          v$scData1$primary.predict -> Idents(v$scData1)
          VlnPlot(v$scData1, input$deg.gene)
        }
        else if(input$deg2 == "newID"){
          v$scData1$newID -> Idents(v$scData1)
          VlnPlot(v$scData1, input$deg.gene)
        }
      })
    }
  })
  
  output$Deg.plot <- renderPlotly({
    plotViolin()
  })
  
  output$download_violn <- downloadHandler(
    filename = function(){"Violin plot.png"}, 
    content = function(fname){
      ggsave(fname,plotViolin(), height = 7, width = 7)
    }
  )
  
  plotFeature <- reactive({
    if(is.null(v$scData1) || is.null(v$isVisdone)){
      plotly_empty()
    }else{
      withProgress(message="Visualizing...", value=0, {
        if(input$deg2 == "seurat_clusters"){
          v$scData1$seurat_clusters -> Idents(v$scData1)
          FeaturePlot(v$scData1, input$deg.gene)
        }
        else if(input$deg2 == "primary.predict"){
          v$scData1$primary.predict -> Idents(v$scData1)
          FeaturePlot(v$scData1, input$deg.gene)
        }
        else if(input$deg2 == "newID"){
          v$scData1$newID -> Idents(v$scData1)
          FeaturePlot(v$scData1, input$deg.gene)
        }
      })
    }
  })
  
  output$Deg1.plot <- renderPlotly({
    plotFeature()
  })
  
  output$download_feature <- downloadHandler(
    filename = function(){"Feature plot.png"}, 
    content = function(fname){
      ggsave(fname,plotFeature(), height = 7, width = 7)
    }
  )
  
  plotRidge <- reactive({
    if(is.null(v$scData1) || is.null(v$isVisdone)){
      plotly_empty()
    }else{
      withProgress(message="Visualizing...", value=0, {
        if(input$deg2 == "seurat_clusters"){
          v$scData1$seurat_clusters -> Idents(v$scData1)
          RidgePlot(v$scData1, features = input$deg.gene)
        }
        else if(input$deg2 == "primary.predict"){
          v$scData1$primary.predict -> Idents(v$scData1)
          RidgePlot(v$scData1, features = input$deg.gene)
        }
        else if(input$deg2 == "newID"){
          v$scData1$newID -> Idents(v$scData1)
          RidgePlot(v$scData1, features = input$deg.gene)
        }
      })
    }
  })
  
  output$Deg2.plot <- renderPlot({
    plotRidge()
  })
  
  output$download_ridge <- downloadHandler(
    filename = function(){"Ridge plot.png"}, 
    content = function(fname){
      ggsave(fname,plotRidge(), height = 7, width = 7)
    }
  )
  
  output$Deg3.plot <- renderPlot({
    if(is.null(v$ips.markers) || is.null(v$isDEGdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        v$ips.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
        DoHeatmap(v$scData1, features = top10$gene, size = 3, angle = 30) + theme(axis.text.y = element_text(size = 4)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
      })
    }
  })
  
  output$Deg.table <- DT::renderDataTable(
    v$ips.markers, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
  
  output$gene1.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      if(input$deg3 == "seurat_clusters"){
      selectInput("gene1", label = "Celltype1",
                  choices = as.vector(v$scData1$seurat_clusters))
      }
      else if(input$deg3 == "primary.predict"){
        selectInput("gene1", label = "Celltype1",
                    choices = as.vector(v$scData1$primary.predict))
      }
      else if(input$deg3 == "newID"){
        selectInput("gene1", label = "Celltype1",
                    choices = as.vector(v$scData1$newID))
      }
    }
  })
  
  output$gene2.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      if(input$deg3 == "seurat_clusters"){
      selectInput("gene2", label = "Celltype2",
                  choices = as.vector(v$scData1$seurat_clusters))
      }
      else if(input$deg3 == "primary.predict"){
        selectInput("gene2", label = "Celltype2",
                    choices = as.vector(v$scData1$primary.predict))
      }
      else if(input$deg3 == "newID"){
        selectInput("gene2", label = "Celltype2",
                    choices = as.vector(v$scData1$newID))
      }
    }
  })
  
  observeEvent(input$doVolcano, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
      withProgress(message="Finding DEGs...", value=0, {
        v$scData1_subset <- subset(v$scData1, subset = primary.predict == input$gene1 | primary.predict == input$gene2)
        if (input$norm1 == "LogNormalize"){ 
          ips.markers_a <- FindAllMarkers(v$scData1_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
          ips.markers_b <- FindMarkers(v$scData1_subset, ident.1 = input$gene1, ident.2 = input$gene2, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
          v$ips.markers_a <- ips.markers_a
          v$ips.markers_b <- ips.markers_b
        }
        if (input$norm1 == "SCTransform"){ 
          ips.markers_a <- FindAllMarkers(v$scData1_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "SCT", test.use = input$test.use_a)
          ips.markers_b <- FindMarkers(v$scData1_subset, ident.1 = input$gene1, ident.2 = input$gene2, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "SCT", test.use = input$test.use_a)
          v$ips.markers_a <- ips.markers_a
          v$ips.markers_b <- ips.markers_b
        }
        shinyalert("Pairwise DEGs done", "Pairwise DEGs done, please run GSEA Analysis", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  output$volcano.plot <- renderPlot({
    if(is.null(v$ips.markers_b)){
      plotly_empty()
    }else{
      withProgress(message="Generating Volcano Plot...", value=0, {
        EnhancedVolcano(toptable = v$ips.markers_b, lab = row.names(v$ips.markers_b), x ="avg_log2FC", y ="p_val_adj", pointSize = 1, labSize = 5, legendLabSize = 12, axisLabSize = 12)
      })
    }
  })
  
  output$dega.plot <- renderPlot({
    if(is.null(v$ips.markers_a)){
      plotly_empty()
    }else{
      withProgress(message="Generating Volcano Plot...", value=0, {
        v$ips.markers_a %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
        DoHeatmap(v$scData1_subset, features = top10$gene, size = 4, angle = 45) + theme(axis.text.y = element_text(size = 6)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
      })
    }
  })
  
  ##---------------GSEA of scRNA-seq-------------------##
  
  output$gsea.ct1.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      selectInput("gsea.ct1", label = "Celltype1",
                  choices = as.vector(v$scData1$primary.predict))
    }
  })
  
  output$gsea.ct2.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      selectInput("gsea.ct2", label = "Celltype2",
                  choices = as.vector(v$scData1$primary.predict))
    }
  })
  
  observeEvent(input$gsea, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message="Generating gene set enrichment analysis...", value=0, {
      if(input$species_gsea == "Homo sapiens" & input$category_gsea == "H"){
        v$msigdbr_hs_go <- msigdbr(species = "Homo sapiens", category = "H")
        print(v$msigdbr_hs_go)
        v$pathways <- split(x = v$msigdbr_hs_go$gene_symbol, f = v$msigdbr_hs_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Mus musculus" & input$category_gsea == "H"){
        v$msigdbr_mm_go <- msigdbr(species = "Mus musculus", category = "H")
        print(v$msigdbr_mm_go)
        v$pathways <- split(x = v$msigdbr_mm_go$gene_symbol, f = v$msigdbr_mm_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Homo sapiens" & input$category_gsea == "C2"){
        v$msigdbr_hs_go <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
        print(v$msigdbr_hs_go)
        v$pathways <- split(x = v$msigdbr_hs_go$gene_symbol, f = v$msigdbr_hs_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Mus musculus" & input$category_gsea == "C2"){
        v$msigdbr_mm_go <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
        print(v$msigdbr_mm_go)
        v$pathways <- split(x = v$msigdbr_mm_go$gene_symbol, f = v$msigdbr_mm_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Homo sapiens" & input$category_gsea == "C5"){
        v$msigdbr_hs_go <- msigdbr(species = "Homo sapiens", category = "C5")
        print(v$msigdbr_hs_go)
        v$pathways <- split(x = v$msigdbr_hs_go$gene_symbol, f = v$msigdbr_hs_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Mus musculus" & input$category_gsea == "C5"){
        v$msigdbr_mm_go <- msigdbr(species = "Mus musculus", category = "C5")
        print(v$msigdbr_mm_go)
        v$pathways <- split(x = v$msigdbr_mm_go$gene_symbol, f = v$msigdbr_mm_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Homo sapiens" & input$category_gsea == "C7"){
        v$msigdbr_hs_go <- msigdbr(species = "Homo sapiens", category = "C7")
        print(v$msigdbr_hs_go)
        v$pathways <- split(x = v$msigdbr_hs_go$gene_symbol, f = v$msigdbr_hs_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Mus musculus" & input$category_gsea == "C7"){
        v$msigdbr_mm_go <- msigdbr(species = "Mus musculus", category = "C7")
        print(v$msigdbr_mm_go)
        v$pathways <- split(x = v$msigdbr_mm_go$gene_symbol, f = v$msigdbr_mm_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Homo sapiens" & input$category_gsea == "C8"){
        v$msigdbr_hs_go <- msigdbr(species = "Homo sapiens", category = "C8")
        print(v$msigdbr_hs_go)
        v$pathways <- split(x = v$msigdbr_hs_go$gene_symbol, f = v$msigdbr_hs_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      if(input$species_gsea == "Mus musculus" & input$category_gsea == "C8"){
        v$msigdbr_mm_go <- msigdbr(species = "Mus musculus", category = "C8")
        print(v$msigdbr_mm_go)
        v$pathways <- split(x = v$msigdbr_mm_go$gene_symbol, f = v$msigdbr_mm_go$gs_name)
        print(v$pathways)
        v$markers <- FindMarkers(v$scData1, ident.1 = input$gsea.ct1, ident.2 = input$gsea.ct2, min.pct = input$min_pct1, logfc.threshold = input$logfc1, test.use = input$test.use1)
        v$markers  <- v$markers %>% arrange(desc(avg_log2FC))
        print(v$markers)
        v$markers.log2FC <- v$markers$avg_log2FC
        names(v$markers.log2FC) <- row.names(v$markers)
        v$markers.log2FC <- sort(na.omit(v$markers.log2FC), decreasing = TRUE)
        print(v$markers.log2FC)
        v$fgseaRes <- fgsea(pathways = v$pathways, stats = v$markers.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
        print(v$fgseaRes)
        v$topPathwaysUp <- v$fgseaRes[ES > 0][head(order(pval), n=10), pathway]
        v$topPathwaysDown <- v$fgseaRes[ES < 0][head(order(pval), n=10), pathway]
        v$topPathways <- c(v$topPathwaysUp, rev(v$topPathwaysDown))
      }
      output$gsea.done <- renderText(paste0("Gene set enrichment done!"))
      v$isGSEAdone <- TRUE
      })
    }
  })
  
  output$gsea.select <- renderUI({
    if(is.null(v$pathways)){
      return(NULL)
    }else{
      selectInput("gsea.pathway", label = "Gene set to visualise",
                  choices = names(v$pathways))
    }
  })
  
  output$gsea_plot <- renderPlot({
    if(is.null(v$pathways) || is.null(v$isGSEAdone)){
      return(NULL)
    }else{
      withProgress(message="Generating GSEA plot...", value=0, {
        plotEnrichment(v$pathways[[input$gsea.pathway]], v$markers.log2FC) + labs(title=input$gsea.pathway)
      })
    }
  })
  
  output$gsea_plot1 <- renderPlot({
    if(is.null(v$pathways) || is.null(v$isGSEAdone)){
      return(NULL)
    }else{
      withProgress(message="Generating GSEA plot...", value=0, {
        plotGseaTable(v$pathways[v$topPathways], v$markers.log2FC, v$fgseaRes, gseaParam=0.5)
      })
    }
  })
  
  output$gsea.table <- DT::renderDataTable(
    v$fgseaRes, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
  
  output$download_gsea.table <- downloadHandler(
    filename = function(){"GSEA Results.csv"}, 
    content = function(fname){
      withProgress(message="Downloading GSEA Results...", value=0, {
        fwrite(v$fgseaRes, fname)
      })
    }
  )
  
  ##---------------Trajectory of scRNA-seq (not working)-------------------
  
  observeEvent(input$doMonocle3, {
    withProgress(message = "Running Monocle3...", value = 0.3, {
      v$scDatad <- DietSeurat(v$scData, graphs = "umap")
      v$scDatatr <- as.cell_data_set(v$scDatad)
      v$scDatatr <- cluster_cells(cds = v$scDatatr, reduction_method = "UMAP")
      v$scDatatr <- learn_graph(v$scDatatr, use_partition = TRUE)
      cell <- WhichCells(v$scData, idents = 0)
      v$scDatatr <- order_cells(v$scDatatr, reduction_method = "UMAP", root_cells = cell)
      output$Trajectory.done <- renderText(paste0("Trajectory done!"))
      v$isTrajectorydone <- TRUE
    })
  })
  
  output$Monocle3_plot <- renderPlotly({
    if(is.null(v$scData) || is.null(v$isTrajectorydone)){
      plotly_empty()
    }else{
      withProgress(message="Generating Trajectory Plot...", value=0, {
        plot_cells(cds = v$scDatatr, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
      })
    }
  })
  
  ##---------------Cell-cell communication of scRNA-seq-------------------
    
  observeEvent(input$doCC, {
    tpmFiles <- v$scData1
    if (is.null(tpmFiles)){
      v$scData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Running cell-cell communication...", value = 0.3, {
      if(input$cc1a == "primary.predict"){
      v$scData1$primary.predict -> Idents(v$scData1)
      v$cellphone <- liana_wrap(v$scData1, method = input$cc_method, resource = input$cc_resource)
      }
      if(input$cc1a == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        v$cellphone <- liana_wrap(v$scData1, method = input$cc_method, resource = input$cc_resource)
      }
      output$cc.done <- renderText(paste0("Cell-cell communication done!"))
      v$isCCdone <- TRUE
      shinyalert("Cell-cell communication done", "Cell-cell communication done", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  output$CC_plot1 <- renderPlot({
    if(is.null(v$isCCdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating cell-cell communication plot...", value=0, {
        heat_freq(v$cellphone, pallette = c("blue", "white", "red"))
        #heatmaps_plot(meta_file = 'metadata.txt', pvalues_file = 'out/pvalues.txt', count_filename = 'heatmap_count.pdf', log_filename = 'heatmap_log_count.pdf', count_network_filename = 'count_network.txt', interaction_count_filename = 'interactions_count.txt', count_network_separator = '\t', interaction_count_separator = '\t', show_rownames = T, show_colnames = T, scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11, fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0, meta_sep='\t', pvalues_sep='\t', pvalue=input$cc_pval)
      })
    }
  })
  
  output$CC.gene.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      selectInput("CC.gene1", label = "Sender",multiple = T,
                  choices = unique(v$scData1$primary.predict), selected = unique(v$scData1$primary.predict)[1])
    }
  })
  
  output$CC.gene1.select <- renderUI({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      selectInput("CC.gene2", label = "Reciever",multiple = T,
                  choices = unique(v$scData1$primary.predict), selected = unique(v$scData1$primary.predict)[2])
    }
  })
  
  output$CC_plot2 <- renderPlot({
    if(is.null(v$isCCdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating cell-cell communication plot...", value=0, {
        chord_freq(v$cellphone, source_groups = input$CC.gene1, target_groups = input$CC.gene2)
      })
    }
  })
  
  output$cc.table <- DT::renderDataTable(
    v$cellphone, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
  
  #output$cc.table1 <- DT::renderDataTable(
  #  v$pvals, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
  
  #observeEvent(input$loadexample_atacseq, {
  # withProgress(message="Loading example Data...", value=0.5, {
  #   tpmFiles_example_atac <- Read10X_h5(filename = "atac_pbmc_5k_v1/filtered_peak_bc_matrix.h5")
  #   meta_example_atac <- read.csv('atac_pbmc_5k_v1/singlecell.csv', header = T, row.names = 1)
  #   chrom_assay_example_atac <- CreateChromatinAssay(counts = tpmFiles_example_atac, sep = c(":", "-"), genome = 'hg19', fragments = 'atac_pbmc_5k_v1/fragments.tsv.gz', min.cells = 10, min.features = 200)
  #   sObj_example_atac <- CreateSeuratObject(chrom_assay_example_atac, assay = "peaks", project = input$projName4, meta.data = meta_example_atac)
  #   print(sObj_example_atac)
  #   annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  #   seqlevelsStyle(annotations_atac) <- 'UCSC'
  #   genome(annotations_atac) <- "hg19"
  #   Annotation(sObj_example_atac) <- annotations_atac
  #   sObj_example_atac <- NucleosomeSignal(object = sObj_example_atac)
  #   sObj_example_atac <- TSSEnrichment(object = sObj_example_atac, fast = FALSE)
  #   sObj_example_atac$pct_reads_in_peaks <- sObj_example_atac$peak_region_fragments / sObj_example_atac$passed_filters * 100
  #   sObj_example_atac$blacklist_ratio <- sObj_example_atac$blacklist_region_fragments / sObj_example_atac$peak_region_fragments
  #   sObj_example_atac$high.tss <- ifelse(sObj_example_atac$TSS.enrichment > 2, 'High', 'Low')
  #   sObj_example_atac$nucleosome_group <- ifelse(sObj_example_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  #   sObj_example_atac <- RunTFIDF(sObj_example_atac)
  #   sObj_example_atac <- FindTopFeatures(sObj_example_atac, min.cutoff = 'q0')
  #   sObj_example_atac <- RunSVD(sObj_example_atac)
  #   sObj_example_atac <- FindNeighbors(object = sObj_example_atac, reduction = 'lsi', dims = 2:30)
  #   sObj_example_atac <- FindClusters(object = sObj_example_atac, resolution = 0.5, verbose = FALSE, algorithm = 3)
  #   sObj_example_atac <- RunUMAP(object = sObj_example_atac, reduction = 'lsi', dims = 2:30)
  #   gene.activities <- GeneActivity(sObj_example_atac, features = VariableFeatures(v$scData1))
  #   sObj_example_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
  #   DefaultAssay(sObj_example_atac) <- "ACTIVITY"
  #   sObj_example_atac <- NormalizeData(sObj_example_atac)
  #   sObj_example_atac <- ScaleData(sObj_example_atac, features = rownames(sObj_example_atac))
  #   transfer.anchors <- FindTransferAnchors(reference = v$scData1, query = sObj_example_atac, features = VariableFeatures(object = v$scData1), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  #   celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = v$scData1$primary.predict, weight.reduction = sObj_example_atac[["lsi"]], dims = 2:30)
  #   sObj_example_atac <- AddMetaData(sObj_example_atac, metadata = celltype.predictions)
  #   v$transfer.anchors <- transfer.anchors
  #   v$atacData <- sObj_example_atac
  # })
  # label1 <- "Example loaded"
  # updateActionButton(inputId = "loadexample_atacseq", label = label1)
  # shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
  #})
  
  observeEvent(input$annotate_atacseq, {
    #tpmFiles_atac <- input$tpmFiles_atac
    tpmFiles_atac <- v$atacData
    if (is.null(tpmFiles_atac)){
      v$atacData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message="Annotating...", value=0.5, {
      transfer.anchors <- FindTransferAnchors(reference = v$scData1, query = v$atacData, reduction = "cca")
      celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = v$scData1$primary.predict, weight.reduction = v$atacData[["lsi"]], dims = 2:30)
      print(celltype.predictions)
      v$atacData <- AddMetaData(v$atacData, metadata = celltype.predictions)
      print("res1")
      print(v$atacData@meta.data)
      v$transfer.anchors <- transfer.anchors
    })
    label1 <- "Annotation done"
    updateActionButton(inputId = "loadexample_atacseq", label = label1)
    shinyalert("Annotation done", "Annotation done.", type = "success", imageWidth = 10, imageHeight = 10)
    v$isCELLAtacsdone <- TRUE
    }
  })
  
  observeEvent(input$annotate_atacseq_a, {
    tpmFiles_atac <- input$tpmFiles_atac
    tpmFiles_atac <- v$atacData
    if (is.null(tpmFiles_atac)){
      v$atacData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
      withProgress(message="Annotating...", value=0.5, {
        transfer.anchors1 <- FindTransferAnchors(reference = v$scData1, query = v$atacData, reduction = "cca")
        celltype.predictions1 <- TransferData(anchorset = transfer.anchors1, refdata = v$scData1$primary.predict, weight.reduction = v$atacData[["lsi"]], dims = 2:30)
        v$atacData <- AddMetaData(v$atacData, metadata = celltype.predictions1)
        v$transfer.anchors1 <- transfer.anchors1
      })
      label1 <- "Annotation done"
      updateActionButton(inputId = "loadexample_atacseq1", label = label1)
      shinyalert("Annotation done", "Annotation done.", type = "success", imageWidth = 10, imageHeight = 10)
    }
  })
  
  observeEvent(input$annotate1_atacseq_a, {
    tpmFiles_atac <- input$tpmFiles_atac
    tpmFiles_atac <- v$atacData
    if (is.null(tpmFiles_atac)){
      v$atacData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
      withProgress(message="Annotating...", value=0.5, {
        transfer.anchors1 <- FindTransferAnchors(reference = v$scData1, query = v$atacData, reduction = "cca")
        celltype.predictions1 <- TransferData(anchorset = transfer.anchors1, refdata = v$scData1$primary.predict, weight.reduction = v$atacData[["lsi"]], dims = 2:30)
        v$atacData <- AddMetaData(v$atacData, metadata = celltype.predictions1)
        v$transfer.anchors1 <- transfer.anchors1
      })
      label1 <- "Annotation done"
      updateActionButton(inputId = "loadexample_atacseq1", label = label1)
      shinyalert("Annotation done", "Annotation done.", type = "success", imageWidth = 10, imageHeight = 10)
    }
  })
  
  #output$annotate_scRNA_ATAC_plot <- renderPlot({
  #  if(is.null(v$scData1)|| is.null(v$atacData)){
  #   return(NULL)
  # }else{
  #   withProgress(message="Generating ATAC UMAP...", value=0, {
  #     DimPlot(v$atacData, group.by = "seurat_clusters", label = T)
  #   })
  # }
  #})
  
  #output$annotate_scRNA_ATAC_plot1 <- renderPlot({
  # if(is.null(v$scData1)|| is.null(v$atacData)){
  #   return(NULL)
  # }else{
  #   withProgress(message="Generating ATAC UMAP...", value=0, {
  #     DimPlot(v$atacData, group.by = "predicted.id", label = T)
  #   })
  # }
  #})
  
  #observe({
  # if(input$process_atacseq > 0){
  #   print('2')
  #   session$sendCustomMessage("myCallbackHandler1", "2")
  # }
  #})
  
  observeEvent(input$annotate1_atacseq, {
    tpmFiles_atac <- input$tpmFiles_atac
    tpmFiles_atac <- v$atacData
    if (is.null(tpmFiles_atac)){
      v$atacData <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message="Annotating...", value=0.5, {
      transfer.anchors <- FindTransferAnchors(reference = v$scData1, query = v$atacData, reduction = "cca")
      celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = v$scData1$primary.predict, weight.reduction = v$atacData[["lsi"]], dims = 2:30)
      v$atacData <- AddMetaData(v$atacData, metadata = celltype.predictions)
      v$transfer.anchors <- transfer.anchors
    })
    label1 <- "Annotation done"
    updateActionButton(inputId = "loadexample_atacseq", label = label1)
    shinyalert("Annotation done", "Annotation done.", type = "success", imageWidth = 10, imageHeight = 10)
    }
  })
  
  ##------------------------Data integration module--------------------------##
  
  output$integration_image <- renderImage({
    list(src = "www/Picture3.png",
         height = 500)
  }, deleteFile = FALSE)
  
  observeEvent(input$loadexample1, {
      withProgress(message="Loading example Data...", value=0.5, {
        tpmFiles1 <- read.table("integration/concatenated_expr_data.txt", header = T, row.names = 1, check.names = F, sep = '\t')
        #scH5 <- input$scH5
        annoFile1 <- read.table("integration/metadata.txt", header = T, row.names = 1, sep = '\t')
      })
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(tpmFiles1)
          #print(tpmFiles$name)
          print(file.exists(paste(tpmFiles1$datapath[1], "/", tpmFiles1$name[1], sep="")))
          v$tpmFiles1 <- tpmFiles1
          v$anno <- annoFile1
          
          sObj1 <- CreateSeuratObject(v$tpmFiles1,
                                      meta.data = v$anno,
                                      project = input$projName1,
                                      min.genes = input$min.genes1,
                                      min.cells = input$min.cells1)
          
          sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
          v$scData2 <- sObj1
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample1", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        })
    }
  })
  
  observeEvent(input$loadexample1a, {
    withProgress(message="Loading example Data...", value=0.5, {
      data1 <- Read10X_h5("integration/filtered_feature_bc_matrix_data1.h5")
      data2 <- Read10X_h5("integration/filtered_feature_bc_matrix_data2.h5")
    })
    if (is.null(data1) || is.null(data2)){
      v$scData1 <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(data1)
        print(data2)
        #print(tpmFiles$name)
        data <- list(data1, data2)
        for (i in 1:length(data)){
          data[[i]] <- CreateSeuratObject(as.matrix(data[[i]]), project = input$projName1, min.genes = input$min.genes1, min.cells = input$min.cells1)
        }
        sObj1 <- merge(data[[1]], data[[2]])
        
        sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
        sObj1[["batch"]] <- ifelse(endsWith(sObj1@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
        v$scData2 <- sObj1
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample1a", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  observeEvent(input$loadButton1, {
    if(input$scInput1 == "Raw Counts Matrix"){
      tpmFiles1 <- input$tpmFiles1
      annoFile1 <- input$cellAnnoFiles1
      names.field <- input$field
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(tpmFiles1$datapath)
          print(tpmFiles1$name)
          print(file.exists(paste(tpmFiles1$datapath[1], "/", tpmFiles1$name[1], sep="")))
          file_list <- tpmFiles1$datapath
          names(file_list) <- basename(file_list)
          exp.data1 <- imap(file_list, function(x, y) {
            cts <- x %>%
              read.table (header = T, row.names = 1, stringsAsFactors = F, check.names = F) %>%
              CreateSeuratObject(project=y)
          })
          print(exp.data1)
          exp.data2 <- merge(exp.data1[[1]], exp.data1[2:length(exp.data1)])
          
          if(!is.null(annoFile1)){
            anno.data1 <- rbindlist(lapply(annoFile1$datapath, fread), use.names = TRUE, fill = TRUE)
            anno.data2 <- data.frame(anno.data1, row.names = 1)
          }
          print(anno.data2)
          
          incProgress(0.5, "Creating Seurat Object")
          
          sObj1 <- AddMetaData(exp.data2, anno.data2)
          
          sObj1$orig.ident <- "Integration"
          Idents(sObj1) <- sObj1$orig.ident
          sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
          print(sObj1@meta.data)
          
          v$scData2 <- sObj1
          print(v$scData2@meta.data)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
    else if(input$scInput1 == "10X cellranger"){
      scH5_1 <- input$scH5_1
      if (is.null(scH5_1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(scH5_1$datapath)
          print(scH5_1$name)
          print(file.exists(paste(scH5_1$datapath[1], "/", scH5_1$name[1], sep="")))
          
          file_list <- scH5_1$datapath
          names(file_list) <- basename(file_list)
          exp.data1 <- imap(file_list, function(x, y) {
            cts <- x %>%
              Read10X_h5 %>%
              CreateSeuratObject(project=y)
          })
          print(exp.data1)
          exp.data2 <- merge(exp.data1[[1]], exp.data1[2:length(exp.data1)])
          print(exp.data2)
          incProgress(0.5, "Creating Seurat Object")
          
          exp.data2$orig.ident <- "Integration"
          Idents(exp.data2) <- exp.data2$orig.ident
          exp.data2[["percent.mt"]] <- PercentageFeatureSet(exp.data2, pattern = "^MT-")
          exp.data2[["batch"]] <- ifelse(endsWith(exp.data2@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
          v$scData2 <- exp.data2
          print(v$scData2@meta.data)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
  })
  
  output$countdataDT1 <- renderDataTable({
    if(!is.null(v$scData2))
    {
      if(ncol(v$scData2) > 20 )
        return(as.matrix(v$scData2@assays$RNA@counts[,1:20]))
    }
  }, server = FALSE)
  
  observeEvent(input$filter_seurat1, {
    if(input$scInput1 == "Raw Counts Matrix"){
      tpmFiles1 <- input$tpmFiles1
      annoFile1 <- input$cellAnnoFiles1
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
    withProgress(message="Loading and Processing Data...", value=0, {
      print(v$scData2)
      v$scData2 <- subset(v$scData2, subset = nFeature_RNA > input$ob1a & nFeature_RNA < input$ob2a & percent.mt < input$ob3a)
      
      #sObj1[["batch"]] <- ifelse(endsWith(sObj1@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
      
      shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  }
    else if(input$scInput1 == "10X cellranger"){
      scH5_1 <- input$scH5_1
      if (is.null(scH5_1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(v$scData2)
        sObj1 <- subset(v$scData2, subset = nFeature_RNA > input$ob1a & nFeature_RNA < input$ob2a & percent.mt < input$ob3a)
        sObj1$orig.ident <- "Integration"
        Idents(sObj1) <- sObj1$orig.ident
        #sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
        v$scData2 <- sObj1
        print(v$scData2@meta.data)
        shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
  })
  
  observeEvent(input$reset_intg, {
    session$reload()
    print("Reset done")
  })
  
  output$nFeature_RNAPlot1 <- renderPlot({
    if(is.null(v$scData2)){
      plotly_empty()
    }else{
      VlnPlot(v$scData2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend()
    }
  })
  
  output$mitoPlot1 <- renderPlot({
    if(is.null(v$scData2)){
      plotly_empty()
    }else{
      VlnPlot(v$scData2, "percent.mt") + NoLegend()
    }
  })
  
  output$nCount_RNAPlot1 <- renderPlot({
    if(is.null(v$scData2)){
      plotly_empty()
    }else{
      VlnPlot(v$scData2, "nCount_RNA") + NoLegend()
    }
  })
  
  output$FeatureScatterPlot1a <- renderPlotly({
    if(is.null(v$scData2)){
      plotly_empty()
    }else{
      print(FeatureScatter(v$scData2, "nCount_RNA", "nFeature_RNA"))
    }
  })
  
  output$FeatureScatterPlot2a <- renderPlotly({
    if(is.null(v$scData2)){
      plotly_empty()
    }else{
      print(FeatureScatter(v$scData2, "nCount_RNA", "percent.mt"))
    }
  })
  
  #observe({if(input$scAnalysis_integ == "Seurat" || input$scAnalysis_integ == "Harmony" || input$scAnalysis_integ == "fastMNN" || input$scAnalysis_integ == "scVI"){
  observeEvent(input$findVarGenes_bef_intg, {
    tpmFiles1 <- v$scData2
    if (is.null(tpmFiles1)){
      v$scData2 <- NULL
      shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    }else{
    withProgress(message = "Finding variable genes...", value = 0, {
        v$scData2 <- NormalizeData(v$scData2)
        v$scData2 <- FindVariableFeatures(v$scData2,
                                          mean.function = ExpMean,
                                          dispersion.function = LogVMR,
                                          nfeatures = input$var.genes_bef_intg,
                                          selection.method = input$selection.method)
        #all.genes <- rownames(v$scData1)
        v$scData2 <- ScaleData(v$scData2)
        print(v$scData2)
        incProgress(0.5)
        #VarGeneText <- paste0("Number of variable genes: ", length(v$scData1@assays$RNA@var.features))
        #output$nVarGenes <- renderText(VarGeneText)
        varGenePlotInput <- function(){
          if(is.null(v$scData2)){
            return(NULL)
          }else{
            withProgress(message="Plotting variable genes...", value=0, {
              top10 <- head(VariableFeatures(v$scData2), 10)
              variable_feature1 <- VariableFeaturePlot(v$scData2)
              variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
              print (variable_feature1)
              print (variable_feature2)
              shinyalert("Highly variable features identified", "Highly variable features identified, please perform PCA (Before Integration) ", type = "success", imageWidth = 10, imageHeight = 10)
              #dev.off()
            })
          }
        }
        output$VarGenes_bef_intg <- renderPlot({
          varGenePlotInput()
        }, height = 500, width = 600)
        observeEvent(input$PDFc, {
          if(!is.null(v$scData2)){
            withProgress(message="Downloading plot PDF files...", value=0, {
              print(getwd())
              pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
              if(!dir.exists(pdfDir)){
                dir.create(pdfDir)
              }
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
              i = 0
              while(file.exists(filename2)){
                filename2 <- paste0(pdfDir, .Platform$file.sep,
                                    "Var_genes_plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
              }
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              plot1 <- VariableFeaturePlot(v$scData2)
              print(plot1)
              dev.off()
              txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
              txtfile <- sub(".pdf", ".txt", txtfile)
              write(v$scData2@assays$RNA@var.features, file = txtfile)
              
              })
            }
          })
        })
      }
    })
    
    observeEvent(input$runPCA_bef_intg, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunPCA(v$scData2, verbose = FALSE)
          print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
          
          v$isPCAdone_bef_intg <- TRUE
          
          PCA_plot1a <- DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
          PCA_plot1b <- DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
          PCA_plot1c <- DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
          print(PCA_plot1a)
          print(PCA_plot1b)
          print(PCA_plot1c)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData2)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          shinyalert("PCA performed", "PCA performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData2 <- RunPCA(v$scData2, verbose = FALSE)
          print(v$scData2@meta.data)
          print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
          v$isPCAdone_bef_intg <- TRUE
          PCA_plot1a <- DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
          PCA_plot1b <- DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
          print(PCA_plot1a)
          print(PCA_plot1b)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData2)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          shinyalert("PCA performed", "PCA performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
        }
      }
    )
    
    output$PCAplot_bef_tpm1 <- renderPlotly({
      if(is.null(v$isPCAdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_tpm2 <- renderPlotly({
      if(is.null(v$isPCAdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_tpm3 <- renderPlotly({
      if(is.null(v$isPCAdone_bef_intg)){
        plotly_empty()
      }else{
       withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
        })
      }
     })
    
    output$PCAplot_bef_h5_1 <- renderPlotly({
      if(is.null(v$isPCAdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_h5_2 <- renderPlotly({
      if(is.null(v$isPCAdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$vizPlot_bef_intg <- renderPlot({
      if(is.null(v$scData2)){
        plotly_empty()
      }else{
        VizDimLoadings(v$scData2, dims = as.numeric(input$select.pc_bef_intg))
      }
    })
    
    output$PCHeatmap_bef_intg <- renderPlot({
      if(is.null(v$scData2)){
        plotly_empty()
      }else{
        DimHeatmap(v$scData2, dims = as.numeric(input$select.pc_bef_intg))
      }
    })
    
    output$PCtable_bef_intg <- DT::renderDataTable({
      if(is.null(v$scData2) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, server = FALSE, options = list(scrollX = TRUE))
    
    output$Elbow_bef_intg <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData2, ndims = 50)
        })
      }
    }, height = 400, width = 450)
    
    observeEvent(input$findCluster_bef_intg, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Finding clusters...", value = 0.3, {
        v$scData2 <- FindNeighbors(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg)
        v$scData2 <- FindClusters(v$scData2, resolution = input$clus.res_bef_intg)
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone1 <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please perform UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Cluster2DPlot_bef_intg <- renderPlotly({
      if(is.null(v$isClusterdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    observeEvent(input$findoptimumCluster1, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Finding optimum resolution...", value = 0.3, {
          #if(input$assays1 == "LogNormalization"){
          #DefaultAssay(v$scData1) <- "RNA"
          reses<-seq(from = input$clus.res_a1, to = input$clus.res_b1, by = 0.1)
          for (res in reses){
            v$scData2<-FindClusters(v$scData2, resolution = res)
            nores<-gsub(pattern = ".", replacement = "", res, fixed = T)
          }
          print(v$scData2@meta.data)
          output$optimumcluster.done <- renderText(paste0("Clustering done!"))
          #}
          #else if(input$assays1 == "SCTransform"){
          #  DefaultAssay(v$scData1) <- "SCT"
          #  v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, assay = "SCT", nn.method = "rann")
          #  v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
          #  output$cluster.done <- renderText(paste0("Clustering done!"))
          #}
          v$isOptimumClusterdone1 <- TRUE
          shinyalert("Number of clusters identified at different resolution", "Number of clusters identified at different resolution, please perform UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    plotOptimumCluster1 <- reactive({
      if(is.null(v$scData2) || is.null(v$isOptimumClusterdone1)){
        plotly_empty()
      }else{
        withProgress(message="Plotting Clustering tree...", value=0, {
          p <- clustree(v$scData2)
        })
      }
    })
    
    output$OptimumCluster2DPlot_a1 <- renderPlot({
      print(plotOptimumCluster1())
    })
    
    output$download_OptimumCluster1 <- downloadHandler(
      filename = function(){"Clustree plot.png"}, 
      content = function(fname){
        ggsave(fname,plotOptimumCluster1(), height = 7, width = 7)
      }
    )
    
    output$download_OptimumClusterTable1 <- downloadHandler(
      
      filename = function(){"Clustering_different_resolution.csv"}, 
      content = function(fname){
        withProgress(message="Downloading Cluster Table at different resolution...", value=0, {
          write.csv(v$scData2@meta.data, fname)
        })
      }
    )
    
    output$subcluster.gene.select1 <- renderUI({
      if(is.null(v$scData2)){
        plotly_empty()
      }else{
        selectInput("cluster_subcluster1", label = "Cluster to subcluster",
                    choices = unique(v$scData2$seurat_clusters), selected = unique(v$scData2$seurat_clusters)[1])
      }
    })
    
    observeEvent(input$subcluster1, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Performing subcluster analysis...", value=0, {
          v$scData2$seurat_clusters -> Idents(v$scData2)
          v$scData2.subset <- FindSubCluster(v$scData2, input$cluster_subcluster1, "RNA_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res1, algorithm = 1)
          print(v$scData2.subset)
          output$Subcluster1.done <- renderText(paste0("Subclustering done!"))
          v$isSubclusterdone1 <- TRUE
        })
      }
    })
    
    plotSubclster1 <- reactive({
      if(is.null(v$scData2.subset) || is.null(v$isSubclusterdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating Celltype similarity plot...", value=0, {
          p <- DimPlot(v$scData2.subset, reduction = "umap", group.by = "sub.cluster", label = T, label.size = 2.5)
        })
      }
    })
    
    output$subcluster_plot1 <- renderPlotly({
      plotSubclster1()
    })
    
    output$download_subcluster1 <- downloadHandler(
      filename = function(){"Subcluster plot.png"}, 
      content = function(fname){
        ggsave(fname,plotSubclster1(), height = 7, width = 7)
      }
    )
    
    observeEvent(input$runUMAP_bef_intg, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg)
          v$isUMAPdone_bef_intg <- TRUE
          UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          shinyalert("UMAP performed", "UMAP performed, please perform tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:30)
          v$isUMAPdone_bef_intg <- TRUE
          UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          shinyalert("UMAP performed", "UMAP performed, please perform tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    }
  })
    
    output$UMAPplot_bef_tpm1 <- renderPlotly({
      if(is.null(v$isUMAPdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_bef_tpm2 <- renderPlotly({
      if(is.null(v$isUMAPdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_bef_tpm3 <- renderPlotly({
      if(is.null(v$isUMAPdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_bef_h5_1 <- renderPlotly({
      if(is.null(v$isUMAPdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_bef_h5_2 <- renderPlotly({
      if(is.null(v$isUMAPdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    observeEvent(input$runTSNE_bef_intg, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg)
          v$isTSNEdone_bef_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          shinyalert("tSNE performed", "tSNE performed, please perform celltype annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg)
          v$isTSNEdone_bef_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          shinyalert("tSNE performed", "tSNE performed, please perform celltype annotation", type = "success", imageWidth = 10, imageHeight = 10)
          }
        })
      }
    })
    
    output$TSNEplot_bef_tpm1 <- renderPlotly({
      if(is.null(v$isTSNEdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_bef_tpm2 <- renderPlotly({
      if(is.null(v$isTSNEdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_bef_tpm3 <- renderPlotly({
      if(is.null(v$isTSNEdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_bef_h5_1 <- renderPlotly({
      if(is.null(v$isTSNEdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_bef_h5_2 <- renderPlotly({
      if(is.null(v$isTSNEdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    observeEvent(input$doCELLiD_bef_intg, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Running CELLiD...", value = 0.3, {
          ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
          v$scData2.rna.data.average = AverageExpression(v$scData2)
          v$scData2.rna.data.average = round(v$scData2.rna.data.average$RNA, 2)
          #v$scData2.rna.data.average = data.frame(v$scData2.rna.data.average$RNA)
          if(input$cellatlas1 == "all"){
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, ref)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "adipose"){
            adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
            adipose1 <- ref[,adipose]
            colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, adipose1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "adrenal_gland"){
            adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
            adrenal_gland1 <- ref[,adrenal_gland]
            colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, adrenal_gland1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "blood"){
            blood <- colnames(ref)[grepl("blood",colnames(ref))] 
            blood1 <- ref[,blood]
            colnames(blood1) <- gsub("--blood","",colnames(blood1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, blood1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "bone_marrow"){
            bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
            bone_marrow1 <- ref[,bone_marrow]
            colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, bone_marrow1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "brain"){
            brain <- colnames(ref)[grepl("brain",colnames(ref))] 
            brain1 <- ref[,brain]
            colnames(brain1) <- gsub("--brain","",colnames(brain1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, brain1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "breast"){
            breast <- colnames(ref)[grepl("breast",colnames(ref))] 
            breast1 <- ref[,breast]
            colnames(breast1) <- gsub("--breast","",colnames(breast1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, breast1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "breast_milk"){
            breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
            breast_milk1 <- ref[,breast_milk]
            colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, breast_milk1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "eye"){
            eye <- colnames(ref)[grepl("eye",colnames(ref))] 
            eye1 <- ref[,eye]
            colnames(eye1) <- gsub("--eye","",colnames(eye1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, eye1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "gut"){
            gut <- colnames(ref)[grepl("gut",colnames(ref))] 
            gut1 <- ref[,gut]
            colnames(gut1) <- gsub("--gut","",colnames(gut1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, gut1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "heart"){
            heart <- colnames(ref)[grepl("heart",colnames(ref))] 
            heart1 <- ref[,heart]
            colnames(heart1) <- gsub("--heart","",colnames(heart1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, heart1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "kidney"){
            kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
            kidney1 <- ref[,kidney]
            colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, kidney1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "liver"){
            liver <- colnames(ref)[grepl("liver",colnames(ref))] 
            liver1 <- ref[,liver]
            colnames(liver1) <- gsub("--liver","",colnames(liver1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, liver1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "lung"){
            lung <- colnames(ref)[grepl("lung",colnames(ref))] 
            lung1 <- ref[,lung]
            colnames(lung1) <- gsub("--lung","",colnames(lung1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, lung1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "pancreas"){
            pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
            pancreas1 <- ref[,pancreas]
            colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, pancreas1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "PDAC"){
            PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
            PDAC1 <- ref[,PDAC]
            colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, PDAC1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "skin"){
            skin <- colnames(ref)[grepl("skin",colnames(ref))] 
            skin1 <- ref[,skin]
            colnames(skin1) <- gsub("--skin","",colnames(skin1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, skin1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "testis"){
            testis <- colnames(ref)[grepl("testis",colnames(ref))] 
            testis1 <- ref[,testis]
            colnames(testis1) <- gsub("--testis","",colnames(testis1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, testis1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "thymus"){
            thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
            thymus1 <- ref[,thymus]
            colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, thymus1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1 == "tonsil"){
            tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
            tonsil1 <- ref[,tonsil]
            colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
            v$res_bef_intg = FastIntegration::CELLiD(v$scData2.rna.data.average, tonsil1)
            print(v$res_bef_intg)
            v$scData2$primary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res_bef_intg[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_bef_intg) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD1.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_bef_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          v$isCELLiDdone_bef_intg <- TRUE
          shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform data integration", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_cellid_bef_intg <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCELLiDdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$Umap_cellid_bef_intg1 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCELLiDdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$ct_bef_intg.table <- DT::renderDataTable(
      v$res_bef_intg, server = FALSE,  options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doCelltypist_bef_intg, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running Celltypist...", value = 0.3, {
        sc <- reticulate::import("scanpy", convert = FALSE)
        ct <- reticulate::import("celltypist", convert = FALSE)
        sceasy::convertFormat(v$scData2, from = "seurat", to = "anndata", outFile = 'ct_scrna.h5ad')
        v$adata = sc$read_h5ad('ct_scrna.h5ad')
        v$res = ct$annotate(filename = 'ct_scrna.h5ad', model = input$celltypistatlas1, majority_voting=T)
        print(v$res)
        v$adata = v$res$to_adata()
        print("fff")
        v$adata$obs$to_csv('celltypist_predict.csv')
        v$meta1 <- read.csv('celltypist_predict.csv', header = T, row.names = 1)
        v$scData2 <- AddMetaData(v$scData2, metadata = v$meta1)
        v$scData2$primary.predict <- v$scData2$majority_voting
        v$scData2$secondary.predict <- v$scData2$predicted_labels
        print(v$scData2@meta.data)
        v$isCelltypistdone_bef_intg <- TRUE
        shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform cell-cell similarity", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_celltypist_bef_intg <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCelltypistdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from Celltypist...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "majority_voting", label = T,  label.size = 3, repel = T) + NoLegend()
        })
      }
    })
    
    output$Umap_celltypist_bef_intg1 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCelltypistdone_bef_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from Celltypist...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "predicted_labels", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$ct_celltypist_bef_intg.table <- DT::renderDataTable(
      v$meta1, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doIntg_seurat, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running Data Integration...", value = 0.3, {
        v$scData1.list <- SplitObject(v$scData2, split.by = "batch")
        print(v$scData2)
        print(v$scData1.list)
        features <- SelectIntegrationFeatures(object.list = v$scData1.list, nfeatures = input$nfeatures_intg_seurat)
        print(features)
        v$scData1.anchors <- FindIntegrationAnchors(object.list = v$scData1.list, anchor.features = features)
        v$scData1.combined <- IntegrateData(anchorset = v$scData1.anchors)
        print(v$scData1.anchors)
        print(v$scData1.combined)
       
        #v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
        #v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_seurat)
        print(v$scData1.combined)
        v$isIntgSeuratdone <- TRUE
        shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    observeEvent(input$runPCA_intg_seurat, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Finding clusters...", value = 0.3, {
          DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- ScaleData(v$scData1.combined, verbose = FALSE)
          v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
          #output$pca1.done <- renderText(paste0("PCA done!"))
          v$isPCAdone_intg_seurat <- TRUE
          shinyalert("PCA done", "PCA done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }  
    })
    
    output$PCAplot_seurat_tpm1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_seurat)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_seurat_tpm2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_seurat)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_seurat_tpm3 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_seurat)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_seurat_h5_1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_seurat)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_seurat_h5_2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_seurat)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$vizPlot_intg_seurat <- renderPlot({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        VizDimLoadings(v$scData1.combined, dims = as.numeric(input$select.pc_intg_seurat))
      }
    })
    
    output$PCHeatmap_intg_seurat <- renderPlot({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        DimHeatmap(v$scData1.combined, dims = as.numeric(input$select.pc_intg_seurat))
      }
    })
    
    output$PCtable_intg_seurat <- DT::renderDataTable({
      if(is.null(v$scData1.combined) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, server = FALSE, options = list(scrollX = TRUE))
    
    output$Elbow_intg_seurat <- renderPlot({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData1.combined, ndims = 50)
        })
      }
    }, height = 400, width = 450)
    
    observeEvent(input$findCluster_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Finding clusters...", value = 0.3, {
        if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Seurat")
        {
        DefaultAssay(v$scData1.combined) <- "integrated"  
        v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg2)
        v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
        output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone_intg <- TRUE
        shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Harmony")
        {
          DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "scVI")
        {
          #DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "scvi", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "fastMNN")
        {
          #DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Seurat")
        {
          DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Harmony")
        {
          DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "scVI")
        {
          #DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "scvi", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "fastMNN")
        {
          #DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg2)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone_intg <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        }
        })
      }  
    })
    
    output$Cluster2DPlot_intg <- renderPlotly({
      if(is.null(v$isClusterdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T) + NoLegend()
        })
      }
    })
    
    observeEvent(input$findoptimumCluster2, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Finding optimum resolution...", value = 0.3, {
          #if(input$assays1 == "LogNormalization"){
          #DefaultAssay(v$scData1) <- "RNA"
          reses<-seq(from = input$clus.res_a2, to = input$clus.res_b2, by = 0.1)
          for (res in reses){
            v$scData1.combined<-FindClusters(v$scData1.combined, resolution = res)
            nores<-gsub(pattern = ".", replacement = "", res, fixed = T)
          }
          print(v$scData1.combined@meta.data)
          output$optimumcluster2.done <- renderText(paste0("Clustering done!"))
          #}
          #else if(input$assays1 == "SCTransform"){
          #  DefaultAssay(v$scData1) <- "SCT"
          #  v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, assay = "SCT", nn.method = "rann")
          #  v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
          #  output$cluster.done <- renderText(paste0("Clustering done!"))
          #}
          v$isOptimumClusterdone2 <- TRUE
          shinyalert("Number of clusters identified at different resolution", "Number of clusters identified at different resolution, please perform UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    plotOptimumCluster2 <- reactive({
      if(is.null(v$scData1.combined) || is.null(v$isOptimumClusterdone2)){
        plotly_empty()
      }else{
        withProgress(message="Plotting Clustering tree...", value=0, {
          p <- clustree(v$scData1.combined)
        })
      }
    })
    
    output$OptimumCluster2DPlot_b1 <- renderPlot({
      print(plotOptimumCluster2())
    })
    
    output$download_OptimumCluster2 <- downloadHandler(
      filename = function(){"Clustree plot.png"}, 
      content = function(fname){
        ggsave(fname,plotOptimumCluster2(), height = 7, width = 7)
      }
    )
    
    output$download_OptimumClusterTable2 <- downloadHandler(
      filename = function(){"Clustering_different_resolution.csv"}, 
      content = function(fname){
        withProgress(message="Downloading Cluster Table at different resolution...", value=0, {
          write.csv(v$scData1.combined@meta.data, fname)
        })
      }
    )
    
    output$subcluster.gene.select2 <- renderUI({
      if(is.null(v$scData1.combined)){
        plotly_empty()
      }else{
        selectInput("cluster_subcluster2", label = "Cluster to subcluster",
                    choices = unique(v$scData1.combined$seurat_clusters), selected = unique(v$scData1.combined$seurat_clusters)[1])
      }
    })
    
    observeEvent(input$subcluster2, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Performing subcluster analysis...", value=0, {
          if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Seurat")
          {
          v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
          v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "integrated_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
          print(v$scData1.combined.subset)
          output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
          v$isSubclusterdone2 <- TRUE
          }
          else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Harmony")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "integrated_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
          else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "scVI")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "RNA_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
          else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "fastMNN")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "RNA_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
          if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Seurat")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "integrated_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Harmony")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "integrated_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "scVI")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "RNA_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "fastMNN")
          {
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.subset <- FindSubCluster(v$scData1.combined, input$cluster_subcluster2, "RNA_snn", subcluster.name = "sub.cluster", resolution = input$subcluster.res2, algorithm = 1)
            print(v$scData1.combined.subset)
            output$Subcluster2.done <- renderText(paste0("Subclustering done!"))
            v$isSubclusterdone2 <- TRUE
          }
        })
      }
    })
    
    plotSubclster2 <- reactive({
      if(is.null(v$scData1.combined.subset) || is.null(v$isSubclusterdone2)){
        plotly_empty()
      }else{
        withProgress(message="Generating Celltype similarity plot...", value=0, {
          p <- DimPlot(v$scData1.combined.subset, reduction = "umap", group.by = "sub.cluster", label = T, label.size = 2.5)
        })
      }
    })
    
    output$subcluster_plot2 <- renderPlotly({
      plotSubclster2()
    })
    
    output$download_subcluster2 <- downloadHandler(
      filename = function(){"Subcluster plot.png"}, 
      content = function(fname){
        ggsave(fname,plotSubclster2(), height = 7, width = 7)
      }
    )
    
    observeEvent(input$runUMAP_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
           if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Seurat")
          {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
           }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Harmony")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "scVI")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "scvi", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "fastMNN")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Seurat")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Harmony")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "scVI")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "scvi", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "fastMNN")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg, spread = 1)
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          v$isUMAPdone_intg <- TRUE
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    }
  })
    
    output$UMAPplot_intg_tpm1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isUMAPdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_intg_tpm2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isUMAPdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_intg_tpm3 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isUMAPdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_intg_h5_1 <- renderPlot({
      if(is.null(v$scData1.combined) || is.null(v$isUMAPdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAPplot_intg_h5_2 <- renderPlot({
      if(is.null(v$scData1.combined) || is.null(v$isUMAPdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$UMAP_lisi_intg <- renderText({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Seurat")
          {
          print(v$scData1.combined@reductions$pca@cell.embeddings)
          print(v$scData1.combined@meta.data)
          v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$pca@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
          print(v$test)
          }
          else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Harmony")
          {
            print(v$scData1.combined@reductions$harmony@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$harmony@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
          else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "scVI")
          {
            print(v$scData1.combined@reductions$scvi@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$scvi@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
          else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "fastMNN")
          {
            print(v$scData1.combined@reductions$mnn@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$mnn@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Seurat")
          {
            print(v$scData1.combined@reductions$pca@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$pca@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Harmony")
          {
            print(v$scData1.combined@reductions$harmony@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$harmony@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "scVI")
          {
            print(v$scData1.combined@reductions$scvi@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$scvi@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
          else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "fastMNN")
          {
            print(v$scData1.combined@reductions$mnn@cell.embeddings)
            print(v$scData1.combined@meta.data)
            v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$mnn@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
            print(v$test)
          }
        })
      }
    })
    
    observeEvent(input$runTSNE_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Seurat")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "Harmony")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "scVI")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "scvi", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "Raw Counts Matrix" & input$scAnalysis_integ == "fastMNN")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Seurat")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "Harmony")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "scVI")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "scvi", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger" & input$scAnalysis_integ == "fastMNN")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg1)
          v$isTSNEdone_intg <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          v$isTSNEdone_intg <- TRUE
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        })
      }
    })
    
    output$TSNEplot_intg_tpm1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isTSNEdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_intg_tpm2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isTSNEdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_intg_tpm3 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isTSNEdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_intg_h5_1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isTSNEdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$TSNEplot_intg_h5_2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isTSNEdone_intg)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    observeEvent(input$doCELLiD_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Running CELLiD...", value = 0.3, {
          ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
          v$scData1.combined.rna.data.average = AverageExpression(v$scData1.combined)
          v$scData1.combined.rna.data.average = round(v$scData1.combined.rna.data.average$RNA, 2)
          if(input$cellatlas1a == "all"){
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, ref)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "adipose"){
            adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
            adipose1 <- ref[,adipose]
            colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adipose1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "adrenal_gland"){
            adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
            adrenal_gland1 <- ref[,adrenal_gland]
            colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adrenal_gland1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "blood"){
            blood <- colnames(ref)[grepl("blood",colnames(ref))] 
            blood1 <- ref[,blood]
            colnames(blood1) <- gsub("--blood","",colnames(blood1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, blood1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "bone_marrow"){
            bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
            bone_marrow1 <- ref[,bone_marrow]
            colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, bone_marrow1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "brain"){
            brain <- colnames(ref)[grepl("brain",colnames(ref))] 
            brain1 <- ref[,brain]
            colnames(brain1) <- gsub("--brain","",colnames(brain1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, brain1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "breast"){
            breast <- colnames(ref)[grepl("breast",colnames(ref))] 
            breast1 <- ref[,breast]
            colnames(breast1) <- gsub("--breast","",colnames(breast1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "breast_milk"){
            breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
            breast_milk1 <- ref[,breast_milk]
            colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast_milk1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "eye"){
            eye <- colnames(ref)[grepl("eye",colnames(ref))] 
            eye1 <- ref[,eye]
            colnames(eye1) <- gsub("--eye","",colnames(eye1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, eye1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "gut"){
            gut <- colnames(ref)[grepl("gut",colnames(ref))] 
            gut1 <- ref[,gut]
            colnames(gut1) <- gsub("--gut","",colnames(gut1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, gut1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "heart"){
            heart <- colnames(ref)[grepl("heart",colnames(ref))] 
            heart1 <- ref[,heart]
            colnames(heart1) <- gsub("--heart","",colnames(heart1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, heart1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "kidney"){
            kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
            kidney1 <- ref[,kidney]
            colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, kidney1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "liver"){
            liver <- colnames(ref)[grepl("liver",colnames(ref))] 
            liver1 <- ref[,liver]
            colnames(liver1) <- gsub("--liver","",colnames(liver1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, liver1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "lung"){
            lung <- colnames(ref)[grepl("lung",colnames(ref))] 
            lung1 <- ref[,lung]
            colnames(lung1) <- gsub("--lung","",colnames(lung1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, lung1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "pancreas"){
            pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
            pancreas1 <- ref[,pancreas]
            colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, pancreas1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "PDAC"){
            PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
            PDAC1 <- ref[,PDAC]
            colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, PDAC1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "skin"){
            skin <- colnames(ref)[grepl("skin",colnames(ref))] 
            skin1 <- ref[,skin]
            colnames(skin1) <- gsub("--skin","",colnames(skin1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, skin1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "testis"){
            testis <- colnames(ref)[grepl("testis",colnames(ref))] 
            testis1 <- ref[,testis]
            colnames(testis1) <- gsub("--testis","",colnames(testis1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, testis1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "thymus"){
            thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
            thymus1 <- ref[,thymus]
            colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, thymus1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas1a == "tonsil"){
            tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
            tonsil1 <- ref[,tonsil]
            colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
            v$res_intg = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, tonsil1)
            print(v$res_intg)
            v$scData1.combined$primary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),1]
            v$scData1.combined$secondary.predict = v$res_intg[as.numeric(v$scData1.combined$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res_intg) <- newheaders
            print(v$scData1.combined@meta.data)
            output$CELLiD2.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res_intg, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          v$isCELLiDdone1a <- TRUE
          shinyalert("Celltype annotation done", "Celltype done, please perform DEG Analysis", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_cellid_intg <- renderPlotly({
      if(is.null(v$isCELLiDdone1a)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "primary.predict", label =  T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$Umap_cellid_intg1 <- renderPlotly({
      if(is.null(v$isCELLiDdone1a)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$ct_intg.table <- DT::renderDataTable(
      v$res_intg, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doCelltypist_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running Celltypist...", value = 0.3, {
        sc <- reticulate::import("scanpy", convert = FALSE)
        ct <- reticulate::import("celltypist", convert = FALSE)
        sceasy::convertFormat(v$scData1.combined, from = "seurat", to = "anndata", outFile = 'ct_scrna.h5ad')
        v$adata = sc$read_h5ad('ct_scrna.h5ad')
        v$res = ct$annotate(filename = 'ct_scrna.h5ad', model = input$celltypistatlas5, majority_voting=T)
        print(v$res)
        v$adata = v$res$to_adata()
        print("fff")
        v$adata$obs$to_csv('celltypist_predict.csv')
        v$meta5 <- read.csv('celltypist_predict.csv', header = T, row.names = 1)
        v$scData1.combined <- AddMetaData(v$scData1.combined, metadata = v$meta5)
        v$scData1.combined$primary.predict <- v$scData1.combined$majority_voting
        v$scData1.combined$secondary.predict <- v$scData1.combined$predicted_labels
        print(v$scData1.combined@meta.data)
        v$isCelltypistdone5 <- TRUE
        shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform cell-cell similarity", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_celltypist_intg <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isCelltypistdone5)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from Celltypist...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "majority_voting", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$Umap_celltypist_intg1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isCelltypistdone5)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from Celltypist...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "predicted_labels", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$ct_celltypist_intg.table <- DT::renderDataTable(
      v$meta5, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doDeg_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          if(input$deg1a == "seurat_clusters"){
          v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
          ips.markers_intg <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg, logfc.threshold = input$logfc_intg, test.use = input$test.use_intg)
          v$ips.markers_intg <- ips.markers_intg
          shinyalert("DEG Analysis done", "DEG Analysis done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
          }
          if(input$deg1a == "Predicted celltype"){
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            ips.markers_intg <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg, logfc.threshold = input$logfc_intg, test.use = input$test.use_intg)
            v$ips.markers_intg <- ips.markers_intg
            shinyalert("DEG Analysis done", "DEG Analysis done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
          }
          if(input$deg1a == "Metadata celltype"){
            v$scData1.combined$celltype -> Idents(v$scData1.combined)
            ips.markers_intg <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg, logfc.threshold = input$logfc_intg, test.use = input$test.use_intg)
            v$ips.markers_intg <- ips.markers_intg
            shinyalert("DEG Analysis done", "DEG Analysis done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
          }
        })
      }
    })
    
    observeEvent(input$Vis_intg, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Visualizing...", value=0, {
          v$isVisIntgdone <- TRUE
        })
      }
    })
    
    output$deg.gene.select_intg <- renderUI({
      if(is.null(v$ips.markers_intg)){
        return(NULL)
      }else{
        selectInput("deg.gene_intg", label = "Gene to visualise",
                    choices = unique(v$ips.markers_intg$gene),  selected = rownames(v$ips.markers_intg$gene)[10])
      }
    })
    
    plotViolin_intg <- reactive({
      if(is.null(v$ips.markers_intg) || is.null(v$isVisIntgdone)){
        plotly_empty()
      }else{
        withProgress(message="Visualizing...", value=0, {
          if(input$deg1b == "seurat_clusters"){
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            VlnPlot(v$scData1.combined, input$deg.gene_intg, group.by = "seurat_clusters")
          }
          else if(input$deg1b == "Predicted celltype"){
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            VlnPlot(v$scData1.combined, input$deg.gene_intg, group.by = "primary.predict")
          }
          else if(input$deg1b == "Metadata celltype"){
            v$scData1.combined$celltype -> Idents(v$scData1.combined)
            VlnPlot(v$scData1.combined, input$deg.gene_intg, group.by = "celltype")
          }
        })
      }
    })
    
    output$Deg.plot_intg <- renderPlotly({
      plotViolin_intg()
    })
    
    output$download_violn_intg <- downloadHandler(
      filename = function(){"Violin plot (Integration).png"}, 
      content = function(fname){
        ggsave(fname,plotViolin_intg(), height = 7, width = 7)
      }
    )
    
    plotFeature_intg <- reactive({
      if(is.null(v$ips.markers_intg) || is.null(v$isVisIntgdone)){
        plotly_empty()
      }else{
        withProgress(message="Visualizing...", value=0, {
          if(input$deg1b == "seurat_clusters"){
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            FeaturePlot(v$scData1.combined, input$deg.gene_intg)
          }
          else if(input$deg1b == "Predicted celltype"){
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            FeaturePlot(v$scData1.combined, input$deg.gene_intg)
          }
          else if(input$deg1b == "Metadata celltype"){
            v$scData1.combined$celltype -> Idents(v$scData1.combined)
            FeaturePlot(v$scData1.combined, input$deg.gene_intg)
          }
        })
      }
    })
    
    output$Deg1.plot_intg <- renderPlotly({
      plotFeature_intg()
    })
    
    output$download_feature_intg <- downloadHandler(
      filename = function(){"Feature plot (Integration).png"}, 
      content = function(fname){
        ggsave(fname,plotFeature_intg(), height = 7, width = 7)
      }
    )
    
    plotRidge_intg <- reactive({
      if(is.null(v$ips.markers_intg) || is.null(v$isVisIntgdone)){
        plotly_empty()
      }else{
        withProgress(message="Visualizing...", value=0, {
          if(input$deg1b == "seurat_clusters"){
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            RidgePlot(v$scData1.combined, input$deg.gene_intg, group.by = "seurat_clusters")
          }
          else if(input$deg1b == "Predicted celltype"){
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            RidgePlot(v$scData1.combined, input$deg.gene_intg, group.by = "primary.predict")
          }
          else if(input$deg1b == "Metadata celltype"){
            v$scData1.combined$celltype -> Idents(v$scData1.combined)
            RidgePlot(v$scData1.combined, input$deg.gene_intg, group.by = "celltype")
          }
        })
      }
    })
    
    output$Deg2.plot_intg <- renderPlot({
      plotRidge_intg()
    })
    
    output$download_ridge_intg <- downloadHandler(
      filename = function(){"Ridge plot (Integration).png"}, 
      content = function(fname){
        ggsave(fname,plotRidge_intg(), height = 7, width = 7)
      }
    )
    
    output$Deg3.plot_intg <- renderPlotly({
      if(is.null(v$ips.markers_intg)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          v$ips.markers_intg %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData1.combined, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
        })
      }
    })
    
    output$Deg.table_intg <- DT::renderDataTable(
      v$ips.markers_intg, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    }
  })
    
    output$gene1a.select <- renderUI({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        if(input$deg3a == "seurat_clusters"){
          selectInput("gene1a", label = "Celltype1",
                      choices = as.vector(v$scData1.combined$seurat_clusters))
        }
        else if(input$deg3a == "Predicted celltype"){
          selectInput("gene1a", label = "Celltype1",
                      choices = as.vector(v$scData1.combined$primary.predict))
        }
        else if(input$deg3a == "Metadata Celltype"){
          selectInput("gene1a", label = "Celltype1",
                      choices = as.vector(v$scData1.combined$celltype))
        }
      }
    })
    
    output$gene2a.select <- renderUI({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        if(input$deg3a == "seurat_clusters"){
          selectInput("gene2a", label = "Celltype2",
                      choices = as.vector(v$scData1.combined$seurat_clusters))
        }
        else if(input$deg3a == "Predicted celltype"){
          selectInput("gene2a", label = "Celltype2",
                      choices = as.vector(v$scData1.combined$primary.predict))
        }
        else if(input$deg3a == "Metadata Celltype"){
          selectInput("gene2a", label = "Celltype2",
                      choices = as.vector(v$scData1.combined$celltype))
        }
      }
    })
    
    observeEvent(input$doVolcano1, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          if (input$deg3a == "seurat_clusters"){ 
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined_subset <- subset(v$scData1.combined, subset = seurat_clusters == input$gene1a | seurat_clusters == input$gene2a)
            ips.markers_a1 <- FindAllMarkers(v$scData1.combined_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
            ips.markers_b1 <- FindMarkers(v$scData1.combined_subset, ident.1 = input$gene1a, ident.2 = input$gene2a, only.pos = F, min.pct = input$min_pct_a1, logfc.threshold = input$logfc_a1, test.use = input$test.use_a1)
            v$ips.markers_a1 <- ips.markers_a1
            v$ips.markers_b1 <- ips.markers_b1
          }
          if (input$deg3a == "Predicted celltype"){ 
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            v$scData1.combined_subset <- subset(v$scData1.combined, subset = primary.predict == input$gene1a | primary.predict == input$gene2a)
            ips.markers_a1 <- FindAllMarkers(v$scData1.combined_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
            ips.markers_b1 <- FindMarkers(v$scData1.combined_subset, ident.1 = input$gene1a, ident.2 = input$gene2a, only.pos = F, min.pct = input$min_pct_a1, logfc.threshold = input$logfc_a1, test.use = input$test.use_a1)
            v$ips.markers_a1 <- ips.markers_a1
            v$ips.markers_b1 <- ips.markers_b1
          }
          if (input$deg3a == "Metadata celltype"){ 
            v$scData1.combined$celltype -> Idents(v$scData1.combined)
            v$scData1.combined_subset <- subset(v$scData1.combined, subset = celltype == input$gene1a | celltype == input$gene2a)
            ips.markers_a1 <- FindAllMarkers(v$scData1.combined_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
            ips.markers_b1 <- FindMarkers(v$scData1.combined_subset, ident.1 = input$gene1a, ident.2 = input$gene2a, only.pos = F, min.pct = input$min_pct_a1, logfc.threshold = input$logfc_a1, test.use = input$test.use_a1)
            v$ips.markers_a1 <- ips.markers_a1
            v$ips.markers_b1 <- ips.markers_b1
          }
          shinyalert("Pairwise DEGs done", "Pairwise DEGs done, please run GSEA Analysis", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$volcano.plot1 <- renderPlot({
      if(is.null(v$ips.markers_b1)){
        return(NULL)
      }else{
        withProgress(message="Generating Volcano Plot...", value=0, {
          EnhancedVolcano(toptable = v$ips.markers_b1, lab = row.names(v$ips.markers_b1), x ="avg_log2FC", y ="p_val_adj", pointSize = 1, labSize = 5, legendLabSize = 12, axisLabSize = 12)
        })
      }
    })
    
    output$dega1.plot <- renderPlot({
      if(is.null(v$ips.markers_a1)){
        return(NULL)
      }else{
        withProgress(message="Generating Volcano Plot...", value=0, {
          v$ips.markers_a1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData1.combined_subset, features = top10$gene, size = 5, angle = 45) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
        })
      }
    })
    
    output$gsea.ct1a.select <- renderUI({
      if(is.null(v$ips.markers_intg)){
        return(NULL)
      }else{
        selectInput("gsea.ct1a", label = "Celltype1",
                    choices = as.vector(v$scData1.combined$primary.predict))
      }
    })
    
    output$gsea.ct2a.select <- renderUI({
      if(is.null(v$ips.markers_intg)){
        return(NULL)
      }else{
        selectInput("gsea.ct2a", label = "Celltype2",
                    choices = as.vector(v$scData1.combined$primary.predict))
      }
    })
    
    observeEvent(input$gsea1, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Generating gene set enrichment analysis...", value=0, {
          if(input$species_gsea1 == "Homo sapiens" & input$category_gsea1 == "H"){
            v$msigdbr_hs_go1 <- msigdbr(species = "Homo sapiens", category = "H")
            print(v$msigdbr_hs_go1)
            v$pathways1 <- split(x = v$msigdbr_hs_go1$gene_symbol, f = v$msigdbr_hs_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Mus musculus" & input$category_gsea1 == "H"){
            v$msigdbr_mm_go1 <- msigdbr(species = "Mus musculus", category = "H")
            print(v$msigdbr_mm_go1)
            v$pathways1 <- split(x = v$msigdbr_mm_go1$gene_symbol, f = v$msigdbr_mm_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Homo sapiens" & input$category_gsea1 == "C2"){
            v$msigdbr_hs_go1 <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
            print(v$msigdbr_hs_go1)
            v$pathways1 <- split(x = v$msigdbr_hs_go1$gene_symbol, f = v$msigdbr_hs_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Mus musculus" & input$category_gsea1 == "C2"){
            v$msigdbr_mm_go1 <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
            print(v$msigdbr_mm_go1)
            v$pathways1 <- split(x = v$msigdbr_mm_go1$gene_symbol, f = v$msigdbr_mm_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Homo sapiens" & input$category_gsea1 == "C5"){
            v$msigdbr_hs_go1 <- msigdbr(species = "Homo sapiens", category = "C5")
            print(v$msigdbr_hs_go1)
            v$pathways1 <- split(x = v$msigdbr_hs_go1$gene_symbol, f = v$msigdbr_hs_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Mus musculus" & input$category_gsea1 == "C5"){
            v$msigdbr_mm_go1 <- msigdbr(species = "Mus musculus", category = "C5")
            print(v$msigdbr_mm_go1)
            v$pathways1 <- split(x = v$msigdbr_mm_go1$gene_symbol, f = v$msigdbr_mm_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Homo sapiens" & input$category_gsea1 == "C7"){
            v$msigdbr_hs_go1 <- msigdbr(species = "Homo sapiens", category = "C7")
            print(v$msigdbr_hs_go1)
            v$pathways1 <- split(x = v$msigdbr_hs_go1$gene_symbol, f = v$msigdbr_hs_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Mus musculus" & input$category_gsea1 == "C7"){
            v$msigdbr_mm_go1 <- msigdbr(species = "Mus musculus", category = "C7")
            print(v$msigdbr_mm_go1)
            v$pathways1 <- split(x = v$msigdbr_mm_go1$gene_symbol, f = v$msigdbr_mm_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Homo sapiens" & input$category_gsea1 == "C8"){
            v$msigdbr_hs_go1 <- msigdbr(species = "Homo sapiens", category = "C8")
            print(v$msigdbr_hs_go1)
            v$pathways1 <- split(x = v$msigdbr_hs_go1$gene_symbol, f = v$msigdbr_hs_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          if(input$species_gsea1 == "Mus musculus" & input$category_gsea1 == "C8"){
            v$msigdbr_mm_go1 <- msigdbr(species = "Mus musculus", category = "C8")
            print(v$msigdbr_mm_go1)
            v$pathways1 <- split(x = v$msigdbr_mm_go1$gene_symbol, f = v$msigdbr_mm_go1$gs_name)
            print(v$pathways1)
            v$markers1 <- FindMarkers(v$scData1.combined, ident.1 = input$gsea.ct1a, ident.2 = input$gsea.ct2a, min.pct = input$min_pct1a, logfc.threshold = input$logfc1a, test.use = input$test.use1a)
            v$markers1  <- v$markers1 %>% arrange(desc(avg_log2FC))
            print(v$markers1)
            v$markers1.log2FC <- v$markers1$avg_log2FC
            names(v$markers1.log2FC) <- row.names(v$markers1)
            v$markers1.log2FC <- sort(na.omit(v$markers1.log2FC), decreasing = TRUE)
            print(v$markers1.log2FC)
            v$fgseaRes1 <- fgsea(pathways = v$pathways1, stats = v$markers1.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes1)
            v$topPathwaysUp1 <- v$fgseaRes1[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown1 <- v$fgseaRes1[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways1 <- c(v$topPathwaysUp1, rev(v$topPathwaysDown1))
          }
          output$gsea.done <- renderText(paste0("Gene set enrichment done!"))
          v$isGSEAdone1 <- TRUE
        })
      }
    })
    
    output$gsea1.select <- renderUI({
      if(is.null(v$pathways1)){
        return(NULL)
      }else{
        selectInput("gsea.pathway1", label = "Gene set to visualise",
                    choices = names(v$pathways1))
      }
    })
    
    output$gsea_plot1a <- renderPlot({
      if(is.null(v$pathways1) || is.null(v$isGSEAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating GSEA plot...", value=0, {
          plotEnrichment(v$pathways1[[input$gsea.pathway1]], v$markers1.log2FC) + labs(title=input$gsea.pathway1)
        })
      }
    })
    
    output$gsea_plot1b <- renderPlot({
      if(is.null(v$pathways1) || is.null(v$isGSEAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating GSEA plot...", value=0, {
          plotGseaTable(v$pathways1[v$topPathways1], v$markers1.log2FC, v$fgseaRes1, gseaParam=0.5)
        })
      }
    })
    
    output$gsea.table1 <- DT::renderDataTable(
      v$fgseaRes1, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
    
    output$download_gsea.table1 <- downloadHandler(
      filename = function(){"GSEA Results.csv"}, 
      content = function(fname){
        withProgress(message="Downloading GSEA Results...", value=0, {
          fwrite(v$fgseaRes1, fname)
        })
      }
    )
    
    observeEvent(input$doCC1, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Running cell-cell communication...", value = 0.3, {
          if(input$cc1b == "primary.predict"){
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            v$cellphone1 <-  liana_wrap(v$scData1.combined, method = input$cc_method1, resource = input$cc_resource1, assay = "RNA")
          }
          if(input$cc1b == "seurat_clusters"){
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$cellphone1 <-  liana_wrap(v$scData1.combined, method = input$cc_method1, resource = input$cc_resource1, assay = "RNA")
          }
          output$cc1.done <- renderText(paste0("Cell-cell communication done!"))
          v$isCC1done <- TRUE
          shinyalert("Cell-cell communication done", "Cell-cell communication done", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$CC.gene.select1 <- renderUI({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        selectInput("CC.gene1a", label = "Source",
                    choices = unique(v$scData1.combined$primary.predict), selected = unique(v$scData1.combined$primary.predict)[1], multiple = T, selectize = T)
      }
    })
    
    output$CC.gene1.select1 <- renderUI({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        selectInput("CC.gene2a", label = "Receptor",
                    choices = unique(v$scData1.combined$primary.predict), selected = unique(v$scData1.combined$primary.predict)[2], multiple = T, selectize = T)
      }
    })
    
    output$CC_plot1a <- renderPlot({
      if(is.null(v$isCC1done)){
        return(NULL)
      }else{
        withProgress(message="Generating cell-cell communication plot...", value=0, {
          heat_freq(v$cellphone1, pallette = c("blue", "white", "red"))
          #heatmaps_plot(meta_file = 'metadata.txt', pvalues_file = 'out/pvalues.txt', count_filename = 'heatmap_count.pdf', log_filename = 'heatmap_log_count.pdf', count_network_filename = 'count_network.txt', interaction_count_filename = 'interactions_count.txt', count_network_separator = '\t', interaction_count_separator = '\t', show_rownames = T, show_colnames = T, scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11, fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0, meta_sep='\t', pvalues_sep='\t', pvalue=input$cc_pval1)
        })
      }
    })
    
    output$CC_plot2a <- renderPlot({
      if(is.null(v$isCC1done)){
        return(NULL)
      }else{
        withProgress(message="Generating cell-cell communication plot...", value=0, {
          chord_freq(v$cellphone1, source_groups = input$CC.gene1a, target_groups = input$CC.gene2a)
                   
        })
      }
    })
    
    output$cc.table1a <- DT::renderDataTable(
      v$cellphone1, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
    
    #observeEvent(input$loadexample_atacseq1, {
    #  withProgress(message="Loading example Data...", value=0.5, {
    #    tpmFiles_example_atac1 <- Read10X_h5(filename = "atac_pbmc_5k_v1/filtered_peak_bc_matrix.h5")
    #    meta_example_atac1 <- read.csv('atac_pbmc_5k_v1/singlecell.csv', header = T, row.names = 1)
    #   chrom_assay_example_atac1 <- CreateChromatinAssay(counts = tpmFiles_example_atac1, sep = c(":", "-"), genome = 'hg19', fragments = 'atac_pbmc_5k_v1/fragments.tsv.gz', min.cells = 10, min.features = 200)
    #   sObj_example_atac1 <- CreateSeuratObject(chrom_assay_example_atac1, assay = "peaks", project = input$projName4, meta.data = meta_example_atac1)
    #   print(sObj_example_atac1)
    #   annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    #   seqlevelsStyle(annotations_atac) <- 'UCSC'
    #   genome(annotations_atac) <- "hg19"
    #   Annotation(sObj_example_atac1) <- annotations_atac
    #   sObj_example_atac1 <- NucleosomeSignal(object = sObj_example_atac1)
    #   sObj_example_atac1 <- TSSEnrichment(object = sObj_example_atac1, fast = FALSE)
    #   sObj_example_atac1$pct_reads_in_peaks <- sObj_example_atac1$peak_region_fragments / sObj_example_atac1$passed_filters * 100
    #   sObj_example_atac1$blacklist_ratio <- sObj_example_atac1$blacklist_region_fragments / sObj_example_atac1$peak_region_fragments
    #   sObj_example_atac1$high.tss <- ifelse(sObj_example_atac1$TSS.enrichment > 2, 'High', 'Low')
    #   sObj_example_atac1$nucleosome_group <- ifelse(sObj_example_atac1$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    #   sObj_example_atac1 <- RunTFIDF(sObj_example_atac1)
    #   sObj_example_atac1 <- FindTopFeatures(sObj_example_atac1, min.cutoff = 'q0')
    #   sObj_example_atac1 <- RunSVD(sObj_example_atac1)
    #   sObj_example_atac1 <- FindNeighbors(object = sObj_example_atac1, reduction = 'lsi', dims = 2:30)
    #   sObj_example_atac1 <- FindClusters(object = sObj_example_atac1, resolution = 0.5, verbose = FALSE, algorithm = 3)
    #   sObj_example_atac1 <- RunUMAP(object = sObj_example_atac1, reduction = 'lsi', dims = 2:30)
    #   gene.activities1 <- GeneActivity(sObj_example_atac1, features = VariableFeatures(v$scData1.combined))
    #   sObj_example_atac1[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities1)
    #   DefaultAssay(sObj_example_atac1) <- "ACTIVITY"
    #   sObj_example_atac1 <- NormalizeData(sObj_example_atac1)
    #   sObj_example_atac1 <- ScaleData(sObj_example_atac1, features = rownames(sObj_example_atac1))
    #   transfer.anchors1 <- FindTransferAnchors(reference = v$scData1.combined, query = sObj_example_atac1, features = VariableFeatures(object = v$scData1.combined), reference.assay = "integrated", query.assay = "ACTIVITY", reduction = "cca")
    #   celltype.predictions1 <- TransferData(anchorset = transfer.anchors1, refdata = v$scData1.combined$primary.predict, weight.reduction = sObj_example_atac1[["lsi"]], dims = 2:30)
    #   sObj_example_atac1 <- AddMetaData(sObj_example_atac1, metadata = celltype.predictions1)
    #   v$transfer.anchors1 <- transfer.anchors1
    #   v$atacData1 <- sObj_example_atac1
    # })
    # label1 <- "Example loaded"
    # updateActionButton(inputId = "loadexample_atacseq", label = label1)
    # shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
    #})
    
    #output$annotate_scRNA_ATAC_plot_a <- renderPlot({
    # if(is.null(v$scData1.combined)|| is.null(v$atacData1)){
    #   return(NULL)
    # }else{
    #   withProgress(message="Generating ATAC UMAP...", value=0, {
    #     DimPlot(v$atacData1, group.by = "seurat_clusters", label = T)
    #   })
    # }
    #})
    
    #output$annotate_scRNA_ATAC_plot_1a <- renderPlot({
    #   if(is.null(v$scData1.combined)|| is.null(v$atacData1)){
    #    return(NULL)
    #  }else{
    #    withProgress(message="Generating ATAC UMAP...", value=0, {
    #      DimPlot(v$atacData1, group.by = "predicted.id", label = T)
    #    })
    #  }
    #})
    
    #observe({
    #  if(input$process_atacseq1 > 0){
    #    print('2')
    #    session$sendCustomMessage("myCallbackHandler1a", "2")
    #  }
    #})
    
    observeEvent(input$annotate1_atacseq, {
      tpmFiles_atac <- input$tpmFiles_atac
      if (is.null(tpmFiles_atac)){
        v$atacData1 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Annotating...", value=0.5, {
          transfer.anchors1 <- FindTransferAnchors(reference = v$scData1.combined, query = v$atacData, reduction = "cca")
          celltype.predictions1 <- TransferData(anchorset = transfer.anchors1, refdata = v$scData1.combined$primary.predict, weight.reduction = v$atacData1[["lsi"]], dims = 2:30)
          v$atacData1 <- AddMetaData(v$atacData1, metadata = celltype.predictions1)
          v$transfer.anchors1 <- transfer.anchors1
        })
        label1 <- "Annotation done"
        updateActionButton(inputId = "loadexample_atacseq1", label = label1)
        shinyalert("Annotation done", "Annotation done.", type = "success", imageWidth = 10, imageHeight = 10)
      }
    })
    
    
    observeEvent(input$doIntg_harmony, {
      tpmFiles1 <- v$scData2
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running Data Integration...", value = 0.3, {
        v$scData1.list <- SplitObject(v$scData2, split.by = "batch")
        print(v$scData2)
        print(v$scData1.list)
        v$scData1.list <- pbmclapply(mc.cores = 20, X = v$scData1.list, FUN = function(x) {
          x <- NormalizeData(x)
          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
        })
        features <- SelectIntegrationFeatures(object.list = v$scData1.list, nfeatures = input$nfeatures_intg_harmony)
        print(features)
        v$scData1.anchors <- FindIntegrationAnchors(object.list = v$scData1.list, anchor.features = features)
        v$scData1.combined <- IntegrateData(anchorset = v$scData1.anchors)
        print(v$scData1.anchors)
        print(v$scData1.combined)
        #v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg)
        #v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
        #print(v$scData1.combined)
        v$isIntgHarmonydone <- TRUE
        shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    observeEvent(input$runPCA_intg_harmony, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Running PCA...", value = 0.3, {
          
          DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- ScaleData(v$scData1.combined, verbose = FALSE)
          v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
          print(v$scData1.combined[["pca"]], dims = 1:5, nfeatures = 5)
          v$scData1.combined <- RunHarmony(v$scData1.combined, "batch", assay.use = "integrated")
          print(v$scData1.combined[["harmony"]], dims = 1:5, nfeatures = 5)
          #output$pca1.done <- renderText(paste0("PCA done!"))
          v$isPCAdone_intg_harmony <- TRUE
          shinyalert("PCA done", "PCA done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }  
    })
    
    output$PCAplot_harmony_tpm1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_harmony)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T) + NoLegend()
        })
      }
    })
    
    output$PCAplot_harmony_tpm2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_harmony)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_harmony_tpm3 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_harmony)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_harmony_h5_1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_harmony)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, label.size = 3) + NoLegend()
        })
      }
    })
    
    output$PCAplot_harmony_h5_2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_harmony)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'batch', label.size = 3) + NoLegend()
        })
      }
    })
    
    output$vizPlot_intg_harmony <- renderPlot({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        VizDimLoadings(v$scData1.combined, dims = as.numeric(input$select.pc_intg_harmony), reduction = "harmony")
      }
    })
    
    output$PCHeatmap_intg_harmony <- renderPlot({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        DimHeatmap(v$scData1.combined, dims = as.numeric(input$select.pc_intg_harmony), reduction = "harmony")
      }
    })
    
    output$PCtable_intg_harmony <- DT::renderDataTable({
      if(is.null(v$scData1.combined) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, server = FALSE, options = list(scrollX = TRUE))
    
    output$Elbow_intg_harmony <- renderPlot({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData1.combined, reduction = "harmony", ndims = 50)
        })
      }
    }, height = 400, width = 450)
    
    observeEvent(input$doDeg_intg_harmony, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          v$scData1.combined$celltype -> Idents(v$scData1.combined)
          ips.markers1 <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg_harmony, logfc.threshold = input$logfc_intg_harmony, test.use = input$test.use_intg_harmony)
          v$ips.markers1 <- ips.markers1
          shinyalert("DEGs estimated", "DEGs estimated, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$deg.gene.select_intg_harmony <- renderUI({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        selectInput("deg.gene_intg_harmony", label = "Gene to visualise",
                    choices = unique(v$ips.markers1$gene))
      }
    })
    
    output$Deg.plot_intg_harmony <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          VlnPlot(v$scData1.combined, input$deg.gene_intg_harmony)
        })
      }
    })
    
    output$Deg1.plot_intg_harmony <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          FeaturePlot(v$scData1.combined, input$deg.gene_intg_harmony)
        })
      }
    })
    
    output$Deg2.plot_intg_harmony <- renderPlot({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          RidgePlot(v$scData1.combined, features = input$deg.gene_intg_harmony)
        })
      }
    })
    
    output$Deg3.plot_intg_harmony <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          v$ips.markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData1.combined, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
        })
      }
    })
    
    output$Deg.table_intg_harmony <- DT::renderDataTable(
      v$ips.markers1, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    #observe({if(input$scAnalysis_integ == "scVI"){
    
      observeEvent(input$doIntg_scvi, {
        tpmFiles1 <- v$scData2
        if (is.null(tpmFiles1)){
          v$scData2 <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running Data Integration...", value = 0.3, {
          sc <- reticulate::import("scanpy", convert = FALSE)
          scvi <- reticulate::import("scvi", convert = FALSE)
          v$adata <- convertFormat(v$scData2, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
          print(v$adata)
          
          scvi$model$SCVI$setup_anndata(v$adata, batch_key = 'batch')
          # create the model
          model = scvi$model$SCVI(v$adata, n_latent = input$scvi_latent)
          print(model)
          # train the model
          model$train()
          print(model)
          v$latent = model$get_latent_representation()
          print(v$latent)
          # put it back in our original Seurat object
          v$latent <- as.matrix(v$latent)
          print(v$latent)
          rownames(v$latent) = colnames(v$scData2)
          #v$scData2[["scvi"]] <- CreateDimReducObject(embeddings = v$latent, key = "scvi_", assay = DefaultAssay(v$scData2))
          #print(v$scData2)
          v$isIntgscVIdone <- TRUE
          shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
      observeEvent(input$runPCA_intg_scvi, {
        #tpmFiles1 <- v$scData1.combined
        tpmFiles1 <- v$scData2
        if (is.null(tpmFiles1)){
          v$scData2 <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message = "Running PCA...", value = 0.3, {
            v$scData2[["scvi"]] <- CreateDimReducObject(embeddings = v$latent, key = "scvi_", assay = DefaultAssay(v$scData2))
            print(v$scData2)
            v$scData2 = v$scData1.combined
            print(v$scData1.combined)
            #output$pca1.done <- renderText(paste0("PCA done!"))
            v$isPCAdone_intg_scvi <- TRUE
            shinyalert("PCA done", "PCA done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }  
      })
  
    
    output$PCAplot_scvi_tpm1 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isPCAdone_intg_scvi)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "scvi", label = T)
        })
      }
    })
    
    output$PCAplot_scvi_tpm2 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isPCAdone_intg_scvi)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "scvi", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$PCAplot_scvi_tpm3 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isPCAdone_intg_scvi)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "scvi", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$PCAplot_scvi_h5_1 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isPCAdone_intg_scvi)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "scvi", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_scvi_h5_2 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isPCAdone_intg_scvi)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "scvi", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$vizPlot_intg_scvi <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        VizDimLoadings(v$scData1.combined, dims = as.numeric(input$select.pc_intg_scvi), reduction = "scvi")
      }
    })
    
    output$PCHeatmap_intg_scvi <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        DimHeatmap(v$scData1.combined, dims = as.numeric(input$select.pc_intg_scvi), reduction = "scvi")
      }
    })
    
    output$PCtable_intg_scvi <- DT::renderDataTable({
      if(is.null(v$scData2) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, server = FALSE, options = list(scrollX = TRUE))
    
    observeEvent(input$doDeg_intg_scvi, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData2 <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
          ips.markers1 <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg_scvi, logfc.threshold = input$logfc_intg_scvi, test.use = input$test.use_intg_scvi)
          v$ips.markers1 <- ips.markers1
          shinyalert("DEGs estimated", "DEGs estimated, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$deg.gene.select_intg_scvi <- renderUI({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        selectInput("deg.gene_intg_scvi", label = "Gene to visualise",
                    choices = unique(v$ips.markers1$gene))
      }
    })
    
    output$Deg.plot_intg_scvi <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          VlnPlot(v$scData1.combined, input$deg.gene_intg_scvi)
        })
      }
    })
    
    output$Deg1.plot_intg_scvi <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          FeaturePlot(v$scData1.combined, input$deg.gene_intg_scvi)
        })
      }
    })
    
    output$Deg2.plot_intg_scvi <- renderPlot({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          RidgePlot(v$scData1.combined, features = input$deg.gene_intg_scvi)
        })
      }
    })
    
    output$Deg3.plot_intg_scvi <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          v$ips.markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData1.combined, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
        })
      }
    })
    
    output$Deg.table_intg_scvi <- DT::renderDataTable(
      v$ips.markers1, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doIntg_fastmnn, {
        tpmFiles1 <- v$scData2
        if (is.null(tpmFiles1)){
          v$scData2 <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running Data Integration...", value = 0.3, {
          v$scData1.list <- SplitObject(v$scData2, split.by = "batch")
          v$scData1.list <- pbmclapply(mc.cores = 20, X = v$scData1.list, FUN = function(x) {
            x <- NormalizeData(x)
            x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
          })
          v$features <- SelectIntegrationFeatures(object.list = v$scData1.list, nfeatures = input$nfeatures_intg_fastmnn)
          print(v$features)
          print(v$scData2)
          print(v$scData1.list)
          #v$scData1.combined <- RunFastMNN(object.list = v$scData1.list, features = features)
          #v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg)
          #v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg)
          #print(v$scData1.combined)
          v$isIntgFastMNNdone <- TRUE
          shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
    observeEvent(input$runPCA_intg_fastmnn, {
      tpmFiles1 <- v$scData1.list
      if (is.null(tpmFiles1)){
        v$scData1.list <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message = "Finding clusters...", value = 0.3, {
          #DefaultAssay(v$scData1.combined) <- "integrated"
          v$scData1.combined <- RunFastMNN(object.list = v$scData1.list, features = v$features)
          #output$pca1.done <- renderText(paste0("PCA done!"))
          v$isPCAdone_intg_fastmnn <- TRUE
          shinyalert("PCA done", "PCA done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }  
      })
    
      output$PCAplot_fastmnn_tpm1 <- renderPlotly({
        if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_fastmnn)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T)
          })
        }
      })
      
      output$PCAplot_fastmnn_tpm2 <- renderPlotly({
        if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_fastmnn)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$PCAplot_fastmnn_tpm3 <- renderPlotly({
        if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_fastmnn)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, group.by = 'celltype', label.size = 3) + NoLegend()
          })
        }
      })
      
      output$PCAplot_fastmnn_h5_1 <- renderPlotly({
        if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_fastmnn)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, label.size = 3)
          })
        }
      })
      
      output$PCAplot_fastmnn_h5_2 <- renderPlotly({
        if(is.null(v$scData1.combined) || is.null(v$isPCAdone_intg_fastmnn)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$vizPlot_intg_fastmnn <- renderPlot({
        if(is.null(v$scData1.combined)){
          return(NULL)
        }else{
          VizDimLoadings(v$scData1.combined, dims = as.numeric(input$select.pc_intg_fastmnn), reduction = "mnn")
        }
      })
      
      output$PCtable_intg_fastmnn <- DT::renderDataTable({
        if(is.null(v$scData1.combined) ){
          return(NULL)
        }else{
          v$pcGenes
        }
      }, server = FALSE, options = list(scrollX = TRUE))
      
      output$Elbow_intg_fastmnn <- renderPlot({
        if(is.null(v$scData1.combined)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData1.combined, reduction = "mnn", ndims = 50)
          })
        }
      }, height = 400, width = 450)
    
    observeEvent(input$doDeg_intg_fastmnn, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          v$scData1.combined$celltype -> Idents(v$scData1.combined)
          ips.markers1 <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg_fastmnn, logfc.threshold = input$logfc_intg_fastmnn, test.use = input$test.use_intg_fastmnn)
          v$ips.markers1 <- ips.markers1
          shinyalert("DEGs estimated", "DEGs estimated, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$deg.gene.select_intg_fastmnn <- renderUI({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        selectInput("deg.gene_intg_fastmnn", label = "Gene to visualise",
                    choices = unique(v$ips.markers1$gene))
      }
    })
    
    output$Deg.plot_intg_fastmnn <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          VlnPlot(v$scData1.combined, input$deg.gene_intg_fastmnn)
        })
      }
    })
    
    output$Deg1.plot_intg_fastmnn <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          FeaturePlot(v$scData1.combined, input$deg.gene_intg_fastmnn)
        })
      }
    })
    
    output$Deg2.plot_intg_fastmnn <- renderPlot({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          RidgePlot(v$scData1.combined, features = input$deg.gene_intg_fastmnn)
        })
      }
    })
    
    output$Deg3.plot_intg_fastmnn <- renderPlotly({
      if(is.null(v$ips.markers1)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          v$ips.markers1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData1.combined, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
        })
      }
    })
    
    output$Deg.table_intg_fastmnn <- DT::renderDataTable(
      v$ips.markers1, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    ##---------------Cell-cell similarity of Data Integration-------------------
    
    observeEvent(input$cell_cell1, {
      tpmFiles1 <- v$scData1.combined
      if (is.null(tpmFiles1)){
        v$scData1.combined <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Generating heatmap...", value=0, {
          if(input$cell1a == "primary.predict"){
            v$scData1.combined$primary.predict -> Idents(v$scData1.combined)
            v$scData1.combined.rna.data.average1 = AverageExpression(v$scData1.combined)
            v$scData1.combined.rna.data.average1 = data.frame(v$scData1.combined.rna.data.average1$RNA)
            print(v$scData1.combined.rna.data.average1)
            v$cor1 <- cor(v$scData1.combined.rna.data.average1, method = input$corr_method1)
            print(v$cor1)
            output$CELL1.done <- renderText(paste0("Celltype similarity done!"))
            v$isCELLdone1 <- TRUE
          }
          else if(input$cell1a == "seurat_clusters"){
            v$scData1.combined$seurat_clusters -> Idents(v$scData1.combined)
            v$scData1.combined.rna.data.average1 = AverageExpression(v$scData1.combined)
            v$scData1.combined.rna.data.average1 = data.frame(v$scData1.combined.rna.data.average1$RNA)
            print(v$scData1.combined.rna.data.average1)
            v$cor1 <- cor(v$scData1.combined.rna.data.average1, method = input$corr_method1)
            rownames(v$cor1) <- substr(rownames(v$cor1),2,nchar(rownames(v$cor1)))
            colnames(v$cor1) <- substr(colnames(v$cor1),2,nchar(colnames(v$cor1)))
            print(v$cor1)
            output$CELL1.done <- renderText(paste0("Celltype similarity done!"))
            v$isCELLdone1 <- TRUE
          }
          if(input$cell1a == "celltype"){
            v$scData1.combined$celltype -> Idents(v$scData1.combined)
            v$scData1.combined.rna.data.average1 = AverageExpression(v$scData1.combined)
            v$scData1.combined.rna.data.average1 = data.frame(v$scData1.combined.rna.data.average1$RNA)
            print(v$scData1.combined.rna.data.average1)
            v$cor1 <- cor(v$scData1.combined.rna.data.average1, method = input$corr_method1)
            print(v$cor1)
            output$CELL1.done <- renderText(paste0("Celltype similarity done!"))
            v$isCELLdone1 <- TRUE
          }
        })
      }
    })
    
    plotCELL1 <- reactive({
      if(is.null(v$scData1.combined.rna.data.average1) || is.null(v$isCELLdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating Celltype similarity plot...", value=0, {
          p <- heatmaply(as.matrix(v$cor1), cexRow = 0.8, cexCol = 0.8, margins = c(10,10), k_col =2, k_row = 2, colors = rev(RColorBrewer::brewer.pal(9, "RdBu")))
        })
      }
    })
    
    output$cell_cell_sim1 <- renderPlotly({
      plotCELL1()
    })
    
    output$download_cell_cell_sim <- downloadHandler(
      filename = "Celltype similarity.png",
      content = function(file) {
        png(file)
        heatmap(as.matrix(v$cor1), col = RColorBrewer::brewer.pal(9, "RdBu"), cexRow = 0.8, cexCol = 0.8, margins = c(10,10))
        dev.off()
      }
    )
    
    output$cor.table1 <- DT::renderDataTable(
      v$cor1, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
    
    output$download_cor.table1 <- downloadHandler(
      filename = function(){"Celltype_similarity.csv"}, 
      content = function(fname){
        withProgress(message="Downloading celltype_similarity...", value=0, {
          write.csv(v$cor1, fname)
        })
      }
    )
    
 
    ##------------------------Single cell multiomics module--------------------------##
    
    output$multiomics_image <- renderImage({
      list(src = "www/multiomics_fig1.png",
           width = 900,
           height = 550)
    }, deleteFile = FALSE)
    
    observeEvent(input$loadexample_cite, {
      withProgress(message="Loading example Data...", value=0.5, {
        cite.data <- Read10X_h5('multiomics/cite/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5')
      })
      if (is.null(cite.data)){
        v$scDatat <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(cite.data)
          #label1 <- "Example loaded"
          #updateActionButton(inputId = "loadexample_cite_seurat", label = label1)
          #shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
          cite.data.rna <- cite.data$`Gene Expression`
          cite.data.ab <- cite.data$`Antibody Capture`
          incProgress(0.5, "Creating Seurat Object")
          sObj3 <- CreateSeuratObject(cite.data.rna, project = input$projName2, min.genes = input$min.genes2, min.cells = input$min.cells2)
          sObj3[["percent.mt"]] <- PercentageFeatureSet(sObj3, pattern = "^MT-")
          v$scDatat <- sObj3
          v$scDatat.rna <- cite.data.rna
          v$scDatat.ab <- cite.data.ab
          print(sObj3@meta.data)
          v$scDatab <- CreateAssayObject(counts = cite.data.ab)
          v$scDatat[["ADT"]] <- v$scDatab
          print(Assays(v$scDatat))
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample_cite_seurat", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    observeEvent(input$loadexample_multiome, {
      withProgress(message="Loading example Data...", value=0.5, {
        multiome.data <- Read10X_h5('multiomics/multiome/pbmc_unsorted_3k_filtered_feature_bc_matrix.h5')
      })
      if (is.null(multiome.data)){
        v$scDatan <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(multiome.data)
          multiome.rna <- multiome.data$`Gene Expression`
          multiome.atac <- multiome.data$Peaks
          incProgress(0.5, "Creating Seurat Object")
          sObj4 <- CreateSeuratObject(multiome.rna, project = input$projName2, min.genes = input$min.genes2, min.cells = input$min.cells2)
          sObj4[["percent.mt"]] <- PercentageFeatureSet(sObj4, pattern = "^MT-")
          print(sObj4)
          print(sObj4@meta.data)
          grange.counts <- StringToGRanges(rownames(multiome.atac), sep = c(":", "-"))
          grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
          multiome.atac <- multiome.atac[as.vector(grange.use), ]
          annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
          genome(annotations) <- "hg19"
          seqlevelsStyle(annotations) <- 'UCSC'
          frag.file <- "multiomics/multiome/pbmc_unsorted_3k_atac_fragments.tsv.gz"
          chrom_assay <- CreateChromatinAssay(
            counts = multiome.atac,
            sep = c(":", "-"),
            genome = 'hg19',
            fragments = frag.file,
            min.cells = 10,
            annotation = annotations)
          print(chrom_assay)
          sObj4[["ATAC"]] <- chrom_assay
          print(sObj4)
          v$scDatan <- sObj4
          v$scDatan.rna <- multiome.rna
          v$scDatan.atac <- multiome.atac
          print(Assays(v$scDatan))
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample_multiome_seurat", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    observeEvent(input$loadButton2, {
      if(input$scAnalysis_type == "CITE-seq"){
        tpmFiles2  <- input$tpmFiles2
        tpmFiles2 <- v$scDatat
        if (is.null(tpmFiles2)){
          v$scDatat <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Loading and Processing Data...", value=0, {
            print(tpmFiles2$datapath)
            print(tpmFiles2$name)
            print(file.exists(paste(tpmFiles2$datapath[1], "/", tpmFiles2$name[1], sep="")))
            exp.data3 <- Read10X_h5(tpmFiles2$datapath)
            exp.data3.rna <- exp.data3$`Gene Expression`
            exp.data3.ab <- exp.data3$`Antibody Capture`
            #additional.ident <- NULL
            incProgress(0.5, "Creating Seurat Object")
            sObj3 <- CreateSeuratObject(exp.data3.rna, project = input$projName2, min.genes = input$min.genes2, min.cells = input$min.cells2)
            #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
            sObj3[["percent.mt"]] <- PercentageFeatureSet(sObj3, pattern = "^MT-")
            v$scDatat <- sObj3
            v$scDatat.rna <- exp.data3.rna
            v$scDatat.ab <- exp.data3.ab
            print(sObj3@meta.data)
            v$scDatab <- CreateAssayObject(counts = exp.data3.ab)
            v$scDatat[["ADT"]] <- v$scDatab
            print(Assays(v$scDatat))
            shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      }
    })
    
    output$countdataDT2a <- renderDataTable({
      if(!is.null(v$scDatat.rna))
      {
        if(ncol(v$scDatat.rna) > 20 )
          return(as.matrix(v$scDatat@assays$RNA@counts[,1:20]))
      }
    }, server = FALSE)
    
    output$countdataDT2b <- renderDataTable({
      if(!is.null(v$scDatat.ab))
      {
        if(ncol(v$scDatat.ab) > 20 )
          return(as.matrix(v$scDatat@assays$ADT@counts[,1:20]))
      }
    }, server = FALSE)
  
    observeEvent(input$filter_seurat2a, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(v$scDatat)
        v$scDatat <- subset(v$scDatat, subset = nFeature_RNA > input$obsa & nFeature_RNA < input$obsa1 & nCount_ADT > input$obsa3 & nCount_ADT < input$obsa4 & percent.mt < input$obsa2)
        print(v$scDatat)
        shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  }
)
    
    observeEvent(input$reset_mult, {
      session$reload()
      print("Reset done")
    })
    
    output$nFeature_RNAPlot2a <- renderPlot({
      if(is.null(v$scDatat)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend()
      }
    })
    
    output$nFeature_RNAPlot2c <- renderPlot({
      if(is.null(v$scDatat)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatat, features = c("nFeature_ADT", "nCount_ADT"), ncol = 2) + NoLegend()
      }
    })
    
    output$nFeature_RNAPlot2b <- renderPlot({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatan, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + NoLegend()
      }
    })
    
    output$nFeature_RNAPlot2d <- renderPlot({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatan, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2) + NoLegend()
      }
    })
    
    #output$FeatureScatterPlot_mult_a <- renderPlotly({
    #  if(is.null(v$scDatat)){
    #    plotly_empty()
    #  }else{
    #    print(FeatureScatter(v$scDatat, "nCount_RNA", "nFeature_RNA"))
    #  }
    #})
    
    #output$FeatureScatterPlot_mult_b <- renderPlotly({
    #  if(is.null(v$scDatat)){
    #    plotly_empty()
    #  }else{
    #    print(FeatureScatter(v$scDatat, "nCount_RNA", "percent.mt"))
    #  }
    #})
    
    #output$FeatureScatterPlot_mult_c <- renderPlotly({
    #  if(is.null(v$scDatan)){
    #    plotly_empty()
    #  }else{
    #    print(FeatureScatter(v$scDatan, "nCount_RNA", "nCount_ATAC"))
    #  }
    #})
    
    #output$FeatureScatterPlot_mult_d <- renderPlotly({
    # if(is.null(v$scDatan)){
    #   plotly_empty()
    # }else{
    #   print(FeatureScatter(v$scDatan, "nCount_RNA", "percent.mt"))
    # }
    #})
    
    observeEvent(input$findVarGenes_mult, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Finding variable genes...", value = 0, {
        DefaultAssay(v$scDatat) <- 'RNA'
        v$scDatat <- NormalizeData(v$scDatat)
        v$scDatat <- FindVariableFeatures(v$scDatat,
                                          mean.function = ExpMean,
                                          dispersion.function = LogVMR,
                                          nfeatures = input$var.genes_mult,
                                          selection.method = input$selection.method)
        #all.genes <- rownames(v$scData1)
        v$scDatat <- ScaleData(v$scDatat)
        incProgress(0.5)
        DefaultAssay(v$scDatat) <- 'ADT'
        VariableFeatures(v$scDatat) <- rownames(v$scDatat[["ADT"]])
        v$scDatat <- NormalizeData(v$scDatat, normalization.method = 'CLR', margin = 2)
        v$scDatat <- FindVariableFeatures(v$scDatat,
                                          mean.function = ExpMean,
                                          dispersion.function = LogVMR,
                                          nfeatures = input$var.genes_mult,
                                          selection.method = input$selection.method)
        #all.genes <- rownames(v$scData1)
        v$scDatat <- ScaleData(v$scDatat)
        DefaultAssay(v$scDatat) <- 'RNA'
        #VarGeneText <- paste0("Number of variable genes: ", length(v$scData1@assays$RNA@var.features))
        #output$nVarGenes <- renderText(VarGeneText)
        varGenePlotInput <- function(){
          if(is.null(v$scDatat)){
            return(NULL)
          }else{
            withProgress(message="Plotting variable genes...", value=0, {
              top10 <- head(VariableFeatures(v$scDatat), 10)
              variable_feature1 <- VariableFeaturePlot(v$scDatat)
              variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
              print (variable_feature1)
              print (variable_feature2)
              shinyalert("Highly variable features identified", "Highly variable features identified, please perform PCA", type = "success", imageWidth = 10, imageHeight = 10)
              #dev.off()
            })
          }
        }
        output$VarGenes_mult <- renderPlot({
          varGenePlotInput()
        }, height = 500, width = 600)
        observeEvent(input$PDFc, {
          if(!is.null(v$scDatat)){
            withProgress(message="Downloading plot PDF files...", value=0, {
              print(getwd())
              pdfDir <- paste0(getwd(), .Platform$file.sep, "Multiomics_results/Generated_reports_", Sys.Date())
              if(!dir.exists(pdfDir)){
                dir.create(pdfDir)
              }
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
              i = 0
              while(file.exists(filename2)){
                filename2 <- paste0(pdfDir, .Platform$file.sep,
                                    "Var_genes_plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
              }
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              plot1 <- VariableFeaturePlot(v$scDatat)
              print(plot1)
              dev.off()
              txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
              txtfile <- sub(".pdf", ".txt", txtfile)
              write(v$scDatat@assays$RNA@var.features, file = txtfile)
              
            })
          }
        })
      })
    }
  })
    
    observeEvent(input$runPCA_mult_seurat, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        DefaultAssay(v$scDatat) <- 'RNA'
        v$scDatat <- RunPCA(v$scDatat, verbose = FALSE)
        print(v$scDatat[["pca"]], dims = 1:5, nfeatures = 5)
        DefaultAssay(v$scDatat) <- 'ADT'
        v$scDatat <- RunPCA(v$scDatat, reduction.name = 'apca')
        DefaultAssay(v$scDatat) <- 'RNA'
        v$isPCAdone1 <- TRUE
        PCA_plot1a <- DimPlot(v$scDatat, reduction = "pca", label = T)
        print(PCA_plot1a)
        
        incProgress(0.4, message = "Getting list of PC genes...")
        pc.table <- list()
        for(i in 1:20){
          pcg <- TopFeatures(v$scDatat)
          pc.table[[i]] <- pcg
        }
        pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
        v$pcGenes <- pc.table
        shinyalert("PCA performed", "PCA performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
          }
        )
      }
    })
    
    output$PCAplot_mult_seurat_h5_1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of multiomics dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "pca", label = T)
        })
      }
    })
    
    output$vizPlot_mult_seurat <- renderPlot({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        VizDimLoadings(v$scDatat, dims = as.numeric(input$select.pc_mult_seurat))
      }
    })
    
    output$PCHeatmap_mult_seurat <- renderPlot({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        DimHeatmap(v$scDatat, dims = as.numeric(input$select.pc_mult_seurat))
      }
    })
    
    output$PCtable_mult_seurat <- DT::renderDataTable({
      if(is.null(v$scDatat) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, server = FALSE, options = list(scrollX = TRUE))
    
    output$Elbow_mult_seurat <- renderPlot({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scDatat, ndims = 50)
        })
      }
    }, height = 400, width = 450)
    
    observeEvent(input$runMOFA_mult, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        mofa <- reticulate::import("mofapy2", convert = FALSE)
        v$mofa <- create_mofa(v$scDatat, assays = c("RNA","ADT"))
        print(v$mofa)
        model_opts <- get_default_model_options(v$mofa)
        model_opts$num_factors <- as.numeric(input$num_factor_mofa)
        v$mofa <- prepare_mofa(v$mofa, model_options = model_opts)
        print(model_opts)
        v$mofa <- run_mofa(v$mofa, use_basilisk = F)
        v$isMOFAdone <- TRUE
        incProgress(0.4, message = "Getting list of PC genes...")
        pc.table <- list()
        shinyalert("MOFA2 performed", "MOFA2 performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
      }
      )}
    })
    
    output$mofaplot_mult_1 <- renderPlot({
      if(is.null(v$isMOFAdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of multiomics dataset...", value=0, {
          plot_data_overview(v$mofa)
        })
      }
    })
    
    output$mofaplot_mult_2 <- renderPlot({
      if(is.null(v$isMOFAdone)){
        return(NULL)
      }else{
        plot_factor_cor(v$mofa)
      }
    })
    
    output$mofaplot_mult_3 <- renderPlot({
      if(is.null(v$isMOFAdone)){
        return(NULL)
      }else{
        plot_variance_explained(v$mofa, max_r2=as.numeric(input$num_factor_mofa))
      }
    })
    
    output$mofaplot_mult_4 <- renderPlot({
      if(is.null(v$isMOFAdone)){
        return(NULL)
      }else{
        plot_factors(v$mofa, factors = c(as.numeric(input$factor1),as.numeric(input$factor2)), color_by = "Factor1")
      }
    })
    
    output$mofaplot_mult_5 <- renderPlot({
      if(is.null(v$isMOFAdone)){
        return(NULL)
      }else{
        plot_data_heatmap(v$mofa, 
                          view = input$selection.view,
                          factor = as.numeric(input$factor1a),  
                          features = as.numeric(input$num.features.mofa),
                          denoise = TRUE,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          show_rownames = TRUE, show_colnames = FALSE,
                          scale = "row")
      }
    })
    
    observeEvent(input$findCluster_mult_seurat, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Finding clusters...", value = 0.3, {
        DefaultAssay(v$scDatat) <- "RNA"
        v$scDatat <- FindNeighbors(v$scDatat, reduction = "pca", dims = 1:input$dim.used_mult_seurat)
        v$scDatat <- FindClusters(v$scDatat, resolution = input$clus.res_mult_seurat)
        DefaultAssay(v$scDatat) <- "ADT"
        #v$scDatat <- FindMultiModalNeighbors(v$scDatat, reduction.list = list("pca", "apca"), dims.list = list(1:input$dim.used_mult_seurat, 1:input$dim.used_mult_seurat), modality.weight.name = "RNA.weight")
        v$scDatat <- FindClusters(v$scDatat, graph.name = "wsnn", algorithm = 3, resolution = input$clus.res_mult_seurat2a, verbose = FALSE)
        #v$scDatat <- RunUMAP(v$scDatat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterMultdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please run CellType Identification", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Cluster2DPlot_mult_seurat_a <- renderPlotly({
      if(is.null(v$isClusterMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scDatat, reduction = "rna.umap", label = T)
        })
      }
    })
    
    output$Cluster2DPlot_mult_seurat_b <- renderPlotly({
      if(is.null(v$isClusterMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scDatat, reduction = "adt.umap", label = T)
        })
      }
    })
    
    output$Cluster2DPlot_mult_seurat_c <- renderPlotly({
      if(is.null(v$isClusterMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scDatat, reduction = "wnn.umap", label = T)
        })
      }
    })
    
    observeEvent(input$findCluster_mofa2, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(v$isMOFAdone) || is.null(v$scDatat)){
        #v$scDatat <- NULL
        plotly_empty()
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Finding clusters...", value = 0.3, {
        
        cluster <- cluster_samples(v$mofa, k=input$num.clus.mofa, factors="all")
        clusters <- cluster$cluster
        v$clusters <- clusters
        samples_metadata(v$mofa) <- as.data.frame(v$clusters) %>% tibble::rownames_to_column("sample") %>% as.data.table
        v$mofa@samples_metadata -> df
        rownames(df) <- df[,1]
        v$scDatat <- MOFA2::add_mofa_factors_to_seurat(mofa_object = v$mofa, seurat_object = v$scDatat, views = "all", factors = "all")
        v$scDatat <- AddMetaData(v$scDatat, metadata = df)
        v$scDatat$v.clusters -> Idents(v$scDatat)
        print(v$scDatat@meta.data)
        v$isClusterMOFAdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Cluster2DPlot_mofa <- renderPlot({
      if(is.null(v$isClusterMOFAdone) || is.null(v$scDatat)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scDatat, reduction = "MOFAUMAP", group.by = "v.clusters", label = F)
        })
      }
    })
    
    observeEvent(input$runUMAP_mult_seurat, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        DefaultAssay(v$scDatat) <- "RNA"
        v$scDatat <- RunUMAP(v$scDatat, reduction = "pca", reduction.name = "rna.umap", dims = 1:input$dim.used_mult_seurat)
        DefaultAssay(v$scDatat) <- "ADT"
        #VariableFeatures(v$scDatat) <- rownames(v$scDatat[["ADT"]])
        #v$scDatat <- NormalizeData(v$scDatat, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% 
        v$scDatat <- RunUMAP(v$scDatat, reduction = "apca", reduction.name = "adt.umap", dims = 1:input$dim.used_mult_seurat)
        v$scDatat <- FindMultiModalNeighbors(v$scDatat, reduction.list = list("pca", "apca"), dims.list = list(1:input$dim.used_mult_seurat, 1:input$dim.used_mult_seurat), modality.weight.name = "RNA.weight")
        #v$scDatat <- FindClusters(v$scDatat, graph.name = "wsnn", algorithm = 3, resolution = input$clus.res_mult_seurat2a, verbose = FALSE)
        v$scDatat <- RunUMAP(v$scDatat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        v$isUMAPdone <- TRUE
        print(v$scDatat)
        UMAP_plot_cite_a <- DimPlot(v$scDatat, reduction = "rna.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        UMAP_plot_cite_b <- DimPlot(v$scDatat, reduction = "adt.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ADT")
        UMAP_plot_cite_c <- DimPlot(v$scDatat, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        print(UMAP_plot_cite_a)
        print(UMAP_plot_cite_b)
        print(UMAP_plot_cite_c)
        shinyalert("UMAP done", "UMAP done, please run Clustering", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$UMAPplot_mult_seurat_1 <- renderPlot({
      if(is.null(v$isUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "rna.umap", label = F)
        })
      }
    })
    
    output$UMAPplot_mult_seurat_2 <- renderPlot({
      if(is.null(v$isUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "adt.umap", label = F)
        })
      }
    })
    
    output$UMAPplot_mult_seurat_3 <- renderPlot({
      if(is.null(v$isUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "wnn.umap", label = F)
        })
      }
    })
    
    observeEvent(input$runUMAP_mofa2, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running UMAP...", value = 0.3, {
        
        v$factors <- 1:get_dimensions(v$mofa)[["K"]]
        v$mofa <- run_umap(v$mofa, factors = v$factors, n_neighbors = input$num.neighbors_mofa2, min_dist = 0.30)
       
        v$isMOFAUMAPdone <- TRUE
        shinyalert("UMAP performed", "Clustering performed, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$UMAPplot_mofa_1 <- renderPlot({
      if(is.null(v$isMOFAUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot...", value=0, {
          plot_dimred(v$mofa, method = "UMAP")
        })
      }
    })
    
    observeEvent(input$runTSNE_mult_seurat, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        DefaultAssay(v$scDatat) <- "RNA"
        v$scDatat <- RunTSNE(v$scDatat, reduction = "pca", reduction.name = "tsne.rna", dims = 1:input$dim.used_mult_seurat)
        DefaultAssay(v$scDatat) <- "ADT"
        v$scDatat <- RunTSNE(v$scDatat, reduction = 'apca', dims = 1:input$dim.used_mult_seurat, reduction.name = "tsne.adt", reduction.key = "adtTSNE_")
        v$scDatat <- RunTSNE(v$scDatat, dims = 1:input$dim.used_mult_seurat, reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_", nn.name = "weighted.nn")
        v$isTSNEdone <- TRUE
        TSNE_plot_cite_a <- DimPlot(v$scDatat, reduction = "tsne.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        TSNE_plot_cite_b <- DimPlot(v$scDatat, reduction = "tsne.adt", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ADT")
        TSNE_plot_cite_c <- DimPlot(v$scDatat, reduction = "wnn.tsne", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        print(TSNE_plot_cite_a)
        print(TSNE_plot_cite_b)
        print(TSNE_plot_cite_c)
        shinyalert("tSNE done", "tSNE done, please run Cell type identification", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$TSNEplot_mult_seurat_1 <- renderPlotly({
      if(is.null(v$isTSNEdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "tsne.rna", label = T)
        })
      }
    })
    
    output$TSNEplot_mult_seurat_2 <- renderPlotly({
      if(is.null(v$isTSNEdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "tsne.adt", label = T)
        })
      }
    })
    
    output$TSNEplot_mult_seurat_3 <- renderPlotly({
      if(is.null(v$isTSNEdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatat, reduction = "wnn.tsne", label = T)
        })
      }
    })
    
    observeEvent(input$runTSNE_mofa2, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running UMAP...", value = 0.3, {
        v$mofa <- run_tsne(v$mofa, factors = "all", groups = "all")
        samples_metadata(v$mofa) <- as.data.frame(v$clusters) %>% tibble::rownames_to_column("sample") %>% as.data.table
        v$mofa@samples_metadata -> df
        rownames(df) <- df[,1]
        v$t1 <- MOFA2::add_mofa_factors_to_seurat(mofa_object = v$mofa, seurat_object = v$scDatat, views = "all", factors = "all")
        v$t1 <- AddMetaData(v$t1, metadata = df)
        print(v$t1@meta.data)
        v$isMOFAdone <- TRUE
        shinyalert("UMAP performed", "Clustering performed, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$TSNEplot_mofa_1 <- renderPlotly({
      if(is.null(v$isMOFAdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$t1, group.by = "v.clusters", label = T, reduction = "MOFATSNE")
        })
      }
    })
    
    observeEvent(input$doCELLiD_mult_cite_seurat, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scDatat.rna.data.average = AverageExpression(v$scDatat)
        v$scDatat.rna.data.average = round(v$scDatat.rna.data.average$RNA, 2)
        if(input$cellatlas_mult_cite_seurat == "all"){
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, ref)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, adipose1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, adrenal_gland1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, blood1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, bone_marrow1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, brain1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, breast1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, breast_milk1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, eye1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, gut1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, heart1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, kidney1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, liver1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, lung1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, pancreas1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, PDAC1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, skin1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, testis1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, thymus1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_seurat == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, tonsil1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDMultSeuratdone <- TRUE
        shinyalert("Cell type identification done", "Cell type identification done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_cellid_mult_cite_seurat <- renderPlot({
      if(is.null(v$scDatat) || is.null(v$isCELLiDMultSeuratdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scDatat, reduction = input$assay_mult_cite_seurat, group.by = "primary.predict", label = F, label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_mult_cite_seurat1 <- renderPlot({
      if(is.null(v$scDatat) || is.null(v$isCELLiDMultSeuratdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scDatat, reduction = input$assay_mult_cite_seurat, group.by = "secondary.predict", label = F, label.size = 3)
        })
      }
    })
    
    output$ct_cite_seurat.table <- DT::renderDataTable(
      v$res, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_cellid_cite_seurat_prediction <- downloadHandler(
      filename = function(){"CELLiD predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$res, fname)
        })
      }
    )
    
    observeEvent(input$doCELLiD_mult_cite_mofa, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scDatat.rna.data.average = AverageExpression(v$scDatat)
        v$scDatat.rna.data.average = round(v$scDatat.rna.data.average$RNA, 2)
        if(input$cellatlas_mult_cite_mofa == "all"){
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, ref)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, adipose1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, adrenal_gland1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, blood1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, bone_marrow1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, brain1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, breast1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, breast_milk1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, eye1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, gut1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, heart1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, kidney1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, liver1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, lung1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, pancreas1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, PDAC1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, skin1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, testis1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, thymus1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res = FastIntegration::CELLiD(v$scDatat.rna.data.average, tonsil1)
          print(v$res)
          v$scDatat$primary.predict = v$res[as.numeric(v$scDatat$v.clusters),1]
          v$scDatat$secondary.predict = v$res[as.numeric(v$scDatat$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatat@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDMultmofadone <- TRUE
        shinyalert("Cell type identification done", "Cell type identification done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_cellid_mult_cite_mofa <- renderPlot({
      if(is.null(v$scDatat) || is.null(v$isCELLiDMultmofadone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scDatat, group.by = "primary.predict", label = F, label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_mult_cite_mofa1 <- renderPlot({
      if(is.null(v$scDatat) || is.null(v$isCELLiDMultmofadone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scDatat, group.by = "secondary.predict", label = F, label.size = 3)
        })
      }
    })
    
    output$ct_cite_mofa.table <- DT::renderDataTable(
      v$res, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_cellid_cite_mofa_prediction <- downloadHandler(
      filename = function(){"CELLiD predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$res, fname)
        })
      }
    )
    
    observeEvent(input$doCelltypist_mult_cite_seurat, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running Celltypist...", value = 0.3, {
        sc <- reticulate::import("scanpy", convert = FALSE)
        ct <- reticulate::import("celltypist", convert = FALSE)
        sceasy::convertFormat(v$scDatat, from = "seurat", to = "anndata", outFile = 'ct_scrna.h5ad')
        v$adata = sc$read_h5ad('ct_scrna.h5ad')
        v$res = ct$annotate(filename = 'ct_scrna.h5ad', model = input$celltypistatlas9, majority_voting=T)
        print(v$res)
        v$adata = v$res$to_adata()
        print("fff")
        v$adata$obs$to_csv('celltypist_predict.csv')
        v$meta9 <- read.csv('celltypist_predict.csv', header = T, row.names = 1)
        v$scDatat <- AddMetaData(v$scDatat, metadata = v$meta9)
        v$scDatat$primary.predict <- v$scDatat$majority_voting
        v$scDatat$secondary.predict <- v$scDatat$predicted_labels
        print(v$scDatat@meta.data)
        v$isCelltypistdone9 <- TRUE
        shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform cell-cell similarity", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_celltypist_mult_cite_seurat <- renderPlotly({
      if(is.null(v$scDatat) || is.null(v$isCelltypistdone9)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from Celltypist...", value=0, {
          DimPlot(v$scDatat, reduction = input$assay_mult1_cite_seurat, group.by = "majority_voting", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$Umap_celltypist_mult_cite_seurat1 <- renderPlotly({
      if(is.null(v$scDatat) || is.null(v$isCelltypistdone9)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from Celltypist...", value=0, {
          DimPlot(v$scDatat, reduction = input$assay_mult1_cite_seurat, group.by = "predicted_labels", label = T,  label.size = 3) + NoLegend()
        })
      }
    })
    
    output$ct_celltypist_mult_cite_seurat.table <- DT::renderDataTable(
      v$meta9, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_celltypist_cite_seurat_prediction <- downloadHandler(
      filename = function(){"Celltypist predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$meta9, fname)
        })
      }
    )
    
    #output$vis.gene.select1 <- renderUI({
    #  if(is.null(v$scDatat)){
    #    return(NULL)
    #  }else{
    #    selectInput("vis.gene1", label = "Gene to visualise",
    #                choices = rownames(v$scDatat[["ADT"]]))
    #  }
    #})
    
    #output$vis1.plot <- renderPlotly({
    #  if(is.null(v$scDatat)){
    #    return(NULL)
    #  }else{
    #    withProgress(message="Generating DEG Plot...", value=0, {
    #      FeaturePlot(v$scDatat, input$vis.gene1, cols = c("lightgrey", "darkgreen"), order = T, reduction = "adt.umap")
    #    })
    #  }
    #})
    
    #observeEvent(input$Vis3, {
    #  tpmFiles2  <- input$tpmFiles2
    #  tpmFiles2 <- v$scDatat
    #  if (is.null(tpmFiles2)){
    #    v$scDatat <- NULL
    #    shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
    #  }else{
    #    withProgress(message = "Visualizing...", value = 0,{
    #      incProgress(0.5, message = "Visualizing...")
    #      DefaultAssay(v$scDatat) <- "ADT"
    #    })
    #  }
    #})
    
    observeEvent(input$Vis3a, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Visualizing...", value = 0,{
        incProgress(0.5, message = "Visualizing...")
        DefaultAssay(v$scDatat) <- "ADT"
        })
      }
    })
    
    output$vis.gene.select_a <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        selectInput("vis.gene_a", label = "Gene to visualise",
                    choices = rownames(v$scDatat[["RNA"]]))
      }
    })
    
    output$vis.plot_a <- renderPlotly({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        withProgress(message="Generating Feature Plot...", value=0, {
          FeaturePlot(v$scDatat, input$vis.gene_a, cols = c("lightgrey", "darkgreen"), order = T, reduction = "MOFAUMAP")
        })
      }
    })
    
    output$vis.gene.select1_a <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        selectInput("vis.gene1a", label = "Gene to visualise",
                    choices = rownames(v$scDatat[["ADT"]]))
      }
    })
    
    output$vis1.plot_a <- renderPlotly({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          FeaturePlot(v$scDatat, input$vis.gene1a, cols = c("lightgrey", "darkgreen"), order = T, reduction = "MOFAUMAP")
        })
      }
    })
    
    observeEvent(input$cell_cell_mult, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Generating heatmap...", value=0, {
          if(input$cell_m1 == "primary.predict"){
            v$scDatat$primary.predict -> Idents(v$scDatat)
            v$scDatat.rna.data.average1 = AverageExpression(v$scDatat)
            v$scDatat.rna.data.average1 = data.frame(v$scDatat.rna.data.average1$RNA)
            print(v$scDatat.rna.data.average1)
            v$cor_mult <- cor(v$scDatat.rna.data.average1, method = input$corr_method2)
            print(v$cor_mult)
            output$CELL_mult.done <- renderText(paste0("Celltype similarity done!"))
            v$isCELLMultdone <- TRUE
          }
          if(input$cell_m1 == "seurat_clusters"){
            v$scDatat$seurat_clusters -> Idents(v$scDatat)
            v$scDatat.rna.data.average1 = AverageExpression(v$scDatat)
            v$scDatat.rna.data.average1 = data.frame(v$scDatat.rna.data.average1$RNA)
            print(v$scDatat.rna.data.average1)
            v$cor_mult <- cor(v$scDatat.rna.data.average1, method = input$corr_method2)
            rownames(v$cor_mult) <- substr(rownames(v$cor_mult),2,nchar(rownames(v$cor_mult)))
            colnames(v$cor_mult) <- substr(colnames(v$cor_mult),2,nchar(colnames(v$cor_mult)))
            print(v$cor_mult)
            output$CELL_mult.done <- renderText(paste0("Celltype similarity done!"))
            v$isCELLMultdone <- TRUE
          }
        })
      }
    })
    
    plotCELLMult <- reactive({
      if(is.null(v$scDatat.rna.data.average1) || is.null(v$isCELLMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating Celltype similarity plot...", value=0, {
          p <- heatmaply(as.matrix(v$cor_mult), cexRow = 0.8, cexCol = 0.8, margins = c(10,10), k_col =2, k_row = 2, colors = rev(RColorBrewer::brewer.pal(9, "RdBu")))
        })
      }
    })
    
    output$cell_cell_mult_sim <- renderPlotly({
      plotCELLMult()
    })
    
    output$download_cell_cell_mult_sim <- downloadHandler(
      filename = "Celltype similarity.png",
      content = function(file) {
        png(file)
        heatmap(as.matrix(v$cor_mult), col = RColorBrewer::brewer.pal(9, "RdBu"), cexRow = 0.8, cexCol = 0.8, margins = c(10,10))
        dev.off()
      }
    )
    
    output$cor_mult.table <- DT::renderDataTable(
      v$cor_mult, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
    
    output$download_cor_mult.table <- downloadHandler(
      filename = function(){"Celltype_similarity.csv"}, 
      content = function(fname){
        withProgress(message="Downloading celltype_similarity...", value=0, {
          write.csv(v$cor_mult, fname)
        })
      }
    )
    
    observeEvent(input$doDeg_mult, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          if(input$deg_mult1 == "seurat_clusters"){
            v$scDatat$seurat_clusters -> Idents(v$scDatat)
              ips.markers_mult <- FindAllMarkers(v$scDatat, only.pos = FALSE, min.pct = input$min_pct_mult, logfc.threshold = input$logfc_mult, assay = "RNA", test.use = input$test.use_mult)
              v$ips.markers_mult <- ips.markers_mult
            shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
          }
          if(input$deg_mult1 == "primary.predict"){
            v$scDatat$primary.predict -> Idents(v$scDatat)
              ips.markers_mult <- FindAllMarkers(v$scDatat, only.pos = FALSE, min.pct = input$min_pct_mult, logfc.threshold = input$logfc_mult, assay = "RNA", test.use = input$test.use_mult)
              v$ips.markers_mult <- ips.markers_mult
            shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
          }
        })
      }
    })
    
    observeEvent(input$Vis_mult, {
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Visualizing...", value=0, {
          v$isVisMultdone <- TRUE
        })
      }
    })
    
    output$vis.gene.select <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        selectInput("vis.gene", label = "Gene to visualise",
                    choices = rownames(v$scDatat[[input$assay_mult_cite_seurat2]]))
      }
    })
    
    plotFeature_mult <- reactive({
      if(is.null(v$scDatat) || is.null(v$isVisMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Visualizing...", value=0, {
          if(input$deg_mult2 == "seurat_clusters"){
            v$scDatat$seurat_clusters -> Idents(v$scDatat)
            FeaturePlot(v$scDatat, input$vis.gene, cols = c("lightgrey", "darkgreen"), order = T, reduction = input$assay_mult_cite_seurat1)
          }
          else if(input$deg_mult2 == "primary.predict"){
            v$scDatat$primary.predict -> Idents(v$scDatat)
            FeaturePlot(v$scDatat, input$vis.gene, cols = c("lightgrey", "darkgreen"), order = T, reduction = input$assay_mult_cite_seurat1)
          }
        })
      }
    })
    
    output$vis.plot <- renderPlotly({
      plotFeature_mult()
    })
    
    output$download_feature_mult <- downloadHandler(
      filename = function(){"Feature plot (Multiomics module).png"}, 
      content = function(fname){
        ggsave(fname,plotFeature_mult(), height = 7, width = 7)
      }
    )
    
    plotViolin_mult <- reactive({
      if(is.null(v$scDatat) || is.null(v$isVisMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Visualizing...", value=0, {
          if(input$deg_mult2 == "seurat_clusters"){
            v$scDatat$seurat_clusters -> Idents(v$scDatat)
            print(v$scDatat)
            VlnPlot(v$scDatat, input$vis.gene, assay = input$assay_mult_cite_seurat2)
          }
          else if(input$deg_mult2 == "primary.predict"){
            v$scDatat$primary.predict -> Idents(v$scDatat)
            VlnPlot(v$scDatat, input$vis.gene, assay = input$assay_mult_cite_seurat2)
          }
        })
      }
    })
    
    output$Deg_mult.plot <- renderPlotly({
      plotViolin_mult()
    })
    
    output$download_violn_mult <- downloadHandler(
      filename = function(){"Violin plot (Multiomics module).png"}, 
      content = function(fname){
        ggsave(fname,plotViolin_mult(), height = 7, width = 7)
      }
    )
    
    plotRidge_mult <- reactive({
      if(is.null(v$scDatat) || is.null(v$isVisMultdone)){
        plotly_empty()
      }else{
        withProgress(message="Visualizing...", value=0, {
          if(input$deg_mult2 == "seurat_clusters"){
            v$scDatat$seurat_clusters -> Idents(v$scDatat)
            RidgePlot(v$scDatat, features = input$vis.gene, assay = input$assay_mult_cite_seurat2)
          }
          else if(input$deg_mult2 == "primary.predict"){
            v$scDatat$primary.predict -> Idents(v$scDatat)
            RidgePlot(v$scDatat, features = input$vis.gene, assay = input$assay_mult_cite_seurat2)
          }
        })
      }
    })
    
    output$Deg_mult2.plot <- renderPlot({
      plotRidge_mult()
    })
    
    output$download_ridge_mult <- downloadHandler(
      filename = function(){"Ridge plot (Multiomics module).png"}, 
      content = function(fname){
        ggsave(fname,plotRidge_mult(), height = 7, width = 7)
      }
    )
    
    output$Deg3_mult.plot <- renderPlot({
      if(is.null(v$ips.markers_mult)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          v$ips.markers_mult %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scDatat, features = top10$gene, assay = "RNA", size = 5, angle = 45) + theme(axis.text.y = element_text(size = 4)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
        })
      }
    })
    
    output$Deg_mult.table <- DT::renderDataTable(
      v$ips.markers_mult, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
    
    output$gene1_mult.select <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        if(input$deg3_mult == "seurat_clusters"){
          selectInput("gene1_mult", label = "Celltype1",
                      choices = as.vector(v$scDatat$seurat_clusters))
        }
        else if(input$deg3_mult == "primary.predict"){
          selectInput("gene1_mult", label = "Celltype1",
                      choices = as.vector(v$scDatat$primary.predict))
        }
      }
    })
    
    output$gene2_mult.select <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        if(input$deg3_mult == "seurat_clusters"){
          selectInput("gene2_mult", label = "Celltype2",
                      choices = as.vector(v$scDatat$seurat_clusters))
        }
        else if(input$deg3_mult == "primary.predict"){
          selectInput("gene2_mult", label = "Celltype2",
                      choices = as.vector(v$scDatat$primary.predict))
        }
      }
    })
    
    observeEvent(input$doVolcano_mult, {
      tpmFiles2  <- input$tpmFiles2
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          v$scDatat_subset <- subset(v$scDatat, subset = primary.predict == input$gene1_mult | primary.predict == input$gene2_mult)
         
            ips.markers_mult_a <- FindAllMarkers(v$scDatat_subset, only.pos = F, min.pct = input$min_pct_a_mult, logfc.threshold = input$logfc_a_mult, assay = "RNA", test.use = input$test.use_a_mult)
            ips.markers_mult_b <- FindMarkers(v$scDatat_subset, ident.1 = input$gene1_mult, ident.2 = input$gene2_mult, only.pos = F, min.pct = input$min_pct_a_mult, logfc.threshold = input$logfc_a_mult, assay = "RNA", test.use = input$test.use_a_mult)
            v$ips.markers_mult_a <- ips.markers_mult_a
            v$ips.markers_mult_b <- ips.markers_mult_b
          shinyalert("Pairwise DEGs done", "Pairwise DEGs done, please run GSEA Analysis", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$volcano_mult.plot <- renderPlot({
      if(is.null(v$ips.markers_mult_b)){
        return(NULL)
      }else{
        withProgress(message="Generating Volcano Plot...", value=0, {
          EnhancedVolcano(toptable = v$ips.markers_mult_b, lab = row.names(v$ips.markers_mult_b), x ="avg_log2FC", y ="p_val_adj", pointSize = 1, labSize = 5, legendLabSize = 12, axisLabSize = 12)
        })
      }
    })
    
    output$dega_mult.plot <- renderPlot({
      if(is.null(v$ips.markers_mult_a)){
        return(NULL)
      }else{
        withProgress(message="Generating Heatmap...", value=0, {
          v$ips.markers_mult_a %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scDatat_subset, features = top10$gene, assay = "RNA", size = 5, angle = 45) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
        })
      }
    })
    
    output$gsea.ct_mult1.select <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        selectInput("gsea_mult.ct1", label = "Celltype1",
                    choices = as.vector(v$scDatat$primary.predict))
      }
    })
    
    output$gsea.ct_mult2.select <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        selectInput("gsea_mult.ct2", label = "Celltype2",
                    choices = as.vector(v$scDatat$primary.predict))
      }
    })
    
    observeEvent(input$gsea_mult, {
      tpmFiles2 <- v$scDatat
      if (is.null(tpmFiles2)){
        v$scDatat <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Generating gene set enrichment analysis...", value=0, {
          if(input$species_gsea_mult == "Homo sapiens" & input$category_gsea_mult == "H"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_hs_go_mult <- msigdbr(species = "Homo sapiens", category = "H")
            print(v$msigdbr_hs_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_hs_go_mult$gene_symbol, f = v$msigdbr_hs_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Mus musculus" & input$category_gsea_mult == "H"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_mm_go_mult <- msigdbr(species = "Mus musculus", category = "H")
            print(v$msigdbr_mm_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_mm_go_mult$gene_symbol, f = v$msigdbr_mm_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Homo sapiens" & input$category_gsea_mult == "C2"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_hs_go_mult <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
            print(v$msigdbr_hs_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_hs_go_mult$gene_symbol, f = v$msigdbr_hs_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Mus musculus" & input$category_gsea_mult == "C2"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_mm_go_mult <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
            print(v$msigdbr_mm_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_mm_go_mult$gene_symbol, f = v$msigdbr_mm_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Homo sapiens" & input$category_gsea_mult == "C5"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_hs_go_mult <- msigdbr(species = "Homo sapiens", category = "C5")
            print(v$msigdbr_hs_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_hs_go_mult$gene_symbol, f = v$msigdbr_hs_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Mus musculus" & input$category_gsea_mult == "C5"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_mm_go_mult <- msigdbr(species = "Mus musculus", category = "C5")
            print(v$msigdbr_mm_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_mm_go_mult$gene_symbol, f = v$msigdbr_mm_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Homo sapiens" & input$category_gsea_mult == "C7"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_hs_go_mult <- msigdbr(species = "Homo sapiens", category = "C7")
            print(v$msigdbr_hs_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_hs_go_mult$gene_symbol, f = v$msigdbr_hs_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Mus musculus" & input$category_gsea_mult == "C7"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_mm_go_mult <- msigdbr(species = "Mus musculus", category = "C7")
            print(v$msigdbr_mm_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_mm_go_mult$gene_symbol, f = v$msigdbr_mm_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Homo sapiens" & input$category_gsea_mult == "C8"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_hs_go_mult <- msigdbr(species = "Homo sapiens", category = "C8")
            print(v$msigdbr_hs_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_hs_go_mult$gene_symbol, f = v$msigdbr_hs_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          if(input$species_gsea_mult == "Mus musculus" & input$category_gsea_mult == "C8"){
            DefaultAssay(v$scDatat) <- "RNA"
            v$msigdbr_mm_go_mult <- msigdbr(species = "Mus musculus", category = "C8")
            print(v$msigdbr_mm_go_mult)
            v$pathways_mult <- split(x = v$msigdbr_mm_go_mult$gene_symbol, f = v$msigdbr_mm_go_mult$gs_name)
            print(v$pathways_mult)
            v$markers_mult <- FindMarkers(v$scDatat, ident.1 = input$gsea_mult.ct1, ident.2 = input$gsea_mult.ct2, min.pct = input$min_pct_mult1, logfc.threshold = input$logfc_mult1, test.use = input$test.use_mult1)
            v$markers_mult  <- v$markers_mult %>% arrange(desc(avg_log2FC))
            print(v$markers_mult)
            v$markers_mult.log2FC <- v$markers_mult$avg_log2FC
            names(v$markers_mult.log2FC) <- row.names(v$markers_mult)
            v$markers_mult.log2FC <- sort(na.omit(v$markers_mult.log2FC), decreasing = TRUE)
            print(v$markers_mult.log2FC)
            v$fgseaRes_mult <- fgsea(pathways = v$pathways_mult, stats = v$markers_mult.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgseaRes_mult)
            v$topPathwaysUp_mult <- v$fgseaRes_mult[ES > 0][head(order(pval), n=10), pathway]
            v$topPathwaysDown_mult <- v$fgseaRes_mult[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_mult <- c(v$topPathwaysUp_mult, rev(v$topPathwaysDown_mult))
          }
          output$gsea_mult.done <- renderText(paste0("Gene set enrichment done!"))
          v$isGSEAmultdone <- TRUE
        })
      }
    })
    
    output$gsea_mult.select <- renderUI({
      if(is.null(v$pathways_mult)){
        return(NULL)
      }else{
        selectInput("gsea_mult.pathway", label = "Gene set to visualise",
                    choices = names(v$pathways_mult))
      }
    })
    
    output$gsea_mult_plot <- renderPlot({
      if(is.null(v$pathways_mult) || is.null(v$isGSEAmultdone)){
        return(NULL)
      }else{
        withProgress(message="Generating GSEA plot...", value=0, {
          plotEnrichment(v$pathways_mult[[input$gsea_mult.pathway]], v$markers_mult.log2FC) + labs(title=input$gsea_mult.pathway)
        })
      }
    })
    
    output$gsea_mult_plot1 <- renderPlot({
      if(is.null(v$pathways_mult) || is.null(v$isGSEAmultdone)){
        return(NULL)
      }else{
        withProgress(message="Generating GSEA plot...", value=0, {
          plotGseaTable(v$pathways_mult[v$topPathways_mult], v$markers_mult.log2FC, v$fgseaRes_mult, gseaParam=0.5)
        })
      }
    })
    
    output$gsea_mult.table <- DT::renderDataTable(
      v$fgseaRes_mult, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
    
    output$download_gsea_mult.table <- downloadHandler(
      filename = function(){"GSEA Results.csv"}, 
      content = function(fname){
        withProgress(message="Downloading GSEA Results...", value=0, {
          fwrite(v$fgseaRes_mult, fname)
        })
      }
    )
    
    observe({if(input$scAnalysis_type == "Multiome"){
      shinyFileChoose(
        input,
        'dir_multi_atac',
        roots = c(home = '.'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw", "gz", "tbi")
      )
      
      dir_multi_atac <- reactive(input$dir_multi_atac)
      output$dir_multi_atac <- renderPrint({  # use renderText instead of renderPrint
        as.character(parseFilePaths(c(home = '.'), dir_multi_atac())$datapath)
      })
    } 
  })  
    
     
    
    output$countdataDT2c <- renderDataTable({
      if(!is.null(v$scDatan.rna))
      {
        if(ncol(v$scDatan.rna) > 20 )
          return(as.matrix(v$scDatan@assays$RNA@counts[,1:20]))
      }
    }, server = FALSE)
    
    output$countdataDT2d <- renderDataTable({
      if(!is.null(v$scDatan.atac))
      {
        if(ncol(v$scDatan.atac) > 20 )
          return(as.matrix(v$scDatan@assays$ATAC@counts[,1:20]))
      }
    }, server = FALSE)
    
    observeEvent(input$filter_seurat2b, {
      tpmFiles3 <- input$tpmFiles3
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message="Loading and Processing Data...", value=0, {
        print(v$scDatan)
        v$scDatan <- subset(v$scDatan, subset = nFeature_RNA > input$obsa & nFeature_RNA < input$obsa1 & nCount_ATAC > input$obsa3 & nCount_ATAC < input$obsa4 & percent.mt < input$obsa2)
        print(v$scDatan)
        shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
      }
    })
    
    observeEvent(input$doSCTransform_multi, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running scTransform...", value = 0,{
        incProgress(0.5, message = "Running scTransform...")
        DefaultAssay(v$scDatan) <- "RNA"
        v$scDatan <- SCTransform(v$scDatan, variable.features.n = input$var.genes_mult, vars.to.regress = "percent.mt", verbose = FALSE, conserve.memory = T)
        incProgress(0.5)
        #VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
        #output$nVarGenes <- renderText(VarGeneText)
        varGenePlotInput <- function(){
          if(is.null(v$scDatan)){
            return(NULL)
          }else{
            withProgress(message="Plotting variable genes...", value=0, {
              top10 <- head(VariableFeatures(v$scDatan), 10)
              variable_feature1 <- VariableFeaturePlot(v$scDatan)
              variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
              print (variable_feature1)
              print (variable_feature2)
              shinyalert("scTransform done", "scTransform done, please perform PCA", type = "success", imageWidth = 10, imageHeight = 10)
            })
          }
        }
        
        output$VarGenes_mult1 <- renderPlot({
          varGenePlotInput()
        }, height = 500, width = 600)
        observeEvent(input$PDFc, {
          if(!is.null(v$scDatan)){
            withProgress(message="Downloading plot PDF files...", value=0, {
              print(getwd())
              pdfDir <- paste0(getwd(), .Platform$file.sep, "Multiomics_results/Generated_reports_", Sys.Date())
              if(!dir.exists(pdfDir)){
                dir.create(pdfDir)
              }
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
              i = 0
              while(file.exists(filename2)){
                filename2 <- paste0(pdfDir, .Platform$file.sep,
                                    "Var_genes_plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
              }
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              plot1 <- VariableFeaturePlot(v$scDatan)
              print(plot1)
              dev.off()
              txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
              txtfile <- sub(".pdf", ".txt", txtfile)
              write(v$scDatan@assays$RNA@var.features, file = txtfile)
              })
            }
          })
        })
      }
    })
    
    observeEvent(input$runPCA_mult_seurat1, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        v$scDatan <- RunPCA(v$scDatan, verbose = FALSE)
        print(v$scDatan[["pca"]], dims = 1:5, nfeatures = 5)
        v$isPCAdone1 <- TRUE
        PCA_plot1a <- DimPlot(v$scDatan, reduction = "pca", label = T)
        print(PCA_plot1a)
        incProgress(0.4, message = "Getting list of PC genes...")
        pc.table <- list()
        for(i in 1:20){
          pcg <- TopFeatures(v$scDatan)
          pc.table[[i]] <- pcg
        }
        pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
        v$pcGenes <- pc.table
        shinyalert("PCA performed", "PCA performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
        }
      )}
    })
    
    output$PCAplot_mult_seurat_h5_2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of multiomics dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "pca", label = T)
        })
      }
    })
    
    output$vizPlot_mult_seurat1 <- renderPlot({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        VizDimLoadings(v$scDatan, dims = as.numeric(input$select.pc_mult_seurat))
      }
    })
    
    output$PCHeatmap_mult_seurat1 <- renderPlot({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        DimHeatmap(v$scDatan, dims = as.numeric(input$select.pc_mult_seurat1))
      }
    })
    
    output$PCtable_mult_seurat1 <- DT::renderDataTable({
      if(is.null(v$scDatan) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, server = FALSE, options = list(scrollX = TRUE))
    
    output$Elbow_mult_seurat1 <- renderPlot({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scDatan)
        })
      }
    }, height = 400, width = 450)
    
    observeEvent(input$findCluster_mult_seurat1, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Finding clusters...", value = 0.3, {
        v$scDatan <- FindNeighbors(v$scDatan, reduction = "pca", dims = 1:input$dim.used_mult_seurat1)
        v$scDatan <- FindClusters(v$scDatan, resolution = input$clus.res_mult_seurat1)
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Cluster2DPlot_mult_seurat1 <- renderPlotly({
      if(is.null(v$isClusterdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scDatan, reduction = "pca", label = T)
        })
      }
    })
    
    observeEvent(input$runUMAP_mult_seurat1, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        v$scDatan <- RunUMAP(v$scDatan, reduction.name = "umap.rna", dims = 1:input$dim.used_mult_seurat1, reduction.key = 'rnaUMAP_')
        print(v$scDatan)
        DefaultAssay(v$scDatan) <- "ATAC"
        v$scDatan <- RunTFIDF(v$scDatan)
        v$scDatan <- FindTopFeatures(v$scDatan)
        v$scDatan <- RunSVD(v$scDatan)
        v$scDatan <- RunUMAP(v$scDatan, reduction = 'lsi', dims = 2:input$dim.used_mult_seurat1, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
        v$scDatan <- FindMultiModalNeighbors(v$scDatan, reduction.list = list("pca", "lsi"), dims.list = list(1:input$dim.used_mult_seurat1, 2:input$dim.used_mult_seurat1))
        v$scDatan <- FindClusters(v$scDatan, resolution = input$clus.res_mult_seurat2b, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
        v$scDatan <- RunUMAP(v$scDatan, dims = 1:input$dim.used_mult_seurat1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        v$isUMAPdone <- TRUE
        print(v$scDatan)
        UMAP_plot1a <- DimPlot(v$scDatan, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        UMAP_plot1b <- DimPlot(v$scDatan, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
        UMAP_plot1c <- DimPlot(v$scDatan, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        print(UMAP_plot1a)
        print(UMAP_plot1b)
        print(UMAP_plot1c)
        shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$UMAPplot1_a <- renderPlotly({
      if(is.null(v$isUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        })
      }
    })
    
    output$UMAPplot1_b <- renderPlotly({
      if(is.null(v$isUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
        })
      }
    })
    
    output$UMAPplot1_c <- renderPlotly({
      if(is.null(v$isUMAPdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        })
      }
    })
    
    observeEvent(input$runTSNE_mult_seurat1, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        v$scDatan <- RunTSNE(v$scDatan, dims = 1:input$dim.used_mult_seurat1, reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')
        DefaultAssay(v$scDatan) <- "ATAC"
        v$scDatan <- RunTSNE(v$scDatan, reduction = 'lsi', dims = 2:input$dim.used_mult_seurat1, reduction.name = "tsne.atac", reduction.key = "atacTSNE_")
        v$scDatan <- RunTSNE(v$scDatan, dims = 1:input$dim.used_mult_seurat1, nn.name = "weighted.nn", reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_")
        v$isTSNEdone <- TRUE
        UMAP_plot1a <- DimPlot(v$scDatan, reduction = "tsne.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        UMAP_plot1b <- DimPlot(v$scDatan, reduction = "tsne.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
        UMAP_plot1c <- DimPlot(v$scDatan, reduction = "wnn.tsne", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        print(UMAP_plot1a)
        print(UMAP_plot1b)
        print(UMAP_plot1c)
        shinyalert("tSNE done", "tSNE done, please run Cell type identification", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$TSNEplot1_a <- renderPlotly({
      if(is.null(v$isTSNEdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating TSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "tsne.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        })
      }
    })
    
    output$TSNEplot1_b <- renderPlotly({
      if(is.null(v$isTSNEdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating TSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "tsne.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
        })
      }
    })
    
    output$TSNEplot1_c <- renderPlotly({
      if(is.null(v$isTSNEdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating TSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scDatan, reduction = "wnn.tsne", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        })
      }
    })
    
    observeEvent(input$doCELLiD_multiome_seurat, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scDatan.rna.data.average = AverageExpression(v$scDatan)
        v$scDatan.rna.data.average = round(v$scDatan.rna.data.average$RNA, 2)
        if(input$cellatlas_mult_multiome_seurat == "all"){
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, ref)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, adipose1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, adrenal_gland1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, blood1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, bone_marrow1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, brain1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, breast1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, breast_milk1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, eye1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, gut1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, heart1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, kidney1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, liver1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, lung1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, pancreas1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, PDAC1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, skin1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, testis1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, thymus1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_multiome_seurat == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res = FastIntegration::CELLiD(v$scDatan.rna.data.average, tonsil1)
          print(v$res)
          v$scDatan$primary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),1]
          v$scDatan$secondary.predict = v$res[as.numeric(v$scDatan$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scDatan@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDMultiomedone <- TRUE
        shinyalert("Cell type identification done", "Cell type identification done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$Umap_cellid_multiome_seurat <- renderPlotly({
      if(is.null(v$scDatan) || is.null(v$isCELLiDMultiomedone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scDatan, reduction = input$assay_mult_multiome_seurat, group.by = "primary.predict", label = T, label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_multiome_seurat1 <- renderPlotly({
      if(is.null(v$scDatan) || is.null(v$isCELLiDMultiomedone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scDatan, reduction = input$assay_mult_multiome_seurat, group.by = "secondary.predict", label = T, label.size = 3)
        })
      }
    })
    
    output$ct_multiome_seurat.table <- DT::renderDataTable(
      v$res, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_cellid_multiome_seurat_prediction <- downloadHandler(
      filename = function(){"CELLiD predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$res, fname)
        })
      }
    )
    
    observeEvent(input$Vis3_multiome, {
      tpmFiles3  <- input$tpmFiles3
      tpmFiles3 <- v$scDatan
      if (is.null(tpmFiles3)){
        v$scDatan <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
      withProgress(message = "Visualizing...", value = 0,{
        incProgress(0.5, message = "Visualizing...")
        DefaultAssay(v$scDatan) <- "RNA"
        })
      }
    })
    
    output$vis.gene.select_multiome <- renderUI({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        selectInput("vis.gene_multiome", label = "Gene to visualise",
                    choices = rownames(v$scDatan[[input$assay_multiome_seurat2a]]))
      }
    })
    
    output$vis.plot_multiome <- renderPlotly({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        withProgress(message="Generating Feature Plot...", value=0, {
          FeaturePlot(v$scDatan, input$vis.gene_multiome, cols = c("lightgrey", "darkgreen"), order = T, reduction = input$assay_multiome_seurat1a)
        })
      }
    })
    
    output$vis.gene.select_multiome1 <- renderUI({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        selectInput("vis.gene_multiome1", label = "Gene to visualise",
                    choices = rownames(v$scDatan[["ATAC"]]))
      }
    })
    
    output$vis1.plot_multiome <- renderPlotly({
      if(is.null(v$scDatan)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          FeaturePlot(v$scDatan, input$vis.gene_multiome1, cols = c("lightgrey", "darkgreen"), order = T, reduction = "umap.atac")
        })
      }
    })
    
      ##---------------Spatial Transcriptomics Analysis-------------------##
      
      output$spatial_image <- renderImage({
        list(src = "www/spatial_fig.png",
             height = 275, width=1000)
      }, deleteFile = FALSE)
      #volumes <- getVolumes()
      shinyDirChoose(
        input,
        'dir',
        roots=c(home = '.'), 
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
      )
      
      dir <- reactive(input$dir)
      output$dir <- renderPrint({  # use renderText instead of renderPrint
        parseDirPath(c(home = '.'), dir())
      })
      
      observeEvent(input$loadexample_seurat_spatial, {
        withProgress(message="Loading example Data...", value=0.5, {
          path <- "spatial/V1_Breast_Cancer_Block_A_Section_1/"
          sp.data <- Load10X_Spatial(data.dir = path, filename = "filtered_feature_bc_matrix.h5")
          sp.data[["percent.mt"]] <- PercentageFeatureSet(sp.data, pattern = "^MT-")
        })
        if (is.null(sp.data)){
          v$scData_spatial <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0, {
            print(sp.data)
            label1 <- "Example loaded"
            updateActionButton(inputId = "loadexample_seurat_spatial", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scData_spatial <- sp.data
          })
        }
      })
      
      observeEvent(input$loadButton3, {
        if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "Seurat"){
          tpmFile_spatial <- input$tpmFile_spatial
          tpmFiles_spatial <- v$scData_spatial
          if (is.null(tpmFile_spatial)){
            v$scData_spatial <- NULL
            shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
          }else{
            withProgress(message="Loading and Processing Data...", value=0, {
              print(tpmFile_spatial$datapath)
              print(tpmFile_spatial$name)
              print(file.exists(paste(tpmFile_spatial$datapath[1], "/", tpmFile_spatial$name[1], sep="")))
              exp.data_spatial <- Load10X_Spatial(parseDirPath(c(home = '.'), dir()), filename = tpmFile_spatial$name, assay = "Spatial")
              additional.ident <- NULL
              incProgress(0.5, "Creating Seurat Object")
              v$scData_spatial <- exp.data_spatial
              shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
            })
          }
        }
      }
    )
      
      observeEvent(input$reset3, {
        session$reload()
        print("Reset done")
      })
      
      observeEvent(input$loadButton3a, {
        if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "GraphST"){
          sc <- reticulate::import("scanpy", convert = FALSE)
          #scvi <- import("scvi", convert = FALSE)
          scipy <- reticulate::import("scipy", convert = FALSE)
          pot <- reticulate::import("ot", convert = FALSE)
          pd <- reticulate::import("pandas", convert = F)
          gt <- reticulate::import("GraphST", convert = F)
          matplotlib <- reticulate::import ("matplotlib", convert = F)
          gt1 <- reticulate::import("GraphST.preprocess", convert = F)
          gt2 <- reticulate::import("GraphST.utils", convert = F)
          #torch <- reticulate::import("torch", convert = FALSE)
          np <- reticulate::import("numpy", convert = F)
          sns <- reticulate::import("seaborn", convert = F)
          skmisc <- reticulate::import("skmisc", convert = F)
          
          #device = torch$device('cuda:0')
          
          tpmFile_graphst <- input$tpmFile_graphst
          if (is.null(tpmFile_graphst)){
            v$scData_graphst <- NULL
            shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
          }else{
            withProgress(message="Loading and Processing Data...", value=0, {
              print(tpmFile_graphst$datapath)
              print(tpmFile_graphst$name)
              print(file.exists(paste(tpmFile_graphst$datapath[1], "/", tpmFile_graphst$name[1], sep="")))
              file_fold <- parseDirPath(c(home = '.'), dir_graphst())
              adata <- sc$read_visium(file_fold, count_file=tpmFile_graphst$name, load_images=TRUE)
              print(adata)
              adata$var_names_make_unique()
              adata$obs_names_make_unique()
              additional.ident <- NULL
              incProgress(0.5, "Creating Seurat Object")
              v$scData_spatial <- exp.data_spatial
              shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
            })
          }
        }
      }
      )
      
      observeEvent(input$reset3a, {
        session$reload()
        print("Reset done")
      })
      
      output$countdataDT_spatial <- renderDataTable({
        if(!is.null(v$scData_spatial))
        {
          if(ncol(v$scData_spatial) > 20 )
            return(as.matrix(v$scData_spatial@assays$Spatial@counts[,1:20]))
        }
      }, server = FALSE)
      
      output$h_e <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          print(SpatialDimPlot(v$scData_spatial, pt.size.factor = 0) + NoLegend())
        }
      })
      
      observeEvent(input$loadexample_graphst_spatial, {
        withProgress(message = "Loading example data...", value = 0.3, {
          if(input$scAnalysis_sp == "GraphST"){
            sc <- reticulate::import("scanpy", convert = FALSE)
            #scvi <- import("scvi", convert = FALSE)
            scipy <- reticulate::import("scipy", convert = FALSE)
            pot <- reticulate::import("ot", convert = FALSE)
            pd <- reticulate::import("pandas", convert = F)
            gt <- reticulate::import("GraphST", convert = F)
            matplotlib <- reticulate::import ("matplotlib", convert = F)
            gt1 <- reticulate::import("GraphST.preprocess", convert = F)
            gt2 <- reticulate::import("GraphST.utils", convert = F)
            np <- reticulate::import("numpy", convert = F)
            sns <- reticulate::import("seaborn", convert = F)
            #torch <- reticulate::import("torch", convert = F)
            skmisc <- reticulate::import("skmisc", convert = F)
            
            #device = torch$device('cuda:0')
            
            withProgress(message="Loading and Processing Data...", value=0, {
              
              adata <- sc$read_visium('spatial/V1_Breast_Cancer_Block_A_Section_1/', count_file='filtered_feature_bc_matrix.h5', load_images=TRUE)
              print(adata)
              adata$var_names_make_unique()
              adata$obs_names_make_unique()
              #sce = zellkonverter::AnnData2SCE(adata = adata, layers = T)
              #counts(sce) <- assay(sce, "X")
              #sce <- logNormCounts(sce)
              #spatial <- as.Seurat(x = sce, counts = "X", data = "logcounts")
              sp.data <- Load10X_Spatial(data.dir = 'spatial/V1_Breast_Cancer_Block_A_Section_1/', filename = "filtered_feature_bc_matrix.h5")
              sp.data[["percent.mt"]] <- PercentageFeatureSet(sp.data, pattern = "^MT-")
              v$scData_spatial <- sp.data
              v$scData_spatial_graphst <- adata
              #print(v$scData_spatial1)
              #print(v$scData_spatial_graphst)
            })
          }
        })
      })
      
      output$h_e_graphst <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          print(SpatialDimPlot(v$scData_spatial, pt.size.factor = 0) + NoLegend())
        }
      })
      
      output$countdataDT_spatial_graphst <- renderDataTable({
        if(!is.null(v$scData_spatial))
        {
          if(ncol(v$scData_spatial) > 20 )
            return(as.matrix(v$scData_spatial@assays$Spatial@counts[,1:20]))
        }
      }, server = FALSE)
      
      output$countdataDT_xenium <- renderDataTable({
        if(!is.null(v$scData_xenium))
        {
          if(ncol(v$scData_xenium) > 20 )
            return(as.matrix(v$scData_xenium@assays$Xenium@counts[,1:20]))
        }
      }, server = FALSE)
      
      output$h_e_xenium <- renderPlot({
        if(is.null(v$scData_xenium)){
          plotly_empty()
        }else{
          print(ImageDimPlot(v$scData_xenium, fov = "fov") + NoLegend())
        }
      })
      
      observeEvent(input$filter_spatial, {
          tpmFile_spatial <- v$scData_spatial
          if (is.null(tpmFile_spatial)){
            v$scData_spatial <- NULL
            shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
          }else{
        withProgress(message = "Filter low quality cells...", value = 0.3, {
          v$scData_spatial <- subset(v$scData_spatial, subset = nCount_Spatial > input$obs & nFeature_Spatial > input$obs1 & percent.mt < input$obs2)
            })
          }
        #}
      })
    
      output$nCount_SpatialPlot <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          VlnPlot(v$scData_spatial, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), pt.size = 0.1) + NoLegend()
        }
      })
      
      output$SpatialFeaturePlot <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          SpatialFeaturePlot(v$scData_spatial, features = c("nCount_Spatial","nFeature_Spatial","percent.mt")) + theme(legend.position = "top")
        }
      })
      
      observeEvent(input$doSCTransform_sp, {
        if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "Seurat" | input$scAnalysis_sp == "GraphST"){
          tpmFile_spatial <- v$scData_spatial
          if (is.null(tpmFile_spatial)){
            v$scData_spatial <- NULL
            shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
          }else{
          withProgress(message = "Running scTransform...", value = 0,{
          incProgress(0.5, message = "Running scTransform...")
          v$scData_spatial <- SCTransform(v$scData_spatial, variable.features.n = input$var.genes_sp, assay = "Spatial", verbose = FALSE, conserve.memory = T)
          incProgress(0.5)
          #VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
          #output$nVarGenes <- renderText(VarGeneText)
          varGenePlotInput <- function(){
            if(is.null(v$scData_spatial)){
              return(NULL)
            }else{
              withProgress(message="Plotting variable genes...", value=0, {
                top10 <- head(VariableFeatures(v$scData_spatial), 10)
                variable_feature1 <- VariableFeaturePlot(v$scData_spatial)
                variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
                print (variable_feature1)
                print (variable_feature2)
                shinyalert("scTransform done", "scTransform done, please visualize gene expression", type = "success", imageWidth = 10, imageHeight = 10)
              })
            }
          }
          
          output$VarGenes_sp <- renderPlot({
            varGenePlotInput()
          }, height = 500, width = 600)
          observeEvent(input$PDFc, {
            if(!is.null(v$scData_spatial)){
              withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Spatial_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                  dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                  filename2 <- paste0(pdfDir, .Platform$file.sep,
                                      "Var_genes_plot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                plot1 <- VariableFeaturePlot(v$scData_spatial)
                print(plot1)
                dev.off()
                txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
                txtfile <- sub(".pdf", ".txt", txtfile)
                write(v$scData_spatial@assays$Spatial@var.features, file = txtfile)
              })
            }
          })
        })
      }
    }
  })
      
      #observeEvent(input$doSCTransform_sp, {
      #  if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "Seurat" | input$scAnalysis_sp == "GraphST"){
      #   withProgress(message = "Running scTransform...", value = 0,{
      #     incProgress(0.5, message = "Running scTransform...")
      #     v$scData_spatial <- SCTransform(v$scData_spatial, variable.features.n = input$var.genes_sp, assay = "Spatial", verbose = FALSE, conserve.memory = T)
      #     incProgress(0.5)
      #     #VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
      #     #output$nVarGenes <- renderText(VarGeneText)
      #      varGenePlotInput <- function(){
      #       if(is.null(v$scData_spatial)){
      #         return(NULL)
      #       }else{
      #         withProgress(message="Plotting variable genes...", value=0, {
      #           top10 <- head(VariableFeatures(v$scData_spatial), 10)
      #           variable_feature1 <- VariableFeaturePlot(v$scData_spatial)
      #           variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
      #           print (variable_feature1)
      #           print (variable_feature2)
      #           shinyalert("scTransform done", "scTransform done, please visualize gene expression", type = "success", imageWidth = 10, imageHeight = 10)
      #         })
      #       }
      #     }
      #     output$VarGenes_sp <- renderPlot({
      #       varGenePlotInput()
      #     }, height = 800, width = 850)
      #     observeEvent(input$PDFc, {
      #       if(!is.null(v$scData_spatial)){
      #         withProgress(message="Downloading plot PDF files...", value=0, {
      #           print(getwd())
      #           pdfDir <- paste0(getwd(), .Platform$file.sep, "Spatial_results/Generated_reports_", Sys.Date())
      #           if(!dir.exists(pdfDir)){
      #             dir.create(pdfDir)
      #           }
      #           filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
      #           i = 0
      #           while(file.exists(filename2)){
      #             filename2 <- paste0(pdfDir, .Platform$file.sep,
      #                                 "Var_genes_plot_",
      #                                 Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
      #             i = i + 1;
      #           }
      #           prePlot()
      #           pdf(filename2,
      #               width=as.numeric(input$pdf_w),
      #               height=as.numeric(input$pdf_h))
      #           plot1 <- VariableFeaturePlot(v$scData_spatial)
      #           print(plot1)
      #           dev.off()
      #           txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
      #           txtfile <- sub(".pdf", ".txt", txtfile)
      #           write(v$scData_spatial@assays$Spatial@var.features, file = txtfile)
      #         })
      #       }
      #     })
      #   })
      # }
      #})
      
      observeEvent(input$doSCTransform_xenium, {
        if(input$scAnalysis_platform == "Xenium"){
           tpmFile_xenium <- v$scData_xenium
          if (is.null(tpmFile_xenium)){
            v$scData_xenium <- NULL
            shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
          }else{
          withProgress(message = "Running scTransform...", value = 0,{
            incProgress(0.5, message = "Running scTransform...")
            v$scData_xenium <- SCTransform(v$scData_xenium, variable.features.n = input$var.genes_xenium, assay = "Xenium", verbose = FALSE, conserve.memory = T)
            incProgress(0.5)
            #VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
            #output$nVarGenes <- renderText(VarGeneText)
            varGenePlotInput <- function(){
              if(is.null(v$scData_xenium)){
                return(NULL)
              }else{
                withProgress(message="Plotting variable genes...", value=0, {
                  top10 <- head(VariableFeatures(v$scData_xenium), 10)
                  variable_feature1 <- VariableFeaturePlot(v$scData_xenium)
                  variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
                  print (variable_feature1)
                  print (variable_feature2)
                  shinyalert("scTransform done", "scTransform done, please visualize gene expression", type = "success", imageWidth = 10, imageHeight = 10)
                })
              }
            }
            
            output$VarGenes_xenium <- renderPlot({
              varGenePlotInput()
            }, height = 500, width = 600)
            observeEvent(input$PDFc, {
              if(!is.null(v$scData_xenium)){
                withProgress(message="Downloading plot PDF files...", value=0, {
                  print(getwd())
                  pdfDir <- paste0(getwd(), .Platform$file.sep, "Spatial_results/Generated_reports_", Sys.Date())
                  if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                  }
                  filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
                  i = 0
                  while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "Var_genes_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                  }
                  prePlot()
                  pdf(filename2,
                      width=as.numeric(input$pdf_w),
                      height=as.numeric(input$pdf_h))
                  plot1 <- VariableFeaturePlot(v$scData_xenium)
                  print(plot1)
                  dev.off()
                  txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
                  txtfile <- sub(".pdf", ".txt", txtfile)
                  write(v$scData_xenium@assays$Xenium@var.features, file = txtfile)
                })
              }
            })
          })
         }
        }
      })
      
      output$sp.gene.select <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("sp.gene", label = "Visualize normalized expression",
                      choices = rownames(v$scData_spatial[["Spatial"]]))
        }
      })
      
      output$sp.plot <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating Feature Plot...", value=0, {
            SpatialFeaturePlot(v$scData_spatial, features = input$sp.gene)
          })
        }
      })
      
      observeEvent(input$runPCA_spatial, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running PCA...", value = 0,{
          incProgress(0.5, message = "Running PCA...")
          v$scData_spatial <- RunPCA(v$scData_spatial, assay = "SCT", verbose = FALSE)
          #DefaultAssay(v$scDatan) <- "SCT"
          print(v$scData_spatial[["pca"]], dims = 1:5, nfeatures = 5)
          v$isPCAdone <- TRUE
          PCA_plot1a <- DimPlot(v$scData_spatial, reduction = "pca", label = T)
          print(PCA_plot1a)
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData_spatial)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
        }
        )}
      })
      
      output$PCAplot_spatial <- renderPlotly({
        if(is.null(v$isPCAdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of spatial dataset...", value=0, {
            DimPlot(v$scData_spatial, reduction = "pca", label = T)
          })
        }
      })
      
      output$vizPlot_spatial <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          VizDimLoadings(v$scData_spatial, dims = as.numeric(input$select.pc_spatial))
        }
      })
      
      output$PCHeatmap_spatial <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          DimHeatmap(v$scData_spatial, dims = as.numeric(input$select.pc_spatial), assays = "SCT")
        }
      })
      
      output$PCtable_spatial <- DT::renderDataTable({
        if(is.null(v$scData_spatial) ){
          return(NULL)
        }else{
          v$pcGenes
        }
      }, server = FALSE, options = list(scrollX = TRUE))
      
      output$Elbow_spatial <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData_spatial, ndims = 50)
          })
        }
      }, height = 400, width = 450)
      
      observeEvent(input$runPCA_xenium, {
        tpmFile_xenium <- v$scData_xenium
        if (is.null(tpmFile_xenium)){
          v$scData_xenium <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running PCA...", value = 0,{
          incProgress(0.5, message = "Running PCA...")
          v$scData_xenium <- RunPCA(v$scData_xenium, npcs = 30, features = rownames(v$scData_xenium))
          #DefaultAssay(v$scDatan) <- "SCT"
          print(v$scData_xenium[["pca"]], dims = 1:5, nfeatures = 5)
          v$isPCAXeniumdone <- TRUE
          PCA_plot1a <- DimPlot(v$scData_xenium, reduction = "pca", label = T)
          print(PCA_plot1a)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData_xenium)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          }
        )}
      })
      
      output$PCAplot_xenium <- renderPlotly({
        if(is.null(v$isPCAXeniumdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of spatial dataset...", value=0, {
            DimPlot(v$scData_xenium, reduction = "pca", label = T)
          })
        }
      })
      
      output$Elbow_xenium <- renderPlot({
        if(is.null(v$scData_xenium)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData_xenium, ndims = 30)
          })
        }
      }, height = 400, width = 450)
      
      shinyDirChoose(
        input,
        'dir_graphst',
        roots = c(home = '.'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
      )
      
      dir_graphst <- reactive(input$dir_graphst)
      output$dir_graphst <- renderPrint({  # use renderText instead of renderPrint
        parseDirPath(c(home = '.'), dir_graphst())
      })
      
      observeEvent(input$findCluster_spatial, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Finding clusters...", value = 0.3, {
          if(input$scAnalysis_sp1 == "Seurat"){
            v$scData_spatial <- FindNeighbors(v$scData_spatial, reduction = "pca", dims = 1:input$dim.used_spatial)
            v$scData_spatial <- FindClusters(v$scData_spatial, resolution = input$clus.res_spatial, graph.name = "SCT_nn")
            #output$cluster1.done <- renderText(paste0("Clustering done!"))
            v$isClusterSpatialdone <- TRUE
            }
          })
        }
      })
      
      output$Cluster2DPlot_spatial <- renderPlotly({
        if(is.null(v$isClusterSpatialdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating 2D Cluster Plot...", value=0, {
            DimPlot(v$scData_spatial, reduction = "umap", label = T)
          })
        }
      })
      
      output$silhouette_seurat <- renderText({
        if(is.null(v$scData_spatial)|| is.null(v$isClusterSpatialdone)){
          return(NULL)
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            print(v$scData_spatial@reductions$pca@cell.embeddings)
            print(v$scData_spatial@meta.data)
            v$test <- paste("Silhouette Score:", mean(silhouette(x = as.numeric(x = as.factor(x = v$scData_spatial$seurat_clusters)), dist = dist(x = v$scData_spatial@reductions$pca@cell.embeddings))[, 3]))
            print(v$test)
          })
        }
      })
      
      observeEvent(input$findCluster_graphst, {
        tpmFile_spatial <- input$tpmFile_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Finding clusters...", value = 0.3, {
          if(input$scAnalysis_sp1 == "GraphST"){
            sc <- reticulate::import("scanpy", convert = FALSE)
            #scvi <- import("scvi", convert = FALSE)
            scipy <- reticulate::import("scipy", convert = FALSE)
            pot <- reticulate::import("ot", convert = FALSE)
            pd <- reticulate::import("pandas", convert = F)
            gt <- reticulate::import("GraphST", convert = F)
            matplotlib <- reticulate::import ("matplotlib", convert = F)
            gt1 <- reticulate::import("GraphST.preprocess", convert = F)
            gt2 <- reticulate::import("GraphST.utils", convert = F)
            np <- reticulate::import("numpy", convert = F)
            sns <- reticulate::import("seaborn", convert = F)
            #torch <- reticulate::import("torch", convert = F)
            skmisc <- reticulate::import("skmisc", convert = F)
            
            #device = torch$device('cuda:0')
            
              withProgress(message="Clustering data...", value=0, {
                model = gt$GraphST$GraphST(adata = v$scData_spatial_graphst)
                v$scData_spatial_graphst  = model$train()
                print(v$scData_spatial_graphst)
                cluster <- gt$clustering(adata = v$scData_spatial_graphst , n_clusters = input$cluster_graphst, method = input$cluster_method_graphst)
                v$isClusterSpatialdone2 <- TRUE
                print(v$scData_spatial_graphst)
                print(cluster)
                v$scData_spatial_graphst_obs = v$scData_spatial_graphst$obs
                v$scData_spatial_graphst_obs$to_csv('adata_obs.csv')
                #v$scData_spatial_graphst$write_csvs(dirname = 'test/')
                print(v$scData_spatial_graphst_obs)
                meta <- read.csv('adata_obs.csv', header = T, row.names = 1)
                v$scData_spatial <- AddMetaData(v$scData_spatial, metadata = meta)
                print(v$scData_spatial)
              })
            }
          })
        }
      }) 
      
      output$Cluster2DPlot_graphst <- renderPlot({
        if(is.null(v$isClusterSpatialdone2)){
          plotly_empty()
        }else{
          withProgress(message="Generating 2D Cluster Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, group.by = "mclust", label = T)
          })
        }
      })
      
      output$DimPlot_GraphST <- renderPlotly({
        if(is.null(v$scData_spatial)|| is.null(v$isClusterSpatialdone2)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP 2D Plot...", value=0, {
            DimPlot(v$scData_spatial, group.by = "mclust", label = T)
          })
        }
      })
      
      output$SpatialDimPlot_GraphST <- renderPlot({
        if(is.null(v$scData_spatial) || is.null(v$isClusterSpatialdone2)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP 2D Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, group.by = "mclust", label = TRUE)
          })
        }
      })
      
      output$silhouette_graphst <- renderText({
        if(is.null(v$scData_spatial)|| is.null(v$isClusterSpatialdone2)){
          return(NULL)
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            print(v$scData_spatial@reductions$pca@cell.embeddings)
            print(v$scData_spatial@meta.data)
            v$test <- paste("Silhouette Score:", mean(silhouette(x = as.numeric(x = as.factor(x = v$scData_spatial$mclust)), dist = dist(x = v$scData_spatial@reductions$pca@cell.embeddings))[, 3]))
            print(v$test)
          })
        }
      })
      
      observeEvent(input$findCluster_xenium, {
        withProgress(message = "Finding clusters...", value = 0.3, {
          v$scData_xenium <- FindNeighbors(v$scData_xenium, reduction = "pca", dims = 1:input$dim.used_xenium)
          v$scData_xenium <- FindClusters(v$scData_xenium, resolution = input$clus.res_xenium)
          #output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterXeniumdone <- TRUE
          
        })
      })
      
      output$Cluster2DPlot_xenium <- renderPlotly({
        if(is.null(v$isClusterXeniumdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating 2D Cluster Plot...", value=0, {
            DimPlot(v$scData_xenium, reduction = "umap", label = T)
          })
        }
      })
      
      observeEvent(input$runUMAP_spatial, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running UMAP...", value = 0.3, {
          v$scData_spatial <- RunUMAP(v$scData_spatial, dims = 1:input$dim.used_spatial, assay = "SCT", spread = 1)
          v$isUMAPSpataildone <- TRUE
          })
        }
      })
      
      output$DimPlot_spatial <- renderPlotly({
        if(is.null(v$scData_spatial) || is.null(v$isUMAPSpataildone)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP 2D Plot...", value=0, {
            DimPlot(v$scData_spatial, reduction = "umap", label = T)
          })
        }
      })
      
      output$SpatialDimPlot <- renderPlot({
        if(is.null(v$scData_spatial) || is.null(v$isClusterSpatialdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP 2D Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
          })
        }
      })
      
      observeEvent(input$runUMAP_xenium, {
        withProgress(message = "Running UMAP...", value = 0.3, {
          v$scData_xenium <- RunUMAP(v$scData_xenium, dims = 1:input$dim.used_xenium, assay = "SCT", spread = 1)
          v$isUMAPXeniumdone <- TRUE
        })
      })
      
      output$umap_xenium <- renderPlotly({
        if(is.null(v$scData_xenium) || is.null(v$isUMAPXeniumdone)){
          plotly_empty()
        }else{
          DimPlot(v$scData_xenium)
        }
      })
      
      output$xenium.gene3.select <- renderUI({
        if(is.null(v$scData_xenium)){
          return(NULL)
        }else{
          selectInput("xenium.gene3", label = "Genes to visualize",
                      choices = rownames(v$scData_xenium@assays$Xenium@counts), selected = rownames(v$scData_xenium@assays$Xenium@counts)[20])
        }
      })
      
      output$xenium_feature.plot <- renderPlotly({
        if(is.null(v$scData_xenium) || is.null(v$isUMAPXeniumdone)){
          plotly_empty()
        }else{
          FeaturePlot(v$scData_xenium, features = input$xenium.gene3)
        }
      })
      
      output$spatialumap_xenium <- renderPlot({
        if(is.null(v$scData_xenium) || is.null(v$isUMAPXeniumdone)){
          plotly_empty()
        }else{
          ImageDimPlot(v$scData_xenium, cols = "polychrome", size = 0.75, molecules = input$xenium.gene3)
        }
      })
      
      observeEvent(input$Vis_sp1, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Visualizing...", value = 0,{
          incProgress(0.5, message = "Visualizing...")
          DefaultAssay(v$scData_spatial) <- "Spatial"
          #Idents(v$scData_spatial) <- "seurat_clusters"
          })
        }
      })
      
      output$sp.cluster.select <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("sp.cluster", label = "Cluster to visualise",
                      choices = unique(v$scData_spatial@meta.data$seurat_clusters))
        }
      })
      
      output$sp1.plot <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating Feature Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, cells.highlight = CellsByIdentities(object = v$scData_spatial, idents = input$sp.cluster), facet.highlight = TRUE)
          })
        }
      })
      
      output$sp.cluster1.select <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("sp.cluster1", label = "Cluster to visualise",
                      choices = unique(v$scData_spatial@meta.data$mclust))
        }
      })
      
      output$sp2.plot <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating Feature Plot...", value=0, {
            v$scData_spatial$mclust -> Idents(v$scData_spatial)
            SpatialDimPlot(v$scData_spatial, cells.highlight = CellsByIdentities(object = v$scData_spatial, idents = input$sp.cluster1), facet.highlight = TRUE)
          })
        }
      })
      
      observeEvent(input$Vis_sp_xenium, {
        withProgress(message = "Visualizing...", value = 0,{
          incProgress(0.5, message = "Visualizing...")
          DefaultAssay(v$scData_xenium) <- "Xenium"
          #Idents(v$scData_spatial) <- "seurat_clusters"
        })
      })
      
      output$xenium.cluster.select <- renderUI({
        if(is.null(v$scData_xenium)){
          return(NULL)
        }else{
          selectInput("xenium.cluster", label = "Cluster to visualise",
                      choices = unique(v$scData_xenium@meta.data$seurat_clusters))
        }
      })
      
      output$xenium1.plot <- renderPlot({
        if(is.null(v$scData_xenium)){
          return(NULL)
        }else{
          withProgress(message="Generating Feature Plot...", value=0, {
            ImageDimPlot(v$scData_xenium, size = 0.75, cells = colnames(subset(v$scData_xenium, seurat_clusters == input$xenium.cluster)))
          })
        }
      })
      
      observeEvent(input$loadexample_ref, {
        withProgress(message="Loading example reference scRNA-seq data...", value=0.5, {
          tpmFiles_ref <- read.csv('gene_counts_reference1.csv', header = T, row.names = 1, check.names = F)
          meta_ref <- read.csv('metadata_reference1.csv', header = T, row.names = 1)
          tpmFiles_ref <- CreateSeuratObject(tpmFiles_ref, meta.data = meta_ref)
          tpmFiles_ref$cell_type = tpmFiles_ref$primary.predict
          tpmFiles_ref$cell_type <- gsub(" ", "_", tpmFiles_ref$cell_type)
          tpmFiles_ref$cell_type <- gsub("/", "_", tpmFiles_ref$cell_type)
          tpmFiles_ref <- Seurat::SCTransform(tpmFiles_ref, ncells = 3000, verbose = FALSE)
          tpmFiles_ref <- Seurat::RunPCA(tpmFiles_ref, verbose = FALSE)
          tpmFiles_ref <- Seurat::RunUMAP(tpmFiles_ref, dims = 1:30, verbose = FALSE)
          #tpmFiles_ref <- Seurat::FindNeighbors(tpmFiles_ref, dims = 1:30, verbose = FALSE)
          #tpmFiles_ref <- Seurat::FindClusters(tpmFiles_ref, verbose = FALSE)
        })
        if (is.null(tpmFiles_ref)){
          v$scRNAData <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0.5, {
            #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
            label1 <- "Example loaded"
            updateActionButton(inputId = "loadexample_ref", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scData1 <- tpmFiles_ref
          })
        }
      })
      
      observe({
        if(input$loaduser_ref > 0){
          print('1')
          session$sendCustomMessage("myCallbackHandler2", "1")
        }
      })
      
      output$scRNAPlot <- renderPlot({
        if(is.null(v$scData1)){
          plotly_empty()
        }else{
          withProgress(message="Generating scRNA-seq Plot...", value=0, {
            DimPlot(v$scData1, group.by = "primary.predict", label = T, label.size = 3, repel = T) + NoLegend()
          })
        }
      })
      
      observeEvent(input$vis_spRNA, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
          incProgress(0.5, message = "Processing...")
          spRNA_plot <- SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
          print (spRNA_plot)
          #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
          v$isVis_spRNAdone <- TRUE
          })
        }
      })
      
      output$spRNAPlot <- renderPlot({
        if(is.null(v$isVis_spRNAdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating SpatialDim Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
          })
        }
      })
      
      observeEvent(input$doDeconv, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Performing deconvolution...", value = 0,{
          incProgress(0.5, message = "Deconvoluting...")
          v$anchors <- FindTransferAnchors(reference = v$scData1, query = v$scData_spatial, normalization.method = "SCT")
          predictions.assay <- TransferData(anchorset = v$anchors, refdata = v$scData1$primary.predict, prediction.assay = TRUE,
                                            weight.reduction = v$scData_spatial[["pca"]], dims = 1:30)
          v$scData_spatial[["predictions"]] <- predictions.assay
          print(v$scData_spatial@meta.data)
          print(v$scData_spatial)
          v$isDeconv_spRNAdone <- TRUE
          #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
          })
        }
      })
      
      output$ct.select <- renderUI({
        if(is.null(v$scData1)){
          return(NULL)
        }else{
          selectInput("ct", label = "Celltype to visualise",
                      choices = unique(v$scData1$primary.predict))
        }
      })
      
      output$DeconvPlot <- renderPlot({
        if(is.null(v$scData_spatial) || is.null(v$isDeconv_spRNAdone)){
          return(NULL)
        }else{
          withProgress(message="Generating SpatialDim Plot...", value=0, {
            DefaultAssay(v$scData_spatial) <- "predictions"
            SpatialFeaturePlot(v$scData_spatial, features = input$ct, pt.size.factor = 1.6, ncol = 2, crop = TRUE)
          })
        }
      })
      
      observeEvent(input$doDeconv_graphst, {
        tpmFile_spatial <- v$scData_spatial
        if (is.null(tpmFile_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Performing deconvolution...", value = 0,{
          incProgress(0.5, message = "Deconvoluting...")
          sc <- reticulate::import("scanpy", convert = FALSE)
          #scvi <- import("scvi", convert = FALSE)
          scipy <- reticulate::import("scipy", convert = FALSE)
          pot <- reticulate::import("ot", convert = FALSE)
          pd <- reticulate::import("pandas", convert = F)
          gt <- reticulate::import("GraphST", convert = F)
          matplotlib <- reticulate::import ("matplotlib", convert = F)
          gt1 <- reticulate::import("GraphST.preprocess", convert = F)
          gt2 <- reticulate::import("GraphST.utils", convert = F)
          np <- reticulate::import("numpy", convert = F)
          sns <- reticulate::import("seaborn", convert = F)
          #torch <- reticulate::import("torch", convert = F)
          skmisc <- reticulate::import("skmisc", convert = F)
          
          #device = torch$device('cuda:0')
          
          sceasy::convertFormat(v$scData1, from = "seurat", to = "anndata", outFile = 'reference.h5ad')
          file_path = 'reference.h5ad' # Please replace 'file_path' with the scRNA download path.
          #Reading reference data
          v$adata_sc = sc$read(file_path)
          v$adata_sc$var_names_make_unique()
          
          #Pre-processing for ST data
          gt$preprocess(v$scData_spatial_graphst)
          gt$construct_interaction(v$scData_spatial_graphst)
          gt$add_contrastive_label(v$scData_spatial_graphst)
          
          #Pre-processing for reference data
          gt$preprocess(v$adata_sc)
          
          #Finding overlap genes between ST and reference data
          v$test = gt1$filter_with_overlap_gene(adata = v$scData_spatial_graphst, adata_sc = v$adata_sc)
          print(v$test[0])
          print(v$test[1])
          gt1$get_feature(v$test[0])
          model = gt$GraphST$GraphST(v$test[0], v$test[1], deconvolution = T)
          v$test = model$train_map()
          print(v$test[0])
          print("test")
          print(v$test[1])
          gt2$project_cell_to_spot(adata = v$test[0], adata_sc = v$test[1], retain_percent = 0.15)
          v$test_obs = v$test[0]$obs
          v$test_obs$to_csv('bdata_obs.csv')
          #v$test[0]$write_csvs('ref1/')
          meta1 <- read.csv('bdata_obs.csv', header = T, row.names = 1)
          v$scData_spatial <- AddMetaData(v$scData_spatial, metadata = meta1)
          })
        }
      })
      
      output$graphst_ct.select <- renderUI({
        if(is.null(v$scData1)){
          return(NULL)
        }else{
          selectInput("ct_graphst", label = "Celltype to visualise",
                      choices = unique(v$scData1$cell_type))
        }
      })
      
      output$DeconvPlot_graphst <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating SpatialDim Plot...", value=0, {
            #DefaultAssay(v$scData_spatial) <- "predictions"
            SpatialFeaturePlot(v$scData_spatial, features = input$ct_graphst, pt.size.factor = 1.6, ncol = 2, crop = TRUE)
          })
        }
      })
      
      shinyDirChoose(
        input,
        'dir_xenium',
        roots = c(home = '.'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
      )
      
      dir_xenium <- reactive(input$dir_xenium)
      output$dir_xenium <- renderPrint({  # use renderText instead of renderPrint
        parseDirPath(c(home = '.'), dir_xenium())
      })
      
      observeEvent(input$loadexample_xenium, {
        withProgress(message="Loading example Data...", value=0.5, {
          path <- "Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs/"
          xenium.obj <- LoadXenium(path, fov = "fov")
          xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
        })
        if (is.null(xenium.obj)){
          v$scData_xenium <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0, {
            print(xenium.obj)
            label_xenium <- "Example loaded"
            updateActionButton(inputId = "loadexample_xenium", label = label_xenium)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scData_xenium <- xenium.obj
          })
        }
      })
      
      observeEvent(input$load_xenium, {
        withProgress(message="Loading and Processing Xenium Data...", value=0, {
          exp.data_xenium <- LoadXenium(data.dir = parseDirPath(c(home = '.'), dir_xenium()), fov = "fov")
          exp.data_xenium <- subset(exp.data_xenium, subset = nCount_Xenium > input$count_xenium)
          incProgress(0.5, "Creating xenium object")
          v$scData_xenium <- exp.data_xenium
          print(v$scData_xenium)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    )
      
      observeEvent(input$reset_xenium, {
        session$reload()
        print("Reset done")
      })
      
    output$nFeature_xenium <- renderPlot({
      if(is.null(v$scData_xenium)){
        plotly_empty()
      }else{
        VlnPlot(v$scData_xenium, features = "nFeature_Xenium", pt.size = 0) + NoLegend()
      }
    }) 
    
    output$nCount_xenium <- renderPlot({
      if(is.null(v$scData_xenium)){
        plotly_empty()
      }else{
        VlnPlot(v$scData_xenium, features = "nCount_Xenium", pt.size = 0) + NoLegend()
      }
    })
    
    output$FeatureScatterPlot_xenium <- renderPlotly({
      if(is.null(v$scData_xenium)){
        plotly_empty()
      }else{
        print(FeatureScatter(v$scData_xenium, "nCount_Xenium", "nFeature_Xenium"))
      }
    })
    
    output$xenium.gene.select <- renderUI({
      if(is.null(v$scData_xenium)){
        return(NULL)
      }else{
        selectInput("xenium.gene", label = "Molecules to visualise",
                    choices = rownames(v$scData_xenium@assays$Xenium@counts), multiple = T, selected = rownames(v$scData_xenium@assays$Xenium@counts)[20])
      }
    })
    
    output$markerplot_xenium <- renderPlot({
      if(is.null(v$scData_xenium)){
        plotly_empty()
      }else{
        ImageDimPlot(v$scData_xenium, fov = "fov", molecules = input$xenium.gene, nmols = 20000)
      }
    })
    
    output$xenium.gene1.select <- renderUI({
      if(is.null(v$scData_xenium)){
        return(NULL)
      }else{
        selectInput("xenium.gene1", label = "Molecules to visualise",
                    choices = rownames(v$scData_xenium@assays$Xenium@counts), selected = rownames(v$scData_xenium@assays$Xenium@counts)[20])
      }
    })
    
    output$featureplot_xenium <- renderPlot({
      if(is.null(v$scData_xenium)){
        plotly_empty()
      }else{
        ImageFeaturePlot(v$scData_xenium, features = input$xenium.gene1, size = 0.75, cols = c("white", "red"))
      }
    })
    
    observeEvent(input$crop_xenium, {
      withProgress(message = "Cropping...", value = 0.3, {
        cropped.coords <- Crop(v$scData_xenium[["fov"]], x = c(input$x1_xenium, input$x2_xenium), y = c(input$y1_xenium, input$y2_xenium), coords = "plot")
        print(cropped.coords)
        v$scData_xenium[["zoom"]] <- cropped.coords
        # visualize cropped area with cell segmentations & selected molecules
        DefaultBoundary(v$scData_xenium[["zoom"]]) <- "segmentation"
        v$isCropXeniumdone <- TRUE
      })
    })
    
    output$xenium.gene2.select <- renderUI({
      if(is.null(v$scData_xenium)){
        return(NULL)
      }else{
        selectInput("xenium.gene2", label = "Molecules to visualize",
                    choices = rownames(v$scData_xenium@assays$Xenium@counts), multiple = T, selected = rownames(v$scData_xenium@assays$Xenium@counts)[20])
      }
    })
    
    output$zoom_xenium <- renderPlot({
      if(is.null(v$scData_xenium) || is.null(v$isCropXeniumdone)){
        plotly_empty()
      }else{
        ImageDimPlot(v$scData_xenium, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome", coord.fixed = FALSE, molecules = input$xenium.gene2, nmols = 10000)
      }
    })
    
    output$gene_sp.select <- renderUI({
      if(is.null(v$scData_spatial)){
        return(NULL)
      }else{
        if(input$deg_sp1 == "seurat_clusters"){
          selectInput("gene_sp1", label = "Celltype1",
                      choices = as.vector(v$scData_spatial$seurat_clusters))
        }
      }
    })
    
    output$gene_sp1.select <- renderUI({
      if(is.null(v$scData_spatial)){
        return(NULL)
      }else{
        if(input$deg_sp1 == "seurat_clusters"){
          selectInput("gene_sp2", label = "Celltype2",
                      choices = as.vector(v$scData_spatial$seurat_clusters))
        }
      }
    })
    
    observeEvent(input$doVolcano_sp, {
      tpmFiles <- v$scData_spatial
      if (is.null(tpmFiles)){
        v$scData_spatial <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
            DefaultAssay(v$scData_spatial) <- "SCT"
            v$scData_spatial_subset <- subset(v$scData_spatial, subset = seurat_clusters == input$gene_sp1 | seurat_clusters == input$gene_sp2)
            ips.markers_sp_a <- FindAllMarkers(v$scData_spatial_subset, only.pos = F, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1, assay = "SCT")
            ips.markers_sp_b <- FindMarkers(v$scData_spatial_subset, ident.1 = input$gene_sp1, ident.2 = input$gene_sp2, only.pos = F, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
            v$ips.markers_sp_a <- ips.markers_sp_a
            v$ips.markers_sp_b <- ips.markers_sp_b
          shinyalert("Pairwise DEGs done", "Pairwise DEGs done, please run GSEA Analysis", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$volcano.plot_sp1 <- renderPlot({
      if(is.null(v$ips.markers_sp_b)){
        return(NULL)
      }else{
        withProgress(message="Generating Volcano Plot...", value=0, {
          EnhancedVolcano(toptable = v$ips.markers_sp_b, lab = row.names(v$ips.markers_sp_b), x ="avg_log2FC", y ="p_val_adj", pointSize = 1, labSize = 5, legendLabSize = 12, axisLabSize = 12)
        })
      }
    })
    
    output$dega_sp.plot <- renderPlot({
      if(is.null(v$scData_spatial_subset)){
        return(NULL)
      }else{
        withProgress(message="Generating Volcano Plot...", value=0, {
          v$ips.markers_sp_a %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData_spatial_subset, features = top10$gene, size = 5, angle = 45) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
        })
      }
    })
      
    observeEvent(input$doDeg_spatial, {
      tpmFile_spatial <- v$scData_spatial
      if (is.null(tpmFile_spatial)){
        v$scData_spatial <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
          withProgress(message="Finding DEGs...", value=0, {
            if(input$deg_spatial == "seurat_clusters"){
            #DefaultAssay(v$scData_spatial) <- "SCT"
            print(v$scData_spatial@meta.data)
            v$scData_spatial$seurat_clusters -> Idents(v$scData_spatial)
            ips.markers_spatial <- FindAllMarkers(v$scData_spatial, only.pos = FALSE, min.pct = input$min_pct_spatial, logfc.threshold = input$logfc_spatial, test.use = input$test.use_spatial, assay = "SCT")
            v$ips.markers_spatial <- ips.markers_spatial
            }
          })
        }
      })
    
    observeEvent(input$Vis_spatial, {
      tpmFile_spatial <- v$scData_spatial
      if (is.null(tpmFile_spatial)){
        v$scData_spatial <- NULL
        shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
      }else{
        withProgress(message="Visualizing...", value=0, {
          v$isVisSpatialdone <- TRUE
        })
      }
    })
      
      output$deg.gene.select_spatial <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("deg.gene_spatial", label = "Gene to visualise",
                      choices = rownames(v$scData_spatial[["SCT"]]))
        }
      })
      
      plotViolin_spatial <- reactive({
        if(is.null(v$scData_spatial) || is.null(v$isVisSpatialdone)){
          plotly_empty()
        }else{
          withProgress(message="Visualizing...", value=0, {
            DefaultAssay(v$scData_spatial) <- "SCT"
            VlnPlot(v$scData_spatial, input$deg.gene_spatial)
          })
        }
      })
      
      output$Deg.plot_spatial <- renderPlotly({
        plotViolin_spatial()
      })
      
      output$download_violn_spatial <- downloadHandler(
        filename = function(){"Violin plot (Spatial module).png"}, 
        content = function(fname){
          ggsave(fname,plotViolin_spatial(), height = 7, width = 7)
        }
      )
      
      plotFeature_spatial <- reactive({
        if(is.null(v$scData_spatial) || is.null(v$isVisSpatialdone)){
          plotly_empty()
        }else{
          withProgress(message="Visualizing...", value=0, {
            DefaultAssay(v$scData_spatial) <- "SCT"
            SpatialFeaturePlot(v$scData_spatial, input$deg.gene_spatial)
          })
        }
      })
      
      output$Deg1.plot_spatial <- renderPlot({
        plotFeature_spatial()
      })
      
      output$download_feature_spatial <- downloadHandler(
        filename = function(){"Feature plot (Spatial module).png"}, 
        content = function(fname){
          ggsave(fname,plotFeature_spatial(), height = 7, width = 7)
        }
      )
      
      plotRidge_spatial <- reactive({
        if(is.null(v$scData_spatial) || is.null(v$isVisSpatialdone)){
          plotly_empty()
        }else{
          withProgress(message="Visualizing...", value=0, {
            DefaultAssay(v$scData_spatial) <- "SCT"
            RidgePlot(v$scData_spatial, features = input$deg.gene_spatial)
          })
        }
      })
      
      output$Deg2.plot_spatial <- renderPlot({
        plotRidge_spatial()
      })
      
      output$download_ridge_spatial <- downloadHandler(
        filename = function(){"Ridge plot (Spatial module).png"}, 
        content = function(fname){
          ggsave(fname,plotRidge_spatial(), height = 7, width = 7)
        }
      )
      
      output$Deg.heatmap_spatial <- renderPlot({
        if(is.null(v$ips.markers_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            v$ips.markers_spatial %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
            DoHeatmap(v$scData_spatial, features = top10$gene, size = 4, assay = "SCT") + theme(axis.text.y = element_text(size = 5)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
          })
        }
      })
      
      output$Deg.table_spatial <- DT::renderDataTable(
        v$ips.markers_spatial, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
      
      output$gsea.ct_sp1.select <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("gsea_sp.ct1", label = "Celltype1",
                      choices = as.vector(v$scData_spatial$seurat_clusters))
        }
      })
      
      output$gsea.ct_sp2.select <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("gsea_sp.ct2", label = "Celltype2",
                      choices = as.vector(v$scData_spatial$seurat_clusters))
        }
      })
      
      observeEvent(input$gsea_sp, {
        tpmFiles_spatial <- v$scData_spatial
        if (is.null(tpmFiles_spatial)){
          v$scData_spatial <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Generating gene set enrichment analysis...", value=0, {
            if(input$species_gsea_sp == "Homo sapiens" & input$category_gsea_sp == "H"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_hs_go_sp <- msigdbr(species = "Homo sapiens", category = "H")
              print(v$msigdbr_hs_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_hs_go_sp$gene_symbol, f = v$msigdbr_hs_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Mus musculus" & input$category_gsea_sp == "H"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_mm_go_sp <- msigdbr(species = "Mus musculus", category = "H")
              print(v$msigdbr_mm_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_mm_go_sp$gene_symbol, f = v$msigdbr_mm_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Homo sapiens" & input$category_gsea_sp == "C2"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_hs_go_sp <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
              print(v$msigdbr_hs_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_hs_go_sp$gene_symbol, f = v$msigdbr_hs_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Mus musculus" & input$category_gsea_sp == "C2"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_mm_go_sp <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
              print(v$msigdbr_mm_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_mm_go_sp$gene_symbol, f = v$msigdbr_mm_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Homo sapiens" & input$category_gsea_sp == "C5"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_hs_go_sp <- msigdbr(species = "Homo sapiens", category = "C5")
              print(v$msigdbr_hs_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_hs_go_sp$gene_symbol, f = v$msigdbr_hs_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Mus musculus" & input$category_gsea_sp == "C5"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_mm_go_sp <- msigdbr(species = "Mus musculus", category = "C5")
              print(v$msigdbr_mm_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_mm_go_sp$gene_symbol, f = v$msigdbr_mm_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Homo sapiens" & input$category_gsea_sp == "C7"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_hs_go_sp <- msigdbr(species = "Homo sapiens", category = "C7")
              print(v$msigdbr_hs_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_hs_go_sp$gene_symbol, f = v$msigdbr_hs_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Mus musculus" & input$category_gsea_sp == "C7"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_mm_go_sp <- msigdbr(species = "Mus musculus", category = "C7")
              print(v$msigdbr_mm_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_mm_go_sp$gene_symbol, f = v$msigdbr_mm_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Homo sapiens" & input$category_gsea_sp == "C8"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_hs_go_sp <- msigdbr(species = "Homo sapiens", category = "C8")
              print(v$msigdbr_hs_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_hs_go_sp$gene_symbol, f = v$msigdbr_hs_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            if(input$species_gsea_sp == "Mus musculus" & input$category_gsea_sp == "C8"){
              DefaultAssay(v$scData_spatial) <- "SCT"
              v$msigdbr_mm_go_sp <- msigdbr(species = "Mus musculus", category = "C8")
              print(v$msigdbr_mm_go_sp)
              v$pathways_sp <- split(x = v$msigdbr_mm_go_sp$gene_symbol, f = v$msigdbr_mm_go_sp$gs_name)
              print(v$pathways_sp)
              v$markers_sp <- FindMarkers(v$scData_spatial, ident.1 = input$gsea_sp.ct1, ident.2 = input$gsea_sp.ct2, min.pct = input$min_pct_sp1, logfc.threshold = input$logfc_sp1, test.use = input$test.use_sp1)
              v$markers_sp  <- v$markers_sp %>% arrange(desc(avg_log2FC))
              print(v$markers_sp)
              v$markers_sp.log2FC <- v$markers_sp$avg_log2FC
              names(v$markers_sp.log2FC) <- row.names(v$markers_sp)
              v$markers_sp.log2FC <- sort(na.omit(v$markers_sp.log2FC), decreasing = TRUE)
              print(v$markers_sp.log2FC)
              v$fgseaRes_sp <- fgsea(pathways = v$pathways_sp, stats = v$markers_sp.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
              print(v$fgseaRes_sp)
              v$topPathwaysUp_sp <- v$fgseaRes_sp[ES > 0][head(order(pval), n=10), pathway]
              v$topPathwaysDown_sp <- v$fgseaRes_sp[ES < 0][head(order(pval), n=10), pathway]
              v$topPathways_sp <- c(v$topPathwaysUp_sp, rev(v$topPathwaysDown_sp))
            }
            output$gsea_sp.done <- renderText(paste0("Gene set enrichment done!"))
            v$isGSEAspdone <- TRUE
          })
        }
      })
      
      output$gsea_sp.select <- renderUI({
        if(is.null(v$pathways_sp)){
          return(NULL)
        }else{
          selectInput("gsea_sp.pathway", label = "Gene set to visualise",
                      choices = names(v$pathways_sp))
        }
      })
      
      output$gsea_sp_plot <- renderPlot({
        if(is.null(v$pathways_sp) || is.null(v$isGSEAspdone)){
          return(NULL)
        }else{
          withProgress(message="Generating GSEA plot...", value=0, {
            plotEnrichment(v$pathways_sp[[input$gsea_sp.pathway]], v$markers_sp.log2FC) + labs(title=input$gsea_sp.pathway)
          })
        }
      })
      
      output$gsea_sp_plot1 <- renderPlot({
        if(is.null(v$pathways_sp) || is.null(v$isGSEAspdone)){
          return(NULL)
        }else{
          withProgress(message="Generating GSEA plot...", value=0, {
            plotGseaTable(v$pathways_sp[v$topPathways_sp], v$markers_sp.log2FC, v$fgseaRes_sp, gseaParam=0.5)
          })
        }
      })
      
      output$gsea_sp.table <- DT::renderDataTable(
        v$fgseaRes_sp, options = list(scrollX = TRUE, scrollY = "400px"), server = FALSE)
      
      output$download_gsea_sp.table <- downloadHandler(
        filename = function(){"GSEA Results.csv"}, 
        content = function(fname){
          withProgress(message="Downloading GSEA Results...", value=0, {
            fwrite(v$fgseaRes_sp, fname)
          })
        }
      )
      
      #---------------------Single cell ATAC-seq-----------------##
      
      output$atac_image <- renderImage({
        list(src = "www/atac_fig.png",
             width = 800,
             height = 450)
      }, deleteFile = FALSE)
      
      observeEvent(input$loadButton_atac, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        meta_atac <- input$meta_atac
        meta_atac <- v$meta_atac
        print(tpmFiles_atac)
        print(meta_atac)
        names.field <- input$field
        
        #tpmFiles_atac <- input$tpmFiles_atac
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Loading and Processing Data...", value=0, {
            print(tpmFiles_atac$datapath)
            print(tpmFiles_atac$name)
            #print(file.exists(paste(tpmFiles_atac$datapath[1], "/", tpmFiles_atac$name[1], sep="")))
            exp.data_atac <- Read10X_h5(tpmFiles_atac$datapath)
            additional.ident <- NULL
            if(!is.null(meta_atac)){
              anno.data_atac <- read.csv(meta_atac$datapath[1], header = T, row.names = 1)
              print(anno.data_atac)
              #to.append <- apply(anno.data, 1, paste, collapse = "_")
              #colnames(exp.data) <- to.append
              #names.field <- match(input$groupby, colnames(anno.data))
              #additional.ident <- data.frame(data.frame(anno.data[,-1], row.names = to.append))
              #additional.ident[] <- lapply(additional.ident, factor)
            }
            incProgress(0.5, "Creating Chromatin Assay")
            print(exp.data_atac)
            #frag <- CreateFragmentObject(path = parseDirPath(c(home = '.'), dir_atac()))
            #print(frag)
            chrom_assay_atac <- CreateChromatinAssay(counts = exp.data_atac, sep = c(":", "-"), genome = 'hg19', fragments = as.character(parseFilePaths(c(home = '.'), dir_atac())$datapath), min.cells = input$min.cells_atac, min.features = input$min.genes_atac)
            print (chrom_assay_atac)
            incProgress(0.5, "Creating Seurat Object")
            sObj_atac <- CreateSeuratObject(chrom_assay_atac,
                                            assay = "peaks",
                                            project = input$projName4,
                                            meta.data = anno.data_atac)
            print(sObj_atac)
            print(sObj_atac[['peaks']])
            print(granges(sObj_atac))
            v$atacData <- sObj_atac
            print(v$atacData)
            shinyalert("Data uploaded", "Data uploaded", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
        dir.create("Seurat_ATACseq_results")
      })
      
      observeEvent(input$loadexample_atac, {
        withProgress(message="Loading example Data...", value=0.5, {
          tpmFiles_atac <- Read10X_h5(filename = "atac_pbmc_500_v1/filtered_peak_bc_matrix.h5")
          meta_atac <- read.csv('atac_pbmc_500_v1/singlecell.csv', header = T, row.names = 1)
          chrom_assay_atac <- CreateChromatinAssay(counts = tpmFiles_atac, sep = c(":", "-"), genome = 'hg19', fragments = 'atac_pbmc_500_v1/fragments.tsv.gz', min.cells = 3, min.features = 200)
          sObj_atac <- CreateSeuratObject(chrom_assay_atac, assay = "peaks", project = input$projName4, meta.data = meta_atac)
          print(sObj_atac)
          annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
          seqlevelsStyle(annotations_atac) <- 'UCSC'
          genome(annotations_atac) <- "hg19"
          print("test")
          Annotation(sObj_atac) <- annotations_atac
          sObj_atac <- NucleosomeSignal(object = sObj_atac)
          print("test2")
          sObj_atac <- TSSEnrichment(object = sObj_atac, fast = FALSE)
          print("test1")
          sObj_atac$pct_reads_in_peaks <- sObj_atac$peak_region_fragments / sObj_atac$passed_filters * 100
          sObj_atac$blacklist_ratio <- sObj_atac$blacklist_region_fragments / sObj_atac$peak_region_fragments
          sObj_atac$high.tss <- ifelse(sObj_atac$TSS.enrichment > 2, 'High', 'Low')
          sObj_atac$nucleosome_group <- ifelse(sObj_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
          v$tpmFiles_atac <- tpmFiles_atac
          v$meta_atac <- meta_atac
          v$atacData <- sObj_atac
        })
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample_atac", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
       
      })
      
      observeEvent(input$reset_atac, {
        session$reload()
        print("Reset done")
      })
      
      shinyFileChoose(
        input,
        'dir_atac',
        roots = c(home = '.'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw", "gz", "tbi")
      )
      
      dir_atac <- reactive(input$dir_atac)
      output$dir_atac <- renderPrint({  # use renderText instead of renderPrint
        as.character(parseFilePaths(c(home = '.'), dir_atac())$datapath)
      })
      
      output$countdataDT_atac <- renderDataTable({
        if(!is.null(v$atacData))
        {
          if(ncol(v$atacData) > 20 )
            return(as.matrix(v$atacData@assays$peaks@counts[1:20,]))
        }
      }, server = FALSE)
      
      observeEvent(input$create_seurat_atac, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(v$atacData)
          chrom_assay_atac <- CreateChromatinAssay(counts = v$tpmFiles_atac, sep = c(":", "-"), genome = 'hg19', fragments = 'atac_pbmc_500_v1/fragments.tsv.gz', min.cells = input$min.cells_atac, min.features = input$min.genes_atac)
          sObj_atac <- CreateSeuratObject(chrom_assay_atac, assay = "peaks", project = input$projName4, meta.data = v$meta_atac)
          print(sObj_atac)
          annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
          seqlevelsStyle(annotations_atac) <- 'UCSC'
          genome(annotations_atac) <- "hg19"
          print("test")
          Annotation(sObj_atac) <- annotations_atac
          sObj_atac <- NucleosomeSignal(object = sObj_atac)
          print("test2")
          sObj_atac <- TSSEnrichment(object = sObj_atac, fast = FALSE)
          print("test1")
          sObj_atac$pct_reads_in_peaks <- sObj_atac$peak_region_fragments / sObj_atac$passed_filters * 100
          sObj_atac$blacklist_ratio <- sObj_atac$blacklist_region_fragments / sObj_atac$peak_region_fragments
          sObj_atac$high.tss <- ifelse(sObj_atac$TSS.enrichment > 2, 'High', 'Low')
          sObj_atac$nucleosome_group <- ifelse(sObj_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
          v$atacData <- sObj_atac
          shinyalert("Seurat object created", "Seurat object created", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
      output$TSSPlot <- renderPlotly({
        if(is.null(v$atacData)){
          plotly_empty()
        }else{
          TSSPlot(v$atacData, group.by = 'high.tss')
        }
      })
      
      output$FragmentHistogram <- renderPlotly({
        if(is.null(v$atacData)){
          plotly_empty()
        }else{
          FragmentHistogram(object = v$atacData, group.by = 'nucleosome_group')
        }
      })
      
      output$Vlnplot_atac_1 <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          VlnPlot(object = v$atacData, features = 'pct_reads_in_peaks')
        }
      })
      
      output$Vlnplot_atac_2 <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          VlnPlot(object = v$atacData, features = 'peak_region_fragments')
        }
      })
      
      output$Vlnplot_atac_3 <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          VlnPlot(object = v$atacData, features = 'TSS.enrichment')
        }
      })
      
      output$Vlnplot_atac_4 <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          VlnPlot(object = v$atacData, features = 'blacklist_ratio')
        }
      })
      
      output$Vlnplot_atac_5 <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          VlnPlot(object = v$atacData, features = 'nucleosome_signal')
        }
      })
      
      observeEvent(input$donorm_ATAC, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Performing Normalization...", value = 0,{
          incProgress(0.5, message = "Running Dimension Reduction...")
          v$atacData <- RunTFIDF(v$atacData)
          v$atacData <- FindTopFeatures(v$atacData, min.cutoff = input$top_features_atac)
          v$atacData <- RunSVD(v$atacData, n = input$npc_atac)
          all_atac.genes <- rownames(v$atacData)
          v$atacData <- ScaleData(v$atacData, features = all_atac.genes)
          output$normalize_atac.done <- renderText(paste0("Normalization done!"))
          v$isNormalizeATACdone <- TRUE
          })
        }
      })
      
      output$DepthCor_ATAC <- renderPlotly({
        if(is.null(v$isNormalizeATACdone)){
          plotly_empty()
        }else{
          DepthCor(v$atacData, n = input$npc_atac)
        }
      })
      
      ##---------------Clustering of ATAC-seq data-------------------
      
      observeEvent(input$doCluster_ATAC, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Performing clustering...", value = 0.3,{
          incProgress(0.5, message = "Running Clustering...")
          v$atacData <- FindNeighbors(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac)
          v$atacData <- FindClusters(object = v$atacData, resolution = input$clus.res_atac, verbose = FALSE, algorithm = 3)
          output$cluster_atac.done <- renderText(paste0("Clustering done!"))
          v$isClusterATACdone <- TRUE
          })
        }
      })
      
      output$cluster_ATAC <- renderPlotly({
        if(is.null(v$isClusterATACdone)){
          plotly_empty()
        }else{
          DimPlot(object = v$atacData, reduction = "umap", label = TRUE)
        }
      })
      
      
      ##---------------UMAP on ATAC-seq data-------------------
      
      observeEvent(input$doUMAP_ATAC, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running UMAP...", value = 0.3,{
          incProgress(0.5, message = "Running UMAP...")
          v$atacData <- RunUMAP(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac)
          output$umap_atac.done <- renderText(paste0("UMAP done!"))
          v$isUMAPATACdone <- TRUE
          })
        }
      })
      
      output$Umap_ATAC <- renderPlotly({
        if(is.null(v$isUMAPATACdone)){
          plotly_empty()
        }else{
          DimPlot(object = v$atacData, label = TRUE)
        }
      })
      
      ##---------------TSNE on ATAC-seq data-------------------
      
      observeEvent(input$doTSNE_ATAC, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message = "Running TSNE...", value = 0.3,{
          incProgress(0.5, message = "Running TSNE...")
          v$atacData <- RunTSNE(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac)
          output$tsne_atac.done <- renderText(paste0("TSNE done!"))
          v$isTSNEATACdone <- TRUE
          })
        }
      })
      
      output$Tsne_ATAC <- renderPlotly({
        if(is.null(v$isTSNEATACdone)){
          plotly_empty()
        }else{
          DimPlot(object = v$atacData, reduction = "tsne", label = TRUE)
        }
      })
      
      observeEvent(input$calc_gene_activity_atac, {
        #tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Calculating gene activity...", value=0, {
            gene.activities <- GeneActivity(v$atacData, assay = 'peaks')
            v$atacData[['RNA']] <- CreateAssayObject(counts = gene.activities)
            v$atacData <- NormalizeData(
              object = v$atacData,
              assay = 'RNA',
              normalization.method = 'LogNormalize',
              scale.factor = median(v$atacData$nCount_RNA)
            )
            DefaultAssay(v$atacData) <- 'RNA'
            v$gene.activities <- gene.activities
            v$isGeneActivitydone <- TRUE
          })
        }
      })
      
      output$gene_activity.select <- renderUI({
        if(is.null(v$atacData) || is.null(v$isGeneActivitydone)){
         return(NULL)
        }else{
          selectInput("gene_activity", label = "Gene activity to visualise",
                      choices = rownames(v$atacData), selected = rownames(v$atacData)[10])
        }
      })
      
      plotGeneActivity <- reactive({
        if(is.null(v$atacData) || is.null(v$isGeneActivitydone)){
          plotly_empty()
        }else{
          withProgress(message="Generating Gene Activity Plot...", value=0, {
            FeaturePlot(object = v$atacData,
              features = input$gene_activity, 
              pt.size = 0.1
            )
          })
        }
      })
      
      output$gene_activity.plot <- renderPlotly({
                                    plotGeneActivity()
                                     })
      
      output$download_gene_activity <- downloadHandler(
        filename = function(){"Gene Activity plot.png"}, 
        content = function(fname){
          ggsave(fname,plotGeneActivity(), height = 7, width = 7)
        }
      )
      
      observeEvent(input$load_eg_scRNA, {
        withProgress(message="Loading example reference scRNA-seq data...", value=0.5, {
          tpmFiles_ex_scRNA <- read.table('scRNA/pbmc2.txt', header = T, row.names = 1, check.names = F)
          tpmFiles_ex_scRNA <- CreateSeuratObject(tpmFiles_ex_scRNA)
          tpmFiles_ex_scRNA[["percent.mt"]] <- PercentageFeatureSet(tpmFiles_ex_scRNA, pattern = "^MT-")
          tpmFiles_ex_scRNA <- NormalizeData(tpmFiles_ex_scRNA)
          tpmFiles_ex_scRNA <- FindVariableFeatures(tpmFiles_ex_scRNA, selection.method = "vst", nfeatures = 2000)
          all.genes <- rownames(tpmFiles_ex_scRNA)
          tpmFiles_ex_scRNA <- ScaleData(tpmFiles_ex_scRNA, features = all.genes)
          tpmFiles_ex_scRNA <- RunPCA(tpmFiles_ex_scRNA, features = VariableFeatures(object = tpmFiles_ex_scRNA))
          tpmFiles_ex_scRNA <- FindNeighbors(tpmFiles_ex_scRNA, dims = 1:30)
          tpmFiles_ex_scRNA <- FindClusters(tpmFiles_ex_scRNA, resolution = 1)
          tpmFiles_ex_scRNA <- RunUMAP(tpmFiles_ex_scRNA, dims = 1:30)
          ref = readRDS('ref.rds')
          tpmFiles_ex_scRNA.data.average = AverageExpression(tpmFiles_ex_scRNA)
          tpmFiles_ex_scRNA.data.average = round(tpmFiles_ex_scRNA.data.average$RNA, 2)
          
          if(input$cellatlas_atac == "all"){
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, ref)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "adipose"){
            adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
            adipose1 <- ref[,adipose]
            colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, adipose1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "adrenal_gland"){
            adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
            adrenal_gland1 <- ref[,adrenal_gland]
            colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, adrenal_gland1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "blood"){
            blood <- colnames(ref)[grepl("blood",colnames(ref))] 
            blood1 <- ref[,blood]
            colnames(blood1) <- gsub("--blood","",colnames(blood1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, blood1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "bone_marrow"){
            bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
            bone_marrow1 <- ref[,bone_marrow]
            colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, bone_marrow1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "brain"){
            brain <- colnames(ref)[grepl("brain",colnames(ref))] 
            brain1 <- ref[,brain]
            colnames(brain1) <- gsub("--brain","",colnames(brain1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, brain1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "breast"){
            breast <- colnames(ref)[grepl("breast",colnames(ref))] 
            breast1 <- ref[,breast]
            colnames(breast1) <- gsub("--breast","",colnames(breast1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, breast1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "breast_milk"){
            breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
            breast_milk1 <- ref[,breast_milk]
            colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, breast_milk1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "eye"){
            eye <- colnames(ref)[grepl("eye",colnames(ref))] 
            eye1 <- ref[,eye]
            colnames(eye1) <- gsub("--eye","",colnames(eye1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, eye1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "gut"){
            gut <- colnames(ref)[grepl("gut",colnames(ref))] 
            gut1 <- ref[,gut]
            colnames(gut1) <- gsub("--gut","",colnames(gut1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, gut1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "heart"){
            heart <- colnames(ref)[grepl("heart",colnames(ref))] 
            heart1 <- ref[,heart]
            colnames(heart1) <- gsub("--heart","",colnames(heart1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, heart1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "kidney"){
            kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
            kidney1 <- ref[,kidney]
            colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, kidney1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "liver"){
            liver <- colnames(ref)[grepl("liver",colnames(ref))] 
            liver1 <- ref[,liver]
            colnames(liver1) <- gsub("--liver","",colnames(liver1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, liver1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "lung"){
            lung <- colnames(ref)[grepl("lung",colnames(ref))] 
            lung1 <- ref[,lung]
            colnames(lung1) <- gsub("--lung","",colnames(lung1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, lung1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "pancreas"){
            pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
            pancreas1 <- ref[,pancreas]
            colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, pancreas1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "PDAC"){
            PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
            PDAC1 <- ref[,PDAC]
            colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, PDAC1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "skin"){
            skin <- colnames(ref)[grepl("skin",colnames(ref))] 
            skin1 <- ref[,skin]
            colnames(skin1) <- gsub("--skin","",colnames(skin1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, skin1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "testis"){
            testis <- colnames(ref)[grepl("testis",colnames(ref))] 
            testis1 <- ref[,testis]
            colnames(testis1) <- gsub("--testis","",colnames(testis1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, testis1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "thymus"){
            thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
            thymus1 <- ref[,thymus]
            colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, thymus1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
          if(input$cellatlas_atac == "tonsil"){
            tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
            tonsil1 <- ref[,tonsil]
            colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
            res = FastIntegration::CELLiD(tpmFiles_ex_scRNA.data.average, tonsil1)
            print(res)
            tpmFiles_ex_scRNA$primary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),1]
            tpmFiles_ex_scRNA$secondary.predict = res[as.numeric(tpmFiles_ex_scRNA$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(res) <- newheaders
            print(tpmFiles_ex_scRNA@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            v$isCELLAtacdone <- TRUE
          }
        })
        if (is.null(tpmFiles_ex_scRNA)){
          v$scData1 <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0.5, {
            #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
            label1 <- "Example loaded"
            updateActionButton(inputId = "load_eg_scRNA", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scData1 <- tpmFiles_ex_scRNA
            v$isCELLAtacdone <- TRUE
            #v$transfer.anchors <- FindTransferAnchors(
            #  reference = tpmFiles_ex_scRNA,
            #  query = v$atacData,
            #  reduction = 'cca'
            #)
            #v$predicted.labels <- TransferData(
            #  anchorset = v$transfer.anchors,
            #  refdata = tpmFiles_ex_scRNA$primary.predict,
            #  weight.reduction = v$atacData[['lsi']],
            #  dims = 2:30
            #)
            #print(v$predicted.labels)
            #v$atacData <- AddMetaData(object = v$atacData, metadata = v$predicted.labels)
            #print(v$atacData@meta.data)
          })
        }
      })
      
      observe({
        if(input$load_user_scRNA > 0){
          print('2')
          session$sendCustomMessage("myCallbackHandler4", "2")
        }
      })
      
      #observeEvent(input$doct_atac, {
      #  if(is.null(v$atacData)){
      #    return(NULL)
      #  }else{
      #    withProgress(message="Annotating celltypes for scATAC-seq data...", value=0, {
      #      v$transfer.anchors <- FindTransferAnchors(
      #        reference = v$scData1,
      #        query = v$atacData,
      #        reduction = 'cca'
      #      )
      #      v$predicted.labels <- TransferData(
      #        anchorset = v$transfer.anchors,
      #        refdata = v$scData1$primary.predict,
      #        weight.reduction = v$atacData[['lsi']],
      #        dims = 2:30
      #      )
      #      v$atacData <- AddMetaData(object = v$atacData, metadata = v$predicted.labels)
      #    })
      #  }
      #})
      
      output$cellanno_atac.plot <- renderPlot({
        if(is.null(v$scData1) || is.null(v$isCELLAtacdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            DimPlot(object = v$scData1, group.by = 'primary.predict', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
          })
        }
      })
      
      plotCellAnnoAtac <- reactive({
        if(is.null(v$atacData) || is.null(v$isCELLAtacsdone)){
          return(NULL)
        }else{
          withProgress(message="Generating Plot...", value=0, {
            DimPlot(object = v$atacData, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
          })
        }
      })
      
      output$cellanno_atac1.plot <- renderPlot({
        plotCellAnnoAtac()
      })
      
      output$download_cellanno_plot <- downloadHandler(
        filename = function(){"Cell annotation plot.png"}, 
        content = function(fname){
          ggsave(fname,plotCellAnnoAtac(), height = 7, width = 7)
        }
      )
      
      #output$cellanno_atac1.plot <- renderPlot({
      #  if(is.null(v$atacData) || is.null(v$isCELLAtacdone)){
      #    plotly_empty()
      #  }else{
      #    withProgress(message="Generating DEG Plot...", value=0, {
      #      DimPlot(object = v$atacData, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
      #    })
      #  }
      #})
      
      observeEvent(input$load_eg1_scRNA, {
        withProgress(message="Loading example reference scRNA-seq data...", value=0.5, {
          tpmFiles_ex_scRNA1 <- read.table('scRNA/pbmc.txt', header = T, row.names = 1, check.names = F)
          tpmFiles_ex_scRNA1 <- CreateSeuratObject(tpmFiles_ex_scRNA1)
          tpmFiles_ex_scRNA1[["percent.mt"]] <- PercentageFeatureSet(tpmFiles_ex_scRNA1, pattern = "^MT-")
          tpmFiles_ex_scRNA1 <- NormalizeData(tpmFiles_ex_scRNA1)
          tpmFiles_ex_scRNA1 <- FindVariableFeatures(tpmFiles_ex_scRNA1, selection.method = "vst", nfeatures = 2000)
          all.genes <- rownames(tpmFiles_ex_scRNA1)
          tpmFiles_ex_scRNA1 <- ScaleData(tpmFiles_ex_scRNA1, features = all.genes)
          tpmFiles_ex_scRNA1 <- RunPCA(tpmFiles_ex_scRNA1, features = VariableFeatures(object = tpmFiles_ex_scRNA1))
          tpmFiles_ex_scRNA1 <- FindNeighbors(tpmFiles_ex_scRNA1, dims = 1:30)
          tpmFiles_ex_scRNA1 <- FindClusters(tpmFiles_ex_scRNA1, resolution = 1)
          tpmFiles_ex_scRNA1 <- RunUMAP(tpmFiles_ex_scRNA1, dims = 1:30)
          sc <- reticulate::import("scanpy", convert = FALSE)
          ct <- reticulate::import("celltypist", convert = FALSE)
          sceasy::convertFormat(tpmFiles_ex_scRNA1, from = "seurat", to = "anndata", outFile = 'ct_scrna.h5ad')
          v$adata = sc$read_h5ad('ct_scrna.h5ad')
          v$res = ct$annotate(filename = 'ct_scrna.h5ad', model = input$celltypistatlas_atac, majority_voting=T)
          print(v$res)
          v$adata = v$res$to_adata()
          print("fff")
          v$adata$obs$to_csv('celltypist_predict.csv')
          v$meta_1 <- read.csv('celltypist_predict.csv', header = T, row.names = 1)
          tpmFiles_ex_scRNA1 <- AddMetaData(tpmFiles_ex_scRNA1, metadata = v$meta_1)
          tpmFiles_ex_scRNA1$primary.predict <- tpmFiles_ex_scRNA1$majority_voting
          tpmFiles_ex_scRNA1$secondary.predict <- tpmFiles_ex_scRNA1$predicted_labels
          print(tpmFiles_ex_scRNA1@meta.data)
        })    
        if (is.null(tpmFiles_ex_scRNA1)){
          v$scData1 <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0.5, {
            #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
            label1 <- "Example loaded"
            updateActionButton(inputId = "load_eg1_scRNA", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scData1 <- tpmFiles_ex_scRNA1
            v$isCELLAtac1done <- TRUE
        })
      }
    })
      
      output$cellanno1_atac.plot <- renderPlot({
        if(is.null(v$scData1) || is.null(v$isCELLAtac1done)){
          plotly_empty()
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            DimPlot(object = v$scData1, group.by = 'primary.predict', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
          })
        }
      })
      
      plotCellAnnoAtac1 <- reactive({
        if(is.null(v$atacData) || is.null(v$isCELLAtac1done)){
          plotly_empty()
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            DimPlot(object = v$atacData, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
          })
        }
      })
      
      output$cellanno1_atac1.plot <- renderPlotly({
        plotCellAnnoAtac1()
      })
      
      output$download_cellanno_plot <- downloadHandler(
        filename = function(){"Cell annotation plot.png"}, 
        content = function(fname){
          ggsave(fname,plotCellAnnoAtac1(), height = 7, width = 7)
        }
      )
          
        observe({
            if(input$load_user1_scRNA > 0){
              print('2')
              session$sendCustomMessage("myCallbackHandler4a", "2")
            }
        })
        
        observeEvent(input$load_eg_scRNA_a, {
          withProgress(message="Loading example reference scRNA-seq data...", value=0.5, {
            tpmFiles_ex_Intg <- read.table("integration/concatenated_expr_data.txt", header = T, row.names = 1, check.names = F, sep = '\t')
            #scH5 <- input$scH5
            annoFile_ex_Intg <- read.table("integration/metadata.txt", header = T, row.names = 1, sep = '\t')
            #tpmFiles_ex_Intg <- read.table('scRNA/pbmc.txt', header = T, row.names = 1, check.names = F)
            tpmFiles_ex_Intg <- CreateSeuratObject(counts = tpmFiles_ex_Intg, meta.data = annoFile_ex_Intg)
            tpmFiles_ex_Intg[["percent.mt"]] <- PercentageFeatureSet(tpmFiles_ex_Intg, pattern = "^MT-")
            tpmFiles_ex_Intg.list <- SplitObject(tpmFiles_ex_Intg, split.by = "batch")
            for (i in 1:length(tpmFiles_ex_Intg.list)) {
              tpmFiles_ex_Intg.list[[i]] <- NormalizeData(tpmFiles_ex_Intg.list[[i]], verbose = FALSE)
              tpmFiles_ex_Intg.list[[i]] <- FindVariableFeatures(tpmFiles_ex_Intg.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
            }
            print(tpmFiles_ex_Intg)
            print(tpmFiles_ex_Intg.list)
            features <- SelectIntegrationFeatures(object.list = tpmFiles_ex_Intg.list, nfeatures = 2000)
            print(features)
            tpmFiles_ex_Intg.anchors <- FindIntegrationAnchors(object.list = tpmFiles_ex_Intg.list, anchor.features = features)
            tpmFiles_ex_Intg <- IntegrateData(anchorset = tpmFiles_ex_Intg.anchors)
            print(tpmFiles_ex_Intg.anchors)
            print(tpmFiles_ex_Intg)
            DefaultAssay(tpmFiles_ex_Intg) <- "integrated"
            tpmFiles_ex_Intg <- ScaleData(tpmFiles_ex_Intg, verbose = FALSE)
            tpmFiles_ex_Intg <- RunPCA(tpmFiles_ex_Intg, verbose = FALSE)
            tpmFiles_ex_Intg <- FindNeighbors(tpmFiles_ex_Intg, reduction = "pca", dims = 1:30)
            tpmFiles_ex_Intg <- FindClusters(tpmFiles_ex_Intg, resolution = 1)
            tpmFiles_ex_Intg <- RunUMAP(tpmFiles_ex_Intg, reduction = "pca", dims = 1:30, spread = 1)
            print(tpmFiles_ex_Intg)
            ref = readRDS('ref.rds')
            tpmFiles_ex_Intg.data.average = AverageExpression(tpmFiles_ex_Intg)
            tpmFiles_ex_Intg.data.average = round(tpmFiles_ex_Intg.data.average$RNA, 2)
            if(input$cellatlas_atac_a == "all"){
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, ref)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "adipose"){
              adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
              adipose1 <- ref[,adipose]
              colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, adipose1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "adrenal_gland"){
              adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
              adrenal_gland1 <- ref[,adrenal_gland]
              colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, adrenal_gland1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "blood"){
              blood <- colnames(ref)[grepl("blood",colnames(ref))] 
              blood1 <- ref[,blood]
              colnames(blood1) <- gsub("--blood","",colnames(blood1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, blood1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "bone_marrow"){
              bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
              bone_marrow1 <- ref[,bone_marrow]
              colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, bone_marrow1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "brain"){
              brain <- colnames(ref)[grepl("brain",colnames(ref))] 
              brain1 <- ref[,brain]
              colnames(brain1) <- gsub("--brain","",colnames(brain1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, brain1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "breast"){
              breast <- colnames(ref)[grepl("breast",colnames(ref))] 
              breast1 <- ref[,breast]
              colnames(breast1) <- gsub("--breast","",colnames(breast1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, breast1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "breast_milk"){
              breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
              breast_milk1 <- ref[,breast_milk]
              colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, breast_milk1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "eye"){
              eye <- colnames(ref)[grepl("eye",colnames(ref))] 
              eye1 <- ref[,eye]
              colnames(eye1) <- gsub("--eye","",colnames(eye1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, eye1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "gut"){
              gut <- colnames(ref)[grepl("gut",colnames(ref))] 
              gut1 <- ref[,gut]
              colnames(gut1) <- gsub("--gut","",colnames(gut1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, gut1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "heart"){
              heart <- colnames(ref)[grepl("heart",colnames(ref))] 
              heart1 <- ref[,heart]
              colnames(heart1) <- gsub("--heart","",colnames(heart1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, heart1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "kidney"){
              kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
              kidney1 <- ref[,kidney]
              colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, kidney1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "liver"){
              liver <- colnames(ref)[grepl("liver",colnames(ref))] 
              liver1 <- ref[,liver]
              colnames(liver1) <- gsub("--liver","",colnames(liver1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, liver1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "lung"){
              lung <- colnames(ref)[grepl("lung",colnames(ref))] 
              lung1 <- ref[,lung]
              colnames(lung1) <- gsub("--lung","",colnames(lung1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, lung1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "pancreas"){
              pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
              pancreas1 <- ref[,pancreas]
              colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, pancreas1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "PDAC"){
              PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
              PDAC1 <- ref[,PDAC]
              colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, PDAC1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "skin"){
              skin <- colnames(ref)[grepl("skin",colnames(ref))] 
              skin1 <- ref[,skin]
              colnames(skin1) <- gsub("--skin","",colnames(skin1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, skin1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "testis"){
              testis <- colnames(ref)[grepl("testis",colnames(ref))] 
              testis1 <- ref[,testis]
              colnames(testis1) <- gsub("--testis","",colnames(testis1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, testis1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "thymus"){
              thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
              thymus1 <- ref[,thymus]
              colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, thymus1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
            if(input$cellatlas_atac_a == "tonsil"){
              tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
              tonsil1 <- ref[,tonsil]
              colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
              res = FastIntegration::CELLiD(tpmFiles_ex_Intg.data.average, tonsil1)
              print(res)
              tpmFiles_ex_Intg$primary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),1]
              tpmFiles_ex_Intg$secondary.predict = res[as.numeric(tpmFiles_ex_Intg$seurat_clusters),2]
              newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
              colnames(res) <- newheaders
              print(tpmFiles_ex_Intg@meta.data)
              output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
              write.table(res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
            }
          })
          if (is.null(tpmFiles_ex_Intg)){
            v$scData1 <- NULL
          }else{
            withProgress(message="Loading and Processing Data...", value=0.5, {
              #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
              label1 <- "Example loaded"
              updateActionButton(inputId = "load_eg_scRNA_a", label = label1)
              shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
              v$scData1 <- tpmFiles_ex_Intg           
            })
          }
        })
        
        observe({
          if(input$load_user_scRNA_a > 0){
            print('2')
            session$sendCustomMessage("myCallbackHandler_4", "2")
          }
        })
        
        output$cellanno_atac.plot_a <- renderPlot({
          if(is.null(v$scData1)){
            return(NULL)
          }else{
            withProgress(message="Generating DEG Plot...", value=0, {
              DimPlot(object = v$scData1, group.by = 'primary.predict', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
            })
          }
        })
        
        output$cellanno_atac1.plot_a <- renderPlot({
          if(is.null(v$atacData)){
            return(NULL)
          }else{
            withProgress(message="Generating DEG Plot...", value=0, {
              DimPlot(object = v$atacData, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
            })
          }
        })
        
        observeEvent(input$load_eg1_scRNA_a, {
          withProgress(message="Loading example reference scRNA-seq data...", value=0.5, {
            tpmFiles_ex_Intg <- read.table("integration/concatenated_expr_data.txt", header = T, row.names = 1, check.names = F, sep = '\t')
            #scH5 <- input$scH5
            annoFile_ex_Intg <- read.table("integration/metadata.txt", header = T, row.names = 1, sep = '\t')
            #tpmFiles_ex_Intg <- read.table('scRNA/pbmc.txt', header = T, row.names = 1, check.names = F)
            tpmFiles_ex_Intg <- CreateSeuratObject(counts = tpmFiles_ex_Intg, meta.data = annoFile_ex_Intg)
            tpmFiles_ex_Intg[["percent.mt"]] <- PercentageFeatureSet(tpmFiles_ex_Intg, pattern = "^MT-")
            tpmFiles_ex_Intg.list <- SplitObject(tpmFiles_ex_Intg, split.by = "batch")
            for (i in 1:length(tpmFiles_ex_Intg.list)) {
              tpmFiles_ex_Intg.list[[i]] <- NormalizeData(tpmFiles_ex_Intg.list[[i]], verbose = FALSE)
              tpmFiles_ex_Intg.list[[i]] <- FindVariableFeatures(tpmFiles_ex_Intg.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
            }
            print(tpmFiles_ex_Intg)
            print(tpmFiles_ex_Intg.list)
            features <- SelectIntegrationFeatures(object.list = tpmFiles_ex_Intg.list, nfeatures = 2000)
            print(features)
            tpmFiles_ex_Intg.anchors <- FindIntegrationAnchors(object.list = tpmFiles_ex_Intg.list, anchor.features = features)
            tpmFiles_ex_Intg <- IntegrateData(anchorset = tpmFiles_ex_Intg.anchors)
            print(tpmFiles_ex_Intg.anchors)
            print(tpmFiles_ex_Intg)
            DefaultAssay(tpmFiles_ex_Intg) <- "integrated"
            tpmFiles_ex_Intg <- ScaleData(tpmFiles_ex_Intg, verbose = FALSE)
            tpmFiles_ex_Intg <- RunPCA(tpmFiles_ex_Intg, verbose = FALSE)
            tpmFiles_ex_Intg <- FindNeighbors(tpmFiles_ex_Intg, reduction = "pca", dims = 1:30)
            tpmFiles_ex_Intg <- FindClusters(tpmFiles_ex_Intg, resolution = 1)
            tpmFiles_ex_Intg <- RunUMAP(tpmFiles_ex_Intg, reduction = "pca", dims = 1:30, spread = 1)
            print(tpmFiles_ex_Intg)
            sc <- reticulate::import("scanpy", convert = FALSE)
            ct <- reticulate::import("celltypist", convert = FALSE)
            sceasy::convertFormat(tpmFiles_ex_Intg, from = "seurat", to = "anndata", outFile = 'ct_intg.h5ad')
            v$adata = sc$read_h5ad('ct_intg.h5ad')
            v$res1 = ct$annotate(filename = 'ct_intg.h5ad', model = input$celltypistatlas_atac_a, majority_voting=T)
            print(v$res1)
            v$adata = v$res1$to_adata()
            print("fff")
            v$adata$obs$to_csv('celltypist_predict.csv')
            v$meta_1a <- read.csv('celltypist_predict.csv', header = T, row.names = 1)
            tpmFiles_ex_Intg <- AddMetaData(tpmFiles_ex_Intg, metadata = v$meta_1a)
            tpmFiles_ex_Intg$primary.predict <- tpmFiles_ex_Intg$majority_voting
            tpmFiles_ex_Intg$secondary.predict <- tpmFiles_ex_Intg$predicted_labels
            print(tpmFiles_ex_Intg@meta.data)
          })    
          if (is.null(tpmFiles_ex_Intg)){
            v$scData1 <- NULL
          }else{
            withProgress(message="Loading and Processing Data...", value=0.5, {
              #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
              label1 <- "Example loaded"
              updateActionButton(inputId = "load_eg1_scRNA_a", label = label1)
              shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
              v$scData1 <- tpmFiles_ex_Intg
            })
          }
        })
        
        observe({
          if(input$load_user1_scRNA_a > 0){
            print('2')
            session$sendCustomMessage("myCallbackHandler_4a", "2")
          }
        })
        
        output$cellanno1_atac.plot_a <- renderPlot({
          if(is.null(v$scData1)){
            return(NULL)
          }else{
            withProgress(message="Generating DEG Plot...", value=0, {
              DimPlot(object = v$scData1, group.by = 'primary.predict', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
            })
          }
        })
        
        output$cellanno1_atac1.plot_a <- renderPlot({
          if(is.null(v$atacData)){
            return(NULL)
          }else{
            withProgress(message="Generating DEG Plot...", value=0, {
              DimPlot(object = v$atacData, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
            })
          }
        })
          
      observeEvent(input$doDeg_atac, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Finding DE Peaks...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            print(v$atacData)
            atac.markers <- FindAllMarkers(v$atacData, only.pos = FALSE, min.pct = input$min_pct_atac, logfc.threshold = input$logfc_atac, test.use = input$test.use_atac)
            v$atac.markers <- atac.markers
          })
        }
      })
      
      output$Deg_atac.table <- DT::renderDataTable(
        v$atac.markers, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
      
      output$Deg_atac.plot <- renderPlotly({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            v$atac.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
            DoHeatmap(v$atacData, features = top10$gene, size = 4, slot = "scale.data",angle = 30) + theme(axis.text.y = element_text(size = 4)) + NoLegend()
          })
        }
      })
      
      observeEvent(input$Vis_atac, {
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Visualizing...", value=0, {
            v$isVisATACdone <- TRUE
          })
        }
      })
      
      output$deg.atac.select <- renderUI({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          selectInput("deg.atac", label = "Peaks to visualise",
                      choices = rownames(v$atac.markers))
        }
      })
      
      plotViolinAtac <- reactive({
        if(is.null(v$atac.markers) || is.null(v$isVisATACdone)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            VlnPlot(v$atacData, input$deg.atac)
          })
        }
      })
      
      output$Deg_atac1.plot <- renderPlotly({
                                            plotViolinAtac()
                                          })
      
      output$download_violn_atac <- downloadHandler(
        filename = function(){"Violin plot for ATAC.png"}, 
        content = function(fname){
          ggsave(fname,plotViolinAtac(), height = 7, width = 7)
        }
      )
      
      plotFeatureAtac <- reactive({
        if(is.null(v$atac.markers) || is.null(v$isVisATACdone)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            FeaturePlot(v$atacData, input$deg.atac)
          })
        }
      })
      
      output$Deg_atac2.plot <- renderPlotly({
                                        plotFeatureAtac()
                                      })
        
      output$download_feature_atac <- downloadHandler(
        filename = function(){"Feature plot for ATAC.png"}, 
        content = function(fname){
          ggsave(fname,plotFeatureAtac(), height = 7, width = 7)
        }
      )
      
      plotRidgeAtac <- reactive({
        if(is.null(v$atac.markers) || is.null(v$isVisATACdone)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            RidgePlot(v$atacData, features = input$deg.atac)
          })
        }
      })
      
      output$Deg_atac3.plot <- renderPlot({
                                    plotRidgeAtac()
                                    })
      
      output$download_ridge_atac <- downloadHandler(
        filename = function(){"Ridge plot for ATAC.png"}, 
        content = function(fname){
          ggsave(fname,plotRidgeAtac(), height = 7, width = 7)
        }
      )
      
      output$peak_gene.atac.select <- renderUI({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          selectInput("gene_peak", label = "Genes to visualise",
                      choices = rownames(v$atacData@assays$RNA@counts), multiple = T, selected = rownames(v$atacData@assays$RNA@counts)[1])
        }
      })
      
      plotCoverageAtac <- reactive({
        if(is.null(v$atacData)){
          plotly_empty()
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            CoveragePlot(
              object = v$atacData,
              region = input$gene_peak,
              features = input$gene_peak,
              extend.upstream = 500,
              extend.downstream = 10000
            )
          })
        }
      })
      
      output$coverage1.plot <- renderPlot({
                                           plotCoverageAtac()
                                          })
      
      observeEvent(input$link_peak_genes, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Linking peaks to genes...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            # first compute the GC content for each peak
            v$atacData <- RegionStats(v$atacData, genome = BSgenome.Hsapiens.UCSC.hg19)
            
            # link peaks to genes
            v$atacData <- LinkPeaks(
              object = v$atacData,
              peak.assay = "peaks",
              expression.assay = "RNA",
              genes.use = input$gene_peak)
            
            v$peak <- as.data.frame(v$atacData@assays$peaks@links)
            v$isLinkPeakdone <- TRUE
          })
        }
      })
      
      output$download_coverage_atac <- downloadHandler(
        filename = function(){"Coverage plot for ATAC.png"}, 
        content = function(fname){
          ggsave(fname,plotCoverageAtac(), height = 7, width = 7)
        }
      )
      
      output$coverage2.plot <- renderPlot({
        if(is.null(v$atacData) || is.null(v$isLinkPeakdone)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            CoveragePlot(
              object = v$atacData,
              region = input$gene_peak,
              features = input$gene_peak,
              extend.upstream = 500,
              extend.downstream = 10000
            )
          })
        }
      })
      
      output$peaks.table <- DT::renderDataTable(
        v$peak, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
      
      ##---------------GSEA of scATAC-seq-------------------##
      
      output$gsea_atac.ct1.select <- renderUI({
        if(is.null(v$atac.markers) || is.null(v$isCELLAtacsdone)){
          return(NULL)
        }else{
          selectInput("gsea_atac.ct1", label = "Celltype1",
                      choices = as.vector(v$atacData$predicted.id))
        }
      })
      
      output$gsea_atac.ct2.select <- renderUI({
        if(is.null(v$atac.markers) || is.null(v$isCELLAtacsdone)){
          return(NULL)
        }else{
          selectInput("gsea_atac.ct2", label = "Celltype2",
                      choices = as.vector(v$atacData$predicted.id))
        }
      })
      
      observeEvent(input$gsea_atac, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
        withProgress(message="Generating gene set enrichment analysis...", value=0, {
          if(input$species_gsea_atac == "Homo sapiens" & input$category_gsea_atac == "H"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_hs_go_atac <- msigdbr(species = "Homo sapiens", category = "H")
            print(v$msigdbr_hs_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_hs_go_atac$gene_symbol, f = v$msigdbr_hs_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Mus musculus" & input$category_gsea_atac == "H"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_mm_go_atac <- msigdbr(species = "Mus musculus", category = "H")
            print(v$msigdbr_mm_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_mm_go_atac$gene_symbol, f = v$msigdbr_mm_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Homo sapiens" & input$category_gsea_atac == "C2"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_hs_go_atac <- msigdbr(species = "Homo sapiens", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
            print(v$msigdbr_hs_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_hs_go_atac$gene_symbol, f = v$msigdbr_hs_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Mus musculus" & input$category_gsea_atac == "C2"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_mm_go_atac <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
            print(v$msigdbr_mm_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_mm_go_atac$gene_symbol, f = v$msigdbr_mm_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Homo sapiens" & input$category_gsea_atac == "C5"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_hs_go_atac <- msigdbr(species = "Homo sapiens", category = "C5")
            print(v$msigdbr_hs_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_hs_go_atac$gene_symbol, f = v$msigdbr_hs_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Mus musculus" & input$category_gsea_atac == "C5"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_mm_go_atac <- msigdbr(species = "Mus musculus", category = "C5")
            print(v$msigdbr_mm_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_mm_go_atac$gene_symbol, f = v$msigdbr_mm_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Homo sapiens" & input$category_gsea_atac == "C7"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_hs_go_atac <- msigdbr(species = "Homo sapiens", category = "C7")
            print(v$msigdbr_hs_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_hs_go_atac$gene_symbol, f = v$msigdbr_hs_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Mus musculus" & input$category_gsea_atac == "C7"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_mm_go_atac <- msigdbr(species = "Mus musculus", category = "C7")
            print(v$msigdbr_mm_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_mm_go_atac$gene_symbol, f = v$msigdbr_mm_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Homo sapiens" & input$category_gsea_atac == "C8"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_hs_go_atac <- msigdbr(species = "Homo sapiens", category = "C8")
            print(v$msigdbr_hs_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_hs_go_atac$gene_symbol, f = v$msigdbr_hs_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          if(input$species_gsea_atac == "Mus musculus" & input$category_gsea_atac == "C8"){
            v$atacData$predicted.id -> Idents(v$atacData)
            v$msigdbr_mm_go_atac <- msigdbr(species = "Mus musculus", category = "C8")
            print(v$msigdbr_mm_go_atac)
            v$pathways_atac <- split(x = v$msigdbr_mm_go_atac$gene_symbol, f = v$msigdbr_mm_go_atac$gs_name)
            print(v$pathways_atac)
            v$markers_atac <- FindMarkers(v$atacData, ident.1 = input$gsea_atac.ct1, ident.2 = input$gsea_atac.ct2, min.pct = input$min_pct_atac1, logfc_atac.threshold = input$logfc_atac1, test.use_atac = input$test.use_atac1, assay = "RNA", min.cells.feature = 0, min.cells.group = 0)
            v$markers_atac  <- v$markers_atac %>% arrange(desc(avg_log2FC))
            print(v$markers_atac)
            v$markers_atac.log2FC <- v$markers_atac$avg_log2FC
            names(v$markers_atac.log2FC) <- row.names(v$markers_atac)
            v$markers_atac.log2FC <- sort(na.omit(v$markers_atac.log2FC), decreasing = TRUE)
            print(v$markers_atac.log2FC)
            v$fgsea_atacRes_atac <- fgsea(pathways = v$pathways_atac, stats = v$markers_atac.log2FC, eps = 0.0, minSize  = 5, maxSize  = 500) %>% arrange((padj))
            print(v$fgsea_atacRes_atac)
            v$topPathways_atacUp_atac <- v$fgsea_atacRes_atac[ES > 0][head(order(pval), n=10), pathway]
            v$topPathways_atacDown_atac <- v$fgsea_atacRes_atac[ES < 0][head(order(pval), n=10), pathway]
            v$topPathways_atac <- c(v$topPathways_atacUp_atac, rev(v$topPathways_atacDown_atac))
          }
          output$gsea_atac.done <- renderText(paste0("Gene set enrichment done!"))
          v$isGSEAatacdone <- TRUE
          })
        }
      })
      
      output$gsea_atac.select <- renderUI({
        if(is.null(v$pathways_atac)){
          return(NULL)
        }else{
          selectInput("gsea_atac.pathway", label = "Gene set to visualise",
                      choices = names(v$pathways_atac))
        }
      })
      
      output$gsea_atac_plot <- renderPlot({
        if(is.null(v$pathways_atac) || is.null(v$isGSEAatacdone)){
          return(NULL)
        }else{
          withProgress(message="Generating gsea_atac plot...", value=0, {
            plotEnrichment(v$pathways_atac[[input$gsea_atac.pathway]], v$markers_atac.log2FC) + labs(title=input$gsea_atac.pathway)
          })
        }
      })
      
      output$gsea_atac_plot1 <- renderPlot({
        if(is.null(v$pathways_atac) || is.null(v$isGSEAatacdone)){
          return(NULL)
        }else{
          withProgress(message="Generating gsea_atac plot...", value=0, {
            plotGseaTable(v$pathways_atac[v$topPathways_atac], v$markers_atac.log2FC, v$fgsea_atacRes_atac, gseaParam=0.2)
          })
        }
      })
      
      output$gsea_atac.table <- DT::renderDataTable(
        v$fgsea_atacRes_atac, server = FALSE, options = list(scrollX = TRUE, scrollY = "400px"))
      
      output$download_gsea_atac.table <- downloadHandler(
        filename = function(){"gsea_atac Results.csv"}, 
        content = function(fname){
          withProgress(message="Downloading gsea_atac Results...", value=0, {
            fwrite(v$fgsea_atacRes_atac, fname)
          })
        }
      )
      
      output$great_atac.ct1.select <- renderUI({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          selectInput("great_atac.ct1", label = "Celltype1",
                      choices = as.vector(v$atacData$predicted.id))
        }
      })
      
      output$great_atac.ct2.select <- renderUI({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          selectInput("great_atac.ct2", label = "Celltype2",
                      choices = as.vector(v$atacData$predicted.id))
        }
      })
      
      observeEvent(input$gsea_atac1, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Generating gene set enrichment analysis...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            v$atacData$predicted.id -> Idents(v$atacData)
              v$da_peaks <- FindMarkers(object = v$atacData, ident.1 = input$great_atac.ct1, ident.2 = input$great_atac.ct2, logfc.threshold = input$logfc_great, min.pct = input$min_pct_great, test.use = input$test.use_atac2, min.cells.feature = 0, min.cells.group = 0)
              v$da_peaks1 <- as.data.frame(rownames(v$da_peaks))
              names(v$da_peaks1)[1] <- "position"
              v$da_peaks1[c('chr', 'start', 'end')] <- str_split_fixed(v$da_peaks1$position, '-', 3)
              v$da_peaks1[2:4]-> v$da_peaks2
              makeGRangesFromDataFrame(v$da_peaks2) -> v$da_peaks3
              granges(v$da_peaks3) -> v$da_peaks4
              v$res_atac = great(v$da_peaks4, input$gene_set_atac, input$tss_atac)
              print(v$res_atac)
              v$tb = getEnrichmentTable(v$res_atac)
              v$tb1 <- v$tb[order(v$tb$fold_enrichment, decreasing = T),][1:50,]
              print(v$tb)
              print(v$tb1)
              output$gsea_atac1.done <- renderText(paste0("Gene set enrichment done!"))
              v$isGSEAatac1done <- TRUE
          })
        }
      })
      
      output$great.select <- renderUI({
        if(is.null(v$res_atac)){
          return(NULL)
        }else{
          selectInput("great.pathway", label = "Pathway to visualise",
                      choices = v$tb$description)
        }
      })
      
      output$great.plot <- renderPlot({
        if(is.null(v$atacData) || is.null(v$isGSEAatac1done)){
          plotly_empty()
        }else{
          withProgress(message="Generating Gene Set Enrichment Plot...", value=0, {
            plotRegionGeneAssociations(v$res_atac, term_id = input$great.pathway)
          })
        }
      })
      
      output$great1.plot <- renderPlotly({
        if(is.null(v$atacData) || is.null(v$isGSEAatac1done)){
          plotly_empty()
        }else{
          withProgress(message="Generating Gene Set Enrichment Plot...", value=0, {
            ggplot(v$tb1, aes(x = fold_enrichment, y = description)) + geom_point(aes(color=p_value_hyper,size=gene_set_size)) + scale_color_gradientn(colours = rainbow(5)) + labs(x='Fold Enrichment', y=NULL, color='P-value',size='Gene number') + theme( axis.title = element_text(face='bold'), axis.text = element_text(face='bold')) + theme_bw()
          })
        }
      })
      
      output$coverage.atac.select <- renderUI({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          selectInput("deg1.atac", label = "Peaks to visualise",
                      choices = rownames(v$atac.markers))
        }
      })
      
      output$coverage.atac_feature.select <- renderUI({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          selectInput("atac_feature", label = "Genes to visualise",
                      choices = rownames(v$atacData@assays$RNA@counts), multiple = T)
        }
      })
      
      output$coverage.plot <- renderPlot({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            CoveragePlot(
              object = v$atacData,
              region = input$deg1.atac,
              features = input$atac_feature,
              extend.upstream = 40000,
              extend.downstream = 20000
            )
          })
        }
      })
      
      observeEvent(input$doDeg_motif, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Finding DEGs...", value=0, {
            pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
            v$atacData <- AddMotifs(object = v$atacData, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = pfm)
            v$da_peaks <- FindAllMarkers(object = v$atacData, only.pos = F, test.use = input$test.use_motif, min.pct = input$min_pct_motif, logfc.threshold = input$logfc_motif)
            v$top.da.peak <- rownames(v$da_peaks[v$da_peaks$p_val < 0.005, ])
            enriched.motifs <- FindMotifs(object = v$atacData, features = v$top.da.peak)
            v$enriched.motifs <- enriched.motifs
            shinyalert("DEG Analysis done", "DEG Analysis done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
      output$motif.select <- renderUI({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          selectInput("motif.atac", label = "Peaks to visualise",
                      choices = rownames(v$enriched.motifs))
        }
      })
      
      output$motif.plot <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            MotifPlot(object = v$atacData, motifs = input$motif.atac)
          })
        }
      })
      
      observeEvent(input$calc_motif_activity, {
        tpmFiles_atac <- input$tpmFiles_atac
        tpmFiles_atac <- v$atacData
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
          shinyalert("Please upload your data", "Please upload your data", type = "warning", imageWidth = 10, imageHeight = 10)
        }else{
          withProgress(message="Finding DEGs...", value=0, {
            v$atacData <- RunChromVAR(object = v$atacData, genome = BSgenome.Hsapiens.UCSC.hg19)
            DefaultAssay(v$atacData) <- 'chromvar'
            differential.activity <- FindAllMarkers(object = v$atacData, only.pos = F, mean.fxn = rowMeans, fc.name = "avg_diff")
            v$differential.activity <- differential.activity
            v$isMotifActivitydone <- TRUE
            shinyalert("Motif analysis done", "Motif analysis done, please perform motif ", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
      output$motif_feature.select <- renderUI({
        if(is.null(v$atacData) || is.null(v$isMotifActivitydone)){
          return(NULL)
        }else{
          selectInput("motif.feature", label = "Peaks to visualise",
                      choices = rownames(v$atacData))
        }
      })
      
      output$motif_feature.plot <- renderPlotly({
        if(is.null(v$atacData) || is.null(v$isMotifActivitydone)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            FeaturePlot(object = v$atacData,
                        features = input$motif.feature,
                        min.cutoff = 'q10',
                        max.cutoff = 'q90',
                        pt.size = 0.1)
          })
        }
      })
      
      output$motif1.select <- renderUI({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          selectInput("motif1.atac", label = "Peaks to visualise",
                      choices = rownames(v$differential.activity))
        }
      })
      
      output$motif1.plot <- renderPlot({
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            MotifPlot(object = v$atacData, motifs = input$motif1.atac, assay = 'peaks')
          })
        }
      })
      
  ##---------------Summary tab
  
  ##------Clean up when ending session----
  session$onSessionEnded(function(){
    prePlot()
  })
})
