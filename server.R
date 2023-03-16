source('ui.R')
source('cellphonedb.R')
source('plot_heatmaps.R')

## max data size
#options(shiny.maxRequestSize = 1024^10)
options(shiny.maxRequestSize = 1024*1024*100*100)
options(shiny.launch.browser = T)
options(bitmapType = 'cairo')
options(timeout=1000000)

shinyServer(function(input, output, session) {
  
  #use_python('/home/george/miniconda3/bin/python', required = NULL)
  
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
    list(src = "www/overview.png",
         width = 500,
         height = 400)
  }, deleteFile = FALSE)
  
  output$demo_image1 <- renderImage({
    list(src = "www/overview2.png",
         width = 350,
         height = 350)
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
      })
      if (is.null(tpmFiles_demo)){
        v$scData <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0.5, {
          print(tpmFiles_demo)
          #print(tpmFiles$name)
          print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample_tpm_human", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
          v$scData <- tpmFiles_demo
        })
      }
  })
  
  observeEvent(input$loadexample_scH5, {
    withProgress(message="Loading example data...", value=0.5, {
      tpmFiles_demo <- Read10X_h5("scRNA/filtered_feature_bc_matrix_human.h5")
    })
    if (is.null(tpmFiles_demo)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0.5, {
        #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample_tpm_human", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        v$scData <- tpmFiles_demo
      })
    }
  })
  
  observeEvent(input$loadexample_rds, {
    withProgress(message="Loading example data...", value=0.5, {
      tpmFiles_demo <- readRDS("scRNA/ERX2756718.rds")
      tpmFiles_demo@assays$RNA@counts = expm1(tpmFiles_demo@assays$RNA@data)
    })
    if (is.null(tpmFiles_demo)){
      v$scData <- NULL
    }else{
      withProgress(message="Loading and Processing Data...", value=0.5, {
        #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample_tpm_human", label = label1)
        shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        v$scData <- tpmFiles_demo@assays$RNA@counts
      })
    }
  })
  
  observeEvent(input$loadButton, {
    if(input$scInput == "Raw Counts Matrix"){
      tpmFiles <- input$tpmFiles
      {
        withProgress(message="Loading and Processing Data...", value=0, {
          print(tpmFiles$datapath)
          print(tpmFiles$name)
          print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
          exp.data <- read.table(tpmFiles$datapath,
                                 sep="\t", header=TRUE, row.names=1, stringsAsFactors = F, check.names = F)
          
          v$scData <- exp.data
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
          v$scData <- exp.data
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
          exp.data@assays$RNA@counts = expm1(exp.data@assays$RNA@data)
          v$scData <- exp.data@assays$RNA@counts
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
    
  observeEvent(input$create_seurat, {
    withProgress(message="Loading and Processing Data...", value=0, {
      print(v$scData)
      sObj <- CreateSeuratObject(v$scData,
                                 project = input$projName,
                                 min.genes = input$min.genes,
                                 min.cells = input$min.cells)
      
      sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
      v$scData1 <- sObj
      print(v$scData1)
      shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
     })
    }
  )
  
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
        return(as.matrix(v$scData[,1:20]))
    }
  })
  
  plotInput <- reactive({
    p <- VlnPlot(v$scData1, "nFeature_RNA") + NoLegend()
  })
  
  output$nFeature_RNAPlot <- renderPlot({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      print(plotInput())
    }
  })
  
  output$download_nFeature_RNA <- downloadHandler(
    filename = function(){"nFeature_RNAPlot.png"}, 
    content = function(fname){
      ggsave(fname,plotInput())
    }
  )
  
  plotInput1 <- reactive({
    p <- VlnPlot(v$scData1, "percent.mt") + NoLegend()
  })
  
  output$mitoPlot <- renderPlot({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      print(plotInput1())
    }
  })
  
  output$download_mito <- downloadHandler(
    filename = function(){"mito_RNAPlot.png"}, 
    content = function(fname){
      ggsave(fname,plotInput1())
    }
  )
  
  plotInput2 <- reactive({
    p <- VlnPlot(v$scData1, "nCount_RNA") + NoLegend()
  })
  
  output$nCount_RNAPlot <- renderPlot({
    if(is.null(v$scData1)){
      plotly_empty()
    }else{
      print(plotInput2())
    }
  })
  
  output$download_nCount_RNA <- downloadHandler(
    filename = function(){"nCount_RNAPlot.png"}, 
    content = function(fname){
      ggsave(fname,plotInput2())
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
        }, height = 800, width = 850)
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
        }, height = 800, width = 850)
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
  })
  
  output$name <- renderPrint({
    s <- event_data("plotly_selected")
    c(s[["key"]], class(s[["key"]]))
  })
  
  
  
  ##---------------PCA of scRNA-seq-------------------
  # PCA plot
  observeEvent(input$doPCA, {
    withProgress(message = "Scaling Data...", value = 0,{
      incProgress(0.5, message = "Running PCA...")
      if(input$assays1 == "LogNormalization"){
      v$scData1 <- RunPCA(v$scData1, features = rownames(v$scData1), assay = "RNA")
      print(v$scData1[["pca"]], dims = 1:5, nfeatures = 5)
      }
      else if(input$assays1 == "SCTransform"){
        v$scData1 <- RunPCA(v$scData1, features = rownames(v$scData1), assay = "SCT")
        print(v$scData1[["pca"]], dims = 1:5, nfeatures = 5)
      }
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
    })
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
    p <- VizDimLoadings(v$scData1, dims = as.numeric(input$select.pc))
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
    p <- DimHeatmap(v$scData1, dims = as.numeric(input$select.pc))
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
  }, options = list(scrollX = TRUE))
  
  output$download_PCTable <- downloadHandler(
    
    filename = function(){"PC_table.csv"}, 
    content = function(fname){
      withProgress(message="Downloading PC Table...", value=0, {
        write.csv(v$pcGenes, fname)
      })
    }
  )
  
  plotElbow <- reactive({
    p <- ElbowPlot(v$scData1, ndims = 50)
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
  
  output$clustUI <- renderUI({
    if(is.null(v$isPCAdone)){
      return(NULL)
    }else{
      tagList(
        fluidRow(
          column(6,
                 numericInput("clus.res",
                              label = "Cluster Resolution",
                              value = 0.6,
                              min = 0.1,
                              step = 0.1)
          ),
          br(),
          selectInput("dim.used",
                      label = "First n PCs",
                      choices = c(2:50)
          ),
          column(6,
                 actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                 textOutput("cluster.done")
          )
        )
      )
    }
  })
  
  observeEvent(input$findCluster, {
    withProgress(message = "Finding clusters...", value = 0.3, {
      if(input$assays1 == "LogNormalization"){
      DefaultAssay(v$scData1) <- "RNA"
      v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, assay = "RNA", nn.method = "rann")
      v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
      output$cluster.done <- renderText(paste0("Clustering done!"))
      }
      else if(input$assays1 == "SCTransform"){
        DefaultAssay(v$scData1) <- "SCT"
        v$scData1 <- FindNeighbors(v$scData1, dims = 1:input$dim.used, assay = "SCT", nn.method = "rann")
        v$scData1 <- FindClusters(v$scData1, resolution = input$clus.res)
        output$cluster.done <- renderText(paste0("Clustering done!"))
      }
      v$isClusterdone <- TRUE
      shinyalert("Clustering performed", "Clustering performed, please perform UMAP", type = "success", imageWidth = 10, imageHeight = 10)
    })
  })
  
  plotCluster <- reactive({
    p <- DimPlot(v$scData1, reduction = "pca", label = T)
  })
  
  output$Cluster2DPlot_1 <- renderPlotly({
    print(plotCluster())
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
  
  ##---------------UMAP of scRNA-seq-------------------
  
  observeEvent(input$doUmap, {
    withProgress(message = "Running UMAP...", value = 0.3, {
      if(input$assays1 == "LogNormalization"){
      v$scData1 <- RunUMAP(v$scData1, dims = 1:input$dim.used, assay = "RNA", spread = 1)
      output$Umap.done <- renderText(paste0("UMAP done!"))
      }
      else if(input$assays1 == "SCTransform"){
        v$scData1 <- RunUMAP(v$scData1, dims = 1:input$dim.used, assay = "SCT", spread = 1)
        output$Umap.done <- renderText(paste0("UMAP done!"))
      }
      v$isUMAPdone <- TRUE
      shinyalert("UMAP performed", "UMAP performed, please perform tSNE", type = "success", imageWidth = 10, imageHeight = 10)
    })
  })
  
  plotUMAP <- reactive({
    if(is.null(v$scData1) || is.null(v$isUMAPdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    p <- DimPlot(v$scData1, reduction = "umap", label = T, label.size = 4)
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
    withProgress(message = "Running tSNE...", value = 0.3, {
      if(input$assays1 == "LogNormalization"){
      v$scData1 <- RunTSNE(v$scData1, dims = 1:input$dim.used, perplexity = input$perplexity, assay = "RNA")
      output$Tsne.done <- renderText(paste0("TSNE done!"))
      }
      else if(input$assays1 == "SCTransform"){
        v$scData1 <- RunTSNE(v$scData1, dims = 1:input$dim.used, perplexity = input$perplexity, assay = "SCT")
        output$Tsne.done <- renderText(paste0("TSNE done!"))
      }
      v$isTSNEdone <- TRUE
      shinyalert("tSNE performed", "tSNE performed, please perform celltype annotation", type = "success", imageWidth = 10, imageHeight = 10)
    })
  })
  
  output$Tsne_2d_plot_1 <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isTSNEdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating TSNE 2D Plot...", value=0, {
        DimPlot(v$scData1, reduction = "tsne", label = T)
      })
    }
  })
  
  plotTsne <- reactive({
    if(is.null(v$scData1) || is.null(v$isTSNEdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "tsne", label = T)
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
    withProgress(message = "Running CELLiD...", value = 0.3, {
      ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
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
  })
  
  output$Umap_cellid <- renderPlotly({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP from CELLiD...", value=0, {
        DimPlot(v$scData1, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3)
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
    v$res, options = list(scrollX = TRUE, scrollY = "400px"))
  
  plotCELLiD <- reactive({
    if(is.null(v$scData1) || is.null(v$isCELLiDdone)){
      plotly_empty()
    }else{
      withProgress(message="Generating UMAP 2D Plot...", value=0, {
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3)
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
        p <- DimPlot(v$scData1, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
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
    
    filename = function(){"CELLiD predictions.csv"}, 
    content = function(fname){
      withProgress(message="Downloading CELLiD predictions...", value=0, {
        write.csv(v$res, fname)
      })
    }
  )
  
  ##---------------Cell-cell similarity of scRNA-seq-------------------
  
  observeEvent(input$cell_cell, {
    withProgress(message="Generating heatmap...", value=0, {
      if(input$cell1 == "primary.predict"){
      v$scData1$primary.predict -> Idents(v$scData1)
      v$scData1.rna.data.average1 = AverageExpression(v$scData1)
      v$scData1.rna.data.average1 = data.frame(v$scData1.rna.data.average1$RNA)
      print(v$scData1.rna.data.average1)
      v$cor <- cor(v$scData1.rna.data.average1, method = "pearson")
      print(v$cor)
      output$CELL.done <- renderText(paste0("Cell-cell similarity done!"))
      v$isCELLdone <- TRUE
      }
      if(input$cell1 == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        v$scData1.rna.data.average1 = AverageExpression(v$scData1)
        v$scData1.rna.data.average1 = data.frame(v$scData1.rna.data.average1$RNA)
        print(v$scData1.rna.data.average1)
        v$cor <- cor(v$scData1.rna.data.average1, method = "pearson")
        rownames(v$cor) <- substr(rownames(v$cor),2,nchar(rownames(v$cor)))
        colnames(v$cor) <- substr(colnames(v$cor),2,nchar(colnames(v$cor)))
        print(v$cor)
        output$CELL.done <- renderText(paste0("Cell-cell similarity done!"))
        v$isCELLdone <- TRUE
      }
    })
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
    v$cor, options = list(scrollX = TRUE, scrollY = "400px"))
  
  output$download_cor.table <- downloadHandler(
    filename = function(){"Celltype_similarity.csv"}, 
    content = function(fname){
      withProgress(message="Downloading celltype_similarity...", value=0, {
        write.csv(v$cor, fname)
      })
    }
  )
  
  ##---------------DEGs of scRNA-seq-------------------
  
  observeEvent(input$doDeg, {
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      withProgress(message="Finding DEGs...", value=0, {
        if(input$deg1 == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        if (input$assays1 == "LogNormalization"){ 
        ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "RNA", test.use = input$test.use)
        v$ips.markers <- ips.markers
        }
        if (input$assays1 == "SCTransform"){ 
        ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "SCT", test.use = input$test.use)
        v$ips.markers <- ips.markers
        }
        shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        }
        if(input$deg1 == "primary.predict"){
          v$scData1$primary.predict -> Idents(v$scData1)
          if (input$assays1 == "LogNormalization"){ 
          ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "RNA", test.use = input$test.use)
          v$ips.markers <- ips.markers
          }
          if (input$assays1 == "SCTransform"){ 
          ips.markers <- FindAllMarkers(v$scData1, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = "SCT", test.use = input$test.use)
          v$ips.markers <- ips.markers
          }
          shinyalert("DEGs estimated", "DEGs estimated, please perform do data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    }
  })
  
  output$deg.gene.select <- renderUI({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
     selectInput("deg.gene", label = "Gene to visualise",
                  choices = rownames(v$scData1@assays$RNA@counts))
     }
  })
  
  output$Deg.plot <- renderPlotly({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
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
      })
    }
  })
  
  output$Deg1.plot <- renderPlotly({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        if(input$deg2 == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        FeaturePlot(v$scData1, input$deg.gene)
        }
        else if(input$deg2 == "primary.predict"){
        v$scData1$primary.predict -> Idents(v$scData1)
        FeaturePlot(v$scData1, input$deg.gene)
        }
      })
    }
  })
  
  output$Deg2.plot <- renderPlot({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        if(input$deg2 == "seurat_clusters"){
        v$scData1$seurat_clusters -> Idents(v$scData1)
        RidgePlot(v$scData1, features = input$deg.gene)
        }
        else if(input$deg2 == "primary.predict"){
        v$scData1$primary.predict -> Idents(v$scData1)
        RidgePlot(v$scData1, features = input$deg.gene)
        }
      })
    }
  })
  
  output$Deg3.plot <- renderPlot({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      withProgress(message="Generating DEG Plot...", value=0, {
        v$ips.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
        DoHeatmap(v$scData1, features = top10$gene, size = 5, angle = 45) + theme(axis.text.y = element_text(size = 4)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
      })
    }
  })
  
  output$Deg.table <- DT::renderDataTable(
    v$ips.markers, options = list(scrollX = TRUE, scrollY = "400px"))
  
  output$gene1.select <- renderUI({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      if(input$deg3 == "seurat_clusters"){
      selectInput("gene1", label = "Celltype1",
                  choices = as.vector(v$scData1$seurat_clusters))
      }
      else if(input$deg3 == "primary.predict"){
        selectInput("gene1", label = "Celltype1",
                    choices = as.vector(v$scData1$primary.predict))
      }
    }
  })
  
  output$gene2.select <- renderUI({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      if(input$deg3 == "seurat_clusters"){
      selectInput("gene2", label = "Celltype2",
                  choices = as.vector(v$scData1$seurat_clusters))
      }
      else if(input$deg3 == "primary.predict"){
        selectInput("gene2", label = "Celltype2",
                    choices = as.vector(v$scData1$primary.predict))
      }
    }
  })
  
  observeEvent(input$doVolcano, {
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      withProgress(message="Finding DEGs...", value=0, {
        v$scData1_subset <- subset(v$scData1, subset = primary.predict == input$gene1 | primary.predict == input$gene2)
        if (input$assays1 == "LogNormalization"){ 
          ips.markers_a <- FindAllMarkers(v$scData1_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
          ips.markers_b <- FindMarkers(v$scData1_subset, ident.1 = input$gene1, ident.2 = input$gene2, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "RNA", test.use = input$test.use_a)
          v$ips.markers_a <- ips.markers_a
          v$ips.markers_b <- ips.markers_b
        }
        if (input$assays1 == "SCTransform"){ 
          ips.markers_a <- FindAllMarkers(v$scData1_subset, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "SCT", test.use = input$test.use_a)
          ips.markers_b <- FindMarkers(v$scData1_subset, ident.1 = input$gene1, ident.2 = input$gene2, only.pos = F, min.pct = input$min_pct_a, logfc.threshold = input$logfc_a, assay = "SCT", test.use = input$test.use_a)
          v$ips.markers_a <- ips.markers_a
          v$ips.markers_b <- ips.markers_b
        }
        shinyalert("Volcano plot done", "Volcano plot done, please run GSEA Analysis", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
  })
  
  output$volcano.plot <- renderPlot({
    if(is.null(v$ips.markers_b)){
      return(NULL)
    }else{
      withProgress(message="Generating Volcano Plot...", value=0, {
        EnhancedVolcano(toptable = v$ips.markers_b, lab = row.names(v$ips.markers_b), x ="avg_log2FC", y ="p_val_adj", pointSize = 1, labSize = 5, legendLabSize = 12, axisLabSize = 12)
      })
    }
  })
  
  output$dega.plot <- renderPlot({
    if(is.null(v$ips.markers_a)){
      return(NULL)
    }else{
      withProgress(message="Generating Volcano Plot...", value=0, {
        v$ips.markers_a %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
        DoHeatmap(v$scData1_subset, features = top10$gene, size = 5, angle = 45) + theme(axis.text.y = element_text(size = 8)) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()
      })
    }
  })
  
  ##---------------GSEA of scRNA-seq-------------------##
  
  output$gsea.ct1.select <- renderUI({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      selectInput("gsea.ct1", label = "Celltype1",
                  choices = as.vector(v$scData1$primary.predict))
    }
  })
  
  output$gsea.ct2.select <- renderUI({
    if(is.null(v$ips.markers)){
      return(NULL)
    }else{
      selectInput("gsea.ct2", label = "Celltype2",
                  choices = as.vector(v$scData1$primary.predict))
    }
  })
  
  observeEvent(input$gsea, {
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
        v$msigdbr_hs_go <- msigdbr(species = "Homo sapiens", category = "C2") %>% filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
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
        v$msigdbr_mm_go <- msigdbr(species = "Mus musculus", category = "C2") %>% filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME"))
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
    v$fgseaRes, options = list(scrollX = TRUE, scrollY = "400px"))
  
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
    withProgress(message = "Running cell-cell communication...", value = 0.3, {
      v$cellphone <- cellphone_for_seurat(v$scData1)
      v$pvals <- read.delim("out/pvalues.txt", check.names = FALSE)
      v$means <- read.delim("out/means.txt", check.names = FALSE)
      v$decon <- read.delim("out/deconvoluted.txt", check.names = FALSE)
      print(v$pvals)
      print(v$means)
      output$CC.done <- renderText(paste0("Cell-cell communication done!"))
      v$isCCdone <- TRUE
      shinyalert("Cell-cell communication done", "Cell-cell communication done", type = "success", imageWidth = 10, imageHeight = 10)
    })
  })
  
  output$CC.gene.select <- renderUI({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      selectInput("CC.gene1", label = "Source",
                  choices = unique(v$scData1$primary.predict))
    }
  })
  
  output$CC.gene1.select <- renderUI({
    if(is.null(v$scData1)){
      return(NULL)
    }else{
      selectInput("CC.gene2", label = "Receptor",
                  choices = unique(v$scData1$primary.predict))
    }
  })
  
  output$CC_plot1 <- renderPlotly({
    if(is.null(v$isCCdone)){
      return(NULL)
    }else{
      withProgress(message="Generating cell-cell communication plot...", value=0, {
        heatmaps_plot(meta_file = 'metadata.txt', pvalues_file = 'out/pvalues.txt', count_filename = 'heatmap_count.pdf', log_filename = 'heatmap_log_count.pdf', count_network_filename = 'count_network.txt', interaction_count_filename = 'interactions_count.txt', count_network_separator = '\t', interaction_count_separator = '\t', show_rownames = T, show_colnames = T, scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11, fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0, meta_sep='\t', pvalues_sep='\t', pvalue=input$cc_pval)
      })
    }
  })
  
  output$CC_plot2 <- renderPlot({
    if(is.null(v$isCCdone)){
      return(NULL)
    }else{
      withProgress(message="Generating cell-cell communication plot...", value=0, {
        plot_cpdb3(cell_type1 = input$CC.gene1, cell_type2 = input$CC.gene2,
                   scdata = as.SingleCellExperiment(v$scData1),
                   idents = 'primary.predict',
                   means = v$means,
                   pvals = v$pvals,
                   deconvoluted = v$decon, # new options from here on specific to plot_cpdb3
                   keep_significant_only = TRUE,
                   standard_scale = TRUE,
                   remove_self = TRUE)
      })
    }
  })
  
  ##------------------------Data integration module--------------------------##
  
  output$integration_image <- renderImage({
    list(src = "www/integrate.png",
         height = 500)
  }, deleteFile = FALSE)
  
  observeEvent(input$loadexample1, {
      withProgress(message="Loading example Data...", value=0.5, {
        tpmFiles1 <- read.table("integration/concatenated_expr_data.txt", header = T, row.names = 1, check.names = F, sep = '\t')
        #scH5 <- input$scH5
        annoFile1 <- read.table("integration/metadata.txt", header = T, row.names = 1, sep = '\t')
      })
      if (is.null(tpmFiles1)){
        v$scData1 <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(tpmFiles1)
          #print(tpmFiles$name)
          print(file.exists(paste(tpmFiles1$datapath[1], "/", tpmFiles1$name[1], sep="")))
          v$tpmFiles1 <- tpmFiles1
          v$anno <- annoFile1
          
          sObj1 <- CreateSeuratObject(v$tpmFiles1,
                                      meta.data = v$anno,
                                      project = input$projName,
                                      min.genes = input$min.genes,
                                      min.cells = input$min.cells)
          
          sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
          v$scData1 <- sObj1
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample1", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        })
    }
  })
  
  observeEvent(input$loadexample1a, {
    withProgress(message="Loading example Data...", value=0.5, {
      data1 <- Read10X_h5("integration/filtered_feature_bc_matrix_X0029_downsample.h5")
      data2 <- Read10X_h5("integration/filtered_feature_bc_matrix_X0030_downsample.h5")
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
          data[[i]] <- CreateSeuratObject(as.matrix(data[[i]]), project = input$projName, min.genes = input$min.genes, min.cells = input$min.cells)
        }
        sObj1 <- merge(data[[1]], data[[2]])
        
        sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
        v$scData1 <- sObj1
        label1 <- "Example loaded"
        updateActionButton(inputId = "loadexample1", label = label1)
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
        v$scData1 <- NULL
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
          # print(exp.data1)
          exp.data2 <- merge(exp.data1[[1]], exp.data1[2:length(exp.data1)])
          
          if(!is.null(annoFile1)){
            anno.data1 <- rbindlist(lapply(annoFile1$datapath, fread), use.names = TRUE, fill = TRUE)
            anno.data2 <- data.frame(anno.data1, row.names = 1)
          }
          print(anno.data2)
          
          incProgress(0.5, "Creating Seurat Object")
          
          sObj1 <- AddMetaData(exp.data2, anno.data2)
          # print(exp.data2)
          v$scData1 <- exp.data2
          v$anno <- anno.data2
          print(v$anno)
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
    else if(input$scInput1 == "10X cellranger"){
      scH5_1 <- input$scH5_1
      
      if (is.null(scH5_1)){
        v$scData1 <- NULL
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
          v$scData1 <- exp.data2
          shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
  })
  
  output$countdataDT1 <- renderDataTable({
    if(!is.null(v$scData1))
    {
      if(ncol(v$scData1) > 20 )
        return(as.matrix(v$scData1@assays$RNA@counts[,1:20]))
    }
  })
  
  observeEvent(input$create_seurat1, {
    if(input$scInput1 == "Raw Counts Matrix"){
    withProgress(message="Loading and Processing Data...", value=0, {
      print(v$scData1)
      print(v$anno)
      
      sObj1 <- CreateSeuratObject(v$scData1@assays$RNA@counts,
                                  meta.data = v$anno,
                                  project = input$projName1,
                                  min.genes = input$min.genes1,
                                  min.cells = input$min.cells1)
      
      
      sObj1$orig.ident <- "Integration"
      Idents(sObj1) <- sObj1$orig.ident
      sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
      #sObj1[["batch"]] <- ifelse(endsWith(sObj1@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
      v$scData2 <- sObj1
      print(v$scData2@meta.data)
      shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
    else if(input$scInput1 == "10X cellranger"){
      withProgress(message="Loading and Processing Data...", value=0, {
        print(v$scData1)
        print(v$anno)
        
        sObj1 <- CreateSeuratObject(v$scData1@assays$RNA@counts,
                                    meta.data = v$anno,
                                    project = input$projName1,
                                    min.genes = input$min.genes1,
                                    min.cells = input$min.cells1)
        
        
        sObj1$orig.ident <- "Integration"
        Idents(sObj1) <- sObj1$orig.ident
        sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^MT-")
        sObj1[["batch"]] <- ifelse(endsWith(sObj1@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
        v$scData2 <- sObj1
        print(v$scData2@meta.data)
        shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
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
      VlnPlot(v$scData2, "nFeature_RNA") + NoLegend()
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
  
  observe({if(input$scAnalysis_integ == "Seurat" || input$scAnalysis_integ == "Harmony" || input$scAnalysis_integ == "fastMNN"){
  observeEvent(input$findVarGenes_bef_intg, {
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
        }, height = 800, width = 850)
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
    })
  }
})  
  
  observe({if(input$scAnalysis_integ == "Seurat"){
    
    observeEvent(input$runPCA_bef_intg_seurat, {
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunPCA(v$scData2, verbose = FALSE)
          print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
          
          v$isPCAdone1 <- TRUE
          
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
          print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
          v$isPCAdone1 <- TRUE
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
    })
    
    output$PCAplot_bef_seurat_tpm1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_seurat_tpm2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_seurat_tpm3 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
       withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
        })
      }
     })
    
    output$PCAplot_bef_seurat_h5_1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_seurat_h5_2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$vizPlot_bef_intg_seurat <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        VizDimLoadings(v$scData2, dims = as.numeric(input$select.pc_bef_intg_seurat))
      }
    })
    
    output$PCHeatmap_bef_intg_seurat <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        DimHeatmap(v$scData2, dims = as.numeric(input$select.pc_bef_intg_seurat))
      }
    })
    
    output$PCtable_bef_intg_seurat <- DT::renderDataTable({
      if(is.null(v$scData2) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, options = list(scrollX = TRUE))
    
    output$Elbow_bef_intg_seurat <- renderPlot({
      if(is.null(v$scData2@reductions$pca@jackstraw)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData2, ndims = 50)
        })
      }
    }, height = 800, width = 850)
    
    observeEvent(input$findCluster_bef_intg_seurat, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        v$scData2 <- FindNeighbors(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_seurat)
        v$scData2 <- FindClusters(v$scData2, resolution = input$clus.res_bef_intg_seurat)
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please perform UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Cluster2DPlot_bef_intg_seurat <- renderPlotly({
      if(is.null(v$isClusterdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    observeEvent(input$runUMAP_bef_intg_seurat, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_seurat)
          v$isUMAPdone1 <- TRUE
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
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          shinyalert("UMAP performed", "UMAP performed, please perform tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$UMAPplot_bef_seurat_tpm1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_seurat_tpm2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_seurat_tpm3 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_seurat_h5_1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_seurat_h5_2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$runTSNE_bef_intg_seurat, {
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          shinyalert("tSNE performed", "tSNE performed, please perform celltype annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          shinyalert("tSNE performed", "tSNE performed, please perform celltype annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$TSNEplot_bef_seurat_tpm1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_seurat_tpm2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_seurat_tpm3 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_seurat_h5_1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_seurat_h5_2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$doCELLiD_bef_intg_seurat, {
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scData2.rna.data.average = AverageExpression(v$scData2)
        v$scData2.rna.data.average = round(v$scData2.rna.data.average$RNA, 2)
        #v$scData2.rna.data.average = data.frame(v$scData2.rna.data.average$RNA)
        if(input$cellatlas1 == "all"){
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, ref)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, adipose1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, adrenal_gland1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, blood1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, bone_marrow1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, brain1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, breast1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, breast_milk1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, eye1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, gut1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, heart1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, kidney1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, liver1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, lung1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, pancreas1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, PDAC1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, skin1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, testis1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, thymus1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1 == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res = FastIntegration::CELLiD(v$scData2.rna.data.average, tonsil1)
          print(v$res)
          v$scData2$primary.predict = v$res[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDdone <- TRUE
        shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform data integration", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Umap_cellid_bef_intg_seurat <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCELLiDdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "primary.predict", label = F)
        })
      }
    })
    
    output$Umap_cellid_bef_intg_seurat1 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCELLiDdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$ct_bef_intg_seurat.table <- DT::renderDataTable(
      v$res, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doIntg_seurat, {
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
        DefaultAssay(v$scData1.combined) <- "integrated"
        v$scData1.combined <- ScaleData(v$scData1.combined, verbose = FALSE)
        v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
        v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
        v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_seurat)
        print(v$scData1.combined)
        shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    observeEvent(input$runPCA_intg_seurat, {
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          PCA_plot1a <- DimPlot(v$scData1.combined, reduction = "pca", label = T, label.size = 3)
          PCA_plot1b <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
          PCA_plot1c <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
          print(PCA_plot1a)
          print(PCA_plot1b)
          print(PCA_plot1c)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData1.combined)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          shinyalert("PCA done", "PCA done, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          PCA_plot1a <- DimPlot(v$scData1.combined, reduction = "pca", label = T, label.size = 3)
          PCA_plot1b <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
          print(PCA_plot1a)
          print(PCA_plot1b)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData1.combined)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          shinyalert("PCA done", "PCA done, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
          }
        }
      )
    })
    
    output$PCAplot_seurat_tpm1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_seurat_tpm2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$PCAplot_seurat_tpm3 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$PCAplot_seurat_h5_1 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_seurat_h5_2 <- renderPlotly({
      if(is.null(v$scData1.combined) || is.null(v$isPCAdone1)){
        return(NULL)
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
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
    }, options = list(scrollX = TRUE))
    
    output$Elbow_intg_seurat <- renderPlot({
      if(is.null(v$scData1.combined@reductions$pca@jackstraw)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData1.combined, ndims = 50)
        })
      }
    }, height = 800, width = 850)
    
    observeEvent(input$findCluster_intg_seurat, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        
        v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
        v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_seurat)
        output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone1 <- TRUE
        shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Cluster2DPlot_intg_seurat <- renderPlotly({
      if(is.null(v$isClusterdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scData1.combined, reduction = "pca", label = T)
        })
      }
    })
    
    observeEvent(input$runUMAP_intg_seurat, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
           if (input$scInput1 == "Raw Counts Matrix")
          {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat, spread = 1)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat, spread = 1)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$UMAPplot_seurat_tpm1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_seurat_tpm2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_seurat_tpm3 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_seurat_h5_1 <- renderPlot({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_seurat_h5_2 <- renderPlot({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAP_lisi_seurat <- renderText({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          print(v$scData1.combined@reductions$pca@cell.embeddings)
          print(v$scData1.combined@meta.data)
          v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$pca@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
          print(v$test)
        })
      }
    })
    
    observeEvent(input$runTSNE_intg_seurat, {
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_seurat)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          shinyalert("tSNE done", "tSNE done, please perform cell type annotation", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$TSNEplot_seurat_tpm1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
        })
      }
    })
    
    output$TSNEplot_seurat_tpm2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_seurat_tpm3 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_seurat_h5_1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, label.size = 3)
        })
      }
    })
    
    output$TSNEplot_seurat_h5_2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$doCELLiD_intg_seurat, {
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scData1.combined.rna.data.average = AverageExpression(v$scData1.combined)
        v$scData1.combined.rna.data.average = round(v$scData1.combined.rna.data.average$RNA, 2)
        if(input$cellatlas1a == "all"){
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, ref)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adipose1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adrenal_gland1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, blood1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, bone_marrow1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, brain1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast_milk1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, eye1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, gut1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, heart1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, kidney1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, liver1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, lung1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, pancreas1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, PDAC1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, skin1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, testis1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, thymus1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1a == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res_intg_seurat = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, tonsil1)
          print(v$res_intg_seurat)
          v$scData1.combined$primary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_seurat[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_seurat) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_seurat, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDdone1 <- TRUE
        shinyalert("Celltype annotation done", "Celltype done, please perform DEG Analysis", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Umap_cellid_intg_seurat <- renderPlotly({
      if(is.null(v$isCELLiDdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "primary.predict", label =  T,  label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_intg_seurat1 <- renderPlotly({
      if(is.null(v$isCELLiDdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$ct_intg_seurat.table <- DT::renderDataTable(
      v$res_intg_seurat, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doDeg_intg_seurat, {
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          v$scData1.combined$celltype -> Idents(v$scData1.combined)
          ips.markers <- FindAllMarkers(v$scData1.combined, only.pos = FALSE, min.pct = input$min_pct_intg_seurat, logfc.threshold = input$logfc_intg_seurat, test.use = input$test.use_intg_seurat)
          v$ips.markers <- ips.markers
          shinyalert("DEG Analysis done", "DEG Analysis done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    })
    
    output$deg.gene.select_intg_seurat <- renderUI({
      if(is.null(v$ips.markers)){
        return(NULL)
      }else{
        selectInput("deg.gene_intg_seurat", label = "Gene to visualise",
                    choices = unique(v$ips.markers$gene))
      }
    })
    
    output$Deg.plot_intg_seurat <- renderPlotly({
      if(is.null(v$ips.markers)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          VlnPlot(v$scData1.combined, input$deg.gene_intg_seurat)
        })
      }
    })
    
    output$Deg1.plot_intg_seurat <- renderPlotly({
      if(is.null(v$ips.markers)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          FeaturePlot(v$scData1.combined, input$deg.gene_intg_seurat)
        })
      }
    })
    
    output$Deg2.plot_intg_seurat <- renderPlot({
      if(is.null(v$ips.markers)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          RidgePlot(v$scData1.combined, features = input$deg.gene_intg_seurat)
        })
      }
    })
    
    output$Deg3.plot_intg_seurat <- renderPlotly({
      if(is.null(v$ips.markers)){
        return(NULL)
      }else{
        withProgress(message="Generating DEG Plot...", value=0, {
          v$ips.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
          DoHeatmap(v$scData1.combined, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
        })
      }
    })
    
    output$Deg.table_intg_seurat <- DT::renderDataTable(
      v$ips.markers, options = list(scrollX = TRUE, scrollY = "400px"))
    
    }
  })    
  
  observe({if(input$scAnalysis_integ == "Harmony"){
    
    observeEvent(input$runPCA_bef_intg_harmony, {
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunPCA(v$scData2, verbose = FALSE)
          print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
          v$isPCAdone1 <- TRUE
          
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
          print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
          v$isPCAdone1 <- TRUE
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
    })
    
    output$PCAplot_bef_harmony_tpm1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_harmony_tpm2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_harmony_tpm3 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
       withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
        })
      }
     })
    
    output$PCAplot_bef_harmony_h5_1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_bef_harmony_h5_2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$vizPlot_bef_intg_harmony <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        VizDimLoadings(v$scData2, dims = as.numeric(input$select.pc_bef_intg_harmony))
      }
    })
    
    output$PCHeatmap_bef_intg_harmony <- renderPlot({
      if(is.null(v$scData2)){
        return(NULL)
      }else{
        DimHeatmap(v$scData2, dims = as.numeric(input$select.pc_bef_intg_harmony))
      }
    })
    
    output$PCtable_bef_intg_harmony <- DT::renderDataTable({
      if(is.null(v$scData2) ){
        return(NULL)
      }else{
        v$pcGenes
      }
    }, options = list(scrollX = TRUE))
    
    output$Elbow_bef_intg_harmony <- renderPlot({
      if(is.null(v$scData2@reductions$pca@jackstraw)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData2, reduction = "pca", ndims = 50)
        })
      }
    }, height = 800, width = 850)
    
    observeEvent(input$findCluster_bef_intg_harmony, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        
        v$scData2 <- FindNeighbors(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_harmony)
        v$scData2 <- FindClusters(v$scData2, resolution = input$clus.res_bef_intg_harmony)
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
        shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Cluster2DPlot_bef_intg_harmony <- renderPlotly({
      if(is.null(v$isClusterdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
        })
      }
    })
    
    observeEvent(input$runUMAP_bef_intg_harmony, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_harmony)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_harmony)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$UMAPplot_bef_harmony_tpm1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_harmony_tpm2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_harmony_tpm3 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_harmony_h5_1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_bef_harmony_h5_2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$runTSNE_bef_intg_harmony, {
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_harmony)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_harmony)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$TSNEplot_bef_harmony_tpm1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_harmony_tpm2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_harmony_tpm3 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_harmony_h5_1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
        })
      }
    })
    
    output$TSNEplot_bef_harmony_h5_2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
          DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$doCELLiD_bef_intg_harmony, {
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scData2.rna.data.average = AverageExpression(v$scData2)
        v$scData2.rna.data.average = round(v$scData2.rna.data.average$RNA, 2)
        if(input$cellatlas2 == "all"){
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, ref)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, adipose1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, adrenal_gland1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, blood1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, bone_marrow1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, brain1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, breast1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, breast_milk1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, eye1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, gut1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, heart1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, kidney1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, liver1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, lung1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, pancreas1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, PDAC1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, skin1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, testis1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, thymus1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas2 == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res1 = FastIntegration::CELLiD(v$scData2.rna.data.average, tonsil1)
          print(v$res1)
          v$scData2$primary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),1]
          v$scData2$secondary.predict = v$res1[as.numeric(v$scData2$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res1) <- newheaders
          print(v$scData2@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res1, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDdone <- TRUE
        shinyalert("Celltype annotation performed", "Celltype annotation performed, please perform data integration", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Umap_cellid_bef_intg_harmony <- renderPlotly({
      if(is.null(v$isCELLiDdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_bef_intg_harmony1 <- renderPlotly({
      if(is.null(v$scData2) || is.null(v$isCELLiDdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData2, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$ct_bef_intg_harmony.table <- DT::renderDataTable(
      v$res, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doIntg_harmony, {
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
        DefaultAssay(v$scData1.combined) <- "integrated"
        v$scData1.combined <- ScaleData(v$scData1.combined, verbose = FALSE)
        v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
        print(v$scData1.combined[["pca"]], dims = 1:5, nfeatures = 5)
        v$scData1.combined <- RunHarmony(v$scData1.combined, "batch", assay.use = "integrated")
        print(v$scData1.combined[["harmony"]], dims = 1:5, nfeatures = 5)
        v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg_harmony)
        v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_harmony)
        print(v$scData1.combined)
        shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    observeEvent(input$runPCA_intg_harmony, {
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
          if (input$scInput1 == "Raw Counts Matrix")
          {
          PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "harmony", label = T)
          PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "harmony", group.by = 'batch', label.size = 3)
          PCA_plot2c <- DimPlot(v$scData1.combined, reduction = "harmony", group.by = 'celltype', label.size = 3)
          print(PCA_plot2a)
          print(PCA_plot2b)
          print(PCA_plot2c)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData1.combined)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          shinyalert("PCA performed", "PCA performed, please run clustering", type = "success", imageWidth = 10, imageHeight = 10)
          }
        else if (input$scInput1 == "10X cellranger")
        {
          PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "harmony", label = T)
          PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "harmony", group.by = 'batch', label.size = 3)
          print(PCA_plot2a)
          print(PCA_plot2b)
          
          incProgress(0.4, message = "Getting list of PC genes...")
          pc.table <- list()
          for(i in 1:20){
            pcg <- TopFeatures(v$scData1.combined)
            pc.table[[i]] <- pcg
          }
          pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
          v$pcGenes <- pc.table
          shinyalert("PCA performed", "PCA performed, please run clustering", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$PCAplot_harmony_tpm1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T)
        })
      }
    })
    
    output$PCAplot_harmony_tpm2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$PCAplot_harmony_tpm3 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$PCAplot_harmony_h5_1 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, label.size = 3)
        })
      }
    })
    
    output$PCAplot_harmony_h5_2 <- renderPlotly({
      if(is.null(v$isPCAdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'batch', label.size = 3)
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
    }, options = list(scrollX = TRUE))
    
    output$Elbow_intg_harmony <- renderPlot({
      if(is.null(v$scData1.combined@reductions$pca@jackstraw)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scData1.combined, reduction = "harmony", ndims = 50)
        })
      }
    }, height = 800, width = 850)
  }
})  
  
  
  observeEvent(input$findCluster_intg_harmony, {
    withProgress(message = "Finding clusters...", value = 0.3, {
      
      v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg_harmony)
      v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_harmony)
      output$cluster2.done <- renderText(paste0("Clustering done!"))
      v$isClusterdone2 <- TRUE
      shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
    })
  })
  
  output$Cluster2DPlot_intg_harmony <- renderPlotly({
    if(is.null(v$isClusterdone2)){
      plotly_empty()
    }else{
      withProgress(message="Generating 2D Cluster Plot...", value=0, {
        DimPlot(v$scData1.combined, reduction = "harmony", label = T, label.size = 3)
      })
    }
  })
    
    observeEvent(input$runUMAP_intg_harmony, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "harmony", dims = 1:input$dim.used_intg_harmony, spread = 1)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_harmony, spread = 1)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$UMAPplot_harmony_tpm1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_harmony_tpm2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_harmony_tpm3 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_harmony_h5_1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_harmony_h5_2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAP_lisi_harmony <- renderText({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          print(v$scData1.combined@reductions$pca@cell.embeddings)
          print(v$scData1.combined@meta.data)
          v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$pca@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
          print(v$test)
        })
      }
    })
    
    observeEvent(input$runTSNE_intg_harmony, {
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_harmony)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:input$dim.used_intg_harmony)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$TSNEplot_harmony_tpm1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T)
        })
      }
    })
    
    output$TSNEplot_harmony_tpm2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_harmony_tpm3 <- renderPlotly({
     if(is.null(v$isTSNEdone1)){
        plotly_empty()
       }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_harmony_h5_1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T)
        })
      }
    })
    
    output$TSNEplot_harmony_h5_2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$doCELLiD_intg_harmony, {
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scData1.combined.rna.data.average = AverageExpression(v$scData1.combined)
        v$scData1.combined.rna.data.average = round(v$scData1.combined.rna.data.average$RNA, 2)
        if(input$cellatlas1b == "all"){
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, ref)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adipose1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adrenal_gland1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, blood1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, bone_marrow1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, brain1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast_milk1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, eye1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, gut1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, heart1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, kidney1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, liver1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, lung1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, pancreas1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, PDAC1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, skin1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, testis1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, thymus1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1b == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res_intg_harmony = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, tonsil1)
          print(v$res_intg_harmony)
          v$scData1.combined$primary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_harmony[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_harmony) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_harmony, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDdone2 <- TRUE
        shinyalert("Celltype annotation done", "Celltype done, please perform DEG Analysis", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Umap_cellid_intg_harmony <- renderPlotly({
      if(is.null(v$isCELLiDdone2)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_intg_harmony1 <- renderPlotly({
      if(is.null(v$isCELLiDdone2)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$ct_intg_harmony.table <- DT::renderDataTable(
      v$res_intg_harmony, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doDeg_intg_harmony, {
      if(is.null(v$scData1.combined)){
        return(NULL)
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
      v$ips.markers1, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observe({if(input$scAnalysis_integ == "fastMNN"){
      
      observeEvent(input$runPCA_bef_intg_fastmnn, {
        withProgress(message = "Running PCA...", value = 0,{
          incProgress(0.5, message = "Running PCA...")
          if (input$scInput1 == "Raw Counts Matrix")
          {
            v$scData2 <- RunPCA(v$scData2, verbose = FALSE)
            print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
            v$isPCAdone1 <- TRUE
            
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
            print(v$scData2[["pca"]], dims = 1:5, nfeatures = 5)
            v$isPCAdone1 <- TRUE
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
      })
      
      output$PCAplot_bef_fastmnn_tpm1 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
          })
        }
      })
      
      output$PCAplot_bef_fastmnn_tpm2 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$PCAplot_bef_fastmnn_tpm3 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'celltype', label.size = 3)
          })
        }
      })
      
      output$PCAplot_bef_fastmnn_h5_1 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
          })
        }
      })
      
      output$PCAplot_bef_fastmnn_h5_2 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "pca", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$vizPlot_bef_intg_fastmnn <- renderPlot({
        if(is.null(v$scData2)){
          return(NULL)
        }else{
          VizDimLoadings(v$scData2, dims = as.numeric(input$select.pc_bef_intg_fastmnn))
        }
      })
      
      output$PCHeatmap_bef_intg_fastmnn <- renderPlot({
        if(is.null(v$scData2)){
          return(NULL)
        }else{
          DimHeatmap(v$scData2, dims = as.numeric(input$select.pc_bef_intg_fastmnn))
        }
      })
      
      output$PCtable_bef_intg_fastmnn <- DT::renderDataTable({
        if(is.null(v$scData2) ){
          return(NULL)
        }else{
          v$pcGenes
        }
      }, options = list(scrollX = TRUE))
      
      output$Elbow_bef_intg_fastmnn <- renderPlot({
        if(is.null(v$scData2@reductions$pca@jackstraw)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData2, ndims = 50)
          })
        }
      }, height = 800, width = 850)
      
      observeEvent(input$findCluster_bef_intg_fastmnn, {
        withProgress(message = "Finding clusters...", value = 0.3, {
          
          v$scData2 <- FindNeighbors(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_fastmnn)
          v$scData2 <- FindClusters(v$scData2, resolution = input$clus.res_bef_intg_fastmnn)
          #output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone <- TRUE
          shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
        })
      })
      
      output$Cluster2DPlot_bef_intg_fastmnn <- renderPlotly({
        if(is.null(v$isClusterdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating 2D Cluster Plot...", value=0, {
            DimPlot(v$scData2, reduction = "pca", label = T, label.size = 3)
          })
        }
      })
      
      observeEvent(input$runUMAP_bef_intg_fastmnn, {
        withProgress(message = "Running UMAP...", value = 0,{
          incProgress(0.5, message = "Running UMAP...")
          if (input$scInput1 == "Raw Counts Matrix")
          {
            v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_fastmnn)
            v$isUMAPdone1 <- TRUE
            UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
            UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
            UMAP_plot1c <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
            print(UMAP_plot1a)
            print(UMAP_plot1b)
            print(UMAP_plot1c)
            shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
          }
          else if (input$scInput1 == "10X cellranger")
          {
            v$scData2 <- RunUMAP(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_fastmnn)
            v$isUMAPdone1 <- TRUE
            UMAP_plot1a <- DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
            UMAP_plot1b <- DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
            print(UMAP_plot1a)
            print(UMAP_plot1b)
            shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
          }
        })
      })
      
      output$UMAPplot_bef_fastmnn_tpm1 <- renderPlotly({
        if(is.null(v$isUMAPdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          })
        }
      })
      
      output$UMAPplot_bef_fastmnn_tpm2 <- renderPlotly({
        if(is.null(v$isUMAPdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$UMAPplot_bef_fastmnn_tpm3 <- renderPlotly({
        if(is.null(v$isUMAPdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          })
        }
      })
      
      output$UMAPplot_bef_fastmnn_h5_1 <- renderPlotly({
        if(is.null(v$isUMAPdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "umap", label = T, label.size = 3)
          })
        }
      })
      
      output$UMAPplot_bef_fastmnn_h5_2 <- renderPlotly({
        if(is.null(v$isUMAPdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      observeEvent(input$runTSNE_bef_intg_fastmnn, {
        withProgress(message = "Running TSNE...", value = 0,{
          incProgress(0.5, message = "Running TSNE...")
          if (input$scInput1 == "Raw Counts Matrix")
          {
            v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_fastmnn)
            v$isTSNEdone1 <- TRUE
            TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
            TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
            TSNE_plot1c <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
            print(TSNE_plot1a)
            print(TSNE_plot1b)
            print(TSNE_plot1c)
            shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
          }
          else if (input$scInput1 == "10X cellranger")
          {
            v$scData2 <- RunTSNE(v$scData2, reduction = "pca", dims = 1:input$dim.used_bef_intg_fastmnn)
            v$isTSNEdone1 <- TRUE
            TSNE_plot1a <- DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
            TSNE_plot1b <- DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
            print(TSNE_plot1a)
            print(TSNE_plot1b)
            shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
          }
        })
      })
      
      output$TSNEplot_bef_fastmnn_tpm1 <- renderPlotly({
        if(is.null(v$isTSNEdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          })
        }
      })
      
      output$TSNEplot_bef_fastmnn_tpm2 <- renderPlotly({
        if(is.null(v$isTSNEdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$TSNEplot_bef_fastmnn_tpm3 <- renderPlotly({
        if(is.null(v$isTSNEdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          })
        }
      })
      
      output$TSNEplot_bef_fastmnn_h5_1 <- renderPlotly({
        if(is.null(v$isTSNEdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "tsne", label = T, label.size = 3)
          })
        }
      })
      
      output$TSNEplot_bef_fastmnn_h5_2 <- renderPlotly({
        if(is.null(v$isTSNEdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating tSNE Plot of integrated dataset...", value=0, {
            DimPlot(v$scData2, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      observeEvent(input$doCELLiD_bef_intg_fastmnn, {
        withProgress(message = "Running CELLiD...", value = 0.3, {
          ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
          v$scData2.rna.data.average = AverageExpression(v$scData2)
          v$scData2.rna.data.average = round(v$scData2.rna.data.average$RNA, 2)
          if(input$cellatlas3 == "all"){
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, ref)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "adipose"){
            adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
            adipose1 <- ref[,adipose]
            colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, adipose1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "adrenal_gland"){
            adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
            adrenal_gland1 <- ref[,adrenal_gland]
            colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, adrenal_gland1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "blood"){
            blood <- colnames(ref)[grepl("blood",colnames(ref))] 
            blood1 <- ref[,blood]
            colnames(blood1) <- gsub("--blood","",colnames(blood1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, blood1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "bone_marrow"){
            bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
            bone_marrow1 <- ref[,bone_marrow]
            colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, bone_marrow1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "brain"){
            brain <- colnames(ref)[grepl("brain",colnames(ref))] 
            brain1 <- ref[,brain]
            colnames(brain1) <- gsub("--brain","",colnames(brain1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, brain1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "breast"){
            breast <- colnames(ref)[grepl("breast",colnames(ref))] 
            breast1 <- ref[,breast]
            colnames(breast1) <- gsub("--breast","",colnames(breast1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, breast1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "breast_milk"){
            breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
            breast_milk1 <- ref[,breast_milk]
            colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, breast_milk1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "eye"){
            eye <- colnames(ref)[grepl("eye",colnames(ref))] 
            eye1 <- ref[,eye]
            colnames(eye1) <- gsub("--eye","",colnames(eye1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, eye1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "gut"){
            gut <- colnames(ref)[grepl("gut",colnames(ref))] 
            gut1 <- ref[,gut]
            colnames(gut1) <- gsub("--gut","",colnames(gut1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, gut1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "heart"){
            heart <- colnames(ref)[grepl("heart",colnames(ref))] 
            heart1 <- ref[,heart]
            colnames(heart1) <- gsub("--heart","",colnames(heart1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, heart1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "kidney"){
            kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
            kidney1 <- ref[,kidney]
            colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, kidney1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "liver"){
            liver <- colnames(ref)[grepl("liver",colnames(ref))] 
            liver1 <- ref[,liver]
            colnames(liver1) <- gsub("--liver","",colnames(liver1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, liver1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "lung"){
            lung <- colnames(ref)[grepl("lung",colnames(ref))] 
            lung1 <- ref[,lung]
            colnames(lung1) <- gsub("--lung","",colnames(lung1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, lung1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "pancreas"){
            pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
            pancreas1 <- ref[,pancreas]
            colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, pancreas1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "PDAC"){
            PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
            PDAC1 <- ref[,PDAC]
            colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, PDAC1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "skin"){
            skin <- colnames(ref)[grepl("skin",colnames(ref))] 
            skin1 <- ref[,skin]
            colnames(skin1) <- gsub("--skin","",colnames(skin1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, skin1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "testis"){
            testis <- colnames(ref)[grepl("testis",colnames(ref))] 
            testis1 <- ref[,testis]
            colnames(testis1) <- gsub("--testis","",colnames(testis1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, testis1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "thymus"){
            thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
            thymus1 <- ref[,thymus]
            colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, thymus1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          if(input$cellatlas3 == "tonsil"){
            tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
            tonsil1 <- ref[,tonsil]
            colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
            v$res2 = FastIntegration::CELLiD(v$scData2.rna.data.average, tonsil1)
            print(v$res2)
            v$scData2$primary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),1]
            v$scData2$secondary.predict = v$res2[as.numeric(v$scData2$seurat_clusters),2]
            newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
            colnames(v$res2) <- newheaders
            print(v$scData2@meta.data)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            write.table(v$res2, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
          }
          v$isCELLiDdone <- TRUE
          shinyalert("Celltype identification done", "Celltype identification done, please run Data Integration", type = "success", imageWidth = 10, imageHeight = 10)
        })
      })
      
      output$Umap_cellid_bef_intg_fastmnn <- renderPlotly({
        if(is.null(v$scData2) || is.null(v$isCELLiDdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP from CELLiD...", value=0, {
            DimPlot(v$scData2, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3)
          })
        }
      })
      
      output$Umap_cellid_bef_intg_fastmnn1 <- renderPlotly({
        if(is.null(v$scData2) || is.null(v$isCELLiDdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP from CELLiD...", value=0, {
            DimPlot(v$scData2, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
          })
        }
      })
      
      output$ct_bef_intg_fastmnn.table <- DT::renderDataTable(
        v$res2, options = list(scrollX = TRUE, scrollY = "400px"))
      
      observeEvent(input$doIntg_fastmnn, {
        withProgress(message = "Running Data Integration...", value = 0.3, {
          v$scData1.list <- SplitObject(v$scData2, split.by = "batch")
          v$scData1.list <- pbmclapply(mc.cores = 20, X = v$scData1.list, FUN = function(x) {
            x <- NormalizeData(x)
            x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
          })
          features <- SelectIntegrationFeatures(object.list = v$scData1.list, nfeatures = input$nfeatures_intg_fastmnn)
          print(features)
          print(v$scData2)
          print(v$scData1.list)
          v$scData1.combined <- RunFastMNN(object.list = v$scData1.list, features = features)
          v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg_fastmnn)
          v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_fastmnn)
          print(v$scData1.combined)
          shinyalert("Data integration done", "Data integration done, please run PCA", type = "success", imageWidth = 10, imageHeight = 10)
        })
      })
      
      observeEvent(input$runPCA_intg_fastmnn, {
        withProgress(message = "Running PCA...", value = 0,{
          incProgress(0.5, message = "Running PCA...")
          if (input$scInput1 == "Raw Counts Matrix")
          {
            print(v$scData1.combined[["mnn"]], dims = 1:5, nfeatures = 5)
            v$isPCAdone1 <- TRUE
            PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "mnn", label = T)
            PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "mnn", group.by = 'batch', label.size = 3)
            PCA_plot2c <- DimPlot(v$scData1.combined, reduction = "mnn", group.by = 'celltype', label.size = 3)
            print(PCA_plot2a)
            print(PCA_plot2b)
            print(PCA_plot2c)
            
            incProgress(0.4, message = "Getting list of PC genes...")
            pc.table <- list()
            for(i in 1:20){
              pcg <- TopFeatures(v$scData1.combined[["mnn"]])
              pc.table[[i]] <- pcg
            }
            pc.table <- as.data.frame(pc.table, col.names = paste0("MNN", 1:20))
            v$pcGenes <- pc.table
            shinyalert("PCA performed", "PCA performed, please run clustering", type = "success", imageWidth = 10, imageHeight = 10)
          }
          else if (input$scInput1 == "10X cellranger")
          {
            print(v$scData1.combined[["mnn"]], dims = 1:5, nfeatures = 5)
            v$isPCAdone1 <- TRUE
            PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "mnn", label = T)
            PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "mnn", group.by = 'batch', label.size = 3)
            print(PCA_plot2a)
            print(PCA_plot2b)
            
            incProgress(0.4, message = "Getting list of PC genes...")
            pc.table <- list()
            for(i in 1:20){
              pcg <- TopFeatures(v$scData1.combined[["mnn"]])
              pc.table[[i]] <- pcg
            }
            pc.table <- as.data.frame(pc.table, col.names = paste0("MNN", 1:20))
            v$pcGenes <- pc.table
            shinyalert("PCA performed", "PCA performed, please run clustering", type = "success", imageWidth = 10, imageHeight = 10)
          }
          
        })
      })
      
      output$PCAplot_fastmnn_tpm1 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T)
          })
        }
      })
      
      output$PCAplot_fastmnn_tpm2 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, group.by = 'batch', label.size = 3)
          })
        }
      })
      
      output$PCAplot_fastmnn_tpm3 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, group.by = 'celltype', label.size = 3)
          })
        }
      })
      
      output$PCAplot_fastmnn_h5_1 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
            DimPlot(v$scData1.combined, reduction = "mnn", label = T, label.size = 3)
          })
        }
      })
      
      output$PCAplot_fastmnn_h5_2 <- renderPlotly({
        if(is.null(v$isPCAdone1)){
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
      }, options = list(scrollX = TRUE))
      
      output$Elbow_intg_fastmnn <- renderPlot({
        if(is.null(v$scData1.combined@reductions$mnn@jackstraw)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData1.combined, reduction = "mnn", ndims = 50)
          })
        }
      }, height = 800, width = 850)
    }
    })
    
    observeEvent(input$findCluster_intg_fastmnn, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg_fastmnn)
        v$scData1.combined <- FindClusters(v$scData1.combined, resolution = input$clus.res_intg_fastmnn)
        #output$cluster2.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone3 <- TRUE
        shinyalert("Clustering done", "Clustering done, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Cluster2DPlot_intg_fastmnn <- renderPlotly({
      if(is.null(v$isClusterdone3)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scData1.combined, reduction = "mnn", label = T, label.size = 3)
        })
      }
    })
    
    observeEvent(input$runUMAP_intg_fastmnn, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg_fastmnn, spread = 1)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          print(UMAP_plot1c)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg_fastmnn, spread = 1)
          v$isUMAPdone1 <- TRUE
          UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
          UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
          print(UMAP_plot1a)
          print(UMAP_plot1b)
          shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$UMAPplot_fastmnn_tpm1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_fastmnn_tpm2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAPplot_fastmnn_tpm3 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = F, group.by = 'celltype')
        })
      }
    })
    
    output$UMAPplot_fastmnn_h5_1 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, label.size = 3)
        })
      }
    })
    
    output$UMAPplot_fastmnn_h5_2 <- renderPlotly({
      if(is.null(v$isUMAPdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$UMAP_lisi_fastmnn <- renderText({
      if(is.null(v$scData1.combined)){
        return(NULL)
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          print(v$scData1.combined@reductions$mnn@cell.embeddings)
          print(v$scData1.combined@meta.data)
          v$test <- paste("LISI Score:", median(lisi::compute_lisi(v$scData1.combined@reductions$mnn@cell.embeddings, v$scData1.combined@meta.data, 'batch')$batch))
          print(v$test)
        })
      }
    })
    
    observeEvent(input$runTSNE_intg_fastmnn, {
      withProgress(message = "Running TSNE...", value = 0,{
        incProgress(0.5, message = "Running TSNE...")
        if (input$scInput1 == "Raw Counts Matrix")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg_harmony)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          print(TSNE_plot1c)
          shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
        }
        else if (input$scInput1 == "10X cellranger")
        {
          v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "mnn", dims = 1:input$dim.used_intg_harmony)
          v$isTSNEdone1 <- TRUE
          TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
          TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
          print(TSNE_plot1a)
          print(TSNE_plot1b)
          shinyalert("tSNE done", "tSNE done, please run Celltype identification", type = "success", imageWidth = 10, imageHeight = 10)
        }
      })
    })
    
    output$TSNEplot_fastmnn_tpm1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T)
        })
      }
    })
    
    output$TSNEplot_fastmnn_tpm2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_fastmnn_tpm3 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype', label.size = 3)
        })
      }
    })
    
    output$TSNEplot_fastmnn_h5_1 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T)
        })
      }
    })
    
    output$TSNEplot_fastmnn_h5_2 <- renderPlotly({
      if(is.null(v$isTSNEdone1)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch', label.size = 3)
        })
      }
    })
    
    observeEvent(input$doCELLiD_intg_fastmnn, {
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$scData1.combined.rna.data.average = AverageExpression(v$scData1.combined)
        v$scData1.combined.rna.data.average = round(v$scData1.combined.rna.data.average$RNA, 2)
        if(input$cellatlas1c == "all"){
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, ref)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adipose1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, adrenal_gland1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, blood1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, bone_marrow1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, brain1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, breast_milk1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, eye1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, gut1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, heart1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, kidney1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, liver1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, lung1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, pancreas1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, PDAC1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, skin1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, testis1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, thymus1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas1c == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res_intg_fastmnn = FastIntegration::CELLiD(v$scData1.combined.rna.data.average, tonsil1)
          print(v$res_intg_fastmnn)
          v$scData1.combined$primary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),1]
          v$scData1.combined$secondary.predict = v$res_intg_fastmnn[as.numeric(v$scData1.combined$seurat_clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res_intg_fastmnn) <- newheaders
          print(v$scData1.combined@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res_intg_fastmnn, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDdone2 <- TRUE
        shinyalert("Celltype identification done", "Celltype identification done, please run DEG Analysis", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Umap_cellid_intg_fastmnn <- renderPlotly({
      if(is.null(v$isCELLiDdone2)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "primary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_intg_fastmnn1 <- renderPlotly({
      if(is.null(v$isCELLiDdone2)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$scData1.combined, reduction = "umap", group.by = "secondary.predict", label = T,  label.size = 3)
        })
      }
    })
    
    output$ct_intg_fastmnn.table <- DT::renderDataTable(
      v$res_intg_fastmnn, options = list(scrollX = TRUE, scrollY = "400px"))
    
    observeEvent(input$doDeg_intg_fastmnn, {
      if(is.null(v$scData1.combined)){
        return(NULL)
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
      v$ips.markers1, options = list(scrollX = TRUE, scrollY = "400px"))
 
    ##------------------------Single cell multiomics module--------------------------##
    
    output$multiomics_image <- renderImage({
      list(src = "www/multiomics_fig.png",
           width = 650,
           height = 700)
    }, deleteFile = FALSE)
    
    observeEvent(input$loadexample_cite_seurat, {
      withProgress(message="Loading example Data...", value=0.5, {
        cite.data <- Read10X_h5('multiomics/cite/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5')
      })
      if (is.null(cite.data)){
        v$scDatat <- NULL
      }else{
        withProgress(message="Loading and Processing Data...", value=0, {
          print(cite.data)
          label1 <- "Example loaded"
          updateActionButton(inputId = "loadexample_cite_seurat", label = label1)
          shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
          cite.data.rna <- cite.data$`Gene Expression`
          cite.data.ab <- cite.data$`Antibody Capture`
          incProgress(0.5, "Creating Seurat Object")
          sObj3 <- CreateSeuratObject(cite.data.rna)
          sObj3[["percent.mt"]] <- PercentageFeatureSet(sObj3, pattern = "^MT-")
          v$scDatat <- sObj3
          v$scDatat.rna <- cite.data.rna
          v$scDatat.ab <- cite.data.ab
          print(sObj3@meta.data)
          v$scDatab <- CreateAssayObject(counts = cite.data.ab)
          v$scDatat[["ADT"]] <- v$scDatab
          print(Assays(v$scDatat))
        })
      }
    })
    
    observeEvent(input$loadexample_multiome_seurat, {
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
          sObj4 <- CreateSeuratObject(multiome.rna)
          sObj4[["percent.mt"]] <- PercentageFeatureSet(sObj4, pattern = "^MT-")
          print(sObj4)
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
        tpmFiles2 <- input$tpmFiles2
        if (is.null(tpmFiles2)){
          v$scDatat <- NULL
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
            sObj3 <- CreateSeuratObject(exp.data3.rna)
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
    })
    
    output$countdataDT2b <- renderDataTable({
      if(!is.null(v$scDatat.ab))
      {
        if(ncol(v$scDatat.ab) > 20 )
          return(as.matrix(v$scDatat@assays$ADT@counts[,1:20]))
      }
    })
    
    observeEvent(input$create_seurat2a, {
      withProgress(message="Loading and Processing Data...", value=0, {
        print(v$scDatat)
        print(Assays(v$scDatat))
        
        sObj2 <- CreateSeuratObject(v$scDatat@assays$RNA@counts,
                                    project = input$projName2,
                                    min.genes = input$min.genes2,
                                    min.cells = input$min.cells2)
        
        sObj2$orig.ident <- "Multiomics"
        Idents(sObj2) <- sObj2$orig.ident
        sObj2[["percent.mt"]] <- PercentageFeatureSet(sObj2, pattern = "^MT-")
        #sObj1[["batch"]] <- ifelse(endsWith(sObj1@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
        v$scDatat <- sObj2
        print(v$scDatat@meta.data)
        v$scDatab <- CreateAssayObject(counts = v$scDatat.ab)
        v$scDatat[["ADT"]] <- v$scDatab
        print(Assays(v$scDatat))
        shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
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
        VlnPlot(v$scDatat, "nFeature_RNA") + NoLegend()
      }
    })
    
    output$nFeature_RNAPlot2b <- renderPlot({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatan, "nCount_ATAC") + NoLegend()
      }
    })
    
    output$mitoPlot2a <- renderPlot({
      if(is.null(v$scDatat)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatat, "percent.mt") + NoLegend()
      }
    })
    
    output$mitoPlot2b <- renderPlot({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatan, "percent.mt") + NoLegend()
      }
    })
    
    output$nCount_RNAPlot2a <- renderPlot({
      if(is.null(v$scDatat)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatat, "nCount_RNA") + NoLegend()
      }
    })
    
    output$nCount_RNAPlot2b <- renderPlot({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        VlnPlot(v$scDatan, "nCount_RNA") + NoLegend()
      }
    })
    
    output$FeatureScatterPlot_mult_a <- renderPlotly({
      if(is.null(v$scDatat)){
        plotly_empty()
      }else{
        print(FeatureScatter(v$scDatat, "nCount_RNA", "nFeature_RNA"))
      }
    })
    
    output$FeatureScatterPlot_mult_b <- renderPlotly({
      if(is.null(v$scDatat)){
        plotly_empty()
      }else{
        print(FeatureScatter(v$scDatat, "nCount_RNA", "percent.mt"))
      }
    })
    
    output$FeatureScatterPlot_mult_c <- renderPlotly({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        print(FeatureScatter(v$scDatan, "nCount_RNA", "nCount_ATAC"))
      }
    })
    
    output$FeatureScatterPlot_mult_d <- renderPlotly({
      if(is.null(v$scDatan)){
        plotly_empty()
      }else{
        print(FeatureScatter(v$scDatan, "nCount_RNA", "percent.mt"))
      }
    })
    
    observeEvent(input$findVarGenes_mult, {
      withProgress(message = "Finding variable genes...", value = 0, {
        
        v$scDatat <- NormalizeData(v$scDatat)
        v$scDatat <- FindVariableFeatures(v$scDatat,
                                          mean.function = ExpMean,
                                          dispersion.function = LogVMR,
                                          nfeatures = input$var.genes_mult,
                                          selection.method = input$selection.method)
        #all.genes <- rownames(v$scData1)
        v$scDatat <- ScaleData(v$scDatat)
        incProgress(0.5)
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
        }, height = 800, width = 850)
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
    })
    
    observeEvent(input$findVarGenes_mult1, {
      withProgress(message = "Finding variable genes...", value = 0, {
        
        v$scDatat <- NormalizeData(v$scDatat, normalization.method = "LogNormalize", assay = "RNA")
        v$scDatat <- ScaleData(v$scDatat, do.center = TRUE, do.scale = FALSE, assay = "RNA")
        v$scDatat <- NormalizeData(v$scDatat, normalization.method = "LogNormalize", assay = "ADT")
        v$scDatat <- ScaleData(v$scDatat, do.center = TRUE, do.scale = FALSE, assay = "ADT")
        v$scDatat <- FindVariableFeatures(v$scDatat, selection.method = input$selection.method_mult1, nfeatures = input$var.genes_mult1, assay = "RNA", verbose = FALSE)
        v$scDatat <- FindVariableFeatures(v$scDatat, selection.method = input$selection.method_mult1, nfeatures = input$var.genes_mult1, assay = "ADT", verbose = FALSE)
        
        incProgress(0.5)
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
        output$VarGenes_mult1aa <- renderPlot({
          varGenePlotInput()
        }, height = 800, width = 850)
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
    })
    
    observeEvent(input$runPCA_mult_seurat, {
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        v$scDatat <- RunPCA(v$scDatat, verbose = FALSE)
        print(v$scDatat[["pca"]], dims = 1:5, nfeatures = 5)
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
    }, options = list(scrollX = TRUE))
    
    output$Elbow_mult_seurat <- renderPlot({
      if(is.null(v$scDatat@reductions$pca@jackstraw)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scDatat, ndims = 50)
        })
      }
    }, height = 800, width = 850)
    
    observeEvent(input$runMOFA_mult, {
      withProgress(message = "Running PCA...", value = 0,{
        incProgress(0.5, message = "Running PCA...")
        v$mofa <- create_mofa(v$scDatat, assays = c("RNA","ADT"))
        print(v$mofa)
        model_opts <- get_default_model_options(v$mofa)
        model_opts$num_factors <- as.numeric(input$num_factor_mofa)
        v$mofa <- prepare_mofa(v$mofa, model_options = model_opts)
        v$mofa <- run_mofa(v$mofa)
        v$isMOFAdone <- TRUE
        incProgress(0.4, message = "Getting list of PC genes...")
        pc.table <- list()
        shinyalert("MOFA2 performed", "MOFA2 performed, please perform clustering", type = "success", imageWidth = 10, imageHeight = 10)
      }
      )
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
      withProgress(message = "Finding clusters...", value = 0.3, {
        
        v$scDatat <- FindNeighbors(v$scDatat, reduction = "pca", dims = 1:input$dim.used_mult_seurat)
        v$scDatat <- FindClusters(v$scDatat, resolution = input$clus.res_mult_seurat)
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Cluster2DPlot_mult_seurat <- renderPlotly({
      if(is.null(v$isClusterdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$scDatat, reduction = "pca", label = T)
        })
      }
    })
    
    observeEvent(input$findCluster_mofa2, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        
        v$factors <- 1:get_dimensions(v$mofa)[["K"]]
        v$factors <- v$factors[!v$factors%in%c(4,7)]
        cluster <- cluster_samples(v$mofa, k=input$num.clus.mofa, factors="all")
        clusters <- cluster$cluster
        v$clusters <- clusters
        samples_metadata(v$mofa) <- as.data.frame(v$clusters) %>% tibble::rownames_to_column("sample") %>% as.data.table
        v$mofa@samples_metadata -> df
        rownames(df) <- df[,1]
        v$t1 <- MOFA2::add_mofa_factors_to_seurat(mofa_object = v$mofa, seurat_object = v$scDatat, views = "all", factors = "all")
        v$t1 <- AddMetaData(v$t1, metadata = df)
        print(v$t1@meta.data)
        v$isMOFAdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Cluster2DPlot_mofa <- renderPlot({
      if(is.null(v$isMOFAdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating 2D Cluster Plot...", value=0, {
          DimPlot(v$t1, reduction = "MOFA", group.by = "v.clusters", label = F)
        })
      }
    })
    
    observeEvent(input$runUMAP_mult_seurat, {
      withProgress(message = "Running UMAP...", value = 0,{
        incProgress(0.5, message = "Running UMAP...")
        v$scDatat <- RunUMAP(v$scDatat, reduction = "pca", reduction.name = "rna.umap", dims = 1:input$dim.used_mult_seurat)
        DefaultAssay(v$scDatat) <- "ADT"
        VariableFeatures(v$scDatat) <- rownames(v$scDatat[["ADT"]])
        v$scDatat <- NormalizeData(v$scDatat, normalization.method = 'CLR', margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca')
        v$scDatat <- RunUMAP(v$scDatat, reduction = "apca", reduction.name = "adt.umap", dims = 1:input$dim.used_mult_seurat)
        v$scDatat <- FindMultiModalNeighbors(v$scDatat, reduction.list = list("pca", "apca"), dims.list = list(1:input$dim.used_mult_seurat, 1:input$dim.used_mult_seurat), modality.weight.name = "RNA.weight")
        v$scDatat <- FindClusters(v$scDatat, graph.name = "wsnn", algorithm = 3, resolution = input$clus.res_mult_seurat2a, verbose = FALSE)
        v$scDatat <- RunUMAP(v$scDatat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        v$isUMAPdone <- TRUE
        print(v$scDatat)
        UMAP_plot_cite_a <- DimPlot(v$scDatat, reduction = "rna.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        UMAP_plot_cite_b <- DimPlot(v$scDatat, reduction = "adt.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ADT")
        UMAP_plot_cite_c <- DimPlot(v$scDatat, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        print(UMAP_plot_cite_a)
        print(UMAP_plot_cite_b)
        print(UMAP_plot_cite_c)
        shinyalert("UMAP done", "UMAP done, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
      })
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
      withProgress(message = "Running UMAP...", value = 0.3, {
        v$mofa <- run_umap(v$mofa, factors = v$factors, n_neighbors = input$num.neighbors_mofa2, min_dist = 0.30)
        samples_metadata(v$mofa) <- as.data.frame(v$clusters) %>% tibble::rownames_to_column("sample") %>% as.data.table
        v$mofa@samples_metadata -> df
        rownames(df) <- df[,1]
        v$t1 <- MOFA2::add_mofa_factors_to_seurat(mofa_object = v$mofa, seurat_object = v$scDatat, views = "all", factors = "all")
        v$t1 <- AddMetaData(v$t1, metadata = df)
        v$t1$v.clusters -> Idents(v$t1)
        print(v$t1@meta.data)
        v$isMOFAdone <- TRUE
        shinyalert("UMAP performed", "Clustering performed, please run tSNE", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$UMAPplot_mofa_1 <- renderPlot({
      if(is.null(v$isMOFAdone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
          DimPlot(v$t1, group.by = "v.clusters", label = F)
        })
      }
    })
    
    observeEvent(input$runTSNE_mult_seurat, {
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
      v$res, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_cellid_cite_seurat_prediction <- downloadHandler(
      filename = function(){"CELLiD predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$res, fname)
        })
      }
    )
    
    observeEvent(input$doCELLiD_mult_cite_mofa, {
      withProgress(message = "Running CELLiD...", value = 0.3, {
        ref = readRDS(url("https://www.immunesinglecell.org/api/vishuo/download/getCellidRef"))
        v$t1.rna.data.average = AverageExpression(v$t1)
        v$t1.rna.data.average = round(v$t1.rna.data.average$RNA, 2)
        if(input$cellatlas_mult_cite_mofa == "all"){
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, ref)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "adipose"){
          adipose <- colnames(ref)[grepl("adipose",colnames(ref))] 
          adipose1 <- ref[,adipose]
          colnames(adipose1) <- gsub("--adipose","",colnames(adipose1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, adipose1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "adrenal_gland"){
          adrenal_gland <- colnames(ref)[grepl("adrenal_gland",colnames(ref))] 
          adrenal_gland1 <- ref[,adrenal_gland]
          colnames(adrenal_gland1) <- gsub("--adrenal_gland","",colnames(adrenal_gland1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, adrenal_gland1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "blood"){
          blood <- colnames(ref)[grepl("blood",colnames(ref))] 
          blood1 <- ref[,blood]
          colnames(blood1) <- gsub("--blood","",colnames(blood1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, blood1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "bone_marrow"){
          bone_marrow <- colnames(ref)[grepl("bone_marrow",colnames(ref))] 
          bone_marrow1 <- ref[,bone_marrow]
          colnames(bone_marrow1) <- gsub("--bone_marrow","",colnames(bone_marrow1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, bone_marrow1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "brain"){
          brain <- colnames(ref)[grepl("brain",colnames(ref))] 
          brain1 <- ref[,brain]
          colnames(brain1) <- gsub("--brain","",colnames(brain1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, brain1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "breast"){
          breast <- colnames(ref)[grepl("breast",colnames(ref))] 
          breast1 <- ref[,breast]
          colnames(breast1) <- gsub("--breast","",colnames(breast1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, breast1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "breast_milk"){
          breast_milk <- colnames(ref)[grepl("breast_milk",colnames(ref))] 
          breast_milk1 <- ref[,breast_milk]
          colnames(breast_milk1) <- gsub("--breast_milk","",colnames(breast_milk1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, breast_milk1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "eye"){
          eye <- colnames(ref)[grepl("eye",colnames(ref))] 
          eye1 <- ref[,eye]
          colnames(eye1) <- gsub("--eye","",colnames(eye1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, eye1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "gut"){
          gut <- colnames(ref)[grepl("gut",colnames(ref))] 
          gut1 <- ref[,gut]
          colnames(gut1) <- gsub("--gut","",colnames(gut1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, gut1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "heart"){
          heart <- colnames(ref)[grepl("heart",colnames(ref))] 
          heart1 <- ref[,heart]
          colnames(heart1) <- gsub("--heart","",colnames(heart1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, heart1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "kidney"){
          kidney <- colnames(ref)[grepl("kidney",colnames(ref))] 
          kidney1 <- ref[,kidney]
          colnames(kidney1) <- gsub("--kidney","",colnames(kidney1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, kidney1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "liver"){
          liver <- colnames(ref)[grepl("liver",colnames(ref))] 
          liver1 <- ref[,liver]
          colnames(liver1) <- gsub("--liver","",colnames(liver1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, liver1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "lung"){
          lung <- colnames(ref)[grepl("lung",colnames(ref))] 
          lung1 <- ref[,lung]
          colnames(lung1) <- gsub("--lung","",colnames(lung1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, lung1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "pancreas"){
          pancreas <- colnames(ref)[grepl("pancreas",colnames(ref))] 
          pancreas1 <- ref[,pancreas]
          colnames(pancreas1) <- gsub("--pancreas","",colnames(pancreas1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, pancreas1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "PDAC"){
          PDAC <- colnames(ref)[grepl("PDAC",colnames(ref))] 
          PDAC1 <- ref[,PDAC]
          colnames(PDAC1) <- gsub("--PDAC","",colnames(PDAC1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, PDAC1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "skin"){
          skin <- colnames(ref)[grepl("skin",colnames(ref))] 
          skin1 <- ref[,skin]
          colnames(skin1) <- gsub("--skin","",colnames(skin1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, skin1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "testis"){
          testis <- colnames(ref)[grepl("testis",colnames(ref))] 
          testis1 <- ref[,testis]
          colnames(testis1) <- gsub("--testis","",colnames(testis1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, testis1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "thymus"){
          thymus <- colnames(ref)[grepl("thymus",colnames(ref))] 
          thymus1 <- ref[,thymus]
          colnames(thymus1) <- gsub("--thymus","",colnames(thymus1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, thymus1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        if(input$cellatlas_mult_cite_mofa == "tonsil"){
          tonsil <- colnames(ref)[grepl("tonsil",colnames(ref))] 
          tonsil1 <- ref[,tonsil]
          colnames(tonsil1) <- gsub("--tonsil","",colnames(tonsil1))
          v$res = FastIntegration::CELLiD(v$t1.rna.data.average, tonsil1)
          print(v$res)
          v$t1$primary.predict = v$res[as.numeric(v$t1$v.clusters),1]
          v$t1$secondary.predict = v$res[as.numeric(v$t1$v.clusters),2]
          newheaders <- c("primary.predict", "secondary.predict", "primary score", "secondary score")
          colnames(v$res) <- newheaders
          print(v$t1@meta.data)
          output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
          write.table(v$res, "output_CELLiD.txt", sep = "\t", quote = F, row.names = F, col.names = F)
        }
        v$isCELLiDMultmofadone <- TRUE
        shinyalert("Cell type identification done", "Cell type identification done, please perform data visualization", type = "success", imageWidth = 10, imageHeight = 10)
      })
    })
    
    output$Umap_cellid_mult_cite_mofa <- renderPlot({
      if(is.null(v$t1) || is.null(v$isCELLiDMultmofadone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$t1, group.by = "primary.predict", label = F, label.size = 3)
        })
      }
    })
    
    output$Umap_cellid_mult_cite_mofa1 <- renderPlot({
      if(is.null(v$t1) || is.null(v$isCELLiDMultmofadone)){
        plotly_empty()
      }else{
        withProgress(message="Generating UMAP from CELLiD...", value=0, {
          DimPlot(v$t1, group.by = "secondary.predict", label = F, label.size = 3)
        })
      }
    })
    
    output$ct_cite_mofa.table <- DT::renderDataTable(
      v$res, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_cellid_cite_mofa_prediction <- downloadHandler(
      filename = function(){"CELLiD predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$res, fname)
        })
      }
    )
    
    observeEvent(input$Vis3, {
      withProgress(message = "Visualizing...", value = 0,{
        incProgress(0.5, message = "Visualizing...")
        #DefaultAssay(v$scDatat) <- "ADT"
      })
    })
    
    output$vis.gene.select <- renderUI({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        selectInput("vis.gene", label = "Gene to visualise",
                    choices = rownames(v$scDatat[[input$assay_mult_cite_seurat2]]))
      }
    })
    
    output$vis.plot <- renderPlotly({
      if(is.null(v$scDatat)){
        return(NULL)
      }else{
        withProgress(message="Generating Feature Plot...", value=0, {
          FeaturePlot(v$scDatat, input$vis.gene, cols = c("lightgrey", "darkgreen"), order = T, reduction = input$assay_mult_cite_seurat1)
        })
      }
    })
    
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
    
    observeEvent(input$Vis3a, {
      withProgress(message = "Visualizing...", value = 0,{
        incProgress(0.5, message = "Visualizing...")
        DefaultAssay(v$scDatat) <- "ADT"
      })
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
          FeaturePlot(v$t1, input$vis.gene_a, cols = c("lightgrey", "darkgreen"), order = T, reduction = "MOFAUMAP")
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
          FeaturePlot(v$t1, input$vis.gene1a, cols = c("lightgrey", "darkgreen"), order = T, reduction = "MOFAUMAP")
        })
      }
    })
    
    observe({if(input$scAnalysis_mult == "Seurat" & input$scAnalysis_type == "Multiome"){
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
      
      observeEvent(input$loadButton2, {
        if(input$scAnalysis_mult == "Seurat" & input$scAnalysis_type == "Multiome"){
          tpmFiles3 <- input$tpmFiles3
          
          if (is.null(tpmFiles3)){
            v$scDatan <- NULL
          }else{
            withProgress(message="Loading and Processing Data...", value=0, {
              print(tpmFiles3$datapath)
              print(tpmFiles3$name)
              print(file.exists(paste(tpmFiles3$datapath[1], "/", tpmFiles3$name[1], sep="")))
              exp.data4 <- Read10X_h5(tpmFiles3$datapath)
              print(exp.data4)
              exp.data4.rna <- exp.data4$`Gene Expression`
              exp.data4.atac <- exp.data4$Peaks
              #additional.ident <- NULL
              incProgress(0.5, "Creating Seurat Object")
              sObj4 <- CreateSeuratObject(exp.data4.rna)
              sObj4[["percent.mt"]] <- PercentageFeatureSet(sObj4, pattern = "^MT-")
              print(sObj4)
              grange.counts <- StringToGRanges(rownames(exp.data4.atac), sep = c(":", "-"))
              grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
              exp.data4.atac <- exp.data4.atac[as.vector(grange.use), ]
              annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
              genome(annotations) <- "hg19"
              seqlevelsStyle(annotations) <- 'UCSC'
              chrom_assay <- CreateChromatinAssay(
                counts = exp.data4.atac,
                sep = c(":", "-"),
                genome = 'hg19',
                fragments = as.character(parseFilePaths(c(home = '.'), dir_multi_atac())$datapath),
                min.cells = 10,
                annotation = annotations)
              print(chrom_assay)
              sObj4[["ATAC"]] <- chrom_assay
              print(sObj4)
              v$scDatan <- sObj4
              v$scDatan.rna <- exp.data4.rna
              v$scDatan.atac <- exp.data4.atac
              print(Assays(v$scDatan))
              shinyalert("Data uploaded", "Data uploaded, please click on Process button", type = "success", imageWidth = 10, imageHeight = 10)
            })
          }
        }
      })
    } 
    })
    
    output$countdataDT2c <- renderDataTable({
      if(!is.null(v$scDatan.rna))
      {
        if(ncol(v$scDatan.rna) > 20 )
          return(as.matrix(v$scDatan@assays$RNA@counts[,1:20]))
      }
    })
    
    output$countdataDT2d <- renderDataTable({
      if(!is.null(v$scDatan.atac))
      {
        if(ncol(v$scDatan.atac) > 20 )
          return(as.matrix(v$scDatan@assays$ATAC@counts[,1:20]))
      }
    })
    
    observeEvent(input$create_seurat2b, {
      withProgress(message="Loading and Processing Data...", value=0, {
        print(v$scDatan)
        print(Assays(v$scDatan))
        print(v$scDatan@assays$RNA@counts)
        
        sObj4 <- CreateSeuratObject(v$scDatan@assays$RNA@counts,
                                    project = input$projName2,
                                    min.genes = input$min.genes2,
                                    min.cells = input$min.cells2)
        
        
        sObj4$orig.ident <- "Multiomics"
        Idents(sObj4) <- sObj4$orig.ident
        sObj4[["percent.mt"]] <- PercentageFeatureSet(sObj4, pattern = "^MT-")
        #sObj1[["batch"]] <- ifelse(endsWith(sObj1@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
        v$scDatan <- sObj4
        print(v$scDatan@meta.data)
        v$scDatatac <- CreateAssayObject(counts = v$scDatan.atac)
        v$scDatan[["ATAC"]] <- v$scDatatac
        print(Assays(v$scDatan))
        print(v$scDatan)
        
        shinyalert("Data processed", "Data processed, please view the Violin Plots", type = "success", imageWidth = 10, imageHeight = 10)
      })
    }
    )
    
    observeEvent(input$doSCTransform_multi, {
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
        }, height = 800, width = 850)
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
    })
    
    observeEvent(input$runPCA_mult_seurat1, {
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
      )
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
    }, options = list(scrollX = TRUE))
    
    output$Elbow_mult_seurat1 <- renderPlot({
      if(is.null(v$scDatan@reductions$pca@jackstraw)){
        return(NULL)
      }else{
        withProgress(message="Generating Elbow Plot...", value=0.5, {
          ElbowPlot(v$scDatan)
        })
      }
    }, height = 800, width = 850)
    
    observeEvent(input$findCluster_mult_seurat1, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        
        v$scDatan <- FindNeighbors(v$scDatan, reduction = "pca", dims = 1:input$dim.used_mult_seurat1)
        v$scDatan <- FindClusters(v$scDatan, resolution = input$clus.res_mult_seurat1)
        #output$cluster1.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
        shinyalert("Clustering performed", "Clustering performed, please run UMAP", type = "success", imageWidth = 10, imageHeight = 10)
      })
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
      v$res, options = list(scrollX = TRUE, scrollY = "400px"))
    
    output$download_cellid_multiome_seurat_prediction <- downloadHandler(
      
      filename = function(){"CELLiD predictions.csv"}, 
      content = function(fname){
        withProgress(message="Downloading CELLiD predictions...", value=0, {
          write.csv(v$res, fname)
        })
      }
    )
    
    observeEvent(input$Vis3_multiome, {
      withProgress(message = "Visualizing...", value = 0,{
        incProgress(0.5, message = "Visualizing...")
        DefaultAssay(v$scDatan) <- "RNA"
      })
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
    
      ##---------------Spatial Transcriptomics Analysis using Seurat-------------------##
      
      output$spatial_image <- renderImage({
        list(src = "www/spatial.png",
             height = 500)
      }, deleteFile = FALSE)
      
      shinyDirChoose(
        input,
        'dir',
        roots = c(home = '.'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
      )
      
      dir <- reactive(input$dir)
      output$dir <- renderPrint({  # use renderText instead of renderPrint
        parseDirPath(c(home = '.'), dir())
      })
      
      observeEvent(input$loadexample_seurat_spatial, {
        withProgress(message="Loading example Data...", value=0.5, {
          sp.data <- Load10X_Spatial(data.dir = "spatial/V1_Mouse_Brain_Sagittal_Anterior/", filename = "filtered_feature_bc_matrix.h5")
        })
        if (is.null(sp.data)){
          v$scData1 <- NULL
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
          if (is.null(tpmFile_spatial)){
            v$scData_spatial <- NULL
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
      
      output$countdataDT_spatial <- renderDataTable({
        if(!is.null(v$scData_spatial))
        {
          if(ncol(v$scData_spatial) > 20 )
            return(as.matrix(v$scData_spatial@assays$Spatial@counts[,1:20]))
        }
      })
      
      output$h_e <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          print(SpatialDimPlot(v$scData_spatial, pt.size.factor = 0) + NoLegend())
        }
      })
      
      output$countdataDT_xenium <- renderDataTable({
        if(!is.null(v$scData_xenium))
        {
          if(ncol(v$scData_xenium) > 20 )
            return(as.matrix(v$scData_xenium@assays$Xenium@counts[,1:20]))
        }
      })
      
      output$h_e_xenium <- renderPlot({
        if(is.null(v$scData_xenium)){
          plotly_empty()
        }else{
          print(ImageDimPlot(v$scData_xenium, fov = "fov") + NoLegend())
        }
      })
      
      observeEvent(input$create_seurat_sp, {
        if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "Seurat" | input$scAnalysis_sp == "GraphST"){
        withProgress(message="Loading and Processing Data...", value=0, {
          print(v$scData_spatial)
          print(Assays(v$scData_spatial))
          v$scData_spatial$orig.ident <- "Spatial"
          Idents(v$scData_spatial) <- v$scData_spatial$orig.ident
         
          print(v$scData_spatial@meta.data)
          
          shinyalert("Data processed", "Data processed, please view the Spatial Plots", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
    }
  )
    
      output$nCount_SpatialPlot <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          VlnPlot(v$scData_spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
        }
      })
      
      output$SpatialFeaturePlot <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          SpatialFeaturePlot(v$scData_spatial, features = "nCount_Spatial") + theme(legend.position = "right")
        }
      })
      
      output$FeatureScatterPlot_sp <- renderPlotly({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          print(FeatureScatter(v$scData_spatial, "nCount_Spatial", "nFeature_Spatial"))
        }
      })
      
      observeEvent(input$doSCTransform_sp, {
        if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "Seurat" | input$scAnalysis_sp == "GraphST"){
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
          }, height = 800, width = 850)
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
    })
      
      observeEvent(input$doSCTransform_sp, {
        if(input$scInput3 == "SpaceRanger output" & input$scAnalysis_sp == "Seurat" | input$scAnalysis_sp == "GraphST"){
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
            }, height = 800, width = 850)
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
      })
      
      observeEvent(input$doSCTransform_xenium, {
        if(input$scAnalysis_platform == "Xenium"){
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
            }, height = 800, width = 850)
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
      })
      
      observeEvent(input$Vis_sp, {
        withProgress(message = "Visualizing...", value = 0,{
          incProgress(0.5, message = "Visualizing...")
          DefaultAssay(v$scData_spatial) <- "Spatial"
        })
      })
      
      output$sp.gene.select <- renderUI({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          selectInput("sp.gene", label = "Gene to visualise",
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
        )
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
      }, options = list(scrollX = TRUE))
      
      output$Elbow_spatial <- renderPlot({
        if(is.null(v$scData_spatial@reductions$pca@jackstraw)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData_spatial, ndims = 50)
          })
        }
      }, height = 800, width = 850)
      
      observeEvent(input$runPCA_xenium, {
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
        )
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
        if(is.null(v$scData_xenium@reductions$pca@jackstraw)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0.5, {
            ElbowPlot(v$scData_xenium, ndims = 30)
          })
        }
      }, height = 800, width = 850)
      
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
        withProgress(message = "Finding clusters...", value = 0.3, {
          if(input$scAnalysis_sp1 == "Seurat"){
          v$scData_spatial <- FindNeighbors(v$scData_spatial, reduction = "pca", dims = 1:input$dim.used_spatial)
          v$scData_spatial <- FindClusters(v$scData_spatial, resolution = input$clus.res_spatial, graph.name = "SCT_nn")
          #output$cluster1.done <- renderText(paste0("Clustering done!"))
          v$isClusterdone <- TRUE
          }
        })
      })
      
      output$Cluster2DPlot_spatial <- renderPlotly({
        if(is.null(v$isClusterdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating 2D Cluster Plot...", value=0, {
            DimPlot(v$scData_spatial, reduction = "pca", label = T)
          })
        }
      })
      
      observeEvent(input$findCluster_graphst, {
        withProgress(message = "Finding clusters...", value = 0.3, {
          if(input$scAnalysis_sp1 == "GraphST"){
            sc <- import("scanpy", convert = FALSE)
            #scvi <- import("scvi", convert = FALSE)
            scipy <- import("scipy", convert = FALSE)
            pot <- import("ot", convert = FALSE)
            pandas <- import("pandas", convert = F)
            gt <- import("GraphST", convert = F)
            matplotlib <- import ("matplotlib", convert = F)
            gt1 <- import("GraphST.preprocess", convert = F)
            gt2 <- import("GraphST.utils", convert = F)
            torch <- import("torch", convert = FALSE)
            np <- import("numpy", convert = F)
            pd <- import("pandas", convert = F)
            sns <- import("seaborn", convert = F)
            tpmFile_graphst <- input$tpmFile_graphst
            if (is.null(tpmFile_graphst)){
              v$scData_graphst <- NULL
            }else{
              withProgress(message="Loading and Processing Data...", value=0, {
                print(tpmFile_graphst$datapath)
                print(tpmFile_graphst$name)
                print(file.exists(paste(tpmFile_graphst$datapath[1], "/", tpmFile_graphst$name[1], sep="")))
                file_fold <- parseDirPath(c(home = '.'), dir_graphst())
                adata <- sc$read_visium(file_fold, count_file=tpmFile_graphst$name, load_images=TRUE)
                print(adata)
                adata$var_names_make_unique()
                model = gt$GraphST$GraphST(adata = adata)
                adata = model$train()
                print(adata)
                cluster <- gt$clustering(adata = adata, n_clusters = input$cluster_graphst, method = input$cluster_method_graphst)
                v$isClusterdone2 <- TRUE
                print(adata)
                print(cluster)
                adata$write_csvs(dirname = 'test/')
                meta <- read.csv('test/obs.csv', header = T, row.names = 1)
                v$scData_spatial <- AddMetaData(v$scData_spatial, metadata = meta)
              })
            }
          }
        })
      }) 
         
      output$Cluster2DPlot_graphst <- renderPlot({
        if(is.null(v$isClusterdone2)){
          plotly_empty()
        }else{
          withProgress(message="Generating 2D Cluster Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, group.by = "mclust", label = T)
          })
        }
      })
      
      output$SpatialDimPlot_GraphST <- renderPlot({
        if(is.null(v$scData_spatial)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP 2D Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, group.by = "mclust", label = TRUE)
          })
        }
      })
      
      output$silhouette_graphst <- renderText({
        if(is.null(v$scData_spatial)){
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
            DimPlot(v$scData_xenium, reduction = "pca", label = T)
          })
        }
      })
      
      observeEvent(input$runUMAP_spatial, {
        withProgress(message = "Running UMAP...", value = 0.3, {
          v$scData_spatial <- RunUMAP(v$scData_spatial, dims = 1:input$dim.used_spatial, assay = "SCT", spread = 1)
          v$isUMAPSpataildone <- TRUE
        })
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
        if(is.null(v$scData_spatial) || is.null(v$isUMAPSpataildone)){
          plotly_empty()
        }else{
          withProgress(message="Generating UMAP 2D Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
          })
        }
      })
      
      output$silhouette_seurat <- renderText({
        if(is.null(v$scData_spatial)){
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
          ImageDimPlot(v$scData_xenium, cols = "polychrome", size = 0.75)
        }
      })
      
      observeEvent(input$Vis_sp1, {
        withProgress(message = "Visualizing...", value = 0,{
          incProgress(0.5, message = "Visualizing...")
          DefaultAssay(v$scData_spatial) <- "Spatial"
          #Idents(v$scData_spatial) <- "seurat_clusters"
        })
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
      
      observeEvent(input$Vis_sp2, {
        withProgress(message = "Visualizing...", value = 0,{
          incProgress(0.5, message = "Visualizing...")
          DefaultAssay(v$scData_spatial) <- "Spatial"
          #Idents(v$scData_spatial) <- "seurat_clusters"
        })
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
          tpmFiles_ref <- readRDS("allen_cortex_down.rds")
        })
        if (is.null(tpmFiles_ref)){
          v$scRNAData <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0.5, {
            #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
            label1 <- "Example loaded"
            updateActionButton(inputId = "loadexample_ref", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scRNAData <- tpmFiles_ref
          })
        }
      })
      
      observeEvent(input$loadButton_process_sp, {
        tpmFiles_scRNA <- input$tpmFiles_scRNA
        if (is.null(tpmFiles_scRNA)){
          v$scRNAData <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0, {
            print(tpmFiles_scRNA$datapath)
            print(tpmFiles_scRNA$name)
            print(file.exists(paste(tpmFiles_scRNA$datapath[1], "/", tpmFiles_scRNA$name[1], sep="")))
            exp.data_scRNA <- readRDS(tpmFiles_scRNA$datapath)
            incProgress(0.5, "Creating Seurat Object")
            v$scRNAData <- exp.data_scRNA
            shinyalert("Single cell reference data loaded", "Single cell reference data loaded, please click on Process reference dataset", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
      observeEvent(input$process_scRNA, {
        withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
          incProgress(0.5, message = "Processing...")
          v$scRNAData <- Seurat::SCTransform(v$scRNAData, ncells = 3000, verbose = FALSE)
          v$scRNAData <- Seurat::RunPCA(v$scRNAData, verbose = FALSE)
          v$scRNAData <- Seurat::RunUMAP(v$scRNAData, dims = 1:30, verbose = FALSE)
          v$scRNAData <- Seurat::FindNeighbors(v$scRNAData, dims = 1:30, verbose = FALSE)
          v$scRNAData <- Seurat::FindClusters(v$scRNAData, verbose = FALSE)
          scRNA_plot <- Seurat::DimPlot(v$scRNAData, group.by = "cell_type")
          print (scRNA_plot)
          #output$process_sc.done <- renderText(paste0("Processing of scRNA-seq data done!"))
          v$isProcess_scRNAdone <- TRUE
        })
      })
      
      output$scRNAPlot <- renderPlot({
        if(is.null(v$isProcess_scRNAdone)){
          plotly_empty()
        }else{
          withProgress(message="Generating scRNA-seq Plot...", value=0, {
            DimPlot(v$scRNAData, group.by = "cell_type", label = T)
          })
        }
      })
      
      observeEvent(input$vis_spRNA, {
        withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
          incProgress(0.5, message = "Processing...")
          spRNA_plot <- SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
          print (spRNA_plot)
          #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
          v$isVis_spRNAdone <- TRUE
        })
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
        withProgress(message = "Performing deconvolution...", value = 0,{
          incProgress(0.5, message = "Deconvoluting...")
          v$anchors <- FindTransferAnchors(reference = v$scRNAData, query = v$scData_spatial, normalization.method = "SCT")
          predictions.assay <- TransferData(anchorset = v$anchors, refdata = v$scRNAData$cell_type, prediction.assay = TRUE,
                                            weight.reduction = v$scData_spatial[["pca"]], dims = 1:input$dim.used_sc)
          v$scData_spatial[["predictions"]] <- predictions.assay
          #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
        })
      })
      
      output$ct.select <- renderUI({
        if(is.null(v$scRNAData)){
          return(NULL)
        }else{
          selectInput("ct", label = "Celltype to visualise",
                      choices = unique(v$scRNAData$cell_type))
        }
      })
      
      output$DeconvPlot <- renderPlot({
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating SpatialDim Plot...", value=0, {
            DefaultAssay(v$scData_spatial) <- "predictions"
            SpatialFeaturePlot(v$scData_spatial, features = input$ct, pt.size.factor = 1.6, ncol = 2, crop = TRUE)
          })
        }
      })
      
      observeEvent(input$loadexample_ref1, {
        withProgress(message="Loading example reference scRNA-seq data...", value=0.5, {
          tpmFiles_ref1 <- readRDS("allen_cortex_down.rds")
        })
        if (is.null(tpmFiles_ref1)){
          v$scRNAData <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0.5, {
            #print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
            label1 <- "Example loaded"
            updateActionButton(inputId = "loadexample_ref1", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            v$scRNAData <- tpmFiles_ref1
          })
        }
      })
      
      observeEvent(input$loadButton_process_sp_graphst1, {
        tpmFiles_scRNA_graphst1 <- input$tpmFiles_scRNA_graphst1
        if (is.null(tpmFiles_scRNA_graphst1)){
          v$scRNAData_graphst1 <- NULL
        }else{
          withProgress(message="Loading and Processing Data...", value=0, {
            print(tpmFiles_scRNA_graphst1$datapath)
            print(tpmFiles_scRNA_graphst1$name)
            print(file.exists(paste(tpmFiles_scRNA_graphst1$datapath[1], "/", tpmFiles_scRNA_graphst1$name[1], sep="")))
            exp.data_scRNA_graphst1 <- readRDS(tpmFiles_scRNA_graphst1$datapath)
            incProgress(0.5, "Creating Seurat Object")
            v$scRNAData_graphst1 <- exp.data_scRNA_graphst1
            shinyalert("Single cell reference data loaded", "Single cell reference data loaded, please click on Process reference dataset", type = "success", imageWidth = 10, imageHeight = 10)
          })
        }
      })
      
      observeEvent(input$process_scRNA_graphst1, {
        withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
          incProgress(0.5, message = "Processing...")
          v$scRNAData_graphst1 <- Seurat::SCTransform(v$scRNAData_graphst1, ncells = 3000, verbose = FALSE)
          v$scRNAData_graphst1 <- Seurat::RunPCA(v$scRNAData_graphst1, verbose = FALSE)
          v$scRNAData_graphst1 <- Seurat::RunUMAP(v$scRNAData_graphst1, dims = 1:30, verbose = FALSE)
          v$scRNAData_graphst1 <- Seurat::FindNeighbors(v$scRNAData_graphst1, dims = 1:30, verbose = FALSE)
          v$scRNAData_graphst1 <- Seurat::FindClusters(v$scRNAData_graphst1, verbose = FALSE)
          scRNA_plot_graphst1 <- Seurat::DimPlot(v$scRNAData_graphst1, group.by = "cell_type")
          print (scRNA_plot_graphst1)
          sceasy::convertFormat(v$scRNAData_graphst1, from = "seurat", to = "anndata", outFile = 'reference.h5ad')
          file_path = 'reference.h5ad' # Please replace 'file_path' with the scRNA download path.
          sc <- import("scanpy", convert = FALSE)
          #scvi <- import("scvi", convert = FALSE)
          scipy <- import("scipy", convert = FALSE)
          pot <- import("ot", convert = FALSE)
          pandas <- import("pandas", convert = F)
          gt <- import("GraphST", convert = F)
          matplotlib <- import ("matplotlib", convert = F)
          gt1 <- import("GraphST.preprocess", convert = F)
          gt2 <- import("GraphST.utils", convert = F)
          torch <- import("torch", convert = FALSE)
          np <- import("numpy", convert = F)
          pd <- import("pandas", convert = F)
          sns <- import("seaborn", convert = F)
          v$adata_sc = sc$read(file_path)
          print(v$adata_sc)
          #output$process_sc.done <- renderText(paste0("Processing of scRNA-seq data done!"))
          v$isProcess_scRNAdone1 <- TRUE
        })
      })
      
      output$scRNAPlot_graphst1 <- renderPlot({
        if(is.null(v$isProcess_scRNAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating scRNA-seq Plot...", value=0, {
            DimPlot(v$scRNAData_graphst1, group.by = "cell_type", label = T)
          })
        }
      })
      
      shinyDirChoose(
        input,
        'dir_graphst1',
        roots = c(home = '.'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
      )
      
      dir_graphst1 <- reactive(input$dir_graphst1)
      output$dir_graphst1 <- renderPrint({  # use renderText instead of renderPrint
        parseDirPath(c(home = '.'), dir_graphst1())
      })
      
      observeEvent(input$vis_spRNA_graphst, {
        withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
          incProgress(0.5, message = "Processing...")
          tpmFile_graphst1 <- input$tpmFile_graphst1
          if (is.null(tpmFile_graphst1)){
            v$scData_graphst1 <- NULL
          }else{
            withProgress(message="Loading and Processing Data...", value=0, {
              print(tpmFile_graphst1$datapath)
              print(tpmFile_graphst1$name)
              print(file.exists(paste(tpmFile_graphst1$datapath[1], "/", tpmFile_graphst1$name[1], sep="")))
              file_fold <- parseDirPath(c(home = '.'), dir_graphst1())
              sc <- import("scanpy", convert = FALSE)
              #scvi <- import("scvi", convert = FALSE)
              scipy <- import("scipy", convert = FALSE)
              pot <- import("ot", convert = FALSE)
              pandas <- import("pandas", convert = F)
              gt <- import("GraphST", convert = F)
              matplotlib <- import ("matplotlib", convert = F)
              gt1 <- import("GraphST.preprocess", convert = F)
              gt2 <- import("GraphST.utils", convert = F)
              torch <- import("torch", convert = FALSE)
              np <- import("numpy", convert = F)
              pd <- import("pandas", convert = F)
              sns <- import("seaborn", convert = F)
              v$adata <- sc$read_visium(file_fold, count_file=tpmFile_graphst1$name, load_images=TRUE)
              print(v$adata)
          #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
          v$isVis_spRNAdone1 <- TRUE
            })
          }
        })
      })
      
      output$spRNAPlot_graphst <- renderPlot({
        if(is.null(v$isVis_spRNAdone1)){
          plotly_empty()
        }else{
          withProgress(message="Generating SpatialDim Plot...", value=0, {
            SpatialDimPlot(v$scData_spatial, group.by = "mclust", label = T)
          })
        }
      })
      
      observeEvent(input$doDeconv_graphst, {
        withProgress(message = "Performing deconvolution...", value = 0,{
          incProgress(0.5, message = "Deconvoluting...")
          sc <- import("scanpy", convert = FALSE)
          #scvi <- import("scvi", convert = FALSE)
          scipy <- import("scipy", convert = FALSE)
          pot <- import("ot", convert = FALSE)
          pandas <- import("pandas", convert = F)
          gt <- import("GraphST", convert = F)
          matplotlib <- import ("matplotlib", convert = F)
          gt1 <- import("GraphST.preprocess", convert = F)
          gt2 <- import("GraphST.utils", convert = F)
          torch <- import("torch", convert = FALSE)
          np <- import("numpy", convert = F)
          pd <- import("pandas", convert = F)
          sns <- import("seaborn", convert = F)
          v$adata$var_names_make_unique()
          gt$preprocess(v$adata)
          gt$construct_interaction(v$adata)
          gt$add_contrastive_label(v$adata)
          v$adata_sc$var_names_make_unique()
          gt$preprocess(v$adata_sc)
          v$test = gt1$filter_with_overlap_gene(adata = v$adata, adata_sc = v$adata_sc)
          gt1$get_feature(v$test[0])
          model = gt$GraphST$GraphST(v$test[0], v$test[1], deconvolution = T)
          v$test = model$train_map()
          gt2$project_cell_to_spot(adata = v$test[0], adata_sc = v$test[1], retain_percent = 0.15)
          v$test[0]$write_csvs('ref1/')
          meta1 <- read.csv('ref1/obs.csv', header = T, row.names = 1)
          v$scData_spatial <- AddMetaData(v$scData_spatial, metadata = meta1)
        })
      })
      
      output$graphst_ct.select <- renderUI({
        if(is.null(v$scRNAData_graphst1)){
          return(NULL)
        }else{
          selectInput("ct_graphst", label = "Celltype to visualise",
                      choices = unique(v$scRNAData_graphst1$cell_type))
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
    
   
      
    observeEvent(input$doDeg_spatial, {
        if(is.null(v$scData_spatial)){
          return(NULL)
        }else{
          withProgress(message="Finding DEGs...", value=0, {
            ips.markers_spatial <- FindAllMarkers(v$scData_spatial, only.pos = FALSE, min.pct = input$min_pct_spatial, logfc.threshold = input$logfc_spatial, test.use = input$test.use_spatial)
            v$ips.markers_spatial <- ips.markers_spatial
          })
        }
      })
      
      output$deg.gene.select_spatial <- renderUI({
        if(is.null(v$ips.markers_spatial)){
          return(NULL)
        }else{
          selectInput("deg.gene_spatial", label = "Gene to visualise",
                      choices = rownames(v$ips.markers_spatial))
        }
      })
      
      output$Deg.plot_spatial <- renderPlotly({
        if(is.null(v$ips.markers_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            VlnPlot(v$scData_spatial, input$deg.gene_spatial)
          })
        }
      })
      
      output$Deg1.plot_spatial <- renderPlotly({
        if(is.null(v$ips.markers_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            FeaturePlot(v$scData_spatial, input$deg.gene_spatial)
          })
        }
      })
      
      output$Deg2.plot_spatial <- renderPlotly({
        if(is.null(v$ips.markers_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            RidgePlot(v$scData_spatial, features = input$deg.gene_spatial)
          })
        }
      })
      
      output$Deg.heatmap_spatial <- renderPlotly({
        if(is.null(v$ips.markers_spatial)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            v$ips.markers_spatial %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
            DoHeatmap(v$scData_spatial, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
          })
        }
      })
      
      output$Deg.table_spatial <- DT::renderDataTable(
        v$ips.markers_spatial, options = list(scrollX = TRUE, scrollY = "400px"))
      
      #---------------------Single cell ATAC-seq-----------------##
      
      output$atac_image <- renderImage({
        list(src = "www/atac.png",
             width = 600,
             height = 450)
      }, deleteFile = FALSE)
      
      observeEvent(input$loadButton_atac, {
        tpmFiles_atac <- input$tpmFiles_atac
        meta_atac <- input$meta_atac
        print(tpmFiles_atac)
        print(meta_atac)
        names.field <- input$field
        
        if (is.null(tpmFiles_atac)){
          v$atacData <- NULL
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
          
          chrom_assay_atac <- CreateChromatinAssay(counts = tpmFiles_atac, sep = c(":", "-"), genome = 'hg19', fragments = 'atac_pbmc_500_v1/fragments.tsv.gz', min.cells = 10, min.features = 200)
          
          sObj_atac <- CreateSeuratObject(chrom_assay_atac, assay = "peaks", meta.data = meta_atac)
          print(sObj_atac)
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
      })
      
      observeEvent(input$create_seurat_atac, {
        withProgress(message="Loading and Processing Data...", value=0, {
          print(v$atacData)
          annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
          seqlevelsStyle(annotations_atac) <- 'UCSC'
          genome(annotations_atac) <- "hg19"
          print("test")
          Annotation(v$atacData) <- annotations_atac
          v$atacData <- NucleosomeSignal(object = v$atacData)
          print("test2")
          v$atacData <- TSSEnrichment(object = v$atacData, fast = FALSE)
          print("test1")
          v$atacData$pct_reads_in_peaks <- v$atacData$peak_region_fragments / v$atacData$passed_filters * 100
          v$atacData$blacklist_ratio <- v$atacData$blacklist_region_fragments / v$atacData$peak_region_fragments
          v$atacData$high.tss <- ifelse(v$atacData$TSS.enrichment > 2, 'High', 'Low')
          v$atacData$nucleosome_group <- ifelse(v$atacData$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
          
          print(v$atacData@meta.data)
          shinyalert("Seurat object created", "Seurat object created", type = "success", imageWidth = 10, imageHeight = 10)
        })
      }
      )
      
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
        withProgress(message = "Performing Normalization...", value = 0,{
          incProgress(0.5, message = "Running Dimension Reduction...")
          v$atacData <- RunTFIDF(v$atacData)
          v$atacData <- FindTopFeatures(v$atacData, min.cutoff = 'q0')
          v$atacData <- RunSVD(v$atacData)
          output$normalize_atac.done <- renderText(paste0("Normalization done!"))
          v$isNormalizeATACdone <- TRUE
        })
      })
      
      output$DepthCor_ATAC <- renderPlotly({
        if(is.null(v$isNormalizeATACdone)){
          plotly_empty()
        }else{
          DepthCor(v$atacData)
        }
      })
      
      ##---------------Clustering of ATAC-seq data-------------------
      
      observeEvent(input$doCluster_ATAC, {
        withProgress(message = "Performing clustering...", value = 0.3,{
          incProgress(0.5, message = "Running Clustering...")
          v$atacData <- FindNeighbors(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac)
          v$atacData <- FindClusters(object = v$atacData, resolution = input$clus.res_atac, verbose = FALSE, algorithm = 3)
          output$cluster_atac.done <- renderText(paste0("Clustering done!"))
          v$isClusterATACdone <- TRUE
        })
      })
      
      output$cluster_ATAC <- renderPlotly({
        if(is.null(v$isClusterATACdone)){
          plotly_empty()
        }else{
          DimPlot(object = v$atacData, reduction = "lsi", label = TRUE)
        }
      })
      
      
      ##---------------UMAP on ATAC-seq data-------------------
      
      observeEvent(input$doUMAP_ATAC, {
        withProgress(message = "Running UMAP...", value = 0.3,{
          incProgress(0.5, message = "Running UMAP...")
          v$atacData <- RunUMAP(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac)
          output$umap_atac.done <- renderText(paste0("UMAP done!"))
          v$isUMAPATACdone <- TRUE
        })
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
        withProgress(message = "Running TSNE...", value = 0.3,{
          incProgress(0.5, message = "Running TSNE...")
          v$atacData <- RunTSNE(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac)
          output$tsne_atac.done <- renderText(paste0("TSNE done!"))
          v$isTSNEATACdone <- TRUE
        })
      })
      
      output$Tsne_ATAC <- renderPlotly({
        if(is.null(v$isTSNEATACdone)){
          plotly_empty()
        }else{
          DimPlot(object = v$atacData, reduction = "tsne", label = TRUE)
        }
      })
      
      observeEvent(input$doDeg_atac, {
        if(is.null(v$atacData)){
          return(NULL)
        }else{
          withProgress(message="Finding DEGs...", value=0, {
            DefaultAssay(v$atacData) <- 'peaks'
            atac.markers <- FindAllMarkers(v$atacData, only.pos = FALSE, min.pct = input$min_pct_atac, logfc.threshold = input$logfc_atac, test.use = input$test.use_atac)
            v$atac.markers <- atac.markers
          })
        }
      })
      
      output$Deg_atac.table <- DT::renderDataTable(
        v$atac.markers, options = list(scrollX = TRUE, scrollY = "400px"))
      
      output$Deg_atac.plot <- renderPlotly({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            v$atac.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
            DoHeatmap(v$atacData, features = top10$gene, size = 4, slot = "data") + theme(axis.text.y = element_text(size = 5)) + NoLegend()
          })
        }
      })
      
      output$deg.atac.select <- renderUI({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          selectInput("deg.atac", label = "Gene to visualise",
                      choices = rownames(v$atac.markers))
        }
      })
      
      output$Deg_atac1.plot <- renderPlotly({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            VlnPlot(v$atacData, input$deg.atac)
          })
        }
      })
      
      output$Deg_atac2.plot <- renderPlotly({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            FeaturePlot(v$atacData, input$deg.atac)
          })
        }
      })
      
      output$Deg_atac3.plot <- renderPlot({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            RidgePlot(v$atacData, features = input$deg.atac)
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
      
      output$coverage.plot <- renderPlot({
        if(is.null(v$atac.markers)){
          return(NULL)
        }else{
          withProgress(message="Generating DEG Plot...", value=0, {
            CoveragePlot(
              object = v$atacData,
              region = input$deg1.atac,
              extend.upstream = 40000,
              extend.downstream = 20000
            )
          })
        }
      })
      
      observeEvent(input$doDeg_motif, {
        if(is.null(v$atacData)){
          return(NULL)
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
        if(is.null(v$atacData)){
          return(NULL)
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
