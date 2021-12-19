source('global.R')
source('ui.R')
source('predict.R')

## max data size
#options(shiny.maxRequestSize = 1024^10)
options(shiny.maxRequestSize = 1024*1024*100*100)
options(shiny.launch.browser = T)
options(bitmapType = 'cairo')

shinyServer(function(input, output, session) {
    v <- reactiveValues(scData = NULL,
                        scData1 = NULL,
                        scData2 = NULL,
                        scData3 = NULL,
                        scData4 = NULL,
                        scData5 = NULL,
                        scData6 = NULL,
                        scData7 = NULL,
                        scData8 = NULL,
                        scData9 = NULL,
                        scData10 = NULL,
                        scData11 = NULL,
                        scData12 = NULL,
                        scDatat = NULL,
                        scDatan = NULL,
                        scDatatr = NULL,
                        sceData = NULL,
                        maskData = NULL,
                        atacData = NULL,
                        idents = NULL,
                        isPCAdone = NULL,
                        isUMAPdone = NULL,
                        isTSNEdone = NULL,
                        ismofadone = NULL,
                        ismofaumapdone = NULL,
                        isTrajectorydone = NULL,
                        isCELLiDdone = NULL,
                        isCCdone = NULL,
                        isPCAdone1 = NULL,
                        isUMAPdone1 = NULL,
                        isTSNEdone1 = NULL,
                        isPCAdone2 = NULL,
                        isUMAPdone2 = NULL,
                        isTSNEdone2 = NULL,
                        isPCAdone3 = NULL,
                        isUMAPdone3 = NULL,
                        isTSNEdone3 = NULL,
                        isUMAPdone4 = NULL,
                        isTSNEdone4 = NULL,
                        isVisdone = NULL,
                        isClusterdone = NULL,
                        isNormalizeATACdone = NULL,
                        isClusterATACdone = NULL,
                        isUMAPATACdone = NULL,
                        isTSNEATACdone = NULL,
                        isDataIntegration = NULL,
                        isSpatialClusterdone = NULL,
                        isSpatialCluster1done = NULL,
                        pcGenes = NULL,
                        plotlySelection = NULL,
                        ips.markers = NULL,
                        markers = NULL,
                        selected_markers = NULL,
                        data = NULL,
                        data1 = NULL,
                        sampleInfo = NULL,
                        sample_selected_index = NULL,
                        sample_choices = NULL,
                        sample_init = FALSE,
                        sample_ready = FALSE,
                        sample_update = FALSE)

    inputs = reactiveValues(
        fcsFiles = NULL
        , markers = NULL
        , l_w = 0.5
        , l_t = 500000
        , l_m = 4.5
        , l_a = 0
    )

    observeEvent(input$merge_method, {
        # browser()
        if(input$merge_method == "ceil" || input$merge_method == "fixed"){
            shinyjs::enable("fix_number")
        } else {
            shinyjs::disable('fix_number')
        }
    })

    observeEvent(input$rawfcs, {
        # browser()
        if(input$Module == "Flow cytometry Analysis" & input$demo_analyse_sc == "Analysis"){
            cur_dir = paste0(dirname(input$rawfcs$datapath)[1], "/")
            file.rename(input$rawfcs$datapath, paste0(cur_dir, input$rawfcs$name))
            inputs$fcsFiles = paste0(cur_dir, input$rawfcs$name)
            # input$rawfcs$datapath = paste0(cur_dir, input$rawfcs$name)

            fcs <- suppressWarnings(read.FCS(inputs$fcsFiles[1]))
            pd <- fcs@parameters@data
            markers <- paste(pd$name, "<", pd$desc, ">", sep = "")
            inputs$markers = markers
            updateSelectInput(session, "markers", choices = markers)
            # updatePickerInput(session, 'pickermarker', choices = markers)
            if (length(markers) == 0) {
                shinyalert(title = "Error", text = "No markers found in the FCS file!", type = "error")
                # stop()
            }
        }
    })

    observeEvent(input$submit, {
        if(input$Module == "Flow cytometry Analysis" & input$demo_analyse_sc == "Analysis"){
            # browser()
            ### need to check data
            withProgress(
                {
                    # browser()
                    cytofkit(fcsFiles = inputs[["fcsFiles"]],
                             markers = input$markers,
                             projectName = input$project_name,
                             mergeMethod = input$merge_method,
                             fixedNum = input$fix_number,
                             transformMethod = input$transform_method,
                             dimReductionMethod = "tsne",
                             clusterMethods = input$cluster_method,
                             visualizationMethods = tolower(input$dr_method),
                             progressionMethod = input$progressionMethods,
                             Rphenograph_k = input$rphenograph_k,
                             FlowSOM_k = input$flowsom_k,
                             seed = input$seed,
                             clusterSampleSize = 500,
                             resultDir = tempdir(),
                             saveResults = TRUE,
                             saveObject = TRUE,
                             l_w = as.numeric(inputs[["l_w"]]),
                             l_t = as.numeric(inputs[["l_t"]]),
                             l_m = as.numeric(inputs[["l_m"]]),
                             l_a = as.numeric(inputs[["l_a"]]))
                }

                , value = 0.5, message = "Analysing")
            if(!is.null(input$sample_anno)){
                res_fn = paste0(input$project_name, ".RData")
                load(res_fn)
                meta_data = read.table(input$sample_anno$datapath, sep = "\t", header = TRUE, row.names = 1)
                analysis_results$meta_data = meta_data
                analysis_results$sampel_name_ori = analysis_results$sampleNames
                # analysis_results$sampleInfo_ori = analysis_results$sampleInfo
                save(analysis_results, file = res_fn)
            }
        }
        else if(input$Module == "Flow cytometry Analysis" & input$demo_analyse_sc == "Demo"){
            withProgress(
                {
                    # browser()
                    cytofkit(fcsFiles = v$data1,
                             markers = v$markers,
                             projectName = input$project_name,
                             mergeMethod = input$merge_method,
                             fixedNum = input$fix_number,
                             transformMethod = input$transform_method,
                             dimReductionMethod = "tsne",
                             clusterMethods = input$cluster_method,
                             visualizationMethods = tolower(input$dr_method),
                             progressionMethod = input$progressionMethods,
                             Rphenograph_k = input$rphenograph_k,
                             FlowSOM_k = input$flowsom_k,
                             seed = input$seed,
                             clusterSampleSize = 500,
                             resultDir = tempdir(),
                             saveResults = TRUE,
                             saveObject = TRUE,
                             l_w = as.numeric(inputs[["l_w"]]),
                             l_t = as.numeric(inputs[["l_t"]]),
                             l_m = as.numeric(inputs[["l_m"]]),
                             l_a = as.numeric(inputs[["l_a"]]))
                }

                , value = 0.5, message = "Analysing")
            if(!is.null(input$sample_anno)){
                res_fn = paste0(input$project_name, ".RData")
                load(res_fn)
                meta_data = read.table(input$sample_anno$datapath, sep = "\t", header = TRUE, row.names = 1)
                analysis_results$meta_data = meta_data
                analysis_results$sampel_name_ori = analysis_results$sampleNames
                # analysis_results$sampleInfo_ori = analysis_results$sampleInfo
                save(analysis_results, file = res_fn)
            }
        }
    })

    output$download_analysis_res = downloadHandler(
        filename = function() {
            paste0(input$project_name, '.zip')
        },
        content = function(file) {
            # browser()
            cur_dir = tempdir()
            files = list.files(path = cur_dir, pattern = input$project_name, full.names = TRUE)
            utils::zip(file, files)
            #stopApp(returnValue = invisible())
        }
    )


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

    ##------------------Reactive Values and Reactive Objects-------------------

    #if?
    c <- reactiveValues(clusterCol = list())
    p <- reactiveValues(progressionCluster = NULL)

    #parseQueryString to use .RData path as analysis results
    output$queryText <- renderText({
        query <- parseQueryString(session$clientData$url_search)
        if(!length(query) == 0){
            load(query[["fcspath"]])
            if(exists("analysis_results")){
                if(!is.null(analysis_results)) {
                    v$data <- analysis_results
                    v$sampleInfo <- data.frame(cellID = row.names(analysis_results$expressionData),
                                               cellSample = factor(sub("_[0-9]*$", "", row.names(analysis_results$expressionData))),
                                               stringsAsFactors = FALSE)
                    p$progressionCluster <- names(analysis_results$clusterRes)[1]
                    paste0("Loaded: ", query[["fcspath"]])
                }
            }
        }else{
            return(NULL)
        }
    })

    ## Scatter plot methods
    visualizationMethods <- reactive({
        if(is.null(v$data) || is.null(v$data$visualizationMethods)){
            return(NULL)
        }else{
            return(v$data$visualizationMethods)
        }
    })

    ## Scatter plot functions
    visualizationFunctions <- reactive({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            return(c(names(v$data$clusterRes),
                     "Sample",
                     "Density",
                     "None"))
        }
    })

    ## cluster methods
    clusterMethods <- reactive({
        if(is.null(v$data))
            return(NULL)
        cMethods <- names(v$data$clusterRes)
        return(cMethods)
    })

    ## progression labs
    progressionLabs <- reactive({
        if(is.null(v$data))
            return(NULL)
        if(is.null(v$data$progressionRes))
            return(NULL)
        progressionLabs <- colnames(v$data$progressionRes[[3]])
        return(progressionLabs)
    })

    ##-------------------Side Panel-------------------

    normMethod <- NULL

    output$name.field <- renderUI({
        if(is.null(input$cellAnnoFiles)){
            numericInput(inputId = "field",
                         label = "Field",
                         value = 1,
                         min = 1)
        }else{
            annoFile <- input$cellAnnoFiles
            anno.data <- read.table(annoFile$datapath[1], header = T,
                                    sep = "\t", stringsAsFactors = FALSE)
            groupings <- colnames(anno.data)
            selectInput("groupby",
                        label = "Group by:",
                        choices = groupings)
        }
    })

    observeEvent(input$loadexample_tpm_human, {
        if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "Raw Counts Matrix" & input$Species == "human"){
            withProgress(message="Loading and Processing Data...", value=0.5, {
                tpmFiles_demo <- read.table("TPM.txt", header = T, row.names = 1, check.names = F)
                #scH5 <- input$scH5
                annoFile_demo <- read.table("metadata_test.txt", header = T, row.names = 1)
            })
            if (is.null(tpmFiles_demo)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0.5, {
                    print(tpmFiles_demo)
                    #print(tpmFiles$name)
                    print(file.exists(paste(tpmFiles_demo$datapath[1], "/", tpmFiles_demo$name[1], sep="")))
                    v$tpmFiles_demo <- tpmFiles_demo
                    v$annoFile_demo <- annoFile_demo
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadexample_tpm_human", label = label1)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                    #withProgress(message="Creating Seurat Object...", value=0.5, {

                    sObj <- CreateSeuratObject(v$tpmFiles_demo,
                                               meta.data = v$annoFile_demo,
                                               project = input$projName,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)


                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
                    v$scData <- sObj
                })
            }
        }
    })


    observeEvent(input$loadexample_scH5_human, {
        if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "H5" & input$Species == "human"){
            withProgress(message="Loading and Processing Data...", value=0.5, {
                scH5_demo <- Read10X_h5('filtered_feature_bc_matrix_human.h5')
            })
            if (is.null(scH5_demo)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(scH5_demo)
                    #print(scH5$name)
                    #print(file.exists(paste(scH5$datapath[1], "/", scH5$name[1], sep="")))
                    v$scH5_demo <- scH5_demo
                    label2 <- "Example loaded"
                    updateActionButton(inputId = "loadexample_scH5_human", label = label2)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)

                    sObj <- CreateSeuratObject(v$scH5_demo,
                                               project = input$projName,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)
                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
                    v$scData <- sObj
                })
            }
        }
    })

    observeEvent(input$loadexample_tpm_mouse, {
        if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "Raw Counts Matrix" & input$Species == "mouse"){
            withProgress(message="Loading and Processing Data...", value=0.5, {
                tpmFiles_demo <- read.table("TPM_mouse_21June2016.txt", header = T, row.names = 1, check.names = F)
                annoFile_demo <- read.table("mouse_sample_25July2016.txt", header = T, row.names = 1)
            })
            #names.field <- input$field
            if (is.null(tpmFiles_demo)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles_demo)
                    #print(tpmFiles$name)
                    #print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
                    v$tpmFiles_demo <- tpmFiles_demo
                    v$annoFile_demo <- annoFile_demo
                    label3 <- "Example loaded"
                    updateActionButton(inputId = "loadexample_tpm_mouse", label = label3)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                    #withProgress(message="Creating Seurat Object...", value=0.5, {

                    sObj <- CreateSeuratObject(v$tpmFiles_demo,
                                               meta.data = v$annoFile_demo,
                                               project = input$projName,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)

                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^mt-")
                    v$scData <- sObj
                })
            }
        }
    })

    observeEvent(input$loadexample_scH5_mouse, {
        if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "H5" & input$Species == "mouse"){
            withProgress(message="Loading and Processing Data...", value=0.5, {
                scH5_demo <- Read10X_h5('filtered_feature_bc_matrix_mouse.h5')
            })
            if (is.null(scH5_demo)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(scH5_demo)
                    #print(scH5$name)
                    #print(file.exists(paste(scH5$datapath[1], "/", scH5$name[1], sep="")))
                    v$scH5_demo <- scH5_demo
                    label4 <- "Example loaded"
                    updateActionButton(inputId = "loadexample_scH5_mouse", label = label4)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                    sObj <- CreateSeuratObject(v$scH5_demo,
                                               project = input$projName,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)
                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^mt-")
                    v$scData <- sObj
                })
            }
        }
    })

    observeEvent(input$loadButton, {
        if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "Raw Counts Matrix" & input$Species == "human"){
            tpmFiles <- input$tpmFiles
            #scH5 <- input$scH5
            annoFile <- input$cellAnnoFiles
            names.field <- input$field
            if (is.null(tpmFiles)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles$datapath)
                    print(tpmFiles$name)
                    print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
                    exp.data <- read.table(tpmFiles$datapath,
                                           sep="\t", header=TRUE, row.names=1, check.names = F, stringsAsFactors = FALSE)

                    if(!is.null(annoFile)){
                        anno.data <- read.table(annoFile$datapath[1], header = T,
                                                sep = "\t", stringsAsFactors = FALSE, row.names=1)

                    }
                    incProgress(0.5, "Creating Seurat Object")

                    sObj <- CreateSeuratObject(exp.data,
                                               meta.data = anno.data,
                                               project = input$projName,
                                               names.field = names.field,
                                               names.delim = input$delim,
                                               is.expr = input$expThres,
                                               normalization.method = normMethod,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)

                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
                    v$scData <- sObj
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
        else if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "H5" & input$Species == "human"){
            scH5 <- input$scH5
            names.field <- input$field
            if (is.null(scH5)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(scH5$datapath)
                    print(scH5$name)
                    print(file.exists(paste(scH5$datapath[1], "/", scH5$name[1], sep="")))
                    exp.data <- Read10X_h5(scH5$datapath)

                    incProgress(0.5, "Creating Seurat Object")

                    sObj <- CreateSeuratObject(exp.data,
                                               project = input$projName,
                                               names.field = names.field,
                                               names.delim = input$delim,
                                               is.expr = input$expThres,
                                               normalization.method = normMethod,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)
                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
                    v$scData <- sObj
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
        else if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "Raw Counts Matrix" & input$Species == "mouse"){
            tpmFiles <- input$tpmFiles
            #scH5 <- input$scH5
            annoFile <- input$cellAnnoFiles
            names.field <- input$field
            if (is.null(tpmFiles)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles$datapath)
                    print(tpmFiles$name)
                    print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
                    exp.data <- read.table(tpmFiles$datapath,
                                           sep="\t", header=TRUE, row.names=1, check.names = F, stringsAsFactors = FALSE)



                    if(!is.null(annoFile)){
                        anno.data <- read.table(annoFile$datapath[1], header = T,
                                                sep = "\t", stringsAsFactors = FALSE, row.names=1)

                    }
                    incProgress(0.5, "Creating Seurat Object")

                    sObj <- CreateSeuratObject(exp.data,
                                               meta.data = anno.data,
                                               project = input$projName,
                                               names.field = names.field,
                                               names.delim = input$delim,
                                               is.expr = input$expThres,
                                               normalization.method = normMethod,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)

                    #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^mt-")
                    #incProgress(0.5, "Adding metadata")
                    #sObj <- AddMetaData(sObj, percent.mt, "percent.mt")
                    #if(!is.null(additional.ident1)){
                    #  sObj1 <- AddMetaData(sObj1, additional.ident1)
                    #}
                    v$scData <- sObj
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }

        else if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "H5" & input$Species == "mouse"){
            scH5 <- input$scH5
            names.field <- input$field
            if (is.null(scH5)){
                v$scData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(scH5$datapath)
                    print(scH5$name)
                    print(file.exists(paste(scH5$datapath[1], "/", scH5$name[1], sep="")))
                    exp.data <- Read10X_h5(scH5$datapath)

                    incProgress(0.5, "Creating Seurat Object")

                    sObj <- CreateSeuratObject(exp.data,
                                               project = input$projName,
                                               names.field = names.field,
                                               names.delim = input$delim,
                                               is.expr = input$expThres,
                                               normalization.method = normMethod,
                                               min.genes = input$min.genes,
                                               min.cells = input$min.cells)
                    sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^mt-")
                    v$scData <- sObj
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }


    })

    dir.create("Seurat_results")
    #})

    observeEvent(input$loadexample1, {
        if(input$Module == "Single cell data integration Analysis" & input$scInput == "Raw Counts Matrix" & input$Species == "human"){
            tpmFiles1 <- read.table("dataset1_sm_uc3.txt", header = T, row.names = 1, check.names = F)
            #scH5 <- input$scH5
            annoFile1 <- read.table("sample_sm_uc3.txt", header = T, row.names = 1)
            if (is.null(tpmFiles1)){
                v$scData1 <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles1)
                    #print(tpmFiles$name)
                    print(file.exists(paste(tpmFiles1$datapath[1], "/", tpmFiles1$name[1], sep="")))
                    v$tpmFiles1 <- tpmFiles1
                    v$annoFile1 <- annoFile1

                    sObj1 <- CreateSeuratObject(v$tpmFiles1,
                                                meta.data = v$annoFile1,
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
        }
    })

    observeEvent(input$loadexample2, {
        if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "Seurat" & input$scInput2 == "H5" & input$scAnalysis_type == "CITE-seq"){
            tpmFiles2 <- Read10X_h5('pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5')
            #scH5 <- input$scH5
            if (is.null(tpmFiles2)){
                v$scDatat <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles2)
                    #print(tpmFiles2$name)
                    print(file.exists(paste(tpmFiles2$datapath[1], "/", tpmFiles2$name[1], sep="")))
                    v$scDatat <- tpmFiles2

                    exp.data3.rna <- v$scDatat$`Gene Expression`
                    exp.data3.ab <- v$scDatat$`Antibody Capture`
                    #additional.ident <- NULL
                    incProgress(0.5, "Creating Seurat Object")
                    sObj3 <- CreateSeuratObject(exp.data3.rna)
                    #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
                    sObj3[["percent.mt"]] <- PercentageFeatureSet(sObj3, pattern = "^MT-")
                    v$scDatat <- sObj3
                    print(sObj3@meta.data)
                    v$scDatab <- CreateAssayObject(counts = exp.data3.ab)
                    v$scDatat[["ADT"]] <- v$scDatab
                    print(Assays(v$scDatat))

                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadexample2", label = label1)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }

        else if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "MOFA" & input$scInput2 == "H5" & input$scAnalysis_type == "CITE-seq"){
            tpmFiles2 <- Read10X_h5('filtered_feature_bc_matrix_mofa.h5')
            meta_mofa <- read.csv('metadata_mofa.csv', header = T, row.names = 1)
            #scH5 <- input$scH5
            if (is.null(tpmFiles2)){
                v$scDatat <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles2)
                    #print(tpmFiles2$name)
                    print(file.exists(paste(tpmFiles2$datapath[1], "/", tpmFiles2$name[1], sep="")))
                    v$scDatat <- tpmFiles2
                    v$meta_mofa <- meta_mofa

                    exp.data3.rna <- v$scDatat$`Gene Expression`
                    exp.data3.ab <- v$scDatat$`Antibody Capture`
                    anno.data_mofa <- v$meta_mofa

                    incProgress(0.5, "Creating Seurat Object")

                    sObj3 <- CreateSeuratObject(exp.data3.rna,
                                                meta.data = anno.data_mofa,
                                                project = input$projName,
                                                min.genes = input$min.genes,
                                                min.cells = input$min.cells)

                    #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
                    sObj3[["percent.mt"]] <- PercentageFeatureSet(sObj3, pattern = "^MT-")
                    print(sObj3@meta.data)
                    #incProgress(0.5, "Adding metadata")
                    #sObj <- AddMetaData(sObj, percent.mt, "percent.mt")
                    #if(!is.null(additional.ident1)){
                    #  sObj1 <- AddMetaData(sObj1, additional.ident1)
                    #}
                    v$scDatat <- sObj3
                    print(v$scDatat@meta.data)
                    v$scDatab <- CreateAssayObject(counts = exp.data3.ab)
                    v$scDatat[["ADT"]] <- v$scDatab
                    print(v$scDatat@meta.data)
                    print(Assays(v$scDatat))
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadexample1", label = label1)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
    })

    observeEvent(input$loadexample3, {
        if(input$Module == "Flow cytometry Analysis" & input$demo_analyse_sc == "Demo"){
            fcs_file1 <- read.FCS('S4_PBMC.fcs')
            fcs_example <- list.files(path = ".", pattern = ".fcs")
            pd <- fcs_file1@parameters@data
            markers <- paste(pd$name, "<", pd$desc, ">", sep = "")
            v$markers = markers[21:27]
            updateSelectInput(session, "markers", choices = markers)
            #scH5 <- input$scH5
            if (is.null(fcs_example)){
                v$data1 <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    v$data1 <- fcs_example
                    print(v$data1)
                    print(v$markers)
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadexample3", label = label1)
                    label2 <- "Markers Loaded"
                    updateSelectInput(inputId = "markers", label = label2, selected = "(Sm150)Di<GranzymeB>")
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
    })

    observeEvent(input$loadexample4, {
        if(input$Module == "Imaging mass cytometry Analysis"){
            sce_file <- readRDS('pancreas_sce.rds')
            mask_file <- readRDS('pancreas_masks.rds')
            if (is.null(sce_file)){
                v$sceData <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(sce_file)
                    v$sceData <- sce_file
                    v$maskData <- mask_file
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadexample4", label = label1)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
    })

    observeEvent(input$loadexample6, {
        if(input$Module == "Single cell ATAC-seq Analysis"){

            tpmFiles_atac <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
            meta_atac <- read.csv('singlecell.csv', header = T)

            chrom_assay_atac <- CreateChromatinAssay(counts = tpmFiles_atac, sep = c(":", "-"), genome = 'hg19', fragments = 'fragments.tsv.gz', min.cells = 10, min.features = 200)

            sObj_atac <- CreateSeuratObject(chrom_assay_atac, assay = "peaks", meta.data = meta_atac)

            annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
            seqlevelsStyle(annotations_atac) <- 'UCSC'
            genome(annotations_atac) <- "hg19"
            Annotation(sObj_atac) <- annotations_atac
            sObj_atac <- NucleosomeSignal(object = sObj_atac)
            sObj_atac <- TSSEnrichment(object = sObj_atac, fast = FALSE)
            sObj_atac$pct_reads_in_peaks <- sObj_atac$peak_region_fragments / sObj_atac$passed_filters * 100
            sObj_atac$blacklist_ratio <- sObj_atac$blacklist_region_fragments / sObj_atac$peak_region_fragments
            sObj_atac$high.tss <- ifelse(sObj_atac$TSS.enrichment > 2, 'High', 'Low')
            sObj_atac$nucleosome_group <- ifelse(sObj_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

            v$atacData <- sObj_atac
            label1 <- "Example loaded"
            updateActionButton(inputId = "loadexample6", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        }
    })

    observeEvent(input$loadexample5, {
        if(input$Module == "Spatial Transcriptomics Analysis"){

            tpmFiles3 <- Load10X_Spatial(data.dir = "stxBrain", filename = "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5")
            v$scData_spatial <- tpmFiles3
            label1 <- "Example loaded"
            updateActionButton(inputId = "loadexample5", label = label1)
            shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
        }
    })

    observeEvent(input$loadexample2_b, {
        if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "Seurat" & input$scInput2 == "H5" & input$scAnalysis_type == "mRNA+scATAC-seq"){
            tpmFiles3b <- Read10X_h5('pbmc_unsorted_3k_filtered_feature_bc_matrix.h5')
            #dir_multi_atac1 <- '/acrc_raman/jinmiao/CJM_lab/Raman/Projects/hyperion_cytofkit2/ShinyApps/ez_single_cell/scATAC-seq_multiomics/pbmc_unsorted_3k_atac_fragments.tsv.gz'
            #print(dir_multi_atac1)
            #scH5 <- input$scH5
            if (is.null(tpmFiles3b)){
                v$scDatan <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles3b)
                    #print(tpmFiles2$name)
                    print(file.exists(paste(tpmFiles3b$datapath[1], "/", tpmFiles3b$name[1], sep="")))
                    v$scDatan <- tpmFiles3b
                    exp.data4.rna <- v$scDatan$`Gene Expression`
                    exp.data4.atac <- v$scDatan$Peaks
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
                        fragments = 'pbmc_unsorted_3k_atac_fragments.tsv.gz')
                    print(chrom_assay)
                    sObj4[["ATAC"]] <- chrom_assay
                    print(sObj4)
                    v$scDatan <- sObj4

                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadexample2_b", label = label1)
                    shinyalert("Loaded", "Example loaded.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
    })


    observeEvent(input$loadButton1, {
        if(input$Module == "Single cell data integration Analysis" & input$scInput1 == "Raw Counts Matrix" & input$Species == "human"){
            tpmFiles1 <- input$tpmFiles1
            #scH5 <- input$scH5
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
                            read.table (header = T, row.names = 1) %>%
                            CreateSeuratObject(project=y)
                    })
                    print(exp.data1)
                    exp.data2 <- merge(exp.data1[[1]], exp.data1[2:length(exp.data1)])
                    print(exp.data2)


                    if(!is.null(annoFile1)){

                        anno.data1 <- rbindlist(lapply(annoFile1$datapath, fread), use.names = TRUE, fill = TRUE)
                        anno.data2 <- data.frame(anno.data1, row.names = 1)

                    }
                    print(anno.data1)

                    incProgress(0.5, "Creating Seurat Object")

                    sObj1 <- AddMetaData(exp.data2, anno.data2)

                    v$scData1 <- sObj1

                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadButton1", label = label1)
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
        else if(input$Module == "Single cell data integration Analysis" & input$scInput1 == "H5" & input$Species == "human"){
            scH5_1 <- input$scH5_1
            #annoFile1 <- input$cellAnnoFiles1
            names.field <- input$field
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
                    exp.data2[["batch"]] <- ifelse(endsWith(exp.data2@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
                    print(exp.data2)

                    incProgress(0.5, "Creating Seurat Object")

                    exp.data2[["percent.mt"]] <- PercentageFeatureSet(exp.data2, pattern = "^MT-")
                    v$scData1 <- exp.data2
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadButton1", label = label1)
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
        else if(input$Module == "Single cell data integration Analysis" & input$scInput1 == "Raw Counts Matrix" & input$Species == "mouse"){
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
                            read.table (header = T, row.names = 1) %>%
                            CreateSeuratObject(project=y)
                    })
                    print(exp.data1)
                    exp.data2 <- merge(exp.data1[[1]], exp.data1[2:length(exp.data1)])
                    print(exp.data2)


                    if(!is.null(annoFile1)){

                        anno.data1 <- rbindlist(lapply(annoFile1$datapath, fread), use.names = TRUE, fill = TRUE)
                        anno.data2 <- data.frame(anno.data1, row.names = 1)

                    }
                    print(anno.data1)

                    incProgress(0.5, "Creating Seurat Object")

                    sObj1 <- AddMetaData(exp.data2, anno.data2)

                    sObj1[["percent.mt"]] <- PercentageFeatureSet(sObj1, pattern = "^mt-")
                    v$scData1 <- sObj1
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadButton1", label = label1)
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }

        else if(input$Module == "Single cell data integration Analysis" & input$scInput1 == "H5" & input$Species == "mouse"){
            scH5_1 <- input$scH5_1
            #annoFile1 <- input$cellAnnoFiles1
            names.field <- input$field
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
                    exp.data2[["batch"]] <- ifelse(endsWith(exp.data2@assays$RNA@data@Dimnames[[2]], "1"), "1", "2")
                    print(exp.data2)

                    incProgress(0.5, "Creating Seurat Object")

                    exp.data2[["percent.mt"]] <- PercentageFeatureSet(exp.data2, pattern = "^mt-")
                    v$scData1 <- exp.data2
                    label1 <- "Example loaded"
                    updateActionButton(inputId = "loadButton1", label = label1)
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
    })

    dir.create("Data_integration_results")

    observeEvent(input$loadButton_atac, {
        tpmFiles_atac <- input$tpmFiles_atac
        meta_atac <- input$meta_atac
        print(tpmFiles_atac)
        print(meta_atac)
        #frag_atac <- input$frag_atac
        #print(frag_atac)
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
                    anno.data_atac <- read.csv(meta_atac$datapath[1], header = T)
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
                chrom_assay_atac <- CreateChromatinAssay(counts = exp.data_atac, sep = c(":", "-"), genome = 'hg19', fragments = as.character(parseFilePaths(c(home = '.'), dir_atac())$datapath), min.cells = input$min.cells_atac, min.features = input$min.features_atac)
                print (chrom_assay_atac)
                incProgress(0.5, "Creating Seurat Object")
                sObj_atac <- CreateSeuratObject(chrom_assay_atac,
                                                assay = "peaks",
                                                meta.data = anno.data_atac)
                print(sObj_atac)
                print(sObj_atac[['peaks']])
                print(granges(sObj_atac))
                #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
                annotations_atac <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
                seqlevelsStyle(annotations_atac) <- 'UCSC'
                genome(annotations_atac) <- "hg19"
                Annotation(sObj_atac) <- annotations_atac
                sObj_atac <- NucleosomeSignal(object = sObj_atac)
                sObj_atac <- TSSEnrichment(object = sObj_atac, fast = FALSE)
                sObj_atac$pct_reads_in_peaks <- sObj_atac$peak_region_fragments / sObj_atac$passed_filters * 100
                sObj_atac$blacklist_ratio <- sObj_atac$blacklist_region_fragments / sObj_atac$peak_region_fragments
                sObj_atac$high.tss <- ifelse(sObj_atac$TSS.enrichment > 2, 'High', 'Low')
                sObj_atac$nucleosome_group <- ifelse(sObj_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
                #sObj <- AddMetaData(sObj, percent.mito, "percent.mito")
                #if(!is.null(additional.ident)){
                #    sObj <- AddMetaData(sObj, additional.ident)
                #}
                v$atacData <- sObj_atac
            })
        }
        dir.create("Seurat_ATACseq_results")
    })

    observeEvent(input$reset1, {
        session$reload()
        print("Reset done")
    })

    observeEvent(input$saveButton, {
        if(!is.null(input$tpmFiles)){
            withProgress(message="Saving Results...", value=0, {
                print(getwd())
                dir.create("Seurat_results")
                resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
                filename <- paste0(resultDir, .Platform$file.sep, v$scData@project.name, "_", Sys.Date())
                sObj <- v$scData
                save(sObj, file= paste0(resultDir, .Platform$file.sep, sObj@project.name, "_", Sys.Date(), ".Robj"))
            })
            ## open the results directory
            opendir(resultDir)
        }
    })

    #output$ident.swap <- renderUI({
    #  if(is.null(v$scData)){
    #    return(NULL)
    #  }else{
    #    groupings1 <- names(v$scData@meta.data[,!names(v$scData@meta.data) %in% c("nFeature_RNA", "nCount_RNA", "percent.mt")])
    #    tagList(
    #      h4("Set current identity:"),
    #      fluidRow(
    #        column(6,
    #               selectInput("active.ident", label = NULL,
    #                           choices = groupings1)
    #        ),
    #        column(6,
    #               actionButton("swap.ident",label = NULL, icon = icon("arrow-right"))
    #        )
    #      )
    #      )
    #    }
    #  })

    # observeEvent(input$swap.ident, {
    #   v$scData <- SetIdent(v$scData, value = as.character(v$scData@meta.data[,input$active.ident]))
    # })

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

    observeEvent(input$OpenDir, {
        resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
        if(!dir.exists(resultDir)){
            dir.create("Seurat_results")
        }
        pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
        if(dir.exists(pdfDir)){
            opendir(pdfDir)
        }else{
            warning("No reports created yet!")
            dir.create(pdfDir)
        }
    })

    ##---------------QC tabset-------------------

    #observe({if(input$Module == "Single Cell RNASeq Analysis" & input$scInput == "Raw Counts Matrix" & input$Species == "human" & input$demo_analyse_sc == "Analysis"){
    output$nFeature_RNAPlot <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            VlnPlot(v$scData, "nFeature_RNA")
        }
    })

    output$mitoPlot <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            VlnPlot(v$scData, "percent.mt")
        }
    })

    output$nCount_RNAPlot <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            VlnPlot(v$scData, "nCount_RNA")
        }
    })

    output$name <- renderPrint({
        s <- event_data("plotly_selected")
        c(s[["key"]], class(s[["key"]]))
    })

    observeEvent(input$PDFa, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_violin_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "QC_violin_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                nG <- VlnPlot(v$scData, "nFeature_RNA")
                pM <- VlnPlot(v$scData, "percent.mt")
                nU <- VlnPlot(v$scData, "nCount_RNA")
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(nG)
                print(pM)
                print(nU)
                dev.off()
            })
        }
    })

    ## FeatureScatter plot

    output$FeatureScatterPlot1 <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            print(FeatureScatter(v$scData, "nCount_RNA", "nFeature_RNA"))
        }
    })

    output$FeatureScatterPlot2 <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            print(FeatureScatter(v$scData, "nCount_RNA", "percent.mt"))
        }
    })

    observeEvent(input$PDFb, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_featurescatter_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "QC_featurescatter_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                feature1 <- FeatureScatter(v$scData, "nCount_RNA", "nFeature_RNA")
                feature2 <- FeatureScatter(v$scData, "nCount_RNA", "percent.mt")
                print(feature1)
                print(feature2)
                dev.off()
            })
        }
    })

    ##---------------Normalization and Variable Genes tabset-------------------

    observeEvent(input$findVarGenes, {
        withProgress(message = "Finding variable genes...", value = 0, {
            if(input$norm1 == "LogNormalize"){
                v$scData <- NormalizeData(v$scData, normalization.method = "LogNormalize")
                v$scData <- FindVariableFeatures(v$scData,
                                                 mean.function = ExpMean,
                                                 dispersion.function = LogVMR,
                                                 nfeatures = input$var.genes,
                                                 selection.method = input$selection.method)
                #all.genes <- rownames(v$scData)
                v$scData <- ScaleData(v$scData)
                incProgress(0.5)
                VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
                output$nVarGenes <- renderText(VarGeneText)
                varGenePlotInput <- function(){
                    if(is.null(v$scData)){
                        return(NULL)
                    }else{
                        withProgress(message="Plotting variable genes...", value=0, {
                            top10 <- head(VariableFeatures(v$scData), 10)
                            variable_feature1 <- VariableFeaturePlot(v$scData)
                            variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
                            print (variable_feature1)
                            print (variable_feature2)
                            #dev.off()
                        })
                    }
                }
                output$VarGenes <- renderPlot({
                    varGenePlotInput()
                }, height = 800, width = 850)
                observeEvent(input$PDFc, {
                    if(!is.null(v$scData)){
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
                            plot1 <- VariableFeaturePlot(v$scData)
                            print(plot1)
                            dev.off()
                            txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
                            txtfile <- sub(".pdf", ".txt", txtfile)
                            write(v$scData@assays$RNA@var.features, file = txtfile)

                        })
                    }
                })
            }

            else if(input$norm1 == "SCTransform"){
                v$scData <- SCTransform(v$scData, variable.features.n = input$var.genes, vars.to.regress = "percent.mt", verbose = FALSE, conserve.memory = T)
                incProgress(0.5)
                #VarGeneText <- paste0("Number of variable genes: ", length(v$scData@assays$RNA@var.features))
                #output$nVarGenes <- renderText(VarGeneText)
                varGenePlotInput <- function(){
                    if(is.null(v$scData)){
                        return(NULL)
                    }else{
                        withProgress(message="Plotting variable genes...", value=0, {
                            top10 <- head(VariableFeatures(v$scData), 10)
                            variable_feature1 <- VariableFeaturePlot(v$scData)
                            variable_feature2 <- LabelPoints(plot = variable_feature1, points = top10, repel = TRUE)
                            print (variable_feature1)
                            print (variable_feature2)
                        })
                    }
                }



                output$VarGenes <- renderPlot({
                    varGenePlotInput()
                }, height = 800, width = 850)
                observeEvent(input$PDFc, {
                    if(!is.null(v$scData)){
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
                            plot1 <- VariableFeaturePlot(v$scData)
                            print(plot1)
                            dev.off()
                            txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
                            txtfile <- sub(".pdf", ".txt", txtfile)
                            write(v$scData@assays$RNA@var.features, file = txtfile)

                        })
                    }
                })
            }
        })
    })


    ##---------------PCA tabset-------------------
    # PCA plot
    observeEvent(input$doPCA, {
        withProgress(message = "Scaling Data...", value = 0,{
            incProgress(0.5, message = "Running PCA...")
            v$scData <- RunPCA(v$scData, features = VariableFeatures(object = v$scData), assay = input$assays1)
            print(v$scData[["pca"]], dims = 1:5, nfeatures = 5)
            v$isPCAdone <- TRUE
            PCA_plot <- DimPlot(v$scData, reduction = "pca", label = T)
            print(PCA_plot)
            incProgress(0.4, message = "Getting list of PC genes...")
            pc.table <- list()
            for(i in 1:20){
                pcg <- TopFeatures(v$scData)
                pc.table[[i]] <- pcg
            }
            pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
            v$pcGenes <- pc.table
        })
    })

    output$PCA2DPlot <- renderPlotly({
        if(is.null(v$isPCAdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating 2D PCA Plot...", value=0, {
                DimPlot(v$scData, reduction = "pca", label = T)
            })
        }
    })

    #output$PCA3DPlot <- renderPlotly({
    #    if(is.null(v$isPCAdone)){
    #        plotly_empty()
    #    }else{
    #        withProgress(message="Generating 3D PCA Plot...", value=0, {
    #            DimPlot(v$scData, reduction = "pca", label = T)
    #        })
    #    }
    #})

    observeEvent(input$PDFd, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"PCA_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "PCA_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pcaplot <- DimPlot(v$scData, reduction = "pca", label = T)
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(pcaplot)
                dev.off()
            })
            withProgress(message="Downloading PCA coordinates...", value=0.5, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"pca_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "pca_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData@reductions$pca@cell.embeddings, file = filename2)
            })
            withProgress(message="Downloading cluster IDs...", value=0.9, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData@active.ident, file = filename2)
            })
        }
    })

    # Viz plot

    output$vizPlot <- renderPlot({
        if(is.null(v$scData)){
            return(NULL)
        }else{
            VizDimLoadings(v$scData, dims = as.numeric(input$select.pc))
        }
    })

    output$PCHeatmap <- renderPlot({
        if(is.null(v$scData)){
            return(NULL)
        }else{
            DimHeatmap(v$scData, dims = as.numeric(input$select.pc))
        }
    })

    output$PCtable <- DT::renderDataTable({
        if(is.null(v$scData) ){
            return(NULL)
        }else{
            v$pcGenes
        }
    }, options = list(scrollX = TRUE))

    observeEvent(input$PDFe, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Viz_Heatmap_plots_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "Viz_Heatmap_plots_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                #isolate({
                vizdim <- VizDimLoadings(v$scData, dims = as.numeric(input$select.pc))
                dimheat <- DimHeatmap(v$scData)
                print (vizdim)
                print (dimheat)
                #})
                dev.off()
                pcGenes <- v$pcGenes
                write.csv(v$pcGenes, file = paste0(pdfDir, .Platform$file.sep,"PC_genes_", Sys.Date(), ".csv"))
            })
        }
    })

    # Elbow
    output$Elbow <- renderPlot({
        if(is.null(v$scData@reductions$pca@jackstraw)){
            return(NULL)
        }else{
            withProgress(message="Generating Elbow Plot...", value=0.5, {
                ElbowPlot(v$scData)
            })
        }
    }, height = 800, width = 850)

    observeEvent(input$PDFh, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Elbow_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"Elbow_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(ElbowPlot(v$scData))
                dev.off()
            })
        }
    })

    ##---------------Clustering-------------------

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
            DefaultAssay(v$scData) <- input$assays1
            v$scData <- FindNeighbors(v$scData, dims = 1:input$dim.used, assay = input$assays1, nn.method = "rann")
            v$scData <- FindClusters(v$scData, resolution = input$clus.res)
            output$cluster.done <- renderText(paste0("Clustering done!"))
            v$isClusterdone <- TRUE
        })
    })

    output$Cluster2DPlot_1 <- renderPlotly({
        if(is.null(v$isClusterdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating 2D Cluster Plot...", value=0, {
                DimPlot(v$scData, reduction = "pca", label = T)
            })
        }
    })

    observeEvent(input$PDFf, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Cluster_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"Cluster_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                clusterplot1 <- DimPlot(v$scData, label = T)
                #umapplot2 <- DimPlot(v$scData, reduction = "umap", label = T, group.by = 'batch')
                #umapplot3 <- DimPlot(v$scData, reduction = "umap", label = T, group.by = 'celltype')
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(clusterplot1)
                #print(umapplot2)
                #print(umapplot3)
                dev.off()
            })
            withProgress(message="Downloading cluster coordinates...", value=0.6, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData@active.ident, file = filename2)
            })
        }
    })

    #output$Cluster2DPlot_2 <- renderPlotly({
    #  if(is.null(v$isClusterdone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating 2D Cluster Plot...", value=0, {
    #      DimPlot(v$scData, reduction = "pca", label = T, group.by = 'batch')
    #    })
    #  }
    #})

    #output$Cluster2DPlot_3 <- renderPlotly({
    #  if(is.null(v$isClusterdone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating 2D Cluster Plot...", value=0, {
    #      DimPlot(v$scData, reduction = "pca", label = T, group.by = 'celltype')
    #    })
    #  }
    #})

    ##---------------UMAP tabset-------------------

    observeEvent(input$doUmap, {
        withProgress(message = "Running UMAP...", value = 0.3, {
            #dims.use <- NULL
            #if(!is.null(v$scData@reductions$pca@jackstraw)){
            #score.df <- JackStraw_pval(v$scData)
            #dims.use <- signif_PCs(score.df)
            #}else{
            #dims.use <- 1:input$dim.used
            #}
            v$scData <- RunUMAP(v$scData, dims = 1:input$dim.used, assay = input$assays1)
            output$Umap.done <- renderText(paste0("UMAP done!"))
            v$isUMAPdone <- TRUE
        })
    })

    output$Umap_2d_plot_1 <- renderPlotly({
        if(is.null(v$scData) || is.null(v$isUMAPdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating UMAP 2D Plot...", value=0, {
                DimPlot(v$scData, reduction = "umap", label = T)
            })
        }
    })

    #  output$Umap_2d_plot_2 <- renderPlotly({
    #    if(is.null(v$scData) || is.null(v$isUMAPdone)){
    #    plotly_empty()
    #   }else{
    #      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    #        DimPlot(v$scData, reduction = "umap", label = T)
    #      })
    #    }
    # })

    # output$Umap_2d_plot_3 <- renderPlotly({
    #  if(is.null(v$scData) || is.null(v$isUMAPdone)){
    #      plotly_empty()
    #    }else{
    #      withProgress(message="Generating UMAP 2D Plot...", value=0, {
    #        DimPlot(v$scData, reduction = "umap", label = T, group.by = 'celltype')
    #      })
    #    }
    #  })

    # output$Umap_3d_plot <- renderPlotly({
    #   if(is.null(v$scData) || is.null(v$isUMAPdone)){
    #     plotly_empty()
    #  }else{
    #     withProgress(message="Generating UMAP 3D Plot...", value=0, {
    #     DimPlot(v$scData, reduction = "umap", label = T)
    #      })
    #    }
    # })

    output$selection.summary1 <- renderText({
        if(is.null(v$plotlySelection)){
            return(NULL)
        }else{
            t1 <- paste0(length(v$plotlySelection), " cells selected")
            t1
        }
    })

    observeEvent(input$create.selection1, {
        ## stash old identity
        if(is.null(v$scData@meta.data$cluster.ident)){
            v$scData <- StashIdent(object = v$scData, save.name = 'cluster.ident')
        }
        v$scData <- Idents (object = v$scData,
                            cells = v$plotlySelection,
                            value = as.character(input$selection.name)
        )
        updateTabsetPanel(session, "tabs", selected = "DEGs")
    })

    observeEvent(input$reset.selection1, {
        v$scData <- SetAllIdent(object = v$scData,  id = 'cluster.ident')
        #event_data("plotly_select") <- NULL
        v$plotlySelection <- NULL
    })

    observeEvent(input$PDFi, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"UMAP_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"UMAP_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                umapplot1 <- DimPlot(v$scData, reduction = "umap", label = T)
                #umapplot2 <- DimPlot(v$scData, reduction = "umap", label = T, group.by = 'batch')
                #umapplot3 <- DimPlot(v$scData, reduction = "umap", label = T, group.by = 'celltype')
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(umapplot1)
                #print(umapplot2)
                #print(umapplot3)
                dev.off()
            })
            withProgress(message="Downloading UMAP coordinates...", value=0.6, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"umap_", Sys.Date(), ".txt")
                j = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    j = j + 1;
                }
                write.csv(v$scData@reductions$umap@cell.embeddings, file = filename2)
            })
        }
    })

    ##---------------TSNE tabset-------------------
    output$perplex.option <- renderUI({
        if(is.null(v$isPCAdone)){
            return(NULL)
        }else{
            ##perplexity test
            n.cells <- isolate(nrow(v$scData@reductions$pca@cell.embeddings))
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
            #dims.use <- NULL
            #if(!is.null(v$scData@reductions$pca@jackstraw)){
            # score.df <- JackStraw_pval(v$scData)
            # dims.use <- signif_PCs(score.df)
            # }else{
            # dims.use <- 1:input$dim.used
            #}
            v$scData <- RunTSNE(v$scData, dims = 1:input$dim.used, perplexity = input$perplexity, assay = input$assays1)
            output$Tsne.done <- renderText(paste0("TSNE done!"))
            v$isTSNEdone <- TRUE
        })
    })

    output$Tsne_2d_plot_1 <- renderPlotly({
        if(is.null(v$scData) || is.null(v$isTSNEdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating TSNE 2D Plot...", value=0, {
                DimPlot(v$scData, reduction = "tsne", label = T)
            })
        }
    })

    # output$Tsne_2d_plot_2 <- renderPlotly({
    #  if(is.null(v$scData) || is.null(v$isTSNEdone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating TSNE 2D Plot...", value=0, {
    #      DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'batch')
    #    })
    #  }
    #})

    #output$Tsne_2d_plot_3 <- renderPlotly({
    #  if(is.null(v$scData) || is.null(v$isTSNEdone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating TSNE 2D Plot...", value=0, {
    #      DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'celltype')
    #    })
    #  }
    #})

    # output$Tsne_3d_plot <- renderPlotly({
    #     if(is.null(v$scData) || is.null(v$isTSNEdone)){
    #         plotly_empty()
    #    }else{
    #          withProgress(message="Generating TSNE 3D Plot...", value=0, {
    #             DimPlot(v$scData, reduction = "tsne", label = T)
    #          })
    #     }
    #  })

    output$selection.summary <- renderText({
        if(is.null(v$plotlySelection)){
            return(NULL)
        }else{
            t <- paste0(length(v$plotlySelection), " cells selected")
            t
        }
    })

    observeEvent(input$create.selection, {
        ## stash old identity
        if(is.null(v$scData@meta.data$cluster.ident)){
            v$scData <- StashIdent(object = v$scData, save.name = 'cluster.ident')
        }
        v$scData <- Idents(object = v$scData,
                           cells = v$plotlySelection,
                           value = as.character(input$selection.name)
        )
        updateTabsetPanel(session, "tabs", selected = "DEGs")
    })

    observeEvent(input$reset.selection, {
        v$scData <- SetAllIdent(object = v$scData,  id = 'cluster.ident')
        #event_data("plotly_select") <- NULL
        v$plotlySelection <- NULL
    })

    observeEvent(input$PDFj, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"TSNE_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"TSNE_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                tsneplot1 <- DimPlot(v$scData, reduction = "tsne", label = T)
                #tsneplot2 <- DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'batch')
                #tsneplot3 <- DimPlot(v$scData, reduction = "tsne", label = T, group.by = 'celltype')
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(tsneplot1)
                #print(tsneplot2)
                #print(tsneplot3)
                dev.off()
            })
            withProgress(message="Downloading tSNE coordinates...", value=0.6, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"tsne_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData@reductions$tsne@cell.embeddings, file = filename2)
            })
        }
    })

    ##---------------DEGs tabset-------------------

    observeEvent(input$doDeg, {
        if(is.null(v$scData)){
            return(NULL)
        }else{
            withProgress(message="Finding DEGs...", value=0, {
                v$scData$seurat_clusters -> Idents(v$scData)
                ips.markers <- FindAllMarkers(v$scData, only.pos = FALSE, min.pct = input$min_pct, logfc.threshold = input$logfc, assay = input$assays1, test.use = input$test.use)
                v$ips.markers <- ips.markers
            })
        }
    })

    output$deg.gene.select <- renderUI({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            selectInput("deg.gene", label = "Gene to visualise",
                        choices = rownames(v$ips.markers))
        }
    })

    output$Deg.plot <- renderPlotly({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            withProgress(message="Generating DEG Plot...", value=0, {
                VlnPlot(v$scData, input$deg.gene)
            })
        }
    })

    output$Deg1.plot <- renderPlotly({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            withProgress(message="Generating DEG Plot...", value=0, {
                FeaturePlot(v$scData, input$deg.gene)
            })
        }
    })

    output$Deg2.plot <- renderPlotly({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            withProgress(message="Generating DEG Plot...", value=0, {
                DotPlot(v$scData, features = input$deg.gene)
            })
        }
    })

    output$Deg3.plot <- renderPlotly({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            withProgress(message="Generating DEG Plot...", value=0, {
                v$ips.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
                DoHeatmap(v$scData, features = top10$gene, size = 4) + theme(axis.text.y = element_text(size = 5)) + NoLegend()
            })
        }
    })

    output$Deg.table <- DT::renderDataTable(
        v$ips.markers, options = list(scrollX = TRUE, scrollY = "400px"))

    observeEvent(input$PDFk, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                degVln <- VlnPlot(v$scData, input$deg.gene)
                degFeature <- FeaturePlot(v$scData, input$deg.gene)
                degDot <- DotPlot(v$scData, features = input$deg.gene)
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(degVln)
                print(degFeature)
                print(degDot)
                dev.off()
                write.csv(v$ips.markers, file = paste0(pdfDir, .Platform$file.sep,"DEG_table_", input$c1, "vs", input$c2, "_", Sys.Date(), ".csv"))
            })
        }
    })

    ##---------------Trajectory Inference-------------------

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

    ##---------------Run CELLiD-------------------

    observeEvent(input$doCELLiD, {
        withProgress(message = "Running CELLiD...", value = 0.3, {
            ref = readRDS("rna_average_v2.rds")
            v$scData.rna.data.average = AverageExpression(v$scData)
            v$scData.rna.data.average = data.frame(v$scData.rna.data.average$RNA)
            write.table(v$scData.rna.data.average, "CELLiD_input.txt", quote = F, col.names = F, row.names = T, sep="\t")
            v$scData.rna.data.average = read.csv("CELLiD_input.txt", sep = "\t", header = F)
            rownames(v$scData.rna.data.average) = v$scData.rna.data.average$V1
            v$scData.rna.data.average = v$scData.rna.data.average[,-1,drop=F]
            v$scData.rna.data.average[,(ncol(v$scData.rna.data.average)+1)] = v$scData.rna.data.average[,1]
            v$res = predict.ct(v$scData.rna.data.average, ref)
            v$res = v$res[-nrow(v$res),,drop=F]
            v$res$original_name -> ot
            names(ot) <- levels(v$scData)
            v$scData <- RenameIdents(v$scData, ot)
            output$CELLiD.done <- renderText(paste0("Cell type identification done!"))
            v$isCELLiDdone <- TRUE
        })
    })

    output$Umap_cellid <- renderPlotly({
        if(is.null(v$scData) || is.null(v$isCELLiDdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating UMAP from CELLiD...", value=0, {
                DimPlot(v$scData, reduction = "umap", label = T)
            })
        }
    })

    observeEvent(input$doCC, {
        withProgress(message = "Running cell-cell communication...", value = 0.3, {
            v$scData_cc <- GetAssayData(v$scData, assay = "RNA", slot = "data")
            v$scData_ccmeta <- v$scData@meta.data
            v$scData_tpm <- createCellChat(object = v$scData_cc, meta = v$scData_ccmeta, group.by = "Celltype")
            CellChatDB <- CellChatDB.human
            v$scData_tpm@DB <- CellChatDB
            v$scData_tpm <- subsetData(v$scData_tpm)
            future::plan("multiprocess", workers = 4) # do parallel
            v$scData_tpm <- identifyOverExpressedGenes(v$scData_tpm)
            v$scData_tpm <- identifyOverExpressedInteractions(v$scData_tpm)
            # project gene expression data onto PPI network (optional)
            v$scData_tpm <- projectData(v$scData_tpm, PPI.human)
            v$scData_tpm <- computeCommunProb(v$scData_tpm)
            # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
            v$scData_tpm <- filterCommunication(v$scData_tpm, min.cells = 10)
            v$scData_tpm <- computeCommunProbPathway(v$scData_tpm)
            v$scData_tpm <- aggregateNet(v$scData_tpm)
            print(v$scData_tpm)
            v$scData_groupSize <- as.numeric(table(v$scData_tpm@idents))
            output$CC.done <- renderText(paste0("Cell-cell communication done!"))
            v$isCCdone <- TRUE
        })
    })

    output$CC_plot1 <- renderPlot({
        if(is.null(v$scData) || is.null(v$isCCdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_circle(v$scData_tpm@net$count, vertex.weight = v$scData_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
            })
        }
    })

    output$CC_plot2 <- renderPlot({
        if(is.null(v$scData) || is.null(v$isCCdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_circle(v$scData_tpm@net$weight, vertex.weight = v$scData_groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
            })
        }
    })

    observeEvent(input$doCC1, {
        withProgress(message = "Running cell-cell communication...", value = 0.3, {
            pathway <- v$scData_tpm@netP$pathways
            v$scData_vertex.receiver = seq(1,4) # a numeric vector.
            output$CC1.done <- renderText(paste0("Cell-cell communication done!"))
            v$isCC1done <- TRUE
            v$scData_pathway <- pathway
            print(v$scData_pathway)
            print(v$scData_vertex.receiver)
        })
    })

    output$CC.gene.select <- renderUI({
        if(is.null(v$scData_pathway)){
            return(NULL)
        }else{
            selectInput("CC.gene", label = "Pathway to visualise",
                        choices = v$scData_pathway)
        }
    })

    output$CC_plot3 <- renderPlot({
        if(is.null(v$scData_pathway)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_aggregate(v$scData_tpm, signaling = input$CC.gene, vertex.receiver = v$scData_vertex.receiver, layout = input$layout_cc)
            })
        }
    })

    output$CC_plot4 <- renderPlot({
        if(is.null(v$scData_pathway)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_heatmap(v$scData_tpm, signaling = input$CC.gene, color.heatmap = "Reds")
            })
        }
    })

    output$CC_plot5 <- renderPlot({
        if(is.null(v$scData_pathway)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netAnalysis_contribution(v$scData_tpm, signaling = input$CC.gene)
            })
        }
    })

    observeEvent(input$doCC2, {
        withProgress(message = "Running cell-cell communication...", value = 0.3, {
            v$scData_pairLR <- extractEnrichedLR(v$scData_tpm, signaling = input$CC.gene, geneLR.return = FALSE)
            v$scData_LR.show <- v$scData_pairLR[,1]
            v$scData_vertex.receiver = seq(1,4) # a numeric vector.
            output$CC2.done <- renderText(paste0("Cell-cell communication done!"))
            v$isCC2done <- TRUE
            print(v$scData_LR.show)
        })
    })

    output$LR.gene.select <- renderUI({
        if(is.null(v$scData_LR.show)){
            return(NULL)
        }else{
            selectInput("LR.gene", label = "LR pair to visualise",
                        choices = v$scData_LR.show)
        }
    })

    output$CC_plot6 <- renderPlot({
        if(is.null(v$scData_pathway)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_individual(v$scData_tpm, signaling = input$CC.gene,  pairLR.use = input$LR.gene, vertex.receiver = v$scData_vertex.receiver, layout = input$layout_cc1)
            })
        }
    })

    output$cell_group.select <- renderUI({
        if(is.null(v$scData)){
            return(NULL)
        }else{
            selectInput("cell_group", label = "Cell group",
                        choices = unique(v$scData$Celltype))
        }
    })

    output$CC_plot7 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_bubble(v$scData_tpm, sources.use = input$cell_group, targets.use = c(1:length(unique(v$scData$Celltype))), remove.isolate = FALSE)
            })
        }
    })

    output$CC_plot8 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netVisual_chord_gene(v$scData_tpm, sources.use = input$cell_group, targets.use = c(1:length(unique(v$scData$Celltype))), lab.cex = 0.5,legend.pos.y = 30)
            })
        }
    })

    output$CC_plot9 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                plotGeneExpression(v$scData_tpm, signaling = input$CC.gene)
            })
        }
    })

    observeEvent(input$doCC3, {
        withProgress(message = "Running cell-cell communication...", value = 0.3, {
            v$scData_tpm <- netAnalysis_computeCentrality(v$scData_tpm, slot.name = "netP")
            output$CC3.done <- renderText(paste0("Cell-cell communication done!"))
            v$isCC3done <- TRUE
        })
    })

    output$CC_plot10 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netAnalysis_signalingRole_network(v$scData_tpm, signaling = input$CC.gene, width = 8, height = 2.5, font.size = 10)
            })
        }
    })

    output$CC_plot11 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netAnalysis_signalingRole_scatter(v$scData_tpm)
            })
        }
    })

    output$CC_plot12 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netAnalysis_signalingRole_heatmap(v$scData_tpm, pattern = "outgoing")
            })
        }
    })

    output$CC_plot13 <- renderPlot({
        if(is.null(v$scData_tpm)){
            return(NULL)
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                netAnalysis_signalingRole_heatmap(v$scData_tpm, pattern = "incoming")
            })
        }
    })
    #}
    #})

    observeEvent(input$do_milo_preprocess, {
        withProgress(message = "Running Milo analysis...", value = 0.3, {
            v$scData_Diet <- DietSeurat(v$scData, graphs = c("umap", "pca"))
            print(v$scData)
            print(v$scData_Diet)
            as.SingleCellExperiment(v$scData_Diet)-> v$scData_milo
            print(v$scData_milo)
            v$ismilo_preprocessdone <- TRUE
        })
    })

    output$vis_milo_plot1 <- renderPlot({
        if(is.null(v$scData_milo) || is.null(v$ismilo_preprocessdone)){
            plotly_empty()
        }else{
            withProgress(message="Visualize the data...", value=0, {
                plotReducedDim(v$scData_milo, colour_by=input$group.by_milo, dimred = "UMAP")
            })
        }
    })

    #output$vis_milo_plot2 <- renderPlot({
    #  if(is.null(v$scData_milo) || is.null(v$ismilo_preprocessdone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating cell-cell communication plot...", value=0, {
    #      plotReducedDim(v$scData_milo, colour_by="Day", dimred = "UMAP")
    #    })
    #  }
    #})

    observeEvent(input$do_milo, {
        withProgress(message = "Running Milo...", value = 0.3, {
            v$scData_milo <- Milo(v$scData_milo)
            print(v$scData_milo)
            v$scData_milo <- buildGraph(v$scData_milo, k = 30, d = 30, reduced.dim = "PCA")
            print(v$scData_milo)
            v$scData_milo <- makeNhoods(v$scData_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")
            print(v$scData_milo)
            v$scData_milo <- countCells(v$scData_milo, meta.data = as.data.frame(colData(v$scData_milo)), samples="Sample")
            print(v$scData_milo)
            v$scData_embryo_design <- data.frame(colData(v$scData_milo))[,c("Sample", "Genotype", "Batch", "Day")]
            v$scData_embryo_design$Batch <- as.factor(v$scData_embryo_design$Batch)
            v$scData_embryo_design <- distinct(v$scData_embryo_design)
            rownames(v$scData_embryo_design) <- v$scData_embryo_design$Sample
            print(v$scData_embryo_design)
            v$scData_milo <- calcNhoodDistance(v$scData_milo, d=input$dim.used_milo, reduced.dim = "PCA")
            v$da_results <- testNhoods(v$scData_milo, design = ~ Batch + Genotype, design.df = v$scData_embryo_design)
            v$scData_milo <- buildNhoodGraph(v$scData_milo)
            v$da_results <- annotateNhoods(v$scData_milo, v$da_results, coldata_col = "Celltype")
            print(v$da_results)
            #v$da_results$Celltype <- ifelse(v$da_results$Celltype_fraction < 0.7, "Mixed", v$da_results$Celltype)
            v$ismilodone <- TRUE
        })
    })

    output$vis_milo_plot3 <- renderPlot({
        if(is.null(v$scData_milo) || is.null(v$ismilodone)){
            plotly_empty()
        }else{
            withProgress(message="Visualize the data...", value=0, {
                plotReducedDim(v$scData_milo, dimred = "UMAP", colour_by="Genotype", text_by = "Celltype",text_size = 3, point_size=0.5) + guides(fill="none")
            })
        }
    })

    output$vis_milo_plot4 <- renderPlot({
        if(is.null(v$scData_milo) || is.null(v$ismilodone)){
            plotly_empty()
        }else{
            withProgress(message="Generating cell-cell communication plot...", value=0, {
                plotNhoodGraphDA(v$scData_milo, v$da_results, layout="UMAP",alpha=0.1)
            })
        }
    })

    #output$vis_milo_plot5 <- renderPlot({
    #  if(is.null(v$scData_milo) || is.null(v$ismilodone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating cell-cell communication plot...", value=0, {
    #      plotDAbeeswarm(v$da_results, group.by = "Celltype")
    #    })
    #  }
    #})

    #output$vis_milo_plot5 <- renderPlot({
    #  if(is.null(v$scData_milo) || is.null(v$ismilodone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Visualize the data...", value=0, {
    #      plotNhoodSizeHist(v$scData_milo)
    #    })
    #  }
    #})

    #output$vis_milo_plot6 <- renderPlot({
    #  if(is.null(v$scData_milo) || is.null(v$ismilodone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating cell-cell communication plot...", value=0, {
    #      ggplot(v$da_results, aes(PValue)) + geom_histogram(bins=50)
    #    })
    #  }
    #})

    #output$vis_milo_plot7 <- renderPlot({
    #  if(is.null(v$scData_milo) || is.null(v$ismilodone)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating cell-cell communication plot...", value=0, {
    #      ggplot(v$da_results, aes(logFC, -log10(SpatialFDR))) + geom_point() + geom_hline(yintercept = 1)
    #    })
    #  }
    #})

    ##---------------Data Integration using Seurat-------------------

    observe({if(input$scAnalysis_integ == "Seurat"){
        #if (input$scAnalysis_integ == "Seurat")

        observeEvent(input$doIntg_seurat, {
            withProgress(message = "Running Data Integration...", value = 0.3, {
                v$scData1.list <- SplitObject(v$scData1, split.by = "batch")
                print(v$scData1)
                print(v$scData1.list)
                v$scData1.list <- pbmclapply(mc.cores = 20, X = v$scData1.list, FUN = function(x) {
                    x <- NormalizeData(x)
                    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                })
                features <- SelectIntegrationFeatures(object.list = v$scData1.list)
                print(features)
                v$scData1.anchors <- FindIntegrationAnchors(object.list = v$scData1.list, anchor.features = features, nn.method = "rann")
                v$scData1.combined <- IntegrateData(anchorset = v$scData1.anchors)
                print(v$scData1.anchors)
                print(v$scData1.combined)

                DefaultAssay(v$scData1.combined) <- "integrated"
                v$scData1.combined <- ScaleData(v$scData1.combined, verbose = FALSE)
                print(v$scData1.combined)
            })
        })

        observeEvent(input$runPCA_intg_seurat, {
            withProgress(message = "Running PCA...", value = 0,{
                incProgress(0.5, message = "Running PCA...")
                if (input$scInput1 == "Raw Counts Matrix")
                {
                    v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
                    print(v$scData1.combined[["pca"]], dims = 1:5, nfeatures = 5)
                    v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:30, nn.method = "rann")
                    v$scData1.combined <- FindClusters(v$scData1.combined, resolution = 0.5)
                    v$isPCAdone1 <- TRUE
                    PCA_plot1a <- DimPlot(v$scData1.combined, reduction = "pca", label = T)
                    PCA_plot1b <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch')
                    PCA_plot1c <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype')
                    print(PCA_plot1a)
                    print(PCA_plot1b)
                    print(PCA_plot1c)
                }
                else if (input$scInput1 == "H5")
                {
                    v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
                    print(v$scData1.combined[["pca"]], dims = 1:5, nfeatures = 5)
                    v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "pca", dims = 1:30, nn.method = "rann")
                    v$scData1.combined <- FindClusters(v$scData1.combined, resolution = 0.5)
                    v$isPCAdone1 <- TRUE
                    PCA_plot1a <- DimPlot(v$scData1.combined, reduction = "pca", label = T)
                    PCA_plot1b <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch')
                    print(PCA_plot1a)
                    print(PCA_plot1b)
                }
            })
        })

        output$PCAplot_seurat_tpm1 <- renderPlotly({
            if(is.null(v$isPCAdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "pca", label = T)
                })
            }
        })

        output$PCAplot_seurat_tpm2 <- renderPlotly({
            if(is.null(v$isPCAdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch')
                })
            }
        })

        output$PCAplot_seurat_tpm3 <- renderPlotly({
            if(is.null(v$isPCAdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype')
                })
            }
        })

        output$PCAplot_seurat_h5_1 <- renderPlotly({
            if(is.null(v$isPCAdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "pca", label = T)
                })
            }
        })

        output$PCAplot_seurat_h5_2 <- renderPlotly({
            if(is.null(v$isPCAdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch')
                })
            }
        })

        observeEvent(input$runUMAP_intg_seurat, {
            withProgress(message = "Running UMAP...", value = 0,{
                incProgress(0.5, message = "Running UMAP...")
                if (input$scInput1 == "Raw Counts Matrix")
                {
                    v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:30)
                    v$isUMAPdone1 <- TRUE
                    UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
                    UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                    UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype')
                    print(UMAP_plot1a)
                    print(UMAP_plot1b)
                    print(UMAP_plot1c)
                }
                else if (input$scInput1 == "H5")
                {
                    v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "pca", dims = 1:30)
                    v$isUMAPdone1 <- TRUE
                    UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
                    UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                    print(UMAP_plot1a)
                    print(UMAP_plot1b)
                }

            })
        })

        output$UMAPplot_seurat_tpm1 <- renderPlotly({
            if(is.null(v$isUMAPdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T)
                })
            }
        })

        output$UMAPplot_seurat_tpm2 <- renderPlotly({
            if(is.null(v$isUMAPdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                })
            }
        })

        output$UMAPplot_seurat_tpm3 <- renderPlotly({
            if(is.null(v$isUMAPdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype')
                })
            }
        })

        output$UMAPplot_seurat_h5_1 <- renderPlotly({
            if(is.null(v$isUMAPdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T)
                })
            }
        })

        output$UMAPplot_seurat_h5_2 <- renderPlotly({
            if(is.null(v$isUMAPdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                })
            }
        })

        observeEvent(input$runTSNE_intg_seurat, {
            withProgress(message = "Running TSNE...", value = 0,{
                incProgress(0.5, message = "Running TSNE...")
                if (input$scInput1 == "Raw Counts Matrix")
                {
                    v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:30)
                    v$isTSNEdone1 <- TRUE
                    TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                    TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                    TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype')
                    print(TSNE_plot1a)
                    print(TSNE_plot1b)
                    print(TSNE_plot1c)
                }
                else if (input$scInput1 == "H5")
                {
                    v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "pca", dims = 1:30)
                    v$isTSNEdone1 <- TRUE
                    TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                    TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                    print(TSNE_plot1a)
                    print(TSNE_plot1b)
                }
            })
        })

        output$TSNEplot_seurat_tpm1 <- renderPlotly({
            if(is.null(v$isTSNEdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                })
            }
        })

        output$TSNEplot_seurat_tpm2 <- renderPlotly({
            if(is.null(v$isTSNEdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                })
            }
        })

        output$TSNEplot_seurat_tpm3 <- renderPlotly({
            if(is.null(v$isTSNEdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype')
                })
            }
        })

        output$TSNEplot_seurat_h5_1 <- renderPlotly({
            if(is.null(v$isTSNEdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                })
            }
        })

        output$TSNEplot_seurat_h5_2 <- renderPlotly({
            if(is.null(v$isTSNEdone1)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                })
            }
        })

        # output$ident.swap1 <- renderUI({
        #  if(is.null(v$scData1.combined)){
        #   return(NULL)
        #   }else{
        #      groupings1 <- names(v$scData1.combined@meta.data[,!names(v$scData1.combined@meta.data) %in% c("nGene", "nUMI", "percent.mito")])
        #     tagList(
        #        h4("Set current identity:"),
        #        fluidRow(
        #          column(3,
        #                selectInput("active.ident1", label = NULL,
        #                             choices = groupings1)
        #         ),
        #          column(3,
        #                actionButton("swap.ident1",label = NULL, icon = icon("arrow-right"))
        #        )
        #     )
        #      )
        #   }
        #  })

        observeEvent(input$PDFl, {
            if(!is.null(v$scData1.combined) ){
                withProgress(message="Downloading plot PDF files...", value=0, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_Seurat_", Sys.Date(), ".pdf")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_Seurat_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                        i = i + 1;
                    }
                    PCA_plot1a <- DimPlot(v$scData1.combined, reduction = "pca", label = T)
                    PCA_plot1b <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch')
                    PCA_plot1c <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype')
                    UMAP_plot1a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
                    UMAP_plot1b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                    UMAP_plot1c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype')
                    TSNE_plot1a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                    TSNE_plot1b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                    TSNE_plot1c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype')
                    prePlot()
                    pdf(filename2,
                        width=as.numeric(input$pdf_w),
                        height=as.numeric(input$pdf_h))
                    print(PCA_plot1a)
                    print(PCA_plot1b)
                    print(PCA_plot1c)
                    print(UMAP_plot1a)
                    print(UMAP_plot1b)
                    print(UMAP_plot1c)
                    print(TSNE_plot1a)
                    print(TSNE_plot1b)
                    print(TSNE_plot1c)
                    dev.off()
                })
                withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2a <- paste0(pdfDir, .Platform$file.sep,"seurat_pca_", Sys.Date(), ".txt")
                    filename2b <- paste0(pdfDir, .Platform$file.sep,"seurat_umap_", Sys.Date(), ".txt")
                    filename2c <- paste0(pdfDir, .Platform$file.sep,"seurat_tsne_", Sys.Date(), ".txt")
                    i = 0
                    while(file.exists(filename2a)){
                        filename2a <- paste0(pdfDir, .Platform$file.sep,"seurat_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    while(file.exists(filename2b)){
                        filename2b <- paste0(pdfDir, .Platform$file.sep,"seurat_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    while(file.exists(filename2c)){
                        filename2c <- paste0(pdfDir, .Platform$file.sep,"seurat_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    write.csv(v$scData1.combined@reductions$pca@cell.embeddings, file = filename2a)
                    write.csv(v$scData1.combined@reductions$umap@cell.embeddings, file = filename2b)
                    write.csv(v$scData1.combined@reductions$tsne@cell.embeddings, file = filename2c)
                })
            }
        })
    }
    })
    ##---------------Data Integration using Harmony-------------------

    observe({if(input$scAnalysis_integ == "Harmony"){

        observeEvent(input$doIntg_harmony, {
            withProgress(message = "Running Data Integration...", value = 0.3, {
                #v$scData <- merge(v$scData4, y = v$scData5, add.cell.ids = c("1", "2"), project = input$projName)
                #v$scData[['type']] <- ifelse(startsWith(v$scData6@assays$RNA@data@Dimnames[[2]], '1'), '1', '2')
                v$scData1.list <- SplitObject(v$scData1, split.by = "batch")
                print(v$scData1)
                print(v$scData1.list)
                v$scData1.list <- pbmclapply(mc.cores = 20, X = v$scData1.list, FUN = function(x) {
                    x <- NormalizeData(x)
                    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                })

                features <- SelectIntegrationFeatures(object.list = v$scData1.list)
                print(features)
                v$scData1.anchors <- FindIntegrationAnchors(object.list = v$scData1.list, anchor.features = features, nn.method = "rann")
                v$scData1.combined <- IntegrateData(anchorset = v$scData1.anchors)
                print(v$scData1.anchors)
                print(v$scData1.combined)

                DefaultAssay(v$scData1.combined) <- "integrated"
                v$scData1.combined <- ScaleData(v$scData1.combined, verbose = FALSE)
                print(v$scData1.combined)
            })
        })

        observeEvent(input$runPCA_intg_harmony, {
            withProgress(message = "Running PCA...", value = 0,{
                incProgress(0.5, message = "Running PCA...")
                if (input$scInput1 == "Raw Counts Matrix")
                {
                    v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
                    print(v$scData1.combined[["pca"]], dims = 1:5, nfeatures = 5)
                    v$isPCAdone2 <- TRUE
                    v$scData1.combined <- RunHarmony(v$scData1.combined, "batch", assay.use = "integrated")
                    v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:30, nn.method = "rann")
                    v$scData1.combined <- FindClusters(v$scData1.combined, resolution = 0.5)
                    PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "harmony", label = T)
                    PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "harmony", group.by = 'batch')
                    PCA_plot2c <- DimPlot(v$scData1.combined, reduction = "harmony", group.by = 'celltype')
                    print(PCA_plot2a)
                    print(PCA_plot2b)
                    print(PCA_plot2c)
                }
                else if (input$scInput1 == "H5")
                {
                    v$scData1.combined <- RunPCA(v$scData1.combined, verbose = FALSE)
                    print(v$scData1.combined[["pca"]], dims = 1:5, nfeatures = 5)
                    v$isPCAdone2 <- TRUE
                    v$scData1.combined <- RunHarmony(v$scData1.combined, "batch", assay.use = "integrated")
                    v$scData1.combined <- FindNeighbors(v$scData1.combined, reduction = "harmony", dims = 1:30, nn.method = "rann")
                    v$scData1.combined <- FindClusters(v$scData1.combined, resolution = 0.5)
                    PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "harmony", label = T)
                    PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "harmony", group.by = 'batch')
                    print(PCA_plot2a)
                    print(PCA_plot2b)
                }

            })
        })

        output$PCAplot_harmony_tpm1 <- renderPlotly({
            if(is.null(v$isPCAdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "harmony", label = T)
                })
            }
        })

        output$PCAplot_harmony_tpm2 <- renderPlotly({
            if(is.null(v$isPCAdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'batch')
                })
            }
        })

        output$PCAplot_harmony_tpm3 <- renderPlotly({
            if(is.null(v$isPCAdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'celltype')
                })
            }
        })

        output$PCAplot_harmony_h5_1 <- renderPlotly({
            if(is.null(v$isPCAdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "harmony", label = T)
                })
            }
        })

        output$PCAplot_harmony_h5_2 <- renderPlotly({
            if(is.null(v$isPCAdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "harmony", label = T, group.by = 'batch')
                })
            }
        })

        observeEvent(input$runUMAP_intg_harmony, {
            withProgress(message = "Running UMAP...", value = 0,{
                incProgress(0.5, message = "Running UMAP...")
                if (input$scInput1 == "Raw Counts Matrix")
                {
                    v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "harmony", dims = 1:30)
                    v$isUMAPdone2 <- TRUE
                    UMAP_plot2a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
                    UMAP_plot2b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                    UMAP_plot2c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype')
                    print(UMAP_plot2a)
                    print(UMAP_plot2b)
                    print(UMAP_plot2c)
                }
                else if (input$scInput1 == "H5")
                {
                    v$scData1.combined <- RunUMAP(v$scData1.combined, reduction = "harmony", dims = 1:30)
                    v$isUMAPdone2 <- TRUE
                    UMAP_plot2a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
                    UMAP_plot2b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                    print(UMAP_plot2a)
                    print(UMAP_plot2b)
                }
            })
        })

        output$UMAPplot_harmony_tpm1 <- renderPlotly({
            if(is.null(v$isUMAPdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T)
                })
            }
        })

        output$UMAPplot_harmony_tpm2 <- renderPlotly({
            if(is.null(v$isUMAPdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                })
            }
        })

        output$UMAPplot_harmony_tpm3 <- renderPlotly({
            if(is.null(v$isUMAPdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype')
                })
            }
        })

        output$UMAPplot_harmony_h5_1 <- renderPlotly({
            if(is.null(v$isUMAPdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T)
                })
            }
        })

        output$UMAPplot_harmony_h5_2 <- renderPlotly({
            if(is.null(v$isUMAPdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                })
            }
        })

        observeEvent(input$runTSNE_intg_harmony, {
            withProgress(message = "Running TSNE...", value = 0,{
                incProgress(0.5, message = "Running TSNE...")
                if (input$scInput1 == "Raw Counts Matrix")
                {
                    v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "harmony", dims = 1:30)
                    v$isTSNEdone2 <- TRUE
                    TSNE_plot2a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                    TSNE_plot2b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                    TSNE_plot2c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype')
                    print(TSNE_plot2a)
                    print(TSNE_plot2b)
                    print(TSNE_plot2c)
                }
                else if (input$scInput1 == "H5")
                {
                    v$scData1.combined <- RunTSNE(v$scData1.combined, reduction = "harmony", dims = 1:30)
                    v$isTSNEdone2 <- TRUE
                    TSNE_plot2a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                    TSNE_plot2b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                    print(TSNE_plot2a)
                    print(TSNE_plot2b)
                }
            })
        })

        output$TSNEplot_harmony_tpm1 <- renderPlotly({
            if(is.null(v$isTSNEdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                })
            }
        })

        output$TSNEplot_harmony_tpm2 <- renderPlotly({
            if(is.null(v$isTSNEdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                })
            }
        })

        output$TSNEplot_harmony_tpm3 <- renderPlotly({
            if(is.null(v$isTSNEdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype')
                })
            }
        })

        output$TSNEplot_harmony_h5_1 <- renderPlotly({
            if(is.null(v$isTSNEdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                })
            }
        })

        output$TSNEplot_harmony_h5_2 <- renderPlotly({
            if(is.null(v$isTSNEdone2)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                })
            }
        })

        #output$ident.swap2 <- renderUI({
        #  if(is.null(v$scData1.combined)){
        #    return(NULL)
        #  }else{
        #    groupings2 <- names(v$scData1.combined@meta.data[,!names(v$scData1.combined@meta.data) %in% c("nGene", "nUMI", "percent.mito")])
        #    tagList(
        #      h4("Set current identity:"),
        #      fluidRow(
        #        column(3,
        #               selectInput("active.ident2", label = NULL,
        #                           choices = groupings2)
        #        ),
        #        column(3,
        #               actionButton("swap.ident2",label = NULL, icon = icon("arrow-right"))
        #        )
        #      )
        #    )
        #  }
        #})

        # observeEvent(input$swap.ident2, {
        #   v$scData1.combined <- SetIdent(v$scData1.combined, ident.use = as.character(v$scData1.combined@meta.data[,input$active.ident2]))
        # })

        observeEvent(input$PDFm, {
            if(!is.null(v$scData1.combined) ){
                withProgress(message="Downloading plot PDF files...", value=0, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_Harmony_", Sys.Date(), ".pdf")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_Harmony_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                        i = i + 1;
                    }
                    PCA_plot2a <- DimPlot(v$scData1.combined, reduction = "pca", label = T)
                    PCA_plot2b <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'batch')
                    PCA_plot2c <- DimPlot(v$scData1.combined, reduction = "pca", label = T, group.by = 'celltype')
                    UMAP_plot2a <- DimPlot(v$scData1.combined, reduction = "umap", label = T)
                    UMAP_plot2b <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'batch')
                    UMAP_plot2c <- DimPlot(v$scData1.combined, reduction = "umap", label = T, group.by = 'celltype')
                    TSNE_plot2a <- DimPlot(v$scData1.combined, reduction = "tsne", label = T)
                    TSNE_plot2b <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'batch')
                    TSNE_plot2c <- DimPlot(v$scData1.combined, reduction = "tsne", label = T, group.by = 'celltype')
                    prePlot()
                    pdf(filename2,
                        width=as.numeric(input$pdf_w),
                        height=as.numeric(input$pdf_h))
                    print(PCA_plot2a)
                    print(PCA_plot2b)
                    print(PCA_plot2c)
                    print(UMAP_plot2a)
                    print(UMAP_plot2b)
                    print(UMAP_plot2c)
                    print(TSNE_plot2a)
                    print(TSNE_plot2b)
                    print(TSNE_plot2c)
                    dev.off()
                })
                withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2a <- paste0(pdfDir, .Platform$file.sep,"harmony_pca_", Sys.Date(), ".txt")
                    filename2b <- paste0(pdfDir, .Platform$file.sep,"harmony_umap_", Sys.Date(), ".txt")
                    filename2c <- paste0(pdfDir, .Platform$file.sep,"harmony_tsne_", Sys.Date(), ".txt")
                    i = 0
                    while(file.exists(filename2a)){
                        filename2a <- paste0(pdfDir, .Platform$file.sep,"harmony_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    while(file.exists(filename2b)){
                        filename2b <- paste0(pdfDir, .Platform$file.sep,"harmony_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    while(file.exists(filename2c)){
                        filename2c <- paste0(pdfDir, .Platform$file.sep,"harmony_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    write.csv(v$scData1.combined@reductions$pca@cell.embeddings, file = filename2a)
                    write.csv(v$scData1.combined@reductions$umap@cell.embeddings, file = filename2b)
                    write.csv(v$scData1.combined@reductions$tsne@cell.embeddings, file = filename2c)
                })
            }
        })
    }
    })


    ##---------------Data Integration using RLiger-------------------

    observeEvent(input$doIntg2, {
        withProgress(message = "Running Data Integration...", value = 0.3, {
            #v$scData16 <- rbind(v$scData30, v$scData31)
            #v$scData9 <- merge(v$scData7, y = v$scData8, add.cell.ids = c("1", "2"), project = input$projName)
            #v$scData9[['type']] <- ifelse(startsWith(v$scData9@assays$RNA@data@Dimnames[[2]], '1'), '1', '2')
            v$scData1.list <- SplitObject(v$scData1, split.by = "batch")
            print(v$scData1)
            print(v$scData1.list)
            #v$scData9.list <- pbmclapply(mc.cores = 20, X = v$scData9.list, FUN = function(x) {
            #  x <-  SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt")
            #})

            v$scData1.list <- pbmclapply(X = v$scData1.list, FUN = SCTransform, mc.cores = 100)

            features <- SelectIntegrationFeatures(object.list = v$scData1.list, nfeatures = 3000)
            print(features)
            v$scData1.list <- PrepSCTIntegration(object.list = v$scData1.list, anchor.features = features)
            v$scData1.anchors <- FindIntegrationAnchors(object.list = v$scData1.list, normalization.method = "SCT", anchor.features = features, nn.method = "rann")
            v$scData1.combined <- IntegrateData(anchorset = v$scData1.anchors, normalization.method = "SCT")
            print(v$scData1.anchors)
            print(v$scData1.combined)

            DefaultAssay(v$scData1.combined) <- "integrated"
            #v$scData9.combined <- ScaleData(v$scData9.combined, verbose = FALSE)
            print(v$scData1.combined)
            v$scData1.combined@meta.data -> v$scData16
            print(v$scData16)
            #v$scData1 <- list(b1 = GetAssayData(subset(v$scData.combined, subset = batch == '1'), slot="counts", assay='SCT'), b2 = GetAssayData(subset(v$scData.combined, subset = batch == '2'), slot="counts", assay='SCT'))
            v$scData2 <- list(b1 = GetAssayData(subset(v$scData1.combined, subset = batch == '1'), slot="counts", assay='SCT'), b2 = GetAssayData(subset(v$scData1.combined, subset = batch == '2'), slot="counts", assay='SCT'), b3 = GetAssayData(subset(v$scData1.combined, subset = batch == '3'), slot="counts", assay='SCT'))
            print(v$scData2)
            #v$scData13 <- list(Tcell = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Tcell'), slot="counts", assay='SCT'), NK = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'NK'), slot="counts", assay='SCT'), Bcell = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Bcell'), slot="counts", assay='SCT'), Monocyte = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Monocyte'), slot="counts", assay='SCT'), Macrophage = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Macrophage'), slot="counts", assay='SCT'), Dendritic = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Dendritic'), slot="counts", assay='SCT'), Endothelial = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Endothelial'), slot="counts", assay='SCT'), SmoothMuscle = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'SmoothMuscle'), slot="counts", assay='SCT'), Neutrophil = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Neutrophil'), slot="counts", assay='SCT'), Stromal = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Stromal'), slot="counts", assay='SCT'), Epithelial = GetAssayData(subset(v$scData9.combined, subset = orig.ident == 'Epithelial'), slot="counts", assay='SCT'))
        })
    })

    #observeEvent(input$runPCA2, {
    # withProgress(message = "Running PCA...", value = 0,{
    #    incProgress(0.5, message = "Running PCA...")
    #    v$scData9.combined <- RunOptimizeALS(v$scData9.combined, k = 30, lambda = 5, split.by = "batch")
    #    v$scData9.combined <- RunQuantileNorm(v$scData9.combined, split.by = "batch")
    #    v$scData9.combined <- RunPCA(v$scData9.combined, verbose = FALSE)
    #    print(v$scData9.combined[["pca"]], dims = 1:5, nfeatures = 5)
    #    v$isPCAdone3 <- TRUE
    #    PCA_plot3a <- DimPlot(v$scData9.combined, reduction = "pca", label = T)
    #    PCA_plot3b <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'batch')
    #    PCA_plot3c <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'orig.ident')
    #    print(PCA_plot3a)
    #    print(PCA_plot3b)
    #    print(PCA_plot3c)
    #  })
    # })

    #output$PCAplot2_a <- renderPlotly({
    #  if(is.null(v$isPCAdone3)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
    #      DimPlot(v$scData9.combined, reduction = "pca", label = T)
    #    })
    #  }
    # })

    #output$PCAplot2_b <- renderPlotly({
    #  if(is.null(v$isPCAdone3)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
    #      DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'batch')
    #    })
    #  }
    #})

    #output$PCAplot2_c <- renderPlotly({
    #  if(is.null(v$isPCAdone3)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating PCA Plot of integrated dataset...", value=0, {
    #      DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'orig.ident')
    #    })
    #  }
    # })


    observeEvent(input$runUMAP2, {
        withProgress(message = "Running UMAP...", value = 0,{
            incProgress(0.5, message = "Running UMAP...")
            v$scData11 <- rliger::createLiger(v$scData2)
            v$scData11 <- rliger::normalize(v$scData11)
            v$scData11 <- rliger::selectGenes(v$scData11, var.thresh = 0.2, do.plot = F)
            v$scData11 <- rliger::scaleNotCenter(v$scData11)
            v$scData11 <- rliger::online_iNMF(v$scData11, k = 20, miniBatch_size = 5000, max.epochs = 5)
            v$scData11 <- rliger::quantile_norm(v$scData11)
            v$scData11 <- rliger::runUMAP(v$scData11)
            v$scData17 <- data.frame(v$scData11@tsne.coords)
            v$scData18 <- rownames(v$scData17)
            #gsub('1_','', v$scData18) -> v$scData19
            #gsub('2_','', v$scData19) -> v$scData20
            v$scData16 <- v$scData16[v$scData18,]
            v$scData11@clusters <- as.factor(v$scData16[,"celltype"])
            names(v$scData11@clusters) <- rownames(v$scData16)
            #plotByDatasetAndCluster(v$scData21, pt.size = 1, text.size = 0)

            #v$scData14 <- rliger::createLiger(v$scData13)
            #v$scData14 <- rliger::normalize(v$scData14)
            #v$scData14 <- rliger::selectGenes(v$scData14, var.thresh = 0.2, do.plot = F)
            #v$scData14 <- rliger::scaleNotCenter(v$scData14)
            #v$scData14 <- rliger::online_iNMF(v$scData14, k = 20, miniBatch_size = 5000, max.epochs = 5)
            #v$scData14 <- rliger::quantile_norm(v$scData14)
            #v$scData14 <- rliger::runUMAP(v$scData14)
            #v$scData11 <- ligerToSeurat(v$scData11, nms = names(v$scData11@H), renormalize = TRUE, use.liger.genes = TRUE, by.dataset = FALSE)
            #v$scData9.combined <- FindNeighbors(v$scData9.combined, reduction = "iNMF", dims = 1:30)
            #v$scData9.combined <- FindClusters(v$scData9.combined, resolution = 0.5)
            #v$scData9.combined <- RunUMAP(v$scData9.combined, dims = 1:30, reduction = "iNMF")
            v$isUMAPdone3 <- TRUE
            UMAP_plot3a <- rliger::plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[1]]
            UMAP_plot3b <- rliger::plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[2]]
            #UMAP_plot3b <- DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'batch')
            #UMAP_plot3c <- DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'orig.ident')
            print(UMAP_plot3a)
            print(UMAP_plot3b)
            #print(UMAP_plot3c)
        })
    })

    output$UMAPplot2_a <- renderPlotly({
        if(is.null(v$isUMAPdone3)){
            plotly_empty()
        }else{
            withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"), return.plots = T)[[1]]
            })
        }
    })

    output$UMAPplot2_b <- renderPlotly({
        if(is.null(v$isUMAPdone3)){
            plotly_empty()
        }else{
            withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"), return.plots = T)[[2]]
            })
        }
    })

    #output$UMAPplot2_c <- renderPlotly({
    #  if(is.null(v$isUMAPdone3)){
    #    plotly_empty()
    #  }else{
    #    withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
    #      DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'orig.ident')
    #    })
    #  }
    # })

    observeEvent(input$runTSNE2, {
        withProgress(message = "Running TSNE...", value = 0,{
            incProgress(0.5, message = "Running TSNE...")
            #v$scData9.combined <- RunTSNE(v$scData9.combined, dims = 1:ncol(v$scData9.combined[["iNMF"]]), reduction = "iNMF", tsne.method = "Rtsne", check_duplicates = FALSE)
            v$scData12 <- rliger::runTSNE(v$scData11)
            v$scData22 <- data.frame(v$scData12@tsne.coords)
            v$scData23 <- rownames(v$scData22)
            #gsub('1_','', v$scData23) -> v$scData24
            #gsub('2_','', v$scData24) -> v$scData25
            v$scData16 <- v$scData16[v$scData23,]
            v$scData12@clusters <- as.factor(v$scData16[,"celltype"])
            names(v$scData12@clusters) <- rownames(v$scData16)
            #plotByDatasetAndCluster(v$scData25, pt.size = 1, text.size = 0)
            #v$scData15 <- rliger::runTSNE(v$scData14)
            v$isTSNEdone3 <- TRUE
            TSNE_plot3a <- rliger::plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[1]]
            TSNE_plot3b <- rliger::plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[2]]
            #TSNE_plot3c <- DimPlot(v$scData9.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
            print(TSNE_plot3a)
            print(TSNE_plot3b)
            #print(TSNE_plot3c)
        })
    })

    output$TSNEplot2_a <- renderPlotly({
        if(is.null(v$isTSNEdone3)){
            plotly_empty()
        }else{
            withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"), return.plots = T)[[1]]
            })
        }
    })

    output$TSNEplot2_b <- renderPlotly({
        if(is.null(v$isTSNEdone3)){
            plotly_empty()
        }else{
            withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"), return.plots = T)[[2]]
            })
        }
    })

    observeEvent(input$PDFn, {
        if(!is.null(v$scData1.combined) ){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Data_integration_RLiger_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"Dimension_reduction_plots_RLiger_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                #PCA_plot3a <- DimPlot(v$scData9.combined, reduction = "pca", label = T)
                #PCA_plot3b <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'batch')
                #PCA_plot3c <- DimPlot(v$scData9.combined, reduction = "pca", label = T, group.by = 'orig.ident')
                UMAP_plot3a <- plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[1]]
                UMAP_plot3b <- plotByDatasetAndCluster(v$scData11, axis.labels = c("UMAP1","UMAP2"))[[2]]
                #UMAP_plot3c <- DimPlot(v$scData9.combined, reduction = "umap", label = T, group.by = 'orig.ident')
                TSNE_plot3a <- plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[1]]
                TSNE_plot3b <- plotByDatasetAndCluster(v$scData12, axis.labels = c("TSNE1","TSNE2"))[[2]]
                #TSNE_plot3c <- DimPlot(v$scData9.combined, reduction = "tsne", label = T, group.by = 'orig.ident')
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                #print(PCA_plot3a)
                #print(PCA_plot3b)
                #print(PCA_plot3c)
                print(UMAP_plot3a)
                print(UMAP_plot3b)
                #print(UMAP_plot3c)
                print(TSNE_plot3a)
                print(TSNE_plot3b)
                #print(TSNE_plot3c)
                dev.off()
            })
            withProgress(message="Downloading Dimension_reduction coordinates...", value=0.6, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2a <- paste0(pdfDir, .Platform$file.sep,"rliger_pca_", Sys.Date(), ".txt")
                filename2b <- paste0(pdfDir, .Platform$file.sep,"rliger_umap_", Sys.Date(), ".txt")
                filename2c <- paste0(pdfDir, .Platform$file.sep,"rliger_tsne_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2a)){
                    filename2a <- paste0(pdfDir, .Platform$file.sep,"rliger_pca_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                while(file.exists(filename2b)){
                    filename2b <- paste0(pdfDir, .Platform$file.sep,"rliger_umap_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                while(file.exists(filename2c)){
                    filename2c <- paste0(pdfDir, .Platform$file.sep,"rliger_tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                    i = i + 1;
                }
                write.csv(v$scData.combined@reductions$pca@cell.embeddings, file = filename2a)
                write.csv(v$scData.combined@reductions$umap@cell.embeddings, file = filename2b)
                write.csv(v$scData.combined@reductions$tsne@cell.embeddings, file = filename2c)
            })
        }
    })

    ##---------------CITE-Seq (Seurat)-------------------

    output$name.field <- renderUI({
        if(is.null(input$meta_mofa)){
            numericInput(inputId = "field",
                         label = "Field",
                         value = 1,
                         min = 1)
        }else{
            annoFile_mofa <- input$meta_mofa
            anno.data_mofa <- read.csv(annoFile_mofa$datapath[1], header = T, row.names = 1)
        }
    })

    observeEvent(input$loadButton2_a, {
        if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "Seurat" & input$scInput2 == "H5" & input$scAnalysis_type == "CITE-seq"){
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
                    print(sObj3@meta.data)
                    v$scDatab <- CreateAssayObject(counts = exp.data3.ab)
                    v$scDatat[["ADT"]] <- v$scDatab
                    print(Assays(v$scDatat))
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }

        else if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "MOFA" & input$scInput2 == "H5" & input$scAnalysis_type == "CITE-seq"){

            tpmFiles2 <- input$tpmFiles2
            #scH5 <- input$scH5
            annoFile_mofa <- input$meta_mofa
            names.field <- input$field
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
                    print(exp.data3)


                    if(!is.null(annoFile_mofa)){
                        anno.data_mofa <- read.csv(annoFile_mofa$datapath[1], header = T,
                                                   stringsAsFactors = FALSE, row.names=1)

                    }
                    incProgress(0.5, "Creating Seurat Object")

                    sObj3 <- CreateSeuratObject(exp.data3.rna,
                                                meta.data = anno.data_mofa,
                                                project = input$projName,
                                                names.field = names.field,
                                                names.delim = input$delim,
                                                is.expr = input$expThres,
                                                normalization.method = normMethod,
                                                min.genes = input$min.genes,
                                                min.cells = input$min.cells)

                    #mito.genes <- grep("^MT-", rownames(sObj@assays$RNA@data), ignore.case = TRUE, value = TRUE)
                    sObj3[["percent.mt"]] <- PercentageFeatureSet(sObj3, pattern = "^MT-")
                    print(sObj3@meta.data)
                    #incProgress(0.5, "Adding metadata")
                    #sObj <- AddMetaData(sObj, percent.mt, "percent.mt")
                    #if(!is.null(additional.ident1)){
                    #  sObj1 <- AddMetaData(sObj1, additional.ident1)
                    #}
                    v$scDatat <- sObj3
                    print(v$scDatat@meta.data)
                    v$scDatab <- CreateAssayObject(counts = exp.data3.ab)
                    v$scDatat[["ADT"]] <- v$scDatab
                    print(v$scDatat@meta.data)
                    print(Assays(v$scDatat))
                    shinyalert("Seurat object created", "Seurat object created for analysis.", type = "success", imageWidth = 10, imageHeight = 10)
                })
            }
        }
    })

    observe({if(input$scAnalysis_mult == "Seurat" & input$scAnalysis_type == "CITE-seq"){

        observeEvent(input$runUMAP3, {
            withProgress(message = "Running UMAP...", value = 0,{
                incProgress(0.5, message = "Running UMAP...")
                DefaultAssay(v$scDatat) <- "RNA"
                v$scDatat <- NormalizeData(v$scDatat)
                v$scDatat <- FindVariableFeatures(v$scDatat)
                v$scDatat <- ScaleData(v$scDatat)
                v$scDatat <- RunPCA(v$scDatat, verbose = FALSE)
                v$scDatat <- FindNeighbors(v$scDatat, dims = 1:input$dim.used1, nn.method = "rann")
                v$scDatat <- FindClusters(v$scDatat, resolution = input$clus.res1, verbose = FALSE)
                v$scDatat <- RunUMAP(v$scDatat, dims = 1:input$dim.used1)
                v$isUMAPdone4 <- TRUE
                UMAP_plot4a <- DimPlot(v$scDatat, label = TRUE, reduction = "umap")
                print(UMAP_plot4a)
            })
        })

        output$UMAPplot4_a <- renderPlotly({
            if(is.null(v$isUMAPdone4)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatat, reduction = "umap", label = T)
                })
            }
        })

        observeEvent(input$runTSNE3, {
            withProgress(message = "Running TSNE...", value = 0,{
                incProgress(0.5, message = "Running TSNE...")
                v$scDatat <- RunTSNE(v$scDatat, dims = 1:input$dim.used1)
                v$isTSNEdone4 <- TRUE
                TSNE_plot4a <- DimPlot(v$scDatat,  label = TRUE, reduction = "tsne")
                print(TSNE_plot4a)
            })
        })

        output$TSNEplot4_a <- renderPlotly({
            if(is.null(v$isTSNEdone4)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatat,  label = TRUE, reduction = "tsne")
                })
            }
        })

        observeEvent(input$Vis3, {
            withProgress(message = "Running TSNE...", value = 0,{
                incProgress(0.5, message = "Running TSNE...")
                DefaultAssay(v$scDatat) <- "ADT"
            })
        })

        output$vis.gene.select <- renderUI({
            if(is.null(v$scDatat)){
                return(NULL)
            }else{
                selectInput("vis.gene", label = "Gene to visualise",
                            choices = rownames(v$scDatat[["RNA"]]))
            }
        })

        output$vis.plot <- renderPlotly({
            if(is.null(v$scDatat)){
                return(NULL)
            }else{
                withProgress(message="Generating Feature Plot...", value=0, {
                    FeaturePlot(v$scDatat, input$vis.gene, cols = c("lightgrey", "darkgreen"))
                })
            }
        })

        output$vis.gene.select1 <- renderUI({
            if(is.null(v$scDatat)){
                return(NULL)
            }else{
                selectInput("vis.gene1", label = "Gene to visualise",
                            choices = rownames(v$scDatat[["ADT"]]))
            }
        })

        output$vis1.plot <- renderPlotly({
            if(is.null(v$scDatat)){
                return(NULL)
            }else{
                withProgress(message="Generating DEG Plot...", value=0, {
                    FeaturePlot(v$scDatat, input$vis.gene1, cols = c("lightgrey", "darkgreen"))
                })
            }
        })
    }
    })

    ##---------------scMultiomics (scRNA + scATAC-seq) using Seurat-------------------

    observe({if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "Seurat" & input$scAnalysis_type == "mRNA+scATAC-seq"){

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

        observeEvent(input$loadButton2_b, {
            if(input$Module == "Single cell multiomics Analysis" & input$scAnalysis_mult == "Seurat" & input$scInput2 == "H5" & input$scAnalysis_type == "mRNA+scATAC-seq"){
                tpmFiles3b <- input$tpmFiles3b
                #fragFile <- input$fragFiles
                #names.field3 <- input$field3
                #if(!is.null(fragFile)){
                #    frag.data1 <- fread(fragFile$datapath)
                #}
                #print(frag.data1)
                if (is.null(tpmFiles3b)){
                    v$scDatan <- NULL
                }else{
                    withProgress(message="Loading and Processing Data...", value=0, {
                        print(tpmFiles3b$datapath)
                        print(tpmFiles3b$name)
                        print(file.exists(paste(tpmFiles3b$datapath[1], "/", tpmFiles3b$name[1], sep="")))
                        exp.data4 <- Read10X_h5(tpmFiles3b$datapath)
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
                    })
                }
            }
        })

        output$nCount_ATAC.plot <- renderPlotly({
            if(is.null(v$scDatan)){
                plotly_empty()
            }else{
                VlnPlot(v$scDatan, "nCount_ATAC")
            }
        })

        output$nCount_RNA.plot <- renderPlotly({
            if(is.null(v$scDatan)){
                plotly_empty()
            }else{
                VlnPlot(v$scDatan, "nCount_RNA")
            }
        })

        output$percent.mt.plot <- renderPlotly({
            if(is.null(v$scDatan)){
                plotly_empty()
            }else{
                VlnPlot(v$scDatan, "percent.mt")
            }
        })

        observeEvent(input$doSCTransform_multi, {
            withProgress(message = "Running scTransform...", value = 0,{
                incProgress(0.5, message = "Running scTransform...")
                DefaultAssay(v$scDatan) <- "RNA"
                v$scDatan <- SCTransform(v$scDatan, verbose = FALSE)
            })
        })

        observeEvent(input$runUMAP4, {
            withProgress(message = "Running UMAP...", value = 0,{
                incProgress(0.5, message = "Running UMAP...")
                v$scDatan <- RunPCA(v$scDatan, verbose = FALSE)
                v$scDatan <- RunUMAP(v$scDatan, dims = 1:input$dim.used1, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

                # ATAC analysis
                # We exclude the first dimension as this is typically correlated with sequencing depth
                DefaultAssay(v$scDatan) <- "ATAC"
                v$scDatan <- RunTFIDF(v$scDatan)
                v$scDatan <- FindTopFeatures(v$scDatan)
                v$scDatan <- RunSVD(v$scDatan)
                v$scDatan <- RunUMAP(v$scDatan, reduction = 'lsi', dims = 2:input$dim.used2, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
                v$scDatan <- FindMultiModalNeighbors(v$scDatan, reduction.list = list("pca", "lsi"), dims.list = list(1:input$dim.used2, 2:input$dim.used2))
                v$scDatan <- FindClusters(v$scDatan, resolution = input$clus.res2, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
                v$scDatan <- RunUMAP(v$scDatan, dims = 1:input$dim.used2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
                v$isUMAPdone5 <- TRUE
                UMAP_plot5a <- DimPlot(v$scDatan, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
                UMAP_plot5b <- DimPlot(v$scDatan, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
                UMAP_plot5c <- DimPlot(v$scDatan, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
                print(UMAP_plot5a)
                print(UMAP_plot5b)
                print(UMAP_plot5c)
            })
        })

        output$UMAPplot5_a <- renderPlotly({
            if(is.null(v$isUMAPdone5)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatan, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
                })
            }
        })

        output$UMAPplot5_b <- renderPlotly({
            if(is.null(v$isUMAPdone5)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatan, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
                })
            }
        })

        output$UMAPplot5_c <- renderPlotly({
            if(is.null(v$isUMAPdone5)){
                plotly_empty()
            }else{
                withProgress(message="Generating UMAP Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatan, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
                })
            }
        })

        observeEvent(input$runTSNE4, {
            withProgress(message = "Running TSNE...", value = 0,{
                incProgress(0.5, message = "Running TSNE...")
                v$scDatan <- RunTSNE(v$scDatan, dims = 1:input$dim.used1, reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')

                # ATAC analysis
                # We exclude the first dimension as this is typically correlated with sequencing depth
                DefaultAssay(v$scDatan) <- "ATAC"
                v$scDatan <- RunTSNE(v$scDatan, reduction = 'lsi', dims = 2:input$dim.used2, reduction.name = "tsne.atac", reduction.key = "atacTSNE_")
                v$scDatan <- RunTSNE(v$scDatan, dims = 1:input$dim.used2, nn.name = "weighted.nn", reduction.name = "wnn.tsne", reduction.key = "wnnTSNE_")
                v$isTSNEdone5 <- TRUE
                UMAP_plot5a <- DimPlot(v$scDatan, reduction = "tsne.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
                UMAP_plot5b <- DimPlot(v$scDatan, reduction = "tsne.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
                UMAP_plot5c <- DimPlot(v$scDatan, reduction = "wnn.tsne", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
                print(UMAP_plot5a)
                print(UMAP_plot5b)
                print(UMAP_plot5c)
            })
        })

        output$TSNEplot5_a <- renderPlotly({
            if(is.null(v$isTSNEdone5)){
                plotly_empty()
            }else{
                withProgress(message="Generating TSNE Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatan, reduction = "tsne.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
                })
            }
        })

        output$TSNEplot5_b <- renderPlotly({
            if(is.null(v$isTSNEdone5)){
                plotly_empty()
            }else{
                withProgress(message="Generating TSNE Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatan, reduction = "tsne.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
                })
            }
        })

        output$TSNEplot5_c <- renderPlotly({
            if(is.null(v$isTSNEdone5)){
                plotly_empty()
            }else{
                withProgress(message="Generating TSNE Plot of integrated dataset...", value=0, {
                    DimPlot(v$scDatan, reduction = "wnn.tsne", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
                })
            }
        })
    }
    })

    ##---------------CITE-Seq (MOFA)-------------------##

    observe({if(input$scAnalysis_mult == "MOFA" & input$scAnalysis_type == "CITE-seq"){

        observeEvent(input$data_overview, {
            withProgress(message = "Plotting data overview...", value = 0,{
                incProgress(0.5, message = "Running UMAP...")
                DefaultAssay(v$scDatat) <- "RNA"
                v$scDatat <- NormalizeData(v$scDatat, normalization.method = "LogNormalize", assay = "RNA")
                v$scDatat <- ScaleData(v$scDatat, do.center = TRUE, do.scale = FALSE)
                DefaultAssay(v$scDatat) <- "ADT"
                v$scDatat <- NormalizeData(v$scDatat, normalization.method = "LogNormalize", assay = "ADT")
                v$scDatat <- ScaleData(v$scDatat, do.center = TRUE, do.scale = FALSE)
                v$scDatat <- FindVariableFeatures(v$scDatat, selection.method = "vst", nfeatures = 5000, assay = "RNA", verbose = FALSE)
                v$scDatat <- FindVariableFeatures(v$scDatat, selection.method = "vst", nfeatures = 5000, assay = "ADT", verbose = FALSE)
                mofa <- create_mofa(v$scDatat, assays = c("RNA","ADT"))
                print(mofa)
                model_opts <- get_default_model_options(mofa)
                model_opts$num_factors <- as.numeric(input$factor)
                mofa <- prepare_mofa(mofa, model_options = model_opts)
                mofa <- run_mofa(mofa)
                print(v$scDatat@meta.data)
                samples_metadata(mofa) <- v$scDatat@meta.data %>%  tibble::rownames_to_column("sample") %>%   as.data.table
                print(mofa)
                v$scDatam <- mofa
                v$ismofadone <- TRUE
            })
        })

        output$mofa_a <- renderPlot({
            if(is.null(v$ismofadone)){
                plotly_empty()
            }else{
                withProgress(message="Plotting variance decomposition...", value=0, {
                    plot_data_overview(v$scDatam)
                })
            }
        })

        output$mofa_b <- renderPlot({
            if(is.null(v$ismofadone)){
                plotly_empty()
            }else{
                withProgress(message="Plotting variance decomposition...", value=0, {
                    plot_factor_cor(v$scDatam)
                })
            }
        })

        # output$mofa_c <- renderPlot({
        #   if(is.null(v$ismofadone)){
        #     plotly_empty()
        #   }else{
        #     withProgress(message="Plotting variance decomposition...", value=0, {
        #       plot_variance_explained(v$scDatam, max_r2=input$factor)
        #     })
        #   }
        #  })


        #output$factor_a <- renderPlot({
        #  if(is.null(v$ismofadone)){
        #    plotly_empty()
        # }else{
        #    withProgress(message="Generating Factor plot...", value=0, {
        #      plot_factor(v$scDatam, factors = input$factor, color_by = input$factors, dodge = TRUE, add_violin = TRUE)
        #     })
        # }
        #  })

        #output$factor_b <- renderPlot({
        # if(is.null(v$ismofadone)){
        #   plotly_empty()
        #  }else{
        #   withProgress(message="Generating Factor plot...", value=0, {
        #      plot_weights(v$scDatam,
        #                  view = input$view,
        #                   factor = input$factor,
        #                   nfeatures = input$nfeatures,     # Top number of features to highlight
        #                  scale = T)           # Scale weights from -1 to 1

        #   })
        #   }
        # })

        observeEvent(input$umap_mofa, {
            withProgress(message = "Plotting UMAP...", value = 0,{
                incProgress(0.5, message = "Running UMAP...")
                v$scDatam <- run_umap(v$scDatam, factors = "all", n_neighbors = as.numeric(input$n_neighbors))
                v$ismofaumapdone <- TRUE
            })
        })

        output$umap_mofa_mult <- renderPlot({
            if(is.null(v$ismofaumapdone)){
                plotly_empty()
            }else{
                withProgress(message="Generating Factor plot...", value=0, {
                    plot_dimred(v$scDatam, method = "UMAP", color_by = "Phenotype", label = TRUE)
                })
            }
        })
    }
    })


    ##---------------Flow cytometry Analysis-------------------##

    ## Load cytofkit RData object
    observeEvent(input$goButton, {
        cytofkitObj <- input$cytofkitObj
        if (is.null(cytofkitObj)){
            v$data <- NULL
        }else{
            # browser()
            message(cytofkitObj$datapath)
            load(cytofkitObj$datapath)
            #print(NROW(analysis_results))
            v$data <- analysis_results
            #print(v$data)

            if(is.null(v$data$projectName)){
                v$data$projectName <- "cytofkit_shinyAPP_output"
            }

            if(!is.null(v$data$progressionRes)){
                ## default the first cluster results are used for progression analysis
                p$progressionCluster <- names(v$data$clusterRes)[1]
            }


            # Need modification later
            # currently doesn't update sampleInfo with v$data$sampleInfo
            v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                                       cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                                       stringsAsFactors = FALSE)
            v$data$sampleInfo <- v$sampleInfo
            v$sample_ready = TRUE
            v$ori_sampleInfo = v$sampleInfo
            v$ori_data = v$data
            updateSelectInput(session, "sample_info_selection", choices = colnames(v$data$meta_data))
        }
    })

    ## For user, set roots option to your server directory
    roots <- c(a=a)
    shinyFileChoose(input, 'serverObj', session = session, roots = roots, filetypes = "RData")

    observeEvent(input$serverObj, {
        inServer <- parseFilePaths(roots= roots, input$serverObj)
        print(inServer$datapath)
        load(as.character(inServer$datapath))
        v$data <- analysis_results
        if(is.null(v$data$projectName)){
            v$data$projectName <- "cytofkit_shinyAPP_output"
        }
        if(!is.null(v$data$progressionRes)){
            ## default the first cluster results are used for progression analysis
            p$progressionCluster <- names(v$data$clusterRes)[1]
        }
        # Need modification later
        # currently doesn't update sampleInfo with v$data$sampleInfo
        v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                                   cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                                   stringsAsFactors = FALSE)
        v$data$sampleInfo <- v$sampleInfo
        #print(v$data$sampleInfo)
    })

    output$rdata_desc <- renderText({
        if(is.null(v$data)){
            paste0("No .RData loaded yet")
        }else if(!length(session$clientData$url_search) == 0){
            return(NULL)
        }else{
            paste0("Loaded: ", v$data$resultDir, "/", v$data$projectName, ".RData")
        }
    })

    observeEvent(input$reset, {
        analysis_results <- NULL
        session$reload()
        print("Reset done")
    })

    output$selectAll <- renderUI({
        if(is.null(v$data) || is.null(v$sampleInfo)){
            return(NULL)
        }else{
            checkboxInput('selectDeselectAll', label = "Select/Deselect All", value = TRUE)
        }
    })

    #### only initial UI once.
    observeEvent(v$sample_ready, {
        output$sampleSelect <- renderUI({
            if(is.null(isolate(v$data)) || is.null(isolate(v$sampleInfo))){
                return(NULL)
            }else{
                # browser()
                sampleNames <- isolate(unique(as.character(v$sampleInfo$cellSample)))
                v$sample_choices = sampleNames
                v$sample_selected_index = 1:length(sampleNames)
                v$sample_init = TRUE
                checkboxGroupInput(inputId = 'samples', label = NULL,
                                   choices = sampleNames, selected = sampleNames)
            }
        })
        if(v$sample_init){
            # browser()
            update_samples()
        }
    })


    observeEvent(unique(as.character(v$sampleInfo$cellSample)), {
        # browser()
        if(v$sample_init){
            # browser()
            sampleNames = unique(as.character(v$sampleInfo$cellSample))
            ### get selected item index
            # selected_index = match(input$samples, v$sample_choices)
            # v$sample_selected_index = selected_index
            v$sample_selected_index = 1:length(sampleNames)
            v$sample_choices = sampleNames
            updateCheckboxGroupInput(session, 'samples', choices = v$sample_choices)
            # browser()
            v$sample_update = !(v$sample_update)
            # updateCheckboxGroupInput(session, 'samples', selected = v$sample_choices[selected_index])
        }
    })

    observeEvent(v$sample_update, {
        # browser()
        if (length(v$sample_selected_index) > length(v$sample_choices)){
            v$sample_selected_index = 1:length(v$sample_choices)
        }
        updateCheckboxGroupInput(session, 'samples', selected = v$sample_choices[v$sample_selected_index])
    })

    update_samples = function(){
        # v$sampleInfo$cellSample = v$data$meta_data[match(v$data$sampleInfo_ori$cellSample, rownames(v$data$meta_data))
        # , input$sample_info_selection]
        # browser()
        allSamp <- input$selectDeselectAll
        sampleNames <- isolate(unique(as.character(v$sampleInfo$cellSample)))
        if(allSamp == TRUE){
            updateCheckboxGroupInput(session, 'samples', selected = sampleNames)
        }else{
            updateCheckboxGroupInput(session, 'samples', selected = character(0))
        }
        # browser()
    }
    observeEvent(input$selectDeselectAll, {
        # browser()
        update_samples()
    })

    observeEvent(input$sample_info_selection, {

        if(!is.null(v$data) && !is.null(input$sample_info_selection)){
            # browser()

            temp_name1 = sapply(1:length(v$data$sampleNames), function(i){
                v$data$sampleNames[[i]][1]
            })
            temp_name2 = sapply(1:length(v$data$sampleNames), function(i){
                v$data$sampleNames[[i]][2]
            })
            new_sample_names = as.character(v$data$meta_data[temp_name1[order(temp_name2)], input$sample_info_selection])
            reset_sample()
            rename_sample(new_sample_names)
            # v$sample_update = !v$sample_update
            # browser()
            # v$data$sampleNames
            # v$data$sampleNames = v$data$mea_data[unlist(v$data$sampel_name_ori), input$sample_info_selection]
            # v$data$sampleInfo$cellSample = v$data$meta_data[match(v$data$sampleInfo_ori$cellSample, rownames(v$data$meta_data))
            # , input$sample_info_selection]
            # update_samples()
            # analysis_results = v$data
            # analysis_results$sampleInfo_ori = analysis_results$sampleInfo
            # save(analysis_results, file = "test2.Rdata")
        }
    })

    observe({
        if(!is.null(v$data) && !is.null(v$sampleInfo) && !is.null(input$samples)){
            x <- input$samples
            sampleNames <- unique(as.character(v$sampleInfo$cellSample))
            if(length(x) == 0){
                x <- character(0)
                updateCheckboxInput(session, 'selectDeselectAll', value = FALSE)
            }
            if(length(x) == length(sampleNames)){
                updateCheckboxInput(session, 'selectDeselectAll', value = TRUE)
            }
        }
    })

    output$summaryText1 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", nrow(v$data[[1]]), " cells x ", ncol(v$data[[1]]), " markers")
    })

    output$summaryText2 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", paste(names(v$data$clusterRes), collapse = " | "))
    })

    output$summaryText3 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", paste(v$data$visualizationMethods, collapse =  " | "))
    })

    output$summaryText4 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", ifelse(is.null(v$data$progressionRes), "NULL",
                             sub("_[0-9]*$", "", colnames(v$data$progressionRes$progressionData)[1])))
    })

    output$summaryText5 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", paste(v$data$dimRedMarkers, collapse =  " | "))
    })

    # ## Save and parse cytofkit RData object
    # observeEvent(input$saveButton, {
    #   if (!is.null(v$data)){
    #     withProgress(message='Saving Results ', value=0, {
    #       ## check results saving path
    #       if(is.null(v$data$resultDir) || !dir.exists(v$data$resultDir)){
    #         v$data$resultDir <- path.expand("~")  ## default save to home if not specified
    #       }
    #       saveToFCS <- input$saveFCS
    #       if(is.null(v$data$rawFCSdir)){
    #         saveToFCS <- FALSE
    #         warning("Path for original FCS files is not provided,
    #                 data cannnot be saved to new copies of FCS files.")
    #       }else if(!dir.exists(v$data$rawFCSdir)){
    #         saveToFCS <- FALSE
    #         warning(paste0("Path for original FCS files doesn't exist,
    #                        data cannnot be saved to new copies of FCS files.",
    #                        "Please check path: ", v$data$rawFCSdir))
    #       }
    #
    #       ## NOTE: if samples are regrouped, then new FCS file cannot be saved
    #       incProgress(1/2, message = paste0("To ", v$data$resultDir))
    #       v$data$sampleInfo <- v$sampleInfo
    #       analysis_results <<- v$data
    #       cytof_writeResults(analysis_results,
    #                          saveToRData = input$saveRData,
    #                          saveToFCS = saveToFCS,
    #                          saveToFiles = input$saveCsv)
    #       incProgress(1/2)
    #       ## open the results directory
    #       opendir(v$data$resultDir)
    #     })
    #     }
    # })

    output$saveButton = downloadHandler(
        filename = function() {
            paste0(input$project_name, '_result.zip')
        },
        content = function(file) {
            # browser()
            withProgress(message='Saving Results ', value=0, {

                # ## check results saving path
                # if(is.null(v$data$resultDir) || !dir.exists(v$data$resultDir)){
                #   v$data$resultDir <- path.expand("~")  ## default save to home if not specified
                # }
                # saveToFCS <- input$saveFCS
                # if(is.null(v$data$rawFCSdir)){
                #   saveToFCS <- FALSE
                #   warning("Path for original FCS files is not provided,
                #           data cannnot be saved to new copies of FCS files.")
                # }else if(!dir.exists(v$data$rawFCSdir)){
                #   saveToFCS <- FALSE
                #   warning(paste0("Path for original FCS files doesn't exist,
                #                  data cannnot be saved to new copies of FCS files.",
                #                  "Please check path: ", v$data$rawFCSdir))
                # }

                # browser()
                res_folder = paste0(input$project_name, "_results")
                if(!dir.exists(res_folder)){
                    dir.create(res_folder, recursive = TRUE)
                }
                ## NOTE: if samples are regrouped, then new FCS file cannot be saved
                incProgress(1/2)
                v$data$sampleInfo <- v$sampleInfo
                analysis_results <<- v$data
                cytof_writeResults(analysis_results,
                                   saveToRData = input$saveRData,
                                   saveToFCS = FALSE,
                                   saveToFiles = input$saveCsv
                                   , resultDir = res_folder)
                incProgress(1/2)
                # ## open the results directory
                # opendir(v$data$resultDir)
                files = list.files(path = res_folder, pattern = ".*", full.names = TRUE)
                utils::zip(file, files)
                #stopApp(returnValue = invisible())
            })

        }
    )

    output$reportButton = downloadHandler(
        filename = function() {
            paste0(input$project_name, '_report.html')
        },
        content = function(file) {
            # browser()
            withProgress(message='Generating report ', value=0, {

                # browser()
                # library(ezTools)
                # library(cytofkit2)

                #### suppress warning
                options(warn=-1)
                analysis_results = v$data

                rownames(analysis_results$sampleInfo) = analysis_results$sampleInfo$cellID
                analysis_results$clusterRes$Rphenograph = as.data.frame(analysis_results$clusterRes$Rphenograph)
                colnames(analysis_results$clusterRes$Rphenograph) = "Rphenograph Cluster"
                #### set output markdown file
                output_file_name = './cytofkit_report.RMD'

                create_script = function(title = "", script = "", values = NULL, fig_width = NULL, fig_height = NULL
                                         , chunk_option = NULL){
                    # browser()
                    param_names = as.list(match.call()$script)
                    param_names = param_names[-1]
                    r_option = "```{r, echo = FALSE"
                    if (!is.null(chunk_option)) {
                        r_option = paste0(r_option, ", ", chunk_option)
                    }
                    if (!is.null(fig_width)) {
                        r_option = paste0(r_option, ", fig.width = ", fig_width)
                    }
                    if (!is.null(fig_height)) {
                        r_option = paste0(r_option, ", fig.height = ", fig_height)
                    }
                    r_option = paste0(r_option, "}")
                    res = paste(title, r_option, paste(param_names, collapse = "\n"), "```", sep = "\n")
                    if (!is.null(values)) {
                        for (i in 1:length(values)) {
                            res = gsub(values[i], paste0("\"", eval_string(values[i], envir = parent.frame()), "\""), res, fixed = TRUE)
                        }
                    }
                    res
                }

                substitute_script = function(script_file = "", script_id = "", scripts = "", output_file = ""){
                    # browser()
                    all_script = read_string(script_file)
                    script_id = paste0("####* ", script_id, " *####")
                    #### substitue script
                    script = paste(scripts, collapse = "\n")
                    script_res = gsub(script_id, script, all_script, fixed = TRUE)
                    write_string(script_res, output_name = output_file)
                }

                create_dr_scripts = function(analysis_results){
                    # browser()
                    dr_names = names(analysis_results$dimReducedRes)
                    scripts = lapply(1:length(dr_names), function(i){
                        if (dr_names[i] == "tsne") {
                            return("")
                            # temp_name = 'tSNE'
                        } else if (dr_names[i] == "umap") {
                            temp_name = 'UMAP'
                        } else {
                            temp_name = dr_names[i]
                        }

                        analysis_results$temp1 <<- str_replace_all(analysis_results$dimRedMarkers, "<", "&lt;")
                        analysis_results$temp2 <<- str_replace_all(analysis_results$temp1, ">", "&gt;")
                        script1 = create_script(paste0("## ", temp_name, " plot color by sample\n", "Based on the markers: "
                                                       , paste0(analysis_results$temp2, collapse = ', ')), {
                                                           plot_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                                                                        , analysis_results$sampleInfo[, "cellSample", drop = FALSE]) + coord_fixed()
                                                       }, values = c("dr_names[i]"))
                        script2 = create_script(paste0("## ", temp_name, " plot color by sample (splitted version)"), {
                            p = plot_split_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                                                   , analysis_results$sampleInfo[, "cellSample", drop = FALSE], ncol = 2, show_legend = FALSE)
                            p[[1]]
                        }, values = c("dr_names[i]"), fig_width = 8, fig_height = 20)
                        forplot3 <<- as.data.frame(cbind(analysis_results$dimReducedRes[[dr_names[i]]], analysis_results$expressionData[,analysis_results$dimRedMarkers, drop = FALSE]))
                        script3 = create_script(paste0("## ", temp_name, " plot color by marker expression"), {
                            for (x in 1:length(analysis_results$dimRedMarkers)){
                                print(ggplot(forplot3, aes(x = umap_1, y = umap_2)) + geom_point(aes(color = get(colnames(forplot3)[x+2]))) +
                                          scale_color_gradientn(colors = c("#00AFBB", "#E7B800", "#FC4E07")) +
                                          theme(legend.position = "right") +
                                          coord_fixed() +
                                          labs(colour = colnames(forplot3)[x+2]))
                            }
                        }, values = c("dr_names[i]"))
                        return(c(script1, script2, script3))
                    })
                    unlist(scripts)
                }

                create_cluster_scripts = function(analysis_results){
                    # browser()
                    dr_names = names(analysis_results$dimReducedRes)
                    scripts = lapply(1:length(dr_names), function(i){
                        if (dr_names[i] == "tsne") {
                            return("")
                            # temp_name = 'tSNE'
                        } else if (dr_names[i] == "umap") {
                            temp_name = 'UMAP'
                        } else {
                            temp_name = dr_names[i]
                        }
                        script1 = create_script(paste0("## ", temp_name, " plot color by cluster"), {
                            plot_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                                         , analysis_results$clusterRes$Rphenograph) + coord_fixed()
                        }, values = c("dr_names[i]"))
                        # fig_height = ceiling(length(unique(analysis_results$clusterRes$Rphenograph[, 1]))/4)*2.5
                        # script2 = create_script(paste0("## ", temp_name, " plot color by cluster (splitted version)"), {
                        #   p = plot_split_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                        #                          , analysis_results$clusterRes$Rphenograph, ncol = 4, show_legend = FALSE)
                        #   p[[1]]
                        # }, values = c("dr_names[i]"), fig_height = fig_height, fig_width = 9)
                        # return(c(script1, script2))
                        return(c(script1))
                    })
                    script2 = create_script(paste0("## ", " Percentage heatmap"), {

                        temp = as.data.frame(analysis_results$clusterRes$Rphenograph)
                        rownames(analysis_results$sampleInfo) = analysis_results$sampleInfo[, 1]
                        temp = ezcbind(temp, analysis_results$sampleInfo[, "cellSample", drop = FALSE])
                        freq = as.data.frame.matrix(table(temp))
                        freq = freq/rowSums(freq)
                        pheatmap(freq, silent = FALSE)
                        # print(p)
                    })
                    unlist(scripts)
                    c(scripts, script2)
                }

                create_markers_script = function(analysis_results){
                    dr_names = names(analysis_results$dimReducedRes)
                    scripts = lapply(1:1, function(i){
                        if (dr_names[i] == "tsne") {
                            return("")
                            # temp_name = 'tSNE'
                        } else if (dr_names[i] == "umap") {
                            temp_name = 'UMAP'
                        } else {
                            temp_name = dr_names[i]
                        }

                        # selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
                        # fig_height = ceiling(length(unique(selected_markers))/6)*2.5
                        script2 = create_script(paste0("## ", temp_name, " plot color by cluster (splitted version)"), {
                            selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
                            marker_list = ez_chunk(selected_markers, ceiling(length(selected_markers)/3/7))
                            lapply(1:length(marker_list), function(k){
                                all_plots = lapply(1:length(marker_list[[k]]), function(j){
                                    plot_scatter(analysis_results$dimReducedRes[[dr_names[i]]]
                                                 , analysis_results$expressionData[, marker_list[[k]][j], drop = FALSE]
                                                 , colors = c("#BEBEBE",brewer.pal(9,"Reds")), color_as_factor = FALSE) +
                                        coord_fixed()
                                })
                                p = plot_grid(plotlist = all_plots, ncol = 3)
                                p
                            })
                        }, values = c("dr_names[i]"), fig_height = 15, fig_width = 9
                        , chunk_option = "results='hide', fig.keep='all', message=FALSE")
                        return(c(script2))
                    })
                    unlist(scripts)
                }

                create_express_heatmap = function(analysis_results){
                    res = create_script(paste0("## Expression heatmap"), {
                        selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
                        expression = as.data.frame(analysis_results$expressionData[, selected_markers, drop = FALSE])
                        expression$clusters = analysis_results$clusterRes$Rphenograph[, 1]
                        cluster_expression = fast_aggr(expression, ncol(expression))
                        dt = seurat_sacle_data(cluster_expression)
                        pheatmap(dt)
                    }, chunk_option = "results='hide', fig.keep='all', message=FALSE")
                    res


                }

                create_expression_histogram = function(analysis_results){
                    script2 = create_script(paste0("## Expression histograms"), {
                        dt = as.data.frame(analysis_results$expressionData)
                        selected_markers = which(!(grepl("NA", colnames(analysis_results$expressionData))))
                        selected_markers = colnames(analysis_results$expressionData)[selected_markers]
                        marker_list = ez_chunk(selected_markers, ceiling(length(selected_markers)/3/6))
                        invisible(lapply(1:length(marker_list), function(j){
                            print(stackDenistyPlot(dt, marker_list[[j]], stackFactor = analysis_results$clusterRes$Rphenograph[, 1]))
                        }))
                    }, fig_width = 9, fig_height = 15, chunk_option = "results='hide', fig.keep='all', message=FALSE")
                    return(c(script2))
                }

                create_abstract_scripts = function(analysis_results){
                    paste0("The project \"", analysis_results$projectName, "\" has ")
                    sample_num = length(unique(analysis_results$sampleInfo$cellSample))
                    cell_num = nrow(analysis_results$sampleInfo)
                    cluster_num = length(unique(analysis_results$clusterRes$Rphenograph$`Rphenograph Cluster`))
                    res = paste0("There are ", sample_num, " sample", ifelse(sample_num > 1, "s", "")
                                 , ", ", cluster_num, " cluster", ifelse(cluster_num > 1, "s", "")
                                 , ", ", cell_num, " cell", ifelse(cluster_num > 1, "s", "")
                                 , " in the project \"", analysis_results$projectName, "\".\n\n")
                    res = paste0(res, "  The dimensionality reduction, clustering and markers expression analysis were conducted on the project data.")
                    res
                }


                dr_scripts = create_dr_scripts(analysis_results)
                substitute_script('./pdf_report_template.Rmd', script_id = "DR analysis"
                                  , scripts = dr_scripts
                                  , output_file = output_file_name)

                abstract_scripts = create_abstract_scripts(analysis_results)
                substitute_script(output_file_name, script_id = "Abstract"
                                  , scripts = abstract_scripts
                                  , output_file = output_file_name)
                # render(output_file_name)

                cluster_scripts = create_cluster_scripts(analysis_results)
                substitute_script(output_file_name, script_id = "Cluster analysis"
                                  , scripts = cluster_scripts
                                  , output_file = output_file_name)

                markers_scripts = create_markers_script(analysis_results)
                substitute_script(output_file_name, script_id = "Expression on DR"
                                  , scripts = markers_scripts
                                  , output_file = output_file_name)

                expression_heatmap_scrip = create_express_heatmap(analysis_results)
                substitute_script(output_file_name, script_id = "Expression heatmap"
                                  , scripts = expression_heatmap_scrip
                                  , output_file = output_file_name)

                expression_histogram_script = create_expression_histogram(analysis_results)
                substitute_script(output_file_name, script_id = "Expression histogram"
                                  , scripts = expression_histogram_script
                                  , output_file = output_file_name)

                render(output_file_name, output_file = file)
                # output_pdf = basename(output_file_name)
                # output_pdf = paste0(get_file_name(output_pdf, with_ext = FALSE), ".pdf")
                # system(paste0("cp ", output_file_name, " ", file))

            })

        }
    )

    #
    # observeEvent(input$OpenDir, {
    #   pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
    #   if(dir.exists(pdfDir)){
    #     opendir(pdfDir)
    #   }else{
    #     stop("PDF not created yet!")
    #   }
    # })

    output$logo <- renderImage({
        return(list(
            src = "vignettes/logo.png",
            contentType = "image/png",
            alt = "Singapore Immunology Network"
        ))
    }, deleteFile = FALSE)

    ##------------------------------Cluster Panel------------------------------

    ##-----cluster plot-----
    output$C_PlotMethod <- renderUI({
        if(is.null(v$data) || is.null(visualizationMethods())){
            return(NULL)
        }else{
            selectInput('c_PlotMethod', 'Visualization Method:', choices = visualizationMethods(),
                        selected = visualizationMethods()[1], width = "100%")
        }
    })

    output$C_PlotFunction <- renderUI({
        if(is.null(v$data) || is.null(visualizationFunctions())){
            return(NULL)
        }else{
            selectInput('c_PlotFunction', 'Cluster By:', choices = visualizationFunctions(),
                        selected = visualizationFunctions()[1], width = "100%")
        }
    })

    output$C_markerSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            markerNames <- colnames(v$data$expressionData)
            markerNames <- markerNames[order(markerNames)]
            checkboxGroupInput('c_markerSelect', strong('Select Markers:'),
                               markerNames, selected = markerNames, inline = TRUE)
        }
    })

    output$C_clusterSelect <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_PlotFunction))
            return(NULL)
        if(input$c_PlotFunction %in% c("Sample", "Density","None")){
            return(NULL)
        }else{
            clusterMethod <- input$c_PlotFunction
            clusterIDs <- sort(unique(v$data$clusterRes[[clusterMethod]]))
            selectizeInput('c_clusterSelect', 'Clusters Filter:',
                           choices = clusterIDs, selected = clusterIDs,
                           multiple = TRUE, width = "100%")
            # checkboxGroupInput('p_clusterSelect', strong('Select Clusters:'),
            #                    clusterIDs, selected = clusterIDs, inline = TRUE)
        }
    })

    ## Complex dependencies here: --> (depends on)
    ## C_ScatterPlotInput --> c_PlotMethod + c_clusterSelect
    ## c_clusterSelect --> c_PlotMethod
    ## carefull checkings are applied to solve concurrency conflicts
    C_ScatterPlotInput <- function(){
        if(is.null(v$data) || is.null(input$c_PlotMethod) ||
           is.null(input$c_PlotFunction) || is.null(input$c_clusterSelect)){
            return(NULL)
        }else if(!all(input$c_clusterSelect %in% v$data$clusterRes[[input$c_PlotFunction]]) &&
                 !(input$c_PlotFunction %in% c("Sample", "Density","None"))){
            return(NULL)
        }else{
            # browser()

            withProgress(message="Generating Cluster Scatter Plot", value=0, {
                if(input$c_PlotFunction %in% c("Sample", "Density", "None")){
                    clusterSelect <- NULL
                    clusterColor <- NULL
                }else{
                    clusterSelect <- input$c_clusterSelect
                    clusterMethod <- input$c_PlotFunction
                    if(!is.null(c$clusterCol[[clusterMethod]])){
                        clusterColor <- c$clusterCol[[clusterMethod]]
                    }else{
                        cluster_num <- length(unique(v$data$clusterRes[[clusterMethod]]))
                        clusterColor <- rainbow(cluster_num)
                    }
                }
                temp_sample_names = lapply(1:length(v$data$sampleNames), function(i){
                    v$data$sampleNames[[i]][length(v$data$sampleNames[[i]])]
                })
                #if(!all(temp_sample_names %in% input$samples)){
                #  return(NULL)
                #}
                # browser()
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$c_PlotMethod,
                                  plotFunction = input$c_PlotFunction,
                                  pointSize = input$C_PointSize,
                                  addLabel = input$C_addLabel,
                                  labelSize = input$C_LabelSize,
                                  sampleLabel = FALSE,
                                  FlowSOM_k = input$C_FlowSOM_k,
                                  selectCluster = clusterSelect,
                                  selectSamples = input$samples,
                                  facetPlot = input$C_facetPlot,
                                  labelRepel = input$C_labelRepel,
                                  removeOutlier = TRUE,
                                  clusterColor = clusterColor)
                incProgress(1/2)
                plot(gp)
                incProgress(1/2)
            })
        }
    }

    output$C_ScatterPlot <- renderPlot({
        C_ScatterPlotInput()
    }, height = 900, width = 950)

    output$PDFClusterPlot = downloadHandler(
        filename = function() {
            filename1 <- paste0(input$project_name, "_shinyAPP_Clusterplot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name, "_shinyAPP_Clusterplot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                # browser()
                withProgress(message="Downloading Clusterplot PDF files...", value=0, {
                    print(getwd())

                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    C_ScatterPlotInput()
                    dev.off()
                })
            }

        }
    )

    observeEvent(input$PDFClusterPlot, {

    })


    ##----- change cluster colour -----
    output$C_colourCluster <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            clusterMethods <- c(names(v$data$clusterRes))
            #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
            selectInput('c_colourCluster', 'Choose Cluster to Change the Colour :',
                        choices = clusterMethods,
                        selected = clusterMethods[1], width = "50%")
        }
    })

    ## currently use 100 as a limit for cluster numbers
    ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ##
    lapply(1:100, function(i) {
        output[[paste0('Cluster_', i, "_col")]] <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_colourCluster)){
                return(NULL)
            }

            clusters <- v$data$clusterRes[[input$c_colourCluster]]
            clusterLabel <- levels(as.factor(clusters))
            if(is.null(c$clusterCol[[input$c_colourCluster]])){
                clusterColor <- rainbow(length(unique(clusters)))
            }else{
                clusterColor <- c$clusterCol[[input$c_colourCluster]]
            }

            if (i <= length(clusterLabel)){
                x <- clusterLabel[i]
                colourpicker::colourInput(inputId=paste0('cluster_', i, '_col'),
                                          label=paste0('Cluster ', x," Colour :"),
                                          value = clusterColor[i], showColour = "both",
                                          palette = "square")
            }
        })
    })

    ## update cluster color
    observeEvent(input$C_updateClusterColor, {
        if(!is.null(v$data) && !is.null(input$c_colourCluster)){
            clusterMethod <- input$c_colourCluster
            clusterVec<- v$data$clusterRes[[clusterMethod]]
            clusters <- levels(as.factor(clusterVec))
            clusterCols <- NULL
            for (i in 1:length(clusters)){
                clusteri <- clusters[i]
                iCol <- input[[paste0('cluster_', i, '_col')]]
                clusterCols <- c(clusterCols, iCol)
            }

            ## update new cluster colours
            c$clusterCol[[clusterMethod]] <- clusterCols

            ## jump to C_tab1
            updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
        }
    })

    ## revert default cluster colors
    observeEvent(input$C_revertClusterColor, {
        if(!is.null(v$data) && !is.null(input$c_colourCluster)){
            clusterMethod <- input$c_colourCluster
            c$clusterCol[[clusterMethod]] <- NULL

            ## jump to C_tab1
            updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
        }
    })


    ## ------annotate clusters-----
    output$C_labelCluster <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            clusterMethods <- c(names(v$data$clusterRes))
            #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
            selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:',
                        choices = clusterMethods,
                        selected = clusterMethods[1], width = "50%")
        }
    })

    output$C_labelCluster_name <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_labelCluster)){
            return(NULL)
        }else{
            textInput("c_labelCluster_name", label = "Type In Your Name for Annotated Cluster",
                      value = paste0("Annotated_", input$c_labelCluster), width = "50%")
        }
    })


    ## currently use 100 as a limit for cluster numbers
    ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ##
    lapply(1:100, function(i) {
        output[[paste0('Cluster', i)]] <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_labelCluster)){
                return(NULL)
            }

            # create new item in RData object
            clusters <- sort(unique(v$data$clusterRes[[input$c_labelCluster]]))
            if (i <= length(clusters)){
                x <- clusters[i]
                textInput(paste0('cluster', i), paste0('Cluster ', x," :"),
                          value = "", width = "30%", placeholder = "Type in the cell type")
            }
        })
    })

    ## update cluster labels
    observeEvent(input$updatelabel, {
        if(!is.null(v$data) && !is.null(input$c_labelCluster) && !is.null(input$c_labelCluster_name)){
            obj <- v$data
            clusterMethod <- input$c_labelCluster
            clusterVec<- obj$clusterRes[[clusterMethod]]
            clusterLabels <- clusterVec
            clusters <- sort(unique(clusterVec))

            for (i in 1:length(clusters)){
                clusteri <- clusters[i]
                ilabel <- input[[paste0('cluster', i)]]
                # if(ilabel == ""){
                #   clusterLabels[clusterLabels==clusteri] <- "Unknown"
                # }else{
                #   clusterLabels[clusterLabels==clusteri] <- ilabel
                # }
                if(ilabel != ""){
                    clusterLabels[clusterLabels==clusteri] <- ilabel
                }
            }

            ## update new cluster results
            labelName <- input$c_labelCluster_name
            obj$clusterRes[[labelName]] <- clusterLabels

            ## update the project name
            obj$projectName <- paste0(obj$projectName, "_annotated_", clusterMethod)

            v$data <- obj

            ## jump to C_tab1
            updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
        }
    })



    ##-----RUN flowSOM-----
    ## result object which will be updated by C_runFlowSOM
    observeEvent(input$C_runFlowSOM, {
        if(!is.null(v$data) && !is.null(input$c_markerSelect)){
            obj <- v$data
            withProgress(message=paste0('Running FlowSOM using k=', input$C_FlowSOM_k), value=0, {
                FlowSOM_cluster <- cytof_cluster(xdata = obj$expressionData[ ,input$c_markerSelect],
                                                 method = "FlowSOM",
                                                 FlowSOM_k = input$C_FlowSOM_k)
                incProgress(1/2)
                ## update FlowSOM cluster results
                obj$clusterRes[["FlowSOM"]] <- FlowSOM_cluster
                ## update the project name
                obj$projectName <- paste0(obj$projectName, "_add_FlowSOM")
                v$data <- obj
                incProgress(1/2)
            })

            ## jump to C_tab1
            updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
        }
    })


    ##------------------------------Marker Panel-------------------------------

    ##-----heat map plot-----
    output$M_plotCluster <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('m_plotCluster', 'Cluster Method:', choices = clusterMethods(),
                        selected = clusterMethods()[1], width = "100%")
        }
    })

    output$M_heatmapmarkerSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            sorted_markerNames <- colnames(v$data$expressionData)
            markerNames <- sorted_markerNames[order(sorted_markerNames)]
            initNum <- ifelse(length(markerNames) >=4, 4, 1)
            selectizeInput('m_heatmapmarkerSelect', 'Select Markers:',
                           choices = markerNames, selected = markerNames[1:initNum],
                           multiple = TRUE, width = "100%")
        }
    })

    observeEvent(input$M_heatmapSelectAll, {
        raw_markers <- colnames(v$data$expressionData)
        markers <- raw_markers[order(raw_markers)]
        updateSelectizeInput(session, "m_heatmapmarkerSelect", selected = markers)
    })

    M_heatmapPlotInput <- reactive({
        if(is.null(v$data) || is.null(input$m_plotCluster) || is.null(input$m_heatmapmarkerSelect))
            return(NULL)
        heatMap(data = v$data,
                clusterMethod = input$m_plotCluster,
                type = input$M_plotMethod,
                dendrogram = input$M_heatmap_dendrogram,
                colPalette = input$M_heatmap_colorPalette,
                selectSamples = input$samples,
                selectMarkers = input$m_heatmapmarkerSelect,
                cex_row_label= input$M_rowLabelSize,
                cex_col_label= input$M_colLabelSize,
                scaleMethod = input$M_scaleMethod)
    })

    output$M_heatmapPlot <- renderPlot({
        M_heatmapPlotInput()
    }, height = 900, width = 950)

    output$PDFHeatmap = downloadHandler(
        filename = function() {
            filename1 <- paste0(input$project_name, "_shinyAPP_Marker_Heatmap_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name,
                                    "_shinyAPP_Marker_Heatmap_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Marker Heatmap PDF files...", value=0, {
                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    heatMap(data = v$data,
                            clusterMethod = input$m_plotCluster,
                            type = input$M_plotMethod,
                            dendrogram = input$M_heatmap_dendrogram,
                            colPalette = input$M_heatmap_colorPalette,
                            selectSamples = input$samples,
                            selectMarkers = input$m_heatmapmarkerSelect,
                            cex_row_label= input$M_rowLabelSize,
                            cex_col_label= input$M_colLabelSize,
                            scaleMethod = input$M_scaleMethod)
                    dev.off()
                })
            }
        }
    )



    session$onSessionEnded(function(){
        # file.remove("cytofkit_shinyAPP_marker_heatmap.pdf")
    })

    ##-----level plot-----
    output$M_PlotMethod <- renderUI({
        if(is.null(v$data) || is.null(visualizationMethods())){
            return(NULL)
        }else{
            selectInput('m_PlotMethod', 'Visualization Method:', choices = visualizationMethods(),
                        selected = visualizationMethods()[1], width = "100%")
        }
    })

    output$M_PlotMarker <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            sorted_markers <- colnames(v$data$expressionData)
            sorted_markers <- sorted_markers[order(sorted_markers)]
            #markers <- c(sorted_markers, "All Markers", "All Markers(scaled)")
            selectizeInput('m_PlotMarker', 'Plot Marker:', choices = sorted_markers,
                           selected = sorted_markers[1], multiple = TRUE, width = "100%")
        }
    })

    observeEvent(input$M_chooseAllMarker, {
        raw_markers <- colnames(v$data$expressionData)
        markers <- raw_markers[order(raw_markers)]
        updateSelectizeInput(session, "m_PlotMarker", selected = markers)
    })

    M_markerExpressionPlotInput <- function(){
        if(is.null(v$data) || is.null(input$m_PlotMethod) || is.null(isolate(input$m_PlotMarker))){
            return(NULL)
        }else{
            withProgress(message="Generating Marker Expression Plot", value=0, {
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$m_PlotMethod,
                                  plotFunction = isolate(input$m_PlotMarker),
                                  pointSize = input$M_PointSize,
                                  alpha = input$M_Alpha,
                                  addLabel = FALSE,
                                  labelSize = input$S_LabelSize,
                                  sampleLabel = FALSE,
                                  FlowSOM_k = input$C_FlowSOM_k,
                                  selectSamples = input$samples,
                                  facetPlot = FALSE,
                                  colorPalette = input$M_colorPalette,
                                  labelRepel = FALSE,
                                  removeOutlier = TRUE,
                                  globalScale = ifelse(input$M_ScaleOptions == "Global", TRUE, FALSE),
                                  centerScale = ifelse(input$M_scaledData == "Centered", TRUE, FALSE))
                incProgress(1/2)
                plot(gp)
                incProgress(1/2)
            })
        }
    }

    output$M_markerExpressionPlot <- renderPlot({
        M_markerExpressionPlotInput()
    }, height = 900, width = 950)

    observeEvent({
        input$M_updateExPlot
        input$m_PlotMethod
        input$M_PointSize
        input$S_LabelSize
        input$M_colorPalette
        input$M_ScaleOptions
        input$M_scaledData
    }, {
        output$M_markerExpressionPlot <- renderPlot({
            M_markerExpressionPlotInput()
        }, height = 900, width = 950)
    })


    output$PDFExpPlot = downloadHandler(
        filename = function() {
            paste0(input$project_name, '_expression_plots.zip')
            filename1 <- paste0(input$project_name, "_shinyAPP_Marker_Expression_Plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name,
                                    "_shinyAPP_Marker_Expression_Plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Marker Expression Plot PDF files...", value=0, {
                    print(getwd())

                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    M_markerExpressionPlotInput()
                    dev.off()
                })
            }
        }
    )

    ##-----histogram plot-----
    output$M_stackFactor <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            stackFactorChoice <- c(names(v$data$clusterRes), "sample")
            selectInput('m_stackFactor', 'Stack Factor:', choices = stackFactorChoice,
                        selected = stackFactorChoice[1], width = "100%")
        }
    })

    output$M_markerSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            sorted_markerNames <- colnames(v$data$expressionData)
            markerNames <- sorted_markerNames[order(sorted_markerNames)]
            initNum <- ifelse(length(markerNames) >=4, 4, 1)
            selectizeInput('m_markerSelect', 'Select Markers:',
                           choices = markerNames, selected = markerNames[1:initNum],
                           multiple = TRUE, width = "100%")
        }
    })

    observeEvent(input$M_histSelectAll, {
        raw_markers <- colnames(v$data$expressionData)
        markers <- raw_markers[order(raw_markers)]
        updateSelectizeInput(session, "m_markerSelect", selected = markers)
    })

    M_stackDensityPlotInput <- function(){
        m_markerSelect <- isolate(input$m_markerSelect)
        if(is.null(v$data) || is.null(input$m_stackFactor) || is.null(m_markerSelect)){
            return(NULL)
        }else{
            withProgress(message="Generating Stack Density Plot", value=0, {
                data <- data.frame(v$data$expressionData, check.names = FALSE)
                samples <- as.character(v$sampleInfo$cellSample)
                mySamples <- samples %in% input$samples
                sfactors <- data.frame(do.call(cbind, v$data$clusterRes),
                                       sample = samples,
                                       stringsAsFactors = FALSE,
                                       check.names = FALSE)
                data <- data[mySamples, ,drop=FALSE]
                stackFactor <- sfactors[mySamples, input$m_stackFactor]

                if(input$m_stackFactor == "sample"){
                    stackFactorColours <- NULL
                }else{
                    clusterMethod <- input$m_stackFactor
                    clusterVec <- v$data$clusterRes[[clusterMethod]]
                    cluster_num <- length(unique(clusterVec))
                    selectColors <- match(levels(as.factor(stackFactor)), levels(as.factor(clusterVec)))
                    if(!is.null(c$clusterCol[[clusterMethod]])){
                        stackFactorColours <- c$clusterCol[[clusterMethod]][selectColors]
                    }else{
                        stackFactorColours <- rainbow(cluster_num)[selectColors]
                    }
                }

                incProgress(1/3)
                gp <- stackDenistyPlot(data = data,
                                       densityCols=m_markerSelect,
                                       stackFactor = stackFactor,
                                       kernel = "gaussian",
                                       bw = "nrd0",
                                       adjust = 1,
                                       stackRotation = 0,
                                       stackSeperation = "auto",
                                       x_text_size = input$M_xlab_size,
                                       strip_text_size = input$M_markerTextSize,
                                       legend_text_size = input$M_legendTextSize,
                                       legendRow = input$M_legendRow,
                                       legend_title = input$m_stackFactor,
                                       stackFactorColours = stackFactorColours)
                incProgress(1/3)
                plot(gp)
                incProgress(1/3)
            })
        }
    }

    observeEvent(input$M_updateDensityPlot, {
        output$M_stackDensityPlot <- renderPlot({
            M_stackDensityPlotInput()
        }, height = 900, width = 950)
    })

    output$PDFHistogram = downloadHandler(
        filename = function() {
            # browser()
            filename1 <- paste0(input$project_name, "_shinyAPP_Stack_Density_Plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name,
                                    "_shinyAPP_Stack_Density_Plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Stack Density Plot PDF files...", value=0, {
                    print(getwd())
                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    M_stackDensityPlotInput()
                    dev.off()
                })
            }
        }
    )

    observeEvent(input$PDFHistogram, {

    })


    ##----- update marker names -----

    ## currently use 100 as a limit for marker number
    ## --- TODO: use reactiveValues to automatically retrive marker numbers --- ##
    lapply(1:100, function(i) {
        output[[paste0('Marker_', i, "_name")]] <- renderUI({
            if(is.null(v$data)){
                return(NULL)
            }
            sorted_markerNames <- colnames(v$data$expressionData)
            markerNames <- sorted_markerNames[order(sorted_markerNames)]

            if (i <= length(markerNames)){
                markeri <- markerNames[i]
                textInput(inputId = paste0('marker_', i, "_name"),
                          label = markeri, value = markeri, width = "30%",
                          placeholder = "Type in your new name for this marker")
            }
        })
    })


    ## update cluster labels
    observeEvent(input$C_updateMarkerNames, {
        if(!is.null(v$data)){
            sorted_markerNames <- colnames(v$data$expressionData)
            markerNames <- sorted_markerNames[order(sorted_markerNames)]
            newMarkerNames <- NULL
            for (i in 1:length(markerNames)){
                iName <- input[[paste0('marker_', i, '_name')]]
                newMarkerNames <- c(newMarkerNames, iName)

            }
            ## update new cluster colours
            mark_pos = which(colnames(v$data$expressionData) %in% markerNames)
            v$data$expressionData[, mark_pos] = v$data$expressionData[, markerNames]
            colnames(v$data$expressionData)[mark_pos] <- newMarkerNames
            ## jump to C_tab1
            updateTabsetPanel(session, "M_markerTabs", selected = "M_tab1")
        }
    })


    ##------------------------------Sample Panel-------------------------------

    ##-----cell percentage heatmap-----
    output$S_plotCluster <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('s_plotCluster', 'Cluster Method:', choices = clusterMethods(),
                        selected = clusterMethods()[1], width = "100%")
        }
    })

    S_heatmapPlotInput <- reactive({
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_plotCluster))
            return(NULL)

        heatMap(data = v$data,
                clusterMethod = input$s_plotCluster,
                type = input$S_plotMethod,
                dendrogram = input$S_heatmap_dendrogram,
                colPalette = input$S_heatmap_colorPalette,
                selectSamples = input$samples,
                cex_row_label= input$S_rowLabelSize,
                cex_col_label= input$S_colLabelSize,
                scaleMethod = input$S_scaleMethod)

    })

    output$S_heatmapPlot <- renderPlot({
        S_heatmapPlotInput()
    }, height = 900, width = 950)

    output$PDFSamHeat = downloadHandler(
        filename = function() {
            filename1 <- paste0(input$project_name, "_shinyAPP_Sample_Heatmap_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name,
                                    "_shinyAPP_Sample_Heatmap_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Sample Heatmap PDF files...", value=0, {
                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    heatMap(data = v$data,
                            clusterMethod = input$s_plotCluster,
                            type = input$S_plotMethod,
                            dendrogram = input$S_heatmap_dendrogram,
                            colPalette = input$S_heatmap_colorPalette,
                            selectSamples = input$samples,
                            cex_row_label= input$S_rowLabelSize,
                            cex_col_label= input$S_colLabelSize,
                            scaleMethod = input$S_scaleMethod)
                    dev.off()
                })
            }
        }
    )

    observeEvent(input$PDFSamHeat, {

    })


    ##-----cell percentage line chart-----
    output$S_clusterMethod2 <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('s_clusterMethod2', 'Cluster Method:', choices = clusterMethods(),
                        selected = clusterMethods()[1], width = "100%")
        }
    })

    output$S_clusterFilter <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2)){
            return(NULL)
        }else{
            clusterIDs <- sort(unique(v$data$clusterRes[[input$s_clusterMethod2]]))
            selectizeInput('s_clusterFilter', 'Filter Clusters:',
                           choices = clusterIDs, selected = clusterIDs,
                           multiple = TRUE, width = "100%")
        }
    })

    S_rateChangePlotInput <- function(){
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2) || is.null(input$s_clusterFilter))
            return(NULL)
        withProgress(message="Generating Rate Change Plot", value=0, {
            ## percentage stat
            data <- data.frame(sample = v$sampleInfo$cellSample,
                               cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod2]]),
                               counts = 1)
            statData1 <- aggregate(counts ~ ., data = data, sum)
            statData2 <- aggregate(counts ~ sample, data = data, sum)
            statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
            statData$percentageInSample <- statData$countsInAll/statData$countsInSample
            incProgress(1/3)
            ## filter clusters
            usedClusters <- input$s_clusterFilter
            clusterCheck <- as.character(statData$cluster) %in% usedClusters
            statData <- statData[clusterCheck, ,drop=FALSE]
            incProgress(1/3)
            gp <- ggplot(data = statData, aes_string(x="sample",
                                                     y="percentageInSample",
                                                     color = "cluster",
                                                     group = "cluster")) +
                geom_point(size = 2) + geom_line(size = 1.5) +
                xlab("Cell Group") + ylab("Percentage of Cells in Group") + theme_bw() +
                theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
            incProgress(1/3)
            plot(gp)
        })
    }

    output$S_rateChangePlot <- renderPlot({
        S_rateChangePlotInput()
    }, height = 500, width = 950)

    output$PDFrateChange = downloadHandler(
        filename = function() {
            filename1 <- paste0(input$project_name, "_shinyAPP_Rate_Change_Plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name,
                                    "_shinyAPP_Rate_Change_Plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Rate Change Plot PDF files...", value=0, {
                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    S_rateChangePlotInput()
                    dev.off()
                })
            }
        }
    )

    observeEvent(input$PDFrateChange, {

    })


    # output$S_clusterTable <- renderTable({
    #     if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2)){
    #         return(NULL)
    #     }else{
    #         data <- data.frame(sample = v$sampleInfo$cellSample,
    #                            cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod2]]),
    #                            counts = 1)
    #
    #         statData1 <- aggregate(counts ~ ., data = data, sum)
    #         statData2 <- aggregate(counts ~ sample, data = data, sum)
    #         statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
    #         if(is.numeric(statData$cluster)) statData$cluster <- as.integer(statData$cluster)
    #         statData$counts <- as.integer(statData$countsInAll)
    #         statData$percentageInAll <- round(statData$countsInAll/nrow(data), 4)
    #         statData$percentageInSample <- round(statData$countsInAll/statData$countsInSample, 2)
    #         statData[, c("sample", "cluster", "counts", "percentageInSample", "percentageInAll")]
    #     }
    # })


    ##-----Regroup samples-----
    output$S_groupSamples <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            clusterMethods <- c(names(v$data$clusterRes))
            #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
            selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:',
                        choices = clusterMethods,
                        selected = clusterMethods[1], width = "30%")
        }
    })

    ## currently use 100 as a limit for sample numbers
    ## --- TODO: use reactiveValues to automatically retrive sample numbers --- ##
    lapply(1:100, function(i) {
        output[[paste0('S_sample', i)]] <- renderUI({
            if(is.null(v$data) || is.null(v$sampleInfo)){
                return(NULL)
            }

            uniqueSampleNames <- sort(unique(v$sampleInfo$cellSample))
            if (i <= length(uniqueSampleNames)){
                x <- uniqueSampleNames[i]
                textInput(paste0('Sample', i), paste0(x," :"),
                          value = "", width = "40%",
                          placeholder = "Type in the group name for this sample")
            }
        })
    })


    rename_sample = function(new_sample_name){
        # browser()
        v$sampleInfo$originalCellSample <- v$sampleInfo$cellSample
        uniqueSampleNames <- sort(unique(v$sampleInfo$originalCellSample))

        temp_names = sapply(1:length(uniqueSampleNames), function(i){
            sample_name_length = length(v$data$sampleNames[[i]])
            v$data$sampleNames[[i]][sample_name_length]
        })
        v$data$sampleNames = v$data$sampleNames[order(temp_names)]

        sampleGroupNames <- NULL
        for(i in 1:length(uniqueSampleNames)){
            sample_name_length = length(v$data$sampleNames[[i]])
            if (new_sample_name[i] != "") {
                sampleGroupNames <- c(sampleGroupNames, new_sample_name[i])
                v$data$sampleNames[[i]] <- c(v$data$sampleNames[[i]][1], new_sample_name[i])
            } else {
                sampleGroupNames <- c(sampleGroupNames, v$data$sampleNames[[i]][sample_name_length])
                # v$data$sampleNames[[i]] <- c(v$data$sampleNames[[i]][sample_name_length], v$data$sampleNames[[i]][sample_name_length])
            }
        }

        groupNameLevels <- strsplit(input$sampleGroupLevels, ";", fixed = TRUE)[[1]]

        if(groupNameLevels != "" && all(sampleGroupNames != "")
           && length(groupNameLevels) == length(unique(sampleGroupNames))
           && all(as.character(groupNameLevels) %in% sampleGroupNames)){
            sampleMatchID <- match(v$sampleInfo$originalCellSample, uniqueSampleNames)
            v$sampleInfo$cellSample <- factor(sampleGroupNames[sampleMatchID],
                                              levels = groupNameLevels)
        }else{
            sampleGroupNames[sampleGroupNames == ""] <- uniqueSampleNames[sampleGroupNames == ""]
            sampleMatchID <- match(v$sampleInfo$originalCellSample, uniqueSampleNames)
            v$sampleInfo$cellSample <- factor(sampleGroupNames[sampleMatchID])
        }

        cellID_number <- do.call(base::c, regmatches(v$sampleInfo$cellID,
                                                     gregexpr("_[0-9]*$", v$sampleInfo$cellID, perl=TRUE)))

        ## update reactive object v$sampleInfo
        ## newCellID = "sampleGroup" + "_cellID" + "globalID" to avoid duplicates
        v$sampleInfo$newCellID <- paste0(as.character(v$sampleInfo$cellSample),
                                         "_",
                                         1:length(cellID_number))


        ## update reactive object v$data
        expressionData <- v$data$expressionData
        row.names(expressionData) <- v$sampleInfo$newCellID
        v$data$expressionData <- expressionData

        ## update the project name
        v$data$projectName <- paste0(v$data$projectName, "_grouped_samples")

        ## update v$data$progressionRes
        if(!is.null(v$data$progressionRes)){
            sampleExpressData <- v$data$progressionRes$sampleData
            row.names(sampleExpressData) <- v$sampleInfo$newCellID[match(row.names(sampleExpressData),
                                                                         v$sampleInfo$cellID)]
            v$data$progressionRes$sampleData <- sampleExpressData
        }

        ## jump to S_tab1
        # updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
    }


    reset_sample = function(){
        # browser()
        v$sampleInfo <- v$ori_sampleInfo
        v$data = v$ori_data
        v$sample_selected_index = 1:length(unique(as.character(v$sampleInfo$cellSample)))

    }

    ## update sample groups
    observeEvent(input$updateSampleGroups, {
        if(!is.null(v$data) && !is.null(v$sampleInfo)){
            # browser()
            v$sampleInfo$originalCellSample <- v$sampleInfo$cellSample
            uniqueSampleNames <- sort(unique(v$sampleInfo$originalCellSample))
            new_sample_names <- NULL
            for(i in 1:length(uniqueSampleNames)){
                new_sample_names <- c(new_sample_names, input[[paste0("Sample", i)]])
            }
            rename_sample(new_sample_names)

            ## jump to S_tab1
            updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
        }
    })

    ## revert old sample names
    observeEvent(input$revertSampleNames, {
        if(!is.null(v$data) && !is.null(v$sampleInfo)){
            if(!is.null(v$sampleInfo$originalCellSample)){
                v$sampleInfo$cellSample <- v$sampleInfo$originalCellSample
                v$sampleInfo$originalCellSample <- NULL

                ## update reactive object v$data
                expressionData <- v$data$expressionData
                row.names(expressionData) <- v$sampleInfo$cellID
                v$data$expressionData <- expressionData

                ## update the project name
                v$data$projectName <- sub("_grouped_samples", "", v$data$projectName)

                ## update reactive object v$sampleInfo
                if(!is.null(v$data$progressionRes)){
                    sampleExpressData <- v$data$progressionRes$sampleData
                    row.names(sampleExpressData) <- v$sampleInfo$cellID[match(row.names(sampleExpressData),
                                                                              v$sampleInfo$newCellID)]
                    v$data$progressionRes$sampleData <- sampleExpressData
                }
            }
            ## jump to S_tab1
            updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
        }
    })

    ##---------------------------Progression Panel------------------------------

    ##-----subset relationship plot-----

    output$P_xlab <- renderUI({
        if(is.null(v$data) || is.null(progressionLabs())){
            return(NULL)
        }else{
            selectInput('p_xlab', 'Plot X:', choices = progressionLabs(),
                        selected = progressionLabs()[1], width = "100%")
        }
    })

    output$P_ylab <- renderUI({
        if(is.null(v$data) || is.null(progressionLabs())){
            return(NULL)
        }else{
            selectInput('p_ylab', 'Plot Y:', choices = progressionLabs(),
                        selected = progressionLabs()[2], width = "100%")
        }
    })

    P_scatterPlotInput <- function(){
        if(is.null(v$data) || is.null(v$data$progressionRes) || is.null(input$p_xlab) || is.null(input$p_ylab)){
            return(NULL)
        }else{
            withProgress(message="Generating Progression Scatter Plot", value=0, {
                obj <- v$data$progressionRes
                data <- data.frame(obj$progressionData,
                                   cluster = obj$sampleCluster,
                                   sample = sub("_[0-9]*$", "", row.names(obj$sampleData)))
                incProgress(1/3)
                data <- data[data$sample %in% input$samples, ,drop=FALSE]

                clusterMethod <- p$progressionCluster
                clusterVec <- v$data$clusterRes[[clusterMethod]]
                cluster_num <- length(unique(clusterVec))
                selectColors <- match(levels(as.factor(data$cluster)), levels(as.factor(clusterVec)))

                if(!is.null(c$clusterCol[[clusterMethod]])){
                    clusterColor <- c$clusterCol[[clusterMethod]][selectColors]
                }else{
                    clusterColor <- rainbow(cluster_num)[selectColors]
                }

                gp <- cytof_clusterPlot(data = data,
                                        xlab = input$p_xlab,
                                        ylab = input$p_ylab,
                                        cluster = "cluster",
                                        sample = "sample",
                                        title = "Subset Relationship",
                                        type = ifelse(input$P_facetPlot, 2, 1),
                                        point_size = input$P_PointSize,
                                        addLabel = input$P_addLabel,
                                        labelSize = input$P_LabelSize,
                                        sampleLabel = FALSE,
                                        labelRepel = input$P_labelRepel,
                                        fixCoord = FALSE,
                                        clusterColor = clusterColor)
                incProgress(1/3)
                plot(gp)
                incProgress(1/3)
            })
        }
    }

    output$P_scatterPlot <- renderPlot({
        P_scatterPlotInput()
    }, height = 900, width = 950)

    output$PDFScatter = downloadHandler(
        filename = function() {
            filename1 <- paste0(input$project_name, "_shinyAPP_Scatterplot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(input$project_name,
                                    "_shinyAPP_Scatterplot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Progression Scatterplot PDF files...", value=0, {
                    print(getwd())
                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    P_scatterPlotInput()
                    dev.off()
                })
            }
        }
    )



    ##-----marker expression profile-----

    output$P_orderBy <- renderUI({
        if(is.null(v$data) || is.null(progressionLabs())){
            return(NULL)
        }else{
            selectInput('p_orderBy', 'Cell Order By:', choices = progressionLabs(),
                        selected = progressionLabs()[1], width = "100%")
        }
    })

    output$P_markerSelect <- renderUI({
        if(is.null(v$data) || is.null(v$data$progressionRes)){
            return(NULL)
        }else{
            sorted_markerNames <- colnames(v$data$progressionRes$sampleData)
            markerNames <- sorted_markerNames[order(sorted_markerNames)]
            initNum <- ifelse(length(markerNames) >=4, 4, 1)
            selectizeInput('p_markerSelect', 'Select Markers:',
                           choices = markerNames, selected = markerNames[1:initNum],
                           multiple = TRUE, width = "100%")
            # checkboxGroupInput('p_markerSelect', strong('Select Markers:'),
            #                    markerNames, selected = markerNames, inline = TRUE)
        }
    })

    output$P_clusterSelect <- renderUI({
        if(is.null(v$data) || is.null(v$data$progressionRes)){
            return(NULL)
        }else{
            clusterIDs <- sort(unique(v$data$progressionRes$sampleCluster))
            selectizeInput('p_clusterSelect', 'Select Clusters:',
                           choices = clusterIDs, selected = clusterIDs,
                           multiple = TRUE, width = "100%")
            # checkboxGroupInput('p_clusterSelect', strong('Select Clusters:'),
            #                    clusterIDs, selected = clusterIDs, inline = TRUE)
        }
    })

    P_markerPlotInput <- function(){
        p_markerSelect <- isolate(input$p_markerSelect)
        p_clusterSelect <- isolate(input$p_clusterSelect)
        if(is.null(v$data) || is.null(v$data$progressionRes) || is.null(p_markerSelect) || is.null(p_clusterSelect) || is.null(input$p_orderBy))
            return(NULL)

        withProgress(message="Generating Marker Expression Profile", value=0, {
            data <- data.frame(v$data$progressionRes$sampleData,
                               cluster = v$data$progressionRes$sampleCluster,
                               v$data$progressionRes$progressionData,
                               check.names = FALSE)

            sampleNames <- sub("_[0-9]*$", "", row.names(v$data$progressionRes$sampleData))
            data <- data[sampleNames %in% input$samples, ,drop=FALSE]
            incProgress(1/3)
            if(input$P_combineTrends){
                pp <- cytof_expressionTrends(data,
                                             markers = p_markerSelect,
                                             clusters = p_clusterSelect,
                                             orderCol = input$p_orderBy,
                                             clusterCol = "cluster",
                                             reverseOrder = input$P_reverseOrder,
                                             addClusterLabel = input$P_addLabel2,
                                             clusterLabelSize = input$P_LabelSize2,
                                             segmentSize = 0.5,
                                             min_expr = NULL)
            }else{
                pp <- cytof_progressionPlot(data,
                                            markers = p_markerSelect,
                                            clusters = p_clusterSelect,
                                            orderCol = input$p_orderBy,
                                            clusterCol = "cluster",
                                            reverseOrder = input$P_reverseOrder,
                                            addClusterLabel = input$P_addLabel2,
                                            clusterLabelSize = input$P_LabelSize2,
                                            segmentSize = 0.5,
                                            min_expr = NULL)
            }
            incProgress(1/3)
            plot(pp)
            incProgress(1/3)
        })
    }

    observeEvent(input$P_updateRegressionPlot, {
        output$P_markerPlot <- renderPlot({
            P_markerPlotInput()
        }, height = 900, width = 950)
    })

    output$PDFmarkerPlot = downloadHandler(
        filename = function() {
            filename1 <- paste0(input$project_name, "_shinyAPP_Marker_Plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename1)){
                filename1 <- paste0(project_name,
                                    "_shinyAPP_Marker_Plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
            }
            filename1
        },
        content = function(file) {
            if(!is.null(v$data)){
                withProgress(message="Downloading Marker Plot PDF files...", value=0, {
                    print(getwd())
                    pdf(file,
                        width=as.integer(input$tab_w),
                        height=as.integer(input$tab_h))
                    P_markerPlotInput()
                    dev.off()
                })
            }
        }
    )



    ##-----Run Diffusionmap-----

    output$P_clusterMethod <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('p_clusterMethod', 'Cluster Method:', choices = clusterMethods(),
                        selected = clusterMethods()[1], width = "100%")
        }
    })

    output$P_clusterTable <- renderTable({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            clusterTable <- t(as.matrix(table(v$data$clusterRes[[input$p_clusterMethod]])))
            out <- as.data.frame(clusterTable, row.names = "Cell Counts")
            colnames(out) <- paste("Cluster", colnames(out))
            out
        }
    })

    output$P_clusterFilter <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            obj <- v$data
            clusterIDs <- sort(unique(obj$clusterRes[[input$p_clusterMethod]]))
            selectizeInput('p_clusterFilter', 'Filter Clusters:',
                           choices = clusterIDs, selected = clusterIDs,
                           multiple = TRUE, width = "100%")
        }
    })


    ## result object which will be updated by P_runDiffusionmap
    observeEvent(input$P_runDiffusionmap, {

        if(!is.null(v$data)){
            obj <- v$data
            usedClusters <- input$p_clusterFilter
            clusterCheck <- obj$clusterRes[[input$p_clusterMethod]] %in% usedClusters
            mdata <- obj$expressionData[clusterCheck, ]
            mcluster <- obj$clusterRes[[input$p_clusterMethod]][clusterCheck]
            withProgress(message="Running Diffusionmap", value=0, {
                diffmapRes <- cytof_progression(data = mdata,
                                                cluster = mcluster,
                                                method = "diffusionmap",
                                                distMethod = input$P_distMethod,
                                                out_dim = input$P_outDim,
                                                clusterSampleMethod = input$P_sampleMethod,
                                                clusterSampleSize = input$P_clusterSampleSize)
                incProgress(1/2)
                ## update progressionRes results
                obj$progressionRes <- diffmapRes

                ## update the project name
                obj$projectName <- paste0(obj$projectName, "_added_diffusionmap")

                v$data <- obj
                incProgress(1/2)
            })
            p$progressionCluster <- input$p_clusterMethod
            ## jump to P_tab1
            updateTabsetPanel(session, "P_progressionTabs", selected = "P_tab1")
        }
    })

    ##---------------Imaging mass cytometry Analysis--------------------##

    observeEvent(input$loadButton_img, {
        sce <- input$sce
        mask <- input$mask
        #images <- input$images

        if (is.null(sce)){
            v$sceData <- NULL
        }else{
            withProgress(message="Loading and Processing Data...", value=0, {
                print(sce$datapath)
                print(sce$name)
                print(file.exists(paste(sce$datapath[1], "/", sce$name[1], sep="")))
                sce1 <- readRDS(sce$datapath)
                #image1 <- readRDS(images$datapath)
                print(sce1)
                v$sceData <- sce1
            })
        }

        if (is.null(mask)){
            v$maskData <- NULL
        }else{
            withProgress(message="Loading and Processing Data...", value=0, {
                print(mask$datapath)
                print(mask$name)
                print(file.exists(paste(mask$datapath[1], "/", mask$name[1], sep="")))
                mask1 <- readRDS(mask$datapath)
                #image1 <- readRDS(images$datapath)
                print(mask1[1])
                v$maskData <- mask1
                print(v$maskData)
            })
        }

        dir.create("Hyperion_results")
    })

    observeEvent(input$reset_img, {
        session$reload()
        print("Reset done")
    })

    observeEvent(input$saveButton, {
        if(!is.null(input$sce) && !is.null(input$mask) && !is.null(input$images)){
            withProgress(message="Saving Results...", value=0, {
                print(getwd())
                resultDir <- paste0(getwd(), .Platform$file.sep, "Hyperion_results")
                #filename <- paste0(resultDir, .Platform$file.sep, sce1, "_", Sys.Date())
                #save(sce1, file= paste0(resultDir, .Platform$file.sep, sce1, "_", Sys.Date(), ".Robj"))
                #save(mask1, file= paste0(resultDir, .Platform$file.sep, mask1, "_", Sys.Date(), ".Robj"))
                #save(image1, file= paste0(resultDir, .Platform$file.sep, imageData, "_", Sys.Date(), ".Robj"))
            })
            ## open the results directory
            opendir(resultDir)
        }
    })

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

    ##---------------Plot cells tabset-------------------

    output$plot_cell1 <- renderPlot({
        if(is.null(v$maskData[1])){
            plotly_empty()
        }else{
            plotCells(v$maskData[1])
        }
    })

    output$plot_cell2 <- renderPlot({
        if(is.null(v$sceData) && is.null(v$maskData[1])){
            plotly_empty()
        }else{
            plotCells(v$maskData[1], v$sceData, cell_id = "CellNumber", img_id = "ImageName", outline_by = "CellType")
        }
    })

    output$plot_cell3 <- renderPlot({
        if(is.null(v$sceData) && is.null(v$maskData[1])){
            plotly_empty()
        }else{
            plotCells(v$maskData[1], v$sceData, cell_id = "CellNumber", img_id = "ImageName", colour_by = "CellType")
        }
    })

    output$sce.select <- renderUI({
        if(is.null(v$sceData)&& is.null(v$maskData[1])){
            return(NULL)
        }else{
            selectInput("sce.gene", label = "Markers",
                        choices = rownames(v$sceData))
        }
    })

    output$plot_cell4 <- renderPlot({
        if(is.null(v$sceData) && is.null(v$maskData[1])){
            plotly_empty()
        }else{
            plotCells(v$maskData[1], v$sceData, cell_id = "CellNumber", img_id = "ImageName", colour_by = input$sce.gene, outline_by = "CellType")
        }
    })

    output$plot_cell5 <- renderPlot({
        if(is.null(v$sceData) && is.null(v$maskData[1])){
            plotly_empty()
        }else{
            plotCells(v$maskData[1], v$sceData, cell_id = "CellNumber", img_id = "ImageName", colour_by = input$sce.gene, exprs_values = "exprs")
        }
    })

    output$name <- renderPrint({
        s <- event_data("plotly_selected")
        c(s[["key"]], class(s[["key"]]))
    })

    observeEvent(input$PDFa, {
        if(!is.null(v$sceData) && !is.null(v$maskData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Hyperion_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Cell_Plot", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "Cell_Plot",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                nA <- plotCells(v$maskData[1])
                nB <- plotCells(v$maskData[1], v$sceData, cell_id = "CellNumber", img_id = "ImageName", colour_by = "CD99", outline_by = "CellType")
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(nA)
                print(nB)
                dev.off()
            })
        }
    })

    ##---------------Select fragments file path -------------------

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

    ##---------------Computing QC Metrics for scATAC-seq-------------------

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

    #output$VlnPlot_atac <- renderPlotly({
    #    if(is.null(v$atacData)){
    #        plotly_empty()
    #    }else{
    #        VlnPlot(object = v$atacData, features = c('pct_reads_in_peaks', 'peak_region_fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1, ncol = 5)
    #    }
    #})

    ##---------------Normalization and linear dimensional reduction for scATAC-seq-------------------

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
            v$atacData <- FindNeighbors(object = v$atacData, reduction = 'lsi', dims = 2:input$dim.used_atac, nn.method = "rann")
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

    output$name <- renderPrint({
        s <- event_data("plotly_selected")
        c(s[["key"]], class(s[["key"]]))
    })

    observeEvent(input$PDFa, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_violin_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "QC_violin_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                nG <- VlnPlot(v$scData, "nFeature_RNA")
                pM <- VlnPlot(v$scData, "percent.mt")
                nU <- VlnPlot(v$scData, "nCount_RNA")
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(nG)
                print(pM)
                print(nU)
                dev.off()
            })
        }
    })

    ##---------------Spatial Transcriptomics Analysis using Seurat-------------------##

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

    observeEvent(input$loadButton3, {
        if(input$Module == "Spatial Transcriptomics Analysis" & input$scInput3 == "H5"){
            tpmFiles3 <- input$tpmFiles3
            if (is.null(tpmFiles3)){
                v$scData_spatial <- NULL
            }else{
                withProgress(message="Loading and Processing Data...", value=0, {
                    print(tpmFiles3$datapath)
                    print(tpmFiles3$name)
                    print(file.exists(paste(tpmFiles3$datapath[1], "/", tpmFiles3$name[1], sep="")))
                    #path = "/acrc_raman/jinmiao/CJM_lab/Raman/Projects/hyperion_cytofkit2/spatial_shiny/stxBrain/"
                    #s <-list.files(path =  parseDirPath(c(home = '~'), dir(), pattern="*.h5"))
                    exp.data_spatial <- Load10X_Spatial(parseDirPath(c(home = '.'), dir()), filename = tpmFiles3$name, assay = "Spatial")
                    additional.ident <- NULL
                    incProgress(0.5, "Creating Seurat Object")
                    v$scData_spatial <- exp.data_spatial
                })
            }
            dir.create("Seurat_spatial_results")
        }
    }
    )

    observeEvent(input$reset3, {
        session$reload()
        print("Reset done")
    })

    ##---------------Spatial QC tabset-------------------

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

        ##---------------PCA tabset-------------------
        # PCA plot
        observeEvent(input$doPCA_spatial, {
            withProgress(message = "Scaling Data...", value = 0,{
                incProgress(0.5, message = "Running Dimension Reduction...")

                v$scData_spatial <- SCTransform(v$scData_spatial, assay = "Spatial", verbose = FALSE)

                output$deg1.gene.select <- renderUI({
                    if(is.null(v$scData_spatial)){
                        return(NULL)
                    }else{
                        selectInput("deg1.gene", label = "Gene to visualise",
                                    choices = rownames(v$scData_spatial))
                    }
                })

                output$Deg1_spatial.plot <- renderPlot({
                    if(is.null(v$scData_spatial)){
                        return(NULL)
                    }else{
                        withProgress(message="Generating DEG Plot...", value=0, {
                            SpatialFeaturePlot(object = v$scData_spatial, features = input$deg1.gene, alpha = c(0.1, 1))
                        })
                    }
                })

                v$scData_spatial <- RunPCA(v$scData_spatial, assay = "SCT", verbose = FALSE)
                v$scData_spatial <- FindNeighbors(v$scData_spatial, reduction = "pca", dims = 1:30, nn.method = "rann")
                v$scData_spatial <- FindClusters(v$scData_spatial, verbose = FALSE)
                v$scData_spatial <- RunUMAP(v$scData_spatial, reduction = "pca", dims = 1:30)
                print(v$scData_spatial[["pca"]], dims = 1:5, nfeatures = 5)
                v$isUMAPdone <- TRUE
                Dim_plot1 <- DimPlot(v$scData_spatial, reduction = "umap", label = TRUE)
                Dim_plot2 <- SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
                print(Dim_plot1)
                print(Dim_plot2)
                incProgress(0.4, message = "Getting list of PC genes...")
                pc.table1 <- list()
                for(i in 1:20){
                    pcg1 <- TopFeatures(v$scData_spatial)
                    pc.table1[[i]] <- pcg1
                }
                pc.table1 <- as.data.frame(pc.table1, col.names = paste0("PC", 1:20))
                v$pcGenes1 <- pc.table1
            })
        })

        output$clustUI <- renderUI({
            if(is.null(v$isUMAPdone)){
                return(NULL)
            }else{
                tagList(
                    fluidRow(
                        column(6,
                               numericInput("clus.res2",
                                            label = "Cluster Resolution",
                                            value = 0.6,
                                            min = 0.1,
                                            step = 0.1)
                        ),
                        #column(6,
                        #       actionButton("findCluster1", "Find Clusters", icon = icon("hand-pointer-o")),
                        #       textOutput("cluster.done")
                        #)
                    )
                )
            }
        })

        # observeEvent(input$findCluster1, {
        #    withProgress(message = "Finding clusters...", value = 0.3, {
        #      v$scData_spatial <- FindNeighbors(v$scData_spatial, dims = 1:10)
        #     v$scData_spatial <- FindClusters(v$scData_spatial, resolution = input$clus.res2)
        #     output$cluster.done <- renderText(paste0("Clustering done!"))
        #     v$isClusterdone <- TRUE
        #   })
        # })

        output$DimPlot_spatial <- renderPlotly({
            if(is.null(v$isUMAPdone)){
                plotly_empty()
            }else{
                withProgress(message="Generating Dim Plot...", value=0, {
                    DimPlot(v$scData_spatial, label = TRUE)
                })
            }
        })

        output$SpatialDimPlot <- renderPlot({
            if(is.null(v$isUMAPdone)){
                plotly_empty()
            }else{
                withProgress(message="Generating SpatialDim Plot...", value=0, {
                    SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
                })
            }
        })

        observeEvent(input$PDFd, {
            if(!is.null(v$scData_spatial)){
                withProgress(message="Downloading plot PDF files...", value=0, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"UMAP_plot_", Sys.Date(), ".pdf")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,
                                            "UMAP_plot_",
                                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                        i = i + 1;
                    }
                    dimplot1 <- DimPlot(v$scData_spatial, reduction = "umap", label = T)
                    dimplot2 <- SpatialDimPlot(v$scData_spatial, label = TRUE, label.size = 3)
                    prePlot()
                    pdf(filename2,
                        width=as.numeric(input$pdf_w),
                        height=as.numeric(input$pdf_h))
                    print(dimplot1)
                    print(dimplot2)
                    dev.off()
                })
                withProgress(message="Downloading UMAP coordinates...", value=0.5, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"umap_", Sys.Date(), ".txt")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,
                                            "umap_",
                                            Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    write.csv(v$scData_spatial@reductions$umap@cell.embeddings, file = filename2)
                })
                withProgress(message="Downloading cluster IDs...", value=0.9, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), ".txt")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), "_", sprintf("%03d", i + 1), ".csv");
                        i = i + 1;
                    }
                    write.csv(v$scData_spatial@active.ident, file = filename2)
                })
            }
        })

        # Viz plot

        output$vizPlot_spatial <- renderPlot({
            if(is.null(v$scData_spatial)){
                return(NULL)
            }else{
                VizDimLoadings(v$scData_spatial, dims = as.numeric(input$select.pc1))
            }
        })

        output$PCHeatmap_spatial <- renderPlot({
            if(is.null(v$scData_spatial)){
                return(NULL)
            }else{
                DimHeatmap(v$scData_spatial, dims = as.numeric(input$select.pc1))
            }
        })

        output$PCtable_spatial <- DT::renderDataTable({
            if(is.null(v$scData_spatial) ){
                return(NULL)
            }else{
                v$pcGenes1
            }
        }, options = list(scrollX = TRUE))

        observeEvent(input$PDFe, {
            if(!is.null(v$scData)){
                withProgress(message="Downloading plot PDF files...", value=0, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"Viz_Heatmap_plots_", Sys.Date(), ".pdf")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,
                                            "Viz_Heatmap_plots_",
                                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                        i = i + 1;
                    }
                    prePlot()
                    pdf(filename2,
                        width=as.numeric(input$pdf_w),
                        height=as.numeric(input$pdf_h))
                    #isolate({
                    vizdim1 <- VizDimLoadings(v$scData, dims = as.numeric(input$select.pc))
                    dimheat1 <- DimHeatmap(v$scData)
                    print (vizdim1)
                    print (dimheat1)
                    #})
                    dev.off()
                    pcGenes1 <- v$pcGenes1
                    write.csv(v$pcGenes1, file = paste0(pdfDir, .Platform$file.sep,"PC_genes_", Sys.Date(), ".csv"))
                })
            }
        })

        ##---------------DEGs tabset-------------------

        observeEvent(input$doDeg_spatial, {
            if(is.null(v$scData_spatial)){
                return(NULL)
            }else{
                withProgress(message="Finding DEGs...", value=0, {
                    sp.markers <- FindAllMarkers(v$scData_spatial, only.pos = FALSE, min.pct = input$min_pct_spatial, logfc.threshold = input$logfc_spatial)
                    v$sp.markers <- sp.markers
                })
            }
        })

        output$deg2.gene.select <- renderUI({
            if(is.null(v$sp.markers)){
                return(NULL)
            }else{
                selectInput("deg2.gene", label = "Gene to visualise",
                            choices = rownames(v$sp.markers))
            }
        })

        output$Deg2_spatial.plot <- renderPlot({
            if(is.null(v$scData_spatial)){
                return(NULL)
            }else{
                withProgress(message="Generating DEG Plot...", value=0, {
                    SpatialFeaturePlot(object = v$scData_spatial, features = input$deg2.gene, alpha = c(0.1, 1))
                })
            }
        })

        # output$Deg1.plot <- renderPlotly({
        #  if(is.null(v$ips.markers)){
        #    return(NULL)
        #  }else{
        #      withProgress(message="Generating DEG Plot...", value=0, {
        #       FeaturePlot(v$scData, input$deg.gene)
        #      })
        #    }
        #  })

        output$Deg_spatial.table <- DT::renderDataTable(
            v$sp.markers, options = list(scrollX = TRUE, scrollY = "400px"))

        observeEvent(input$doDegn, {
            if(is.null(v$scData_spatial)){
                return(NULL)
            }else{
                withProgress(message="Finding DEGs...", value=0, {
                    v$scData_spatial <- FindSpatiallyVariableFeatures(v$scData_spatial, assay = "SCT", features = VariableFeatures(v$scData_spatial)[1:1000], selection.method = "markvariogram")
                    sp.spatialfeatures <- SpatiallyVariableFeatures(v$scData_spatial, selection.method = "markvariogram")
                    v$sp.spatialfeatures <- sp.spatialfeatures
                })
            }
        })

        output$degn.gene.select <- renderUI({
            if(is.null(v$sp.spatialfeatures)){
                return(NULL)
            }else{
                selectInput("degn.gene", label = "Gene to visualise",
                            choices = v$sp.spatialfeatures)
            }
        })

        output$Degn.plot <- renderPlot({
            if(is.null(v$scData_spatial)){
                return(NULL)
            }else{
                withProgress(message="Generating DEG Plot...", value=0, {
                    SpatialFeaturePlot(object = v$scData_spatial, features = input$degn.gene, alpha = c(0.1, 1))
                })
            }
        })

        observeEvent(input$PDFk, {
            if(!is.null(v$scData_spatial)){
                withProgress(message="Downloading plot PDF files...", value=0, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                    if(!dir.exists(pdfDir)){
                        dir.create(pdfDir)
                    }
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), ".pdf")
                    i = 0
                    while(file.exists(filename2)){
                        filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                        i = i + 1;
                    }
                    degSpatial <- SpatialFeaturePlot(object = v$scData_spatial, features = input$degn.gene, alpha = c(0.1, 1))
                    #degFeature <- FeaturePlot(v$scData, input$deg.gene)
                    prePlot()
                    pdf(filename2,
                        width=as.numeric(input$pdf_w),
                        height=as.numeric(input$pdf_h))
                    print(degSpatial)
                    #print(degFeature)
                    dev.off()
                    write.csv(v$sp.markers, file = paste0(pdfDir, .Platform$file.sep,"DEG_table_", input$c1, "vs", input$c2, "_", Sys.Date(), ".csv"))
                })
            }
        })

        ##---------------Deconvolution-------------------

        observeEvent(input$loadexample_deconv, {
            if(input$Module == "Spatial Transcriptomics Analysis"){

                tpmFiles_scRNA <- readRDS('allen_cortex_dwn.rds')
                v$scRNAData <- tpmFiles_scRNA
                label1 <- "example scRNA-seq Reference dataset loaded"
                updateActionButton(inputId = "loadexample_deconv", label = label1)
                shinyalert("Loaded", "example scRNA-seq Reference dataset loaded.", type = "success", imageWidth = 10, imageHeight = 10)
            }
        })

        observeEvent(input$loadButton4, {
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
                })
            }
        })

        observeEvent(input$process_scRNA, {
            withProgress(message = "Processing scRNA-seq dataset...", value = 0,{
                incProgress(0.5, message = "Processing...")
                v$scRNAData <- Seurat::SCTransform(v$scRNAData, verbose = FALSE)
                v$scRNAData <- Seurat::RunPCA(v$scRNAData, verbose = FALSE)
                v$scRNAData <- Seurat::RunUMAP(v$scRNAData, dims = 1:30, verbose = FALSE)
                v$scRNAData <- Seurat::FindNeighbors(v$scRNAData, dims = 1:30, verbose = FALSE, nn.method = "rann")
                v$scRNAData <- Seurat::FindClusters(v$scRNAData, verbose = FALSE)
                scRNA_plot <- Seurat::DimPlot(v$scRNAData, group.by = "subclass")
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
                    DimPlot(v$scRNAData, group.by = "subclass", label = T)
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

        observeEvent(input$get_markers, {
            withProgress(message = "Getting marker genes...", value = 0,{
                incProgress(0.5, message = "Processing...")
                Seurat::Idents(object = v$scRNAData) <- v$scRNAData@meta.data$subclass
                ips.markers1 <- Seurat::FindAllMarkers(object = v$scRNAData,
                                                       assay = "SCT",
                                                       slot = "data",
                                                       verbose = TRUE,
                                                       only.pos = TRUE,
                                                       logfc.threshold = input$logfc,
                                                       min.pct = input$min_pct)
                v$ips.markers1 <- ips.markers1
                #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
            })
        })

        observeEvent(input$doDeconv, {
            withProgress(message = "Performing deconvolution...", value = 0,{
                incProgress(0.5, message = "Deconvoluting...")
                spotlight_ls <- spotlight_deconvolution(v$scRNAData,
                                                        counts_spatial = v$scData_spatial@assays$Spatial@counts,
                                                        clust_vr = "subclass",
                                                        cluster_markers = v$ips.markers1,
                                                        cl_n = 50,
                                                        hvg = 3000,
                                                        ntop = NULL,
                                                        transf = "uv",
                                                        method = "nsNMF",
                                                        min_cont = 0.09)
                v$spotlight_ls <- spotlight_ls
                v$isDeconvdone <- TRUE
                decon_mtrx <- v$spotlight_ls[[2]]
                v$decon_mtrx <- decon_mtrx
                cell_types_all <- colnames(v$decon_mtrx )[which(colnames(v$decon_mtrx ) != "res_ss")]
                v$cell_types_all <- cell_types_all
                v$scData_spatial@meta.data <- cbind(v$scData_spatial@meta.data, v$decon_mtrx)
                img_path = "tissue_lowres_image.png"
                v$img_path <- img_path
                deconv_plot <- SPOTlight::spatial_scatterpie(se_obj = v$scData_spatial,
                                                             cell_types_all = v$cell_types_all,
                                                             img_path = v$img_path,
                                                             pie_scale = 0.4)
                print(deconv_plot)
                #output$vis_sp.done <- renderText(paste0("Processing of scRNA-seq data done!"))
            })
        })

        output$DeconvPlot <- renderPlot({
            if(is.null(v$isDeconvdone)){
                plotly_empty()
            }else{
                withProgress(message="Generating SpatialDim Plot...", value=0, {
                    SPOTlight::spatial_scatterpie(se_obj = v$scData_spatial,
                                                  cell_types_all = v$cell_types_all,
                                                  img_path = v$img_path,
                                                  pie_scale = 0.4)
                })
            }
        })


    ##---------------Summary tab

    ##------Clean up when ending session----
    session$onSessionEnded(function(){
        prePlot()
    })
})
