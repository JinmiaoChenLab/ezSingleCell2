## loading package
require(cytofkit2)
require(ggplot2)
require(reshape2)
require(plyr)
require(VGAM)
require(colourpicker)
require(gplots)
require(rmarkdown)
require(RColorBrewer)
require(cowplot)
require(pheatmap)

forplot3 = NULL

## Main function for scatter plot
scatterPlot <- function(obj, plotMethod, plotFunction, pointSize=1, alpha = 1,
                        addLabel=TRUE, labelSize=1, sampleLabel = TRUE,
                        FlowSOM_k = 40, selectCluster=NULL, selectSamples, 
                        facetPlot = FALSE, colorPalette = "bluered", labelRepel = FALSE, 
                        removeOutlier = TRUE, clusterColor, globalScale = TRUE, centerScale = FALSE){
  # browser()
  data <- data.frame(obj$expressionData, 
                     obj$dimReducedRes[[plotMethod]], 
                     do.call(cbind, obj$clusterRes), 
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
  
  Markers <- obj$allMarkers
  
  xlab <- colnames(obj$dimReducedRes[[plotMethod]])[1]
  ylab <- colnames(obj$dimReducedRes[[plotMethod]])[2]
  row.names(data) <- row.names(obj$expressionData)
  
  clusterMethods <- names(obj$clusterRes)
  samples <- sub("_[0-9]*$", "", row.names(obj$expressionData))
  data <- data[samples %in% selectSamples, ,drop=FALSE]
  nsamples <- samples[samples %in% selectSamples]
  data$sample <- nsamples
  sample_num <- length(unique(nsamples))
  
  if(plotFunction == "Density"){
    colPalette <- colorRampPalette(c("blue", "turquoise", "green", 
                                     "yellow", "orange", "red"))
    densCol <- densCols(data[, c(xlab, ylab)], colramp = colPalette)
    data$densCol <- densCol
    gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
      geom_point(colour=densCol, size = pointSize) + ggtitle("Density Plot") +
      theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
  }else if(plotFunction == "None"){
    gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
      geom_point(size = pointSize) + ggtitle("Dot Plot") +
      xlab(xlab) + ylab(ylab) + theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
  }else if(plotFunction == "Sample"){
    size_legend_row <- ceiling(sample_num/4)
    sample <- "sample"
    gp <- ggplot(data, aes_string(x=xlab, y=ylab, colour = sample)) +
      geom_point(size = pointSize) + ggtitle("Color By Sample") +
      xlab(xlab) + ylab(ylab) + theme_bw() + theme(legend.position = "bottom") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
      guides(colour = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
  }else if(plotFunction == "All Markers"){
    gp <- cytof_wrap_colorPlot(data = data, 
                               xlab = xlab, 
                               ylab = ylab, 
                               markers = colnames(obj$expressionData), 
                               colorPalette = colorPalette,
                               limits = NULL,
                               pointSize = pointSize, 
                               removeOutlier = TRUE)
    
  }else if(plotFunction == "All Markers(scaled)"){
    gp <- cytof_wrap_colorPlot(data = data, 
                               xlab = xlab, 
                               ylab = ylab, 
                               markers = colnames(obj$expressionData), 
                               scaleMarker = TRUE,
                               colorPalette = colorPalette,
                               limits = NULL,
                               pointSize = pointSize, 
                               removeOutlier = TRUE)
    
  }else if(plotFunction %in% clusterMethods){
    
    if(!is.null(selectCluster)){
      clusterIDs <- as.character(data[,plotFunction])
      selectCluster <- as.character(selectCluster)
      data <- data[clusterIDs %in% selectCluster, ,drop=FALSE]
    }
    clusterVec <- obj$clusterRes[[plotFunction]]
    ## make sure they are not factors before transforming to factors
    selectColors <- match(levels(as.factor(data[,plotFunction])), levels(as.factor(clusterVec)))
    clusterColor <- clusterColor[selectColors]
    
    gp <- cytof_clusterPlot(data = data, 
                            xlab = xlab, 
                            ylab = ylab, 
                            cluster = plotFunction, 
                            sample = "sample",
                            title = plotFunction, 
                            type = ifelse(facetPlot, 2, 1),
                            point_size = pointSize, 
                            addLabel = addLabel, 
                            labelSize = labelSize, 
                            sampleLabel = sampleLabel,
                            labelRepel = labelRepel,
                            fixCoord = FALSE,
                            clusterColor = clusterColor)
  }else{
    limits <- NULL
    if(globalScale){
      exprData <- obj$expressionData
      markers <- colnames(exprData)
      glimits <- quantile(exprData, probs=c(.02, .98), na.rm = TRUE)
      local.bounds <- as.data.frame(lapply(markers, function(x) quantile(exprData[,x], probs=c(.02, .98), na.rm = TRUE)), col.names = markers)
      gmax <- ifelse(max(local.bounds[2,]) < glimits[2], glimits[2], max(local.bounds[2,]))
      gmin <- ifelse(min(local.bounds[1,]) > glimits[1],min(local.bounds[1,]), glimits[1])
      limits <- c(gmin, gmax)
    }
    if(length(plotFunction) > 1){
      gp <- cytof_wrap_colorPlot(data = data, 
                                 xlab = xlab, 
                                 ylab = ylab, 
                                 markers = plotFunction, 
                                 colorPalette = colorPalette,
                                 limits = limits,
                                 scaleMarker = centerScale,
                                 pointSize = pointSize,
                                 alpha = alpha,
                                 removeOutlier = TRUE)
    }else{
      gp <- cytof_colorPlot(data = data, 
                            xlab = xlab, 
                            ylab = ylab, 
                            zlab = plotFunction, 
                            colorPalette = colorPalette,
                            limits = limits,
                            pointSize = pointSize,
                            alpha = alpha,
                            removeOutlier = TRUE)
    }
  }
  
  return(gp)
}

## Facet wrap plot of marker expression
cytof_wrap_colorPlot <- function(data, xlab, ylab, markers, scaleMarker = FALSE,
                                 colorPalette = c("bluered", "spectral1", "spectral2", "heat"),
                                 limits = NA,
                                 pointSize=1,
                                 alpha = 1,
                                 removeOutlier = TRUE){
  
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.02, .98), na.rm = na.rm, ...)
    x[x <= qnt[1]] <- qnt[1]
    x[x >= qnt[2]] <- qnt[2]
    x
  }
  
  data <- as.data.frame(data)
  title <- "Marker Expression Level Plot"
  data <- data[,c(xlab, ylab, markers)]
  
  if(removeOutlier){
    for(m in markers){
      data[[m]] <- remove_outliers(data[ ,m])
    }
  }
  
  if(scaleMarker){
    data[ ,markers] <- scale(data[ ,markers], center = TRUE, scale = TRUE)
    ev <- "ScaledExpression"
    data <- melt(data, id.vars = c(xlab, ylab), 
                 measure.vars = markers,
                 variable.name = "markers", 
                 value.name = ev)
  }else{
    ev <- "Expression"
    data <- melt(data, id.vars = c(xlab, ylab), 
                 measure.vars = markers,
                 variable.name = "markers", 
                 value.name = ev)
  }
  
  
  colorPalette <- match.arg(colorPalette)
  switch(colorPalette,
         bluered = {
           myPalette <- colorRampPalette(c("blue", "white", "red"))
         },
         spectral1 = {
           myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                           "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                           "#F46D43", "#D53E4F", "#9E0142"))
         },
         spectral2 = {
           myPalette <- colorRampPalette(rev(c("#7F0000","red","#FF7F00","yellow","white", 
                                               "cyan", "#007FFF", "blue","#00007F")))
         },
         heat = {
           myPalette <- colorRampPalette(heat.colors(50))
         }
  )
  zlength <- nrow(data)
  grid_row_num <- round(sqrt(length(markers)))
  gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) + 
    facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
    scale_colour_gradientn(limits = limits, name = ev, colours = myPalette(zlength * 2)) +
    geom_point(size = pointSize, alpha = alpha) + theme_bw() + coord_fixed() +
    theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
  
  return(gp)
}

## Heat Map
heatMap <- function(data, clusterMethod = "DensVM", type = "mean", 
                    dendrogram = "both", colPalette = "bluered", selectSamples, selectMarkers = NULL,
                    cex_row_label = 1, cex_col_label = 1, scaleMethod = "none") {
  exprs <- data$expressionData
  samples <- sub("_[0-9]*$", "", row.names(exprs))
  if(!(is.null(selectMarkers))) {
    marker_id <- selectMarkers
  }else{
    marker_id <- colnames(exprs)
  }
  markers <- colnames(exprs)
  mySamples <- samples %in% selectSamples
  myMarkers <- markers %in% marker_id
  exprs <- exprs[mySamples, , drop = FALSE]
  exprs <- exprs[, myMarkers, drop = FALSE]
  dataj <- data$clusterRes[[clusterMethod]][mySamples]
  exprs_cluster <- data.frame(exprs, cluster = dataj, check.names = FALSE )
  
  cluster_stat <- cytof_clusterStat(data = exprs_cluster,
                                    cluster = "cluster", 
                                    statMethod = type)
  
  cytof_heatmap(data = as.matrix(cluster_stat), 
                baseName = paste(clusterMethod, type), 
                scaleMethod = scaleMethod, 
                dendrogram = dendrogram,
                colPalette = colPalette,
                cex_row_label = cex_row_label, 
                cex_col_label = cex_col_label,
                margins = c(8, 8), 
                keysize = 1, 
                key.par=list(mgp=c(1.5, 0.5, 0), mar=c(3, 2.5, 3.5, 1))) 
}

## density plot

#' @param densData Data frame.
#' @param stackRotation Rotation degree of density plot to the right side, range (0-90).
#' @param stackSeperation Control factor for stack seperation interval, numeric value from 0-1, or auto.
#'
#' @importFrom plyr ldply
#' @importFrom reshape2 melt
#' @import ggplot2
stackDenistyPlot <- function(data, densityCols, stackFactor,
                             kernel = c("gaussian", "epanechnikov", "rectangular",
                                        "triangular", "biweight",
                                        "cosine", "optcosine"),
                             bw = "nrd0", adjust = 1,
                             reomoveOutliers = FALSE, 
                             stackRotation = 0, 
                             stackSeperation = "auto",
                             x_text_size = 2, 
                             strip_text_size = 7,
                             legend_text_size = 0.5, 
                             legendRow = 1,
                             legend_title = "stackName",
                             stackFactorColours = NULL){
  
  if(!is.numeric(stackRotation)){
    stop("stackRotation must be a numeric number")
  }else if(stackRotation < 0 || stackRotation > 90){
    stop("stackRotation must be a numeric number in range 0-90")
  }
  
  if(missing(densityCols)){
    densityCols <- colnames(data)
  }else if(any(!(densityCols %in% colnames(data)))){
    stop("Unmatch densityCols found:", paste(densityCols[!(densityCols %in% colnames(data))], collapse = " "))
  }
  
  if(missing(stackFactor)){
    warning("no stackFactor was provided!")
    stackFactor <- rep("stack", length = nrow(data))
  }else if(length(stackFactor) != nrow(data)){
    stop("Length of stackFactor unequal row number of input data")
  }
  kernel <- match.arg(kernel)
  
  stackCount <- length(unique(stackFactor))
  densityCount <- length(densityCols)
  
  if(missing(stackFactorColours) || is.null(stackFactorColours)){
    stackFactorColours <- rainbow(stackCount)
  }else if(length(stackFactorColours) == 0 || length(stackFactorColours) != stackCount){
    stackFactorColours <- rainbow(stackCount)
  }
  
  data <- data.frame(data[ ,densityCols, drop=FALSE], stackFactor = stackFactor, check.names = FALSE)
  
  densData <- .densityCal(data, kernel = kernel, bw = bw, adjust = adjust, reomoveOutliers = reomoveOutliers)
  ## dataframe densData contains {stackName, x , y , densityName}
  xStat <- aggregate(x ~ stackName + densityName, densData, max)
  yStat <- aggregate(y ~ stackName + densityName, densData, max)
  
  if(stackSeperation == "auto"){
    stackIntervals <- aggregate(y ~ densityName, yStat, function(x){0.8*median(x) * (1-(stackRotation/90)^0.2)^2})
  }else if(stackSeperation < 0 || stackSeperation > 1){
    stop("stackSeperation must be value in range 0-1")
  }else{
    stackIntervals <- aggregate(y ~ densityName, yStat, function(x){median(x)*stackSeperation})
  }
  
  stackShifts <- aggregate(x ~ densityName, xStat, function(x){max(x) * (stackRotation/90)})
  
  densData$stack_x <- densData$x + (as.numeric(densData$stackName)-1) * stackShifts$x[match(densData$densityName, stackShifts$densityName)]
  densData$stack_y <- densData$y + (as.numeric(densData$stackName)-1) * stackIntervals$y[match(densData$densityName, stackIntervals$densityName)]
  
  ## segment lines, x tick, x label
  alignSegments <- ldply(split(densData$x, densData$densityName),
                         function(x){seq(min(x), max(x), length.out=5)},
                         .id = "densityName")
  alignSegments <- melt(alignSegments, id.vars="densityName", variable.name="x_tick", value.name = "x")
  alignSegments$y <- min(densData$y)
  alignSegments$xend <- alignSegments$x + (length(unique(densData$stackName))-1) * stackShifts$x[match(alignSegments$densityName, stackShifts$densityName)]
  alignSegments$yend <- min(densData$y) + (length(unique(densData$stackName))-1) * stackIntervals$y[match(alignSegments$densityName, stackIntervals$densityName)]
  
  densityHeights <- aggregate(y ~ densityName, yStat, max)
  alignSegments$tickXend <- alignSegments$x
  alignSegments$tickYend <- alignSegments$y - densityHeights$y[match(alignSegments$densityName, densityHeights$densityName)] * 0.01
  alignSegments$tickText <- format(alignSegments$x,scientific=TRUE, digits=3)
  alignSegments$textY <- alignSegments$y - densityHeights$y[match(alignSegments$densityName, densityHeights$densityName)] * 0.03
  
  message(" Plotting ...\n")
  stackDensityPlot_theme <- theme(legend.position = "top",
                                  legend.title = element_text(size = rel(1)),
                                  legend.text = element_text(size = rel(legend_text_size)),
                                  strip.text = element_text(size=strip_text_size, lineheight=1, hjust = 0.5, vjust = 0.5),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_blank(),
                                  strip.background=element_rect(fill = "grey90", colour = NA))
  
  gp <- ggplot(densData, aes(x=stack_x, y=stack_y)) +
    geom_segment(data = alignSegments,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey80", size=0.3) +
    geom_segment(data = alignSegments,
                 aes(x = x, y = y, xend = tickXend, yend = tickYend),
                 color = "grey20", size=0.3) +
    geom_text(data = alignSegments, aes(x = x, y = textY, label = tickText),
              hjust = 0.3, vjust = 1.1, size = x_text_size) +
    geom_polygon(aes(fill=stackName, color=stackName), alpha = 0.15) + 
    scale_colour_manual(values = stackFactorColours) + 
    scale_fill_manual(values = stackFactorColours) +
    facet_wrap(~densityName, scale = "free") +
    xlab("") + ylab("") +
    guides(col = guide_legend(title = legend_title, nrow = legendRow, byrow = TRUE),
           fill = guide_legend(title = legend_title, nrow = legendRow, byrow = TRUE)) +
    theme_bw() + stackDensityPlot_theme
  
  gp
}


#' Internal density calculation function serves for \code{stackDenistyPlot}
#'
#' Output data frame with columns: stackName, x , y , densityName
.densityCal <- function(data, kernel, bw, adjust, reomoveOutliers = FALSE){
  message("  Calculating Density for each stack column...\n")
  print(table(data$stackFactor))
  dataBystackFactor <- split(subset(data, select = -stackFactor), data$stackFactor)
  densityWrap <- function(d, ...){
    resOut <- NULL
    for(i in colnames(d)){
      x <- d[,i]
      if(reomoveOutliers){
        message("  Remove outliers...\n")
        x_IQR <- IQR(x)
        x_lowLimit <- quantile(x, 0.25) - 1.5*x_IQR
        x_highLimit <- quantile(x, 0.75) + 1.5*x_IQR
        x <- x[x >= x_lowLimit && x <= x_highLimit]
      }
      dens <- density(x, ...)
      densOut <- data.frame(x=dens$x, y=dens$y, densityName = i)
      resOut <- rbind(resOut, densOut)
    }
    return(resOut)
  }
  
  r <- ldply(dataBystackFactor, densityWrap,
             kernel = kernel, bw = bw, adjust = adjust,
             .progress = "text",
             .id = "stackName")
  return(r)
}


## Combined marker expression trend
cytof_expressionTrends <- function(data, markers, clusters, 
                                   orderCol="isomap_1", 
                                   clusterCol = "cluster", 
                                   reverseOrder = FALSE,
                                   addClusterLabel = TRUE,
                                   clusterLabelSize = 5,
                                   segmentSize = 0.5,
                                   min_expr = NULL, 
                                   trend_formula="expression ~ sm.ns(Pseudotime, df=3)"){
  
  if(!is.data.frame(data)) data <- data.frame(data, check.names = FALSE)
  if(!all(markers %in% colnames(data))) stop("Unmatching markers found!")
  if(!(length(orderCol)==1 && orderCol %in% colnames(data)))
    stop("Can not find orderCol in data!")
  if(!(length(clusterCol)==1 && clusterCol %in% colnames(data)))
    stop("Can not find clusterCol in data!")
  if(!missing(clusters)){
    if(!all(clusters %in% data[[clusterCol]]))
      stop("Wrong clusters selected!")
    data <- data[data[[clusterCol]] %in% clusters, , drop=FALSE]
  }
  
  if(reverseOrder){
    newOrderCol <- paste0(orderCol, "(reverse)")
    data[[newOrderCol]] <- -data[[orderCol]]
    orderCol <- newOrderCol
  }
  orderValue <- data[[orderCol]]
  data <- data[order(orderValue), c(markers, clusterCol)]
  data$Pseudotime <- sort(orderValue)
  
  mdata <- melt(data, id.vars = c("Pseudotime", clusterCol), 
                variable.name = "markers", value.name= "expression")
  colnames(mdata) <- c("Pseudotime", clusterCol, "markers", "expression")
  mdata$markers <- factor(mdata$markers)
  mdata[[clusterCol]] <- factor(mdata[[clusterCol]])
  min_expr <- min(mdata$expression)
  
  ## tobit regression
  vgamPredict <- ddply(mdata, .(markers), function(x) { 
    fit_res <- tryCatch({
      vg <- suppressWarnings(vgam(formula = as.formula(trend_formula), 
                                  family = VGAM::tobit(Lower = min_expr, lmu = "identitylink"), 
                                  data = x, maxit=30, checkwz=FALSE))
      res <- VGAM::predict(vg, type="response")
      res[res < min_expr] <- min_expr
      res
    }
    ,error = function(e) {
      print("Error!")
      print(e)
      res <- rep(NA, nrow(x))
      res
    }
    )
    expectation = fit_res
    data.frame(Pseudotime=x[["Pseudotime"]], expectation=expectation)
  })
  
  color_by <- clusterCol
  plot_cols <- round(sqrt(length(markers)))
  cell_size <- 1
  x_lab <- orderCol
  y_lab <- "Expression"
  legend_title <- "Cluster"
  
  ## copied from monocle package
  monocle_theme_opts <- function(){
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
      #theme(panel.border = element_blank(), axis.line = element_line()) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
      theme(panel.background = element_rect(fill='white')) +
      theme(legend.position = "right") +
      theme(axis.title = element_text(size = 15)) +
      theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))}
  
  q <- ggplot(data=vgamPredict, aes_string(x="Pseudotime", y="expectation", col="markers")) + geom_line(size = 1.5)
  q <- q + ylab(y_lab) + xlab(x_lab) + theme_bw()
  q <- q + guides(colour = guide_legend(title = legend_title, override.aes = list(size = cell_size*3)))
  q <- q + monocle_theme_opts() 
  
  # if(addClusterLabel){
  #     # edata <- data[ ,c("Pseudotime", clusterCol)]
  #     # colnames(edata) <- c('x', "z")
  #     # center <- aggregate(x ~ z, data = edata, median)
  #     # center$y <- -0.5 ## add to the botom
  #     # q <- q + geom_text_repel(data=center, aes(x=x, y=y, label=z), parse=TRUE)
  #     mdata$cluster <- mdata[[clusterCol]]
  #     center <- aggregate(cbind(Pseudotime, expression) ~ cluster + markers, data = mdata, median)
  #     q <- q + geom_text_repel(data=center, aes(x=Pseudotime, y=expression, label=cluster),
  #                              size = clusterLabelSize, fontface = 'bold',
  #                              box.padding = unit(0.5, 'lines'),
  #                              point.padding = unit(1.6, 'lines'),
  #                              segment.color = '#555555',
  #                              segment.size = segmentSize,
  #                              arrow = arrow(length = unit(0.02, 'npc')))
  # }
  
  q
}

## function for opening the results directory
opendir <- function(dir = getwd()){
  if (.Platform['OS.type'] == "windows"){
    shell.exec(dir)
  } else {
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}

read_string = function (file_name) 
{
  content = readChar(file_name, file.info(file_name)$size)
  return(content)
}

write_string = function (content, output_name) 
{
  f <- file(output_name, "wb")
  tryCatch({
    writeBin(charToRaw(content), f)
  }, finally = {
    close(f)
  })
}

ezcbind = function (..., type = c("sync", "keep_all", "all", "outer", "keep_left", 
                                  "left", "keep_right", "right", "keep_common", "common", "inner", 
                                  "cross")) 
{
  type = tolower(type[1])
  if (type == "sync") {
    res = sync_cbind(...)
  }
  else {
    if (type %in% c("keep_common", "common", "inner")) {
      argl = list(...)
      if (length(argl) == 1 && class(argl[[1]]) == "list") {
        argl = argl[[1]]
      }
      if (all(is.null(argl))) {
        return(NULL)
      }
      res = NULL
      argl = argl[!sapply(argl, is.null)]
      if (length(argl) > 0) {
        common_genes = rownames(argl[[1]])
        if (length(argl) > 1) {
          for (i in 2:length(argl)) {
            common_genes = intersect(common_genes, rownames(argl[[i]]))
          }
          res = argl[[1]][common_genes, , drop = FALSE]
          for (i in 2:length(argl)) {
            res = cbind(res, argl[[i]][common_genes, 
                                       , drop = FALSE])
          }
        }
      }
    }
  }
  return(res)
  argl = list(...)
  if (length(argl) == 1 && class(argl[[1]]) == "list") {
    argl = argl[[1]]
  }
  if (all(is.null(argl))) {
    return(NULL)
  }
  argl = argl[!sapply(argl, is.null)]
  if (length(argl) == 1) {
    return(as.data.frame(argl[[1]]))
  }
  else {
    res = argl[[1]]
    if (class(res) != "data.frame") {
      res = as.data.frame(res)
    }
    res$row_name = rownames(res)
    res = setkey(as.data.table(res), "row_name")
    for (i in 2:length(argl)) {
      print(i)
      temp1 = res
      temp2 = argl[[i]]
      if (class(temp2) != "data.frame") {
        temp2 = as.data.frame(temp2)
      }
      temp2$row_name = rownames(temp2)
      temp2 = setkey(as.data.table(temp2), "row_name")
      if (type == "keep_all" || type == "all" || type == 
          "outer") {
        res = merge(temp1, temp2, by = "row_name", all = TRUE, 
                    allow.cartesian = TRUE)
      }
      else if (type == "keep_left" || type == "left") {
        res = temp2[temp1, allow.cartesian = TRUE]
      }
      else if (type == "keep_right" || type == "right") {
        res = temp1[temp2, allow.cartesian = TRUE]
      }
      else if (type == "keep_common" || type == "common" || 
               type == "inner") {
        res = temp1[temp2, nomatch = 0L, allow.cartesian = TRUE]
      }
      else if (type == "cross") {
        res = merge(res, argl[[i]], by = NULL)
      }
      res = as.data.frame(res)
      rownames(res) = res[, 1]
      res = res[, -1, drop = FALSE]
    }
    return(res)
  }
}

sync_cbind = function (...) 
{
  argl = list(...)
  if (length(argl) == 1 && class(argl[[1]]) == "list") {
    argl = argl[[1]]
  }
  if (all(is.null(argl))) {
    return(NULL)
  }
  argl = argl[!sapply(argl, is.null)]
  if (length(argl) == 1) {
    return(as.data.frame(argl[[1]]))
  }
  else {
    res = argl[[1]]
    if (class(res) != "data.frame") {
      res = as.data.frame(res)
    }
    for (i in 2:length(argl)) {
      temp2 = argl[[i]]
      if (class(temp2) != "data.frame") {
        temp2 = as.data.frame(temp2)
      }
      temp2 = sync_variable(res, temp2)
      res = cbind(res, temp2)
      res = as.data.frame(res)
    }
    return(res)
  }
}

sync_variable = function (base_order, to_be_sorted, base_dim = 1, to_be_sorted_dim = 1) 
{
  name1 = if (base_dim == 1) 
    rownames(base_order)
  else colnames(base_order)
  name2 = if (to_be_sorted_dim == 1) 
    rownames(to_be_sorted)
  else colnames(to_be_sorted)
  if (all(name1 %in% name2) || all(name2 %in% name1)) {
    name_i = intersect(name1, name2)
    sorted = if (to_be_sorted_dim == 1) 
      to_be_sorted[name_i, , drop = FALSE]
    else to_be_sorted[, name_i, drop = FALSE]
  }
  else {
    stop("Can not synchronize two variables!")
  }
}

fast_aggr = function (mat, col_id, func = sum) 
{
  dt = data.table(mat)
  col_name = colnames(mat)[col_id]
  sum_index = 1:dim(mat)[2]
  sum_index = setdiff(sum_index, col_id)
  sum_item = colnames(dt)[sum_index]
  b = dt[, lapply(.SD, func), by = col_name, .SDcols = sum_item]
  b = as.data.frame(b)
  b = remove_all_NA(b)
  rownames(b) = b[, 1]
  b = b[, -1, drop = FALSE]
  return(b)
}

remove_all_NA = function (raw_data) 
{
  raw_data = raw_data[complete.cases(raw_data), , drop = FALSE]
}

seurat_sacle_data = function (data.use, genes.use = NULL, do.scale = TRUE, do.center = TRUE, 
                              scale.max = 10) 
{
  genes.use <- rownames(data.use)
  scale.data <- matrix(data = NA, nrow = length(x = genes.use), 
                       ncol = ncol(x = data.use))
  dimnames(x = scale.data) <- dimnames(x = data.use)
  if (do.scale | do.center) {
    bin.size <- 1000
    max.bin <- floor(length(genes.use)/bin.size) + 1
    message("Scaling data matrix")
    pb <- txtProgressBar(min = 0, max = max.bin, style = 3, 
                         file = stderr())
    for (i in 1:max.bin) {
      my.inds <- ((bin.size * (i - 1)):(bin.size * i - 
                                          1)) + 1
      my.inds <- my.inds[my.inds <= length(x = genes.use)]
      new.data <- t(x = scale(x = t(x = as.matrix(x = data.use[genes.use[my.inds], 
                                                               ])), center = do.center, scale = do.scale))
      new.data[new.data > scale.max] <- scale.max
      scale.data[genes.use[my.inds], ] <- new.data
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  scale.data[is.na(scale.data)] <- 0
  return(scale.data)
}

eval_string = function (string, envir = parent.frame(), ...) 
{
  eval(parse(text = string), envir = envir, ...)
}

ggsettings = new.env()
ggsettings$pointsize = 1.5
ggsettings$alpha = 0.7
ggsettings$title_size = 10
ggsettings$axis_title_size = 9
ggsettings$axis_label_size = 8
ggsettings$legend_title_size = 9
ggsettings$legend_label_size = 8
ggsettings$page_width = 4
ggsettings$page_height = 3
ggsettings$stroke = 0.2
ggsettings$dpi = 300
plot_scatter = function (pos, color = NULL, colors = NULL, colors_lims = NULL, 
                         shape = NULL, global_point_size = get_ggplot_settings("pointsize"), 
                         size = NULL, size_breaks = NULL, size_labels = NULL, color_as_factor = "auto", 
                         title = NULL, return_type = "ggplot") 
{
  if (is.null(color)) {
    color = "Undefined"
  }
  pos = as.data.frame(pos)
  if (is.character(color)) {
    color = create_group(list(rownames(pos)), c(color))
  }
  color = as.data.frame(color)
  color = color[ez_order(color[, 1]), , drop = FALSE]
  pos = sync_variable(color, pos)
  color = sync_variable(pos, color)
  if (color_as_factor == "auto") {
    color_as_factor = ifelse((nrow(color)/length(unique(color[, 
                                                              1]))) > 4, yes = TRUE, no = FALSE)
  }
  if (color_as_factor) {
    color[, 1] = factor(color[, 1], levels = ez_sort(unique(color[, 
                                                                  1])))
  }
  aes_param = list(x = pos[, 1], y = pos[, 2], color = color[, 
                                                             1], shape = shape[, 1], text = rownames(pos), key = rownames(pos))
  main_param = list(data = pos, alpha = ggsettings$alpha, pch = 19, 
                    stroke = ggsettings$stroke)
  if (!is.null(size)) {
    size = sync_variable(color, size)
    aes_param = ez_param(aes_param, size = size[, 1])
  }
  else {
    size = global_point_size
    main_param = ez_param(main_param, size = size)
  }
  main_param = ez_param(main_param, mapping = do.call(aes, 
                                                      aes_param))
  p = ggplot() + do.call(geom_point, main_param) + theme_classic() + 
    labs(color = colnames(color)[1], x = colnames(pos)[1], 
         y = colnames(pos)[2]) + theme(title = element_text(size = ggsettings$title_size), 
                                       plot.title = element_text(hjust = 0.5), axis.text = element_text(size = ggsettings$axis_title_size), 
                                       axis.title = element_text(size = ggsettings$axis_label_size), 
                                       legend.title = element_text(size = ggsettings$legend_label_size), 
                                       legend.text = element_text(size = ggsettings$legend_label_size)) + 
    labs(title = title)
  if (!is.null(size)) {
    p = p + labs(size = colnames(size)[1])
  }
  size_param = list()
  if (!is.null(size_breaks) && is.numeric(size[, 1])) {
    size_param[["breaks"]] = size_breaks
  }
  if (!is.null(size_labels)) {
    size_param[["labels"]] = size_labels
  }
  if (length(size_param) > 0) {
    p = p + do.call(scale_size_continuous, args = ez_param(size_param))
  }
  if (is.numeric(color[, 1])) {
    if (is.null(colors)) {
      colors = "red"
    }
    colors = tolower(colors)
    if (colors == "red") {
      colours = brewer.pal(10, "Reds")
    }
    else if (colors == "blue") {
      colours = rev(brewer.pal(11, "RdYlBu"))
    }
    else if (colors == "rainbow") {
      colours = rev(rainbow(10, end = 0.7))
    }
    else {
      colours = colors
    }
    p = p + scale_color_gradientn(colours = colours, limits = colors_lims)
  }
  else if (is.factor(color[, 1]) | is.character(color[, 1])) {
    if (!is.null(colors)) {
      p = p + scale_color_manual(values = colors)
    }
  }
  if (return_type == "ggplot") {
    return(p)
  }
  else if (return_type == "plotly") {
    return(ggplotly(p))
  }
}


ez_order = function (dt, order_nchar_first = TRUE) 
{
  if (is.character(dt) && order_nchar_first) {
    res = order(nchar(dt), dt)
  }
  else {
    res = order(dt)
  }
  res
}

ez_sort = function (dt, sort_nchar_first = TRUE) 
{
  res = dt[ez_order(dt, sort_nchar_first)]
  res
}

get_ggplot_settings = function (name) 
{
  get(name, envir = ggsettings)
}

ez_param = function (...) 
{
  params = list(...)
  params = params[!sapply(params, is.null)]
  params = ez_flattern_list(params)
  res = list()
  for (i in seq_along(params)) {
    res[[names(params)[i]]] = params[[i]]
  }
  return(res)
}

ez_flattern_list = function (x, recursive = TRUE, use.name = TRUE) 
{
  temp = x
  res = list()
  num = 1
  finished = FALSE
  for (i in seq_along(temp)) {
    if (class(temp[[i]])[1] == "list") {
      for (j in seq_along(temp[[i]])) {
        if (is.null(temp[[i]][[j]])) {
          res[num] = list(NULL)
        }
        else {
          res[[num]] = temp[[i]][[j]]
        }
        if (use.name) {
          names(res)[num] = names(temp[[i]])[j]
        }
        num = num + 1
      }
    }
    else {
      res[[num]] = temp[[i]]
      if (use.name) {
        if (!is.null(names(temp)) && !is.null(names(temp)[i])) {
          names(res)[num] = names(temp)[i]
        }
      }
      num = num + 1
    }
  }
  if (recursive) {
    while (!finished) {
      finished = TRUE
      temp = res
      res = list()
      num = 1
      for (i in seq_along(temp)) {
        if (class(temp[[i]])[1] == "list") {
          for (j in seq_along(temp[[i]])) {
            if (is.null(temp[[i]][[j]])) {
              res[num] = list(NULL)
            }
            else {
              res[[num]] = temp[[i]][[j]]
            }
            if (use.name) {
              names(res)[num] = names(temp[[i]])[j]
            }
            num = num + 1
          }
          finished = FALSE
        }
        else {
          if (is.null(temp[[i]])) {
            res[num] = list(NULL)
          }
          else {
            res[[num]] = temp[[i]]
          }
          if (use.name) {
            names(res)[num] = names(temp)[i]
          }
          num = num + 1
        }
      }
    }
  }
  return(res)
}

plot_split_scatter = function (pos, color = NULL, colors = NULL, shape = NULL, size = NULL, 
                               color_as_factor = FALSE, order = NULL, sort_color = TRUE, title = NULL, 
                               return_type = "ggplot", ncol = NULL, show_legend = TRUE, coor_fix = TRUE, 
                               ...) 
{
  if (is.null(order)) {
    if (is.factor(color[, 1])) {
      sub_type = levels(color[, 1])
    }
    else {
      sub_type = unique(color[, 1])
    }
    if (sort_color) {
      sub_type = ez_sort(sub_type)
    }
  }
  else {
    sub_type = order
  }
  if (is.null(colors)) {
    colors = gg_color_hue(length(sub_type))
  }
  pos = sync_variable(color, pos)
  color = sync_variable(pos, color)
  res = lapply(1:length(sub_type), function(i) {
    sub_cells = get_selected_clusters_cell_name(color, sub_type[i])
    p = plot_underline_points(pos, sub_cells, colors = colors[i], 
                              selected_points_name = sub_type[i], title = sub_type[i])
    if (!show_legend) {
      p = p + gg_remove_all_legend()
    }
    if (coor_fix) {
      p = p + coord_fixed()
    }
    p
  })
  names(res) = sub_type
  p_all_type = plot_scatter(pos = pos, color = color, colors = colors, 
                            shape = shape, color_as_factor = color_as_factor, title = title, 
                            return_type = return_type, ...)
  if (is.null(ncol)) {
    ncol = ceiling(sqrt(length(res)))
  }
  p = plot_grid(plotlist = res, ncol = ncol)
  final_res = list(one_figure = p, all_figures = res, all_type = p_all_type)
  return(final_res)
}

ez_chunk = function (x, n) 
{
  split(x, cut(seq_along(x), n, labels = FALSE))
}

gg_color_hue = function (n) 
{
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

get_selected_clusters_cell_name = function (cell_type, selected_clusters) 
{
  selected_cell_name = rownames(cell_type)[cell_type[, 1] %in% 
                                             selected_clusters]
  return(selected_cell_name)
}

plot_underline_points = function (pos, selected_points, color = NULL, colors = "red", 
                                  background_color = "gray", global_point_size = get_ggplot_settings("pointsize"), 
                                  selected_points_name = "selected", background_name = "Others", 
                                  ...) 
{
  cell_type = pos[, 1, drop = FALSE]
  colnames(cell_type) = ""
  cell_type[, 1] = background_name
  cell_type[selected_points, 1] = selected_points_name
  cell_type = as.data.frame(cell_type)
  cell_type = dataframe_reorder(cell_type, n_col = 1, order = c(background_name, 
                                                                selected_points_name))
  if (!is.null(color)) {
    common_cells = intersect(selected_points, rownames(color))
  }
  plot_scatter(pos, color = cell_type, colors = c(background_color, 
                                                  colors), global_point_size = global_point_size, ...) + 
    labs(color = "")
}

dataframe_reorder = function (ori_data, n_col = 1, order = NULL, order_item = TRUE) 
{
  if (is.null(order)) {
    print(paste0("c(\"", paste(unique(ori_data[, n_col]), 
                               collapse = "\", \""), "\")"), quote = FALSE)
    return()
  }
  else {
    ori_data[, n_col] = factor(ori_data[, n_col], levels = order)
  }
  if (order_item) {
    ori_data = ori_data[order(ori_data[, n_col]), , drop = FALSE]
  }
  ori_data
}

gg_remove_all_legend = function () 
{
  theme(legend.position = "none")
}