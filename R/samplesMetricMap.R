#' Plot a Grid of Sample-wise Predictive Metrics
#' 
#' A grid of coloured tiles is drawn. There is one column for each sample and
#' one row for each cross-validation result.
#' 
#' The names of \code{results} determine the row names that will be in the
#' plot. The length of \code{metricColours} determines how many bins the metric
#' values will be discretised to.
#' 
#' @aliases samplesMetricMap samplesMetricMap,list-method
#' samplesMetricMap,matrix-method
#' @param results A list of \code{\link{ClassifyResult}} objects. Could also be
#' a \code{matrix} of pre-calculated metrics, for backwards compatibility.
#' @param classes If \code{results} is a matrix, this is a factor vector of the
#' same length as the number of columns that \code{results} has.
#' @param comparison Default: \code{"auto"}. The aspect of the experimental
#' design to compare. Can be any characteristic that all results share.
#' @param metric Default: \code{"auto"}. The name of the
#' performance measure or "auto". If the results are classification then
#' sample accuracy will be displayed. Otherwise, the results would be survival risk
#' predictions and then a sample C-index will be displayed. Valid values are \code{"Sample Error"},
#' \code{"Sample Error"} or \code{"Sample C-index"}. If the metric is not stored in the
#' results list, the performance metric will be calculated automatically.
#' @param featureValues If not NULL, can be a named factor or named numeric
#' vector specifying some variable of interest to plot above the heatmap.
#' @param featureName A label describing the information in
#' \code{featureValues}. It must be specified if \code{featureValues} is.
#' @param metricColours If the outcome is categorical, a list of vectors of colours
#' for metric levels for each class. If the outcome is numeric, such as a risk score,
#' then a single vector of colours for the metric levels for all samples.
#' @param classColours Either a vector of colours for class levels if both
#' classes should have same colour, or a list of length 2, with each component
#' being a vector of the same length. The vector has the colour gradient for
#' each class.
#' @param groupColours A vector of colours for group levels. Only useful if
#' \code{featureValues} is not NULL.
#' @param fontSizes A vector of length 5. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values. The fourth number is the size of the
#' legends' titles. The fifth number is the font size of the legend labels.
#' @param mapHeight Height of the map, relative to the height of the class
#' colour bar.
#' @param title The title to place above the plot.
#' @param showLegends Logical. IF FALSE, the legend is not drawn.
#' @param xAxisLabel The name plotted for the x-axis. NULL suppresses label.
#' @param showXtickLabels Logical. IF FALSE, the x-axis labels are hidden.
#' @param showYtickLabels Logical. IF FALSE, the y-axis labels are hidden.
#' @param yAxisLabel Default: \code{"auto"} for list of \code{ClassifyResult}s and
#' \code{"Analysis"} for a \code{matrix}. The axis name plotted for the y-axis.
#' If \code{"auto"}, automatically set depending on value of \code{comparison}.
#' \code{NULL} suppresses the label.
#' @param legendSize The size of the boxes in the legends.
#' @param ... Parameters not used by the \code{ClassifyResult} method that does
#' list-packaging but used by the main \code{list} method.
#' @return A grob is returned that can be drawn on a graphics device.
#' @author Dario Strbenac
#' @examples
#' 
#'   predicted <- DataFrame(sample = LETTERS[sample(10, 100, replace = TRUE)],
#'                           class = rep(c("Healthy", "Cancer"), each = 50))
#'   actual <- factor(rep(c("Healthy", "Cancer"), each = 5), levels = c("Healthy", "Cancer"))
#'   features <- sapply(1:100, function(index) paste(sample(LETTERS, 3), collapse = ''))
#'   result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Example", "t-test", "Differential Expression", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], features, list(1:100), list(sample(10, 10)),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   predicted[, "class"] <- sample(predicted[, "class"])
#'   result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Example", "Bartlett Test", "Differential Variability", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], features, list(1:100), list(sample(10, 10)),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   result1 <- calcCVperformance(result1)
#'   result2 <- calcCVperformance(result2)
#'   groups <- factor(rep(c("Male", "Female"), length.out = 10))
#'   names(groups) <- LETTERS[1:10]
#'   cholesterol <- c(4.0, 5.5, 3.9, 4.9, 5.7, 7.1, 7.9, 8.0, 8.5, 7.2)
#'   names(cholesterol) <- LETTERS[1:10]
#'   
#'   wholePlot <- samplesMetricMap(list(Gene = result1, Protein = result2))
#'   wholePlot <- samplesMetricMap(list(Gene = result1, Protein = result2),
#'                                 featureValues = groups, featureName = "Gender")
#'   wholePlot <- samplesMetricMap(list(Gene = result1, Protein = result2),
#'                                 featureValues = cholesterol, featureName = "Cholesterol")                                
#'
#' @export
#' @usage NULL
setGeneric("samplesMetricMap", function(results, ...)
standardGeneric("samplesMetricMap"))

#' @rdname samplesMetricMap
#' @export
setMethod("samplesMetricMap", "ClassifyResult", function(results, ...) {
    samplesMetricMap(list(assay = results), ...)
})

#' @rdname samplesMetricMap
#' @export
setMethod("samplesMetricMap", "list", 
          function(results,
                   comparison = "auto",
                   metric = "auto",
                   featureValues = NULL, featureName = NULL,
                   metricColours = list(c("#FFFFFF", "#CFD1F2", "#9FA3E5", "#6F75D8", "#3F48CC"),
                                        c("#FFFFFF", "#E1BFC4", "#C37F8A", "#A53F4F", "#880015")),
                   classColours = c("#3F48CC", "#880015"), groupColours = c("darkgreen", "yellow2"),
                   fontSizes = c(24, 16, 12, 12, 12),
                   mapHeight = 4, title = "auto",
                   showLegends = TRUE, xAxisLabel = "Sample Name", showXtickLabels = TRUE,
                   yAxisLabel = "auto", showYtickLabels = TRUE, legendSize = grid::unit(1, "lines"))
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")  
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")       
  if(!requireNamespace("gtable", quietly = TRUE))
    stop("The package 'gtable' could not be found. Please install it.")
  
  characteristicsCounts <- table(unlist(lapply(results, function(result) result@characteristics[["characteristic"]])))
  if(comparison == "auto")
  {
    if(max(characteristicsCounts) == length(results))
    { # Choose a characteristic which varies the most across the results.
      candidates <- names(characteristicsCounts)[characteristicsCounts == length(results)]
      allCharacteristics <- do.call(rbind, lapply(results, function(result) result@characteristics))
      distinctValues <- by(allCharacteristics[, "value"], allCharacteristics[, "characteristic"], function(values) length(unique(values)))
      comparison <- names(distinctValues)[which.max(distinctValues)][1]
    } else {
      stop("No characteristic is present for all results but must be.")
    }
  }
  if(yAxisLabel == "auto")
    yAxisLabel <- comparison
  isSurvival <- "risk" %in% colnames(results[[1]]@predictions)
  validMetrics <- c("Sample Error", "Sample Accuracy", "Sample C-index")
  if(metric == "auto")
    metric <- ifelse(isSurvival, "Sample C-index", "Sample Accuracy")
  else
    if(!metric %in% validMetrics) stop("metric must be one of ", validMetrics, " but is ", metric, '.')   
  if(title == "auto") title <- switch(metric, `Sample Error` = "Error Comparison", `Sample Accuracy` = "Accuracy Comparison", `Sample C-index` = "Risk Score Comparison")
  if(isSurvival && is.list(metricColours)) metricColours <- metricColours[[1]]
  metricText <- gsub("Sample ", '', metric) # For legend labelling.
  if(showXtickLabels == FALSE && xAxisLabel == "Sample Name") xAxisLabel <- "Sample"
  
  resultsWithComparison <- sum(sapply(results, function(result) any(result@characteristics[, "characteristic"] == comparison)))
  if(resultsWithComparison < length(results))
    stop("Not all results have comparison characteristic ", comparison, ' but need to.')
  
  if(comparison == "Cross-validation")
    compareFactor <- sapply(results, function(result) .validationText(result))  
  else
    compareFactor <- sapply(results, function(result) {
                     useRow <- result@characteristics[, "characteristic"] == comparison
                     result@characteristics[useRow, "value"]
                    })
  
  metrics <- unlist(lapply(results, function(result)
    if(!is.null(result@performance)) names(result@performance)))
  namesCounts <- table(metrics)
  commonNames <- names(namesCounts)[namesCounts == length(results)]
  if(!metric %in% commonNames)
  {
    warning(paste(metric, "not found in all elements of results. Calculating it now."))
    results <- lapply(results, function(result) calcCVperformance(result, metric))
  }
  if(!is.null(featureValues) && is.null(featureName))
    stop("featureValues is specified by featureNames isn't. Specify both.")
  if(!is.null(featureValues) && is.null(names(featureValues)))
    stop("featureValues vector must be named with sample IDs.")
  comparisonValuesCounts <- table(compareFactor)
  if(any(comparisonValuesCounts > 1))
    stop("Some classification results have same the comparison value. Check that each
  clasification result is distinctive for the comparison type specified by
  'comparison'.")
  
  if(length(unique(lengths(lapply(results, function(result) unique(predictions(result)[, "sample"]))))) > 1)
      stop("Cross-validation results contain different sets of samples.")
  
  nColours <- if(is.list(metricColours)) length(metricColours[[1]]) else length(metricColours)
  metricBinEnds <- seq(0, 1, 1/nColours)

  metricValues <- lapply(results, function(result)
  {
    sampleMetricValues <- result@performance[[metric]]
    cut(sampleMetricValues, metricBinEnds, include.lowest = TRUE)
  })
  
  if(metric != "Sample C-index")
  {
    knownClasses <- actualOutcome(results[[1]])
    classedMetricValues <- lapply(metricValues, function(metricSet)
    {
      metricSet <- factor(paste(knownClasses, metricSet, sep = ','),
                         levels = c(t(outer(levels(knownClasses), levels(metricSet), paste, sep = ','))))
    })
  } else {knownClasses <- NULL}

  metricMatrix <- do.call(rbind, metricValues)
  meanSample <- colMeans(metricMatrix, na.rm = TRUE)
  rowOrder <- hclust(dist(metricMatrix, "manhattan"))[["order"]]
  metricMatrixSamplewise <- t(metricMatrix)
  uninformativeSamples <- which(apply(metricMatrixSamplewise, 1, function(sampleMetrics) all(is.na(sampleMetrics))))
  if(length(uninformativeSamples) > 0)
    metricMatrixSamplewise[uninformativeSamples, ] <- 0
  
  if(is.null(featureValues))
  {
    if(metric != "Sample C-index") # Sort within each class.
    {
      colOrder <- unlist(lapply(levels(knownClasses), function(class)
                         {
                           classIndices <- which(as.character(knownClasses) == class)
                           classMatrix <- metricMatrixSamplewise[classIndices, ]
                           classIndices[hclust(dist(classMatrix, "manhattan"))[["order"]]]
                         }))
    } else { # Sort all samples together.
        colOrder <- hclust(dist(metricMatrixSamplewise, "manhattan"))[["order"]]
    }
  } else { # Primarily sort by some sample characteristic and only experimental values if there are characteristic ties.
    featureValues <- featureValues[match(sampleNames(results[[1]]), names(featureValues))]
    if(metric != "Sample C-index") # Sort within each class.
      colOrder <- order(knownClasses, featureValues, meanSample)
    else # Sort all samples together.
      colOrder <- order(featureValues, meanSample)
  }
  if(metric != "Sample C-index")
    knownClasses <- knownClasses[colOrder]
  if(!is.null(featureValues))
    featureValues <- featureValues[colOrder]
  
  metricValues <- lapply(metricValues, function(resultMetricValues) resultMetricValues[colOrder])
  if(metric != "Sample C-index") classedMetricValues <- lapply(classedMetricValues, function(resultmetricValues) resultmetricValues[colOrder])
  
  plotData <- data.frame(name = factor(rep(sampleNames(results[[1]])[colOrder], length(results)), levels = sampleNames(results[[1]])[colOrder]),
                         type = factor(rep(compareFactor, sapply(metricValues, length)), levels = compareFactor[rowOrder]),
                         Metric = unlist(metricValues))
  
  if(metric != "Sample C-index") plotData <- cbind(plotData, class = rep(knownClasses, length(results)))

  originalLegends <- showLegends                       
  originalmetricColours <- metricColours
  showLegends <- FALSE
  if(is.list(metricColours))
  {
    metricColours <- unlist(metricColours)
    plotData[, "Metric"] <- unlist(classedMetricValues)
  }                         
  
  if(metric != "Sample C-index")
  {
    classData <- data.frame(Class = knownClasses)
    classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, 0.01, 0), "npc"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
  }
  
  if(!is.null(featureValues))
  {
    if(is.factor(featureValues))
    {
      featureValuesData <- data.frame(Group = featureValues)
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, ggplot2::aes(1:length(featureValues), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(name = featureName, values = groupColours) + ggplot2::geom_tile(ggplot2::aes(fill = Group)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "npc"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
    } else # Numeric data about the samples.
    {
      featureValuesData <- data.frame(measurements = featureValues, Class = 1)
      if(metric != "Sample C-index") featureValuesData[, "Class"] <- knownClasses
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, environment = environment()) +
      ggplot2::geom_point(ggplot2::aes(x = 1:length(featureValues), y = measurements, colour = Class)) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0.5), breaks = NULL, limits = c(1, length(featureValues))) +
      ggplot2::scale_y_continuous(breaks = c(min(featureValues), max(featureValues))) +
      ggplot2::labs(x = featureName, y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 1, 0), "lines"),
                                                     axis.text = ggplot2::element_text(colour = "black"),
                                                     axis.title.x = ggplot2::element_text(size = fontSizes[4]),          
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize, panel.background = ggplot2::element_blank(),
                                                     panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
                                                     )
      if(metric != "Sample C-index")
          featureValuesPlot <- featureValuesPlot + ggplot2::scale_colour_manual(values = classColours)
    }
  }
                                                   
  metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric), show.legend = TRUE) +
    ggplot2::scale_fill_manual(values = metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                   axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                   axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                   plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                   legend.position = ifelse(showLegends, "right", "none"),
                   legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)

  if(metric != "Sample C-index")
    classGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))
  metricGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
  if(!is.null(featureValues))
    featureValuesGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(featureValuesPlot))
    
  if(!is.null(featureValues))
  {
    if(metric != "Sample C-index")      
      commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]], featureValuesGrob[["widths"]])
    else
      commonWidth <- grid::unit.pmax(metricGrob[["widths"]], featureValuesGrob[["widths"]])
  } else {
    if(metric != "Sample C-index")  
      commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]])
    else
      commonWidth <- metricGrob[["widths"]]
  }
  if(metric != "Sample C-index") classGrob[["widths"]] <- commonWidth
  metricGrob[["widths"]] <- commonWidth
  if(!is.null(featureValues))
    featureValuesGrob[["widths"]] <- commonWidth
  
  if(originalLegends == TRUE)
  {
    showLegends <- TRUE
    metricColours <- originalmetricColours
    if(metric != "Sample C-index")
    {  classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
        ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
        ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
        ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, 0.01, 0), "lines"),
                                                       legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                       legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                       legend.position = ifelse(showLegends, "right", "none"),
                                                       legend.key.size = legendSize)
      classGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot)) 
    }
    
    if(!is.null(featureValues) && is.factor(featureValues))
    {
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, ggplot2::aes(1:length(featureValues), factor(1)), environment = environment()) +
        ggplot2::scale_fill_manual(name = featureName, values = groupColours) + ggplot2::geom_tile(ggplot2::aes(fill = Group)) +
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(1, length(featureValues))) +
        ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
        ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"),
                                                       legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                       legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                       legend.position = ifelse(showLegends, "right", "none"),
                                                       legend.key.size = legendSize)
      featureValuesGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(featureValuesPlot))     
    } else {featureValuesGrobUnused <- grid::grob()}

    plotData[, "Metric"] <- unlist(metricValues)
    classLegend <- NULL
    if(is.list(metricColours) && metric != "Sample C-index")
      classLegend <- paste(levels(knownClasses)[1], NULL)
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric), show.legend = TRUE) +
      ggplot2::scale_fill_manual(name = paste(classLegend, metricText, sep = ''),
                                 values = if(is.list(metricColours)) metricColours[[1]] else metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 0, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
    if(!is.null(featureValues) && is.factor(featureValues))
    {
      if(metric != "Sample C-index")         
        commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]], featureValuesGrobUnused[["widths"]])                   
      else
        commonWidth <- grid::unit.pmax(metricGrobUnused[["widths"]], featureValuesGrobUnused[["widths"]])                   
    } else {
      if(metric != "Sample C-index")    
        commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]])                   
      else
        commonWidth <- metricGrobUnused[["widths"]]
    }
    metricGrobUnused[["widths"]] <- commonWidth
    if(metric != "Sample C-index")   
      classGrobUnused[["widths"]] <- commonWidth
    if(!is.null(featureValues) && is.factor(featureValues))
      featureValuesGrobUnused[["widths"]] <- commonWidth
    if(metric != "Sample C-index") 
      classLegend <- classGrobUnused[["grobs"]][[which(sapply(classGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    if(!is.null(featureValues) && is.factor(featureValues))
      featureValuesLegend <- featureValuesGrobUnused[["grobs"]][[which(sapply(featureValuesGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    else
      featureValuesLegend <- grid::grob()
    if(showLegends == TRUE)    
      firstLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric), show.legend = TRUE) +
      ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[2], metricText),
                                 values = if(is.list(metricColours)) metricColours[[2]] else metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot)) 
    metricGrobUnused[["widths"]] <- commonWidth  
    secondLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
  }
  
  if(showLegends == TRUE)
  {
    annosHeight <- grid::unit(1 / (mapHeight + 1), "npc")
    if(!is.null(featureValues) && is.factor(featureValues))
    {
      if(metric != "Sample C-index") 
        legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]), sum(featureValuesLegend[["widths"]]))
      else
        legendWidth <- max(sum(firstLegend[["widths"]]), sum(featureValuesLegend[["widths"]]))          
    } else {
      if(metric != "Sample C-index")         
        legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]))
      else
        legendWidth <- sum(firstLegend[["widths"]])
    }
    featureValuesHeight <- grid::unit(0, "cm")
    if(!is.null(featureValues))
    {
      if(is.factor(featureValues))
        featureValuesHeight <- featureValuesLegend[["heights"]][3]
      else
        featureValuesHeight <- annosHeight
    }
    
    if(is.list(metricColours))
      legendHeight <- (grid::unit(1, "npc") - featureValuesHeight) * (mapHeight / (mapHeight + 1) / 2)
    else
      legendHeight <- (grid::unit(1, "npc") - featureValuesHeight) * (mapHeight / (mapHeight + 1))
    
    widths <- grid::unit.c(unit(1, "npc") - legendWidth, legendWidth)
    if(metric != "Sample C-index") classesHeight <-  unit(1 / (mapHeight + 1), "npc") else classesHeight <-  unit(0, "npc")
    heights <- grid::unit.c(classesHeight, featureValuesHeight, legendHeight)
    if(is.list(metricColours))
    {
      heights <- grid::unit.c(heights, legendHeight)
    } else # Greyscale legend.
    {
      heights <- grid::unit.c(heights, unit(1, "npc") - legendHeight - featureValuesHeight - grid::unit(1 / (mapHeight + 1), "npc"))
    }
  }
  else
  {
    widths <- grid::unit(1, "npc")
    heights <- grid::unit.c(grid::unit(1 / (mapHeight + 1), "npc"), featureValuesHeight, grid::unit(1, "npc") - featureValuesHeight - grid::unit(1 / (mapHeight + 1), "npc"))
  }
  
  if(metric != "Sample C-index")
    classLegend[["vp"]][["valid.just"]] <- c(0.7, 0.5)
  if(!is.null(featureValues) && is.factor(featureValues))
    featureValuesLegend[["vp"]][["valid.just"]] <- c(0.7, 0.33)

  grobTable <- gtable::gtable(widths, heights)
  if(metric != "Sample C-index")
    grobTable <- gtable::gtable_add_grob(grobTable, classGrob, 1, 1)
  else
    grobTable <- gtable::gtable_add_grob(grobTable, grid::grob(), 1, 1)
  if(!is.null(featureValues))
    grobTable <- gtable::gtable_add_grob(grobTable, featureValuesGrob, 2, 1)
  else
    grobTable <- gtable::gtable_add_grob(grobTable, grid::grob(), 2, 1)
  if(showLegends == TRUE)
  {
    if(is.list(metricColours))
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1, 4, 1)
    else
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1, 3, 1)
  } else
  {
    grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1)
  }
  if(showLegends == TRUE)
  {
    if(metric != "Sample C-index")
      grobTable <- gtable::gtable_add_grob(grobTable, classLegend, 1, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, featureValuesLegend, 2, 2)
  }
  if(showLegends == TRUE && !is.list(metricColours))
  {
    firstLegend[["vp"]][["valid.just"]] <- c(0.62, 0.5)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 3, 2)
  }
  if(showLegends == TRUE && is.list(metricColours))
  {
    firstLegend[["vp"]][["valid.just"]] <- c(0.62, 0.5)
    secondLegend[["vp"]][["valid.just"]] <- c(0.62, 0.4)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 3, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, secondLegend, 4, 2)
  }
  wholePlot <- gridExtra::arrangeGrob(grobTable, top = grid::textGrob(title, vjust = 0.5, gp = grid::gpar(fontsize = fontSizes[1])))
  grid::grid.draw(wholePlot)
  wholePlot
})

#' @rdname samplesMetricMap
#' @export
setMethod("samplesMetricMap", "matrix", 
          function(results, classes,
                   metric = c("Sample Error", "Sample Accuracy"),
                   featureValues = NULL, featureName = NULL,
                   metricColours = list(c("#3F48CC", "#6F75D8", "#9FA3E5", "#CFD1F2", "#FFFFFF"),
                                        c("#880015", "#A53F4F", "#C37F8A", "#E1BFC4", "#FFFFFF")),
                   classColours = c("#3F48CC", "#880015"), groupColours = c("darkgreen", "yellow2"),
                   fontSizes = c(24, 16, 12, 12, 12),
                   mapHeight = 4, title = "Error Comparison", showLegends = TRUE, xAxisLabel = "Sample Name", showXtickLabels = TRUE,
                   yAxisLabel = "Analysis", showYtickLabels = TRUE, legendSize = grid::unit(1, "lines"))
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")  
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")       
  if(!requireNamespace("gtable", quietly = TRUE))
    stop("The package 'gtable' could not be found. Please install it.")
       
  metric <- match.arg(metric)
  metricText <- switch(metric, `Sample Error` = "Error", `Sample Accuracy` = "Accuracy")

  if(!is.null(featureValues) && is.null(featureName))
    stop("featureValues is specified by featureNames isn't. Specify both.")
  
  nColours <- if(is.list(metricColours)) length(metricColours[[1]]) else length(metricColours)
  metricBinEnds <- seq(0, 1, 1/nColours)
  knownClasses <- classes

  sampleIDs <- colnames(results)
  characteristic <- rownames(results)
  results <- as.list(as.data.frame(t(results)))
  
  metricValues <- lapply(results, function(result)
  {
    cut(result, metricBinEnds, include.lowest = TRUE)
  })
  
  classedMetricValues <- lapply(metricValues, function(metricSet)
  {
    metricSet <- factor(paste(knownClasses, metricSet, sep = ','),
                       levels = c(t(outer(levels(knownClasses), levels(metricSet), paste, sep = ','))))
  })

  rowOrder <- hclust(dist(results, "manhattan"))[["order"]]
  meanSample <- colMeans(do.call(rbind, metricValues))

  metricMatrixSamplewise <- t(results)
  uninformativeSamples <- which(apply(metricMatrixSamplewise, 1, function(sampleMetrics) all(is.na(sampleMetrics))))
  if(length(uninformativeSamples) > 0)
      metricMatrixSamplewise[uninformativeSamples, ] <- 0
  
  if(is.null(featureValues))
  {
    if(metric != "Sample C-index") # Sort within each class.
    {
      colOrder <- unlist(lapply(levels(knownClasses), function(class)
      {
        classIndices <- which(as.character(knownClasses) == class)
        classMatrix <- metricMatrixSamplewise[classIndices, ]
        classIndices[hclust(dist(classMatrix, "manhattan"))[["order"]]]
      }))
    } else { # Sort all samples together.
      colOrder <- hclust(dist(metricMatrixSamplewise, "manhattan"))[["order"]]
    }
  }  else {
    featureValues <- featureValues[match(colnames(results), names(featureValues))]
    colOrder <- order(knownClasses, featureValues, meanSample)
  }
  
  knownClasses <- knownClasses[colOrder]
  if(!is.null(featureValues))
    featureValues <- featureValues[colOrder]
  
  metricValues <- lapply(metricValues, function(resultmetricValues) resultmetricValues[colOrder])
  classedMetricValues <- lapply(classedMetricValues, function(resultmetricValues) resultmetricValues[colOrder])
  sampleIDs <- sampleIDs[colOrder]
  plotData <- data.frame(name = factor(rep(sampleIDs, length(characteristic)), levels = sampleIDs),
                         type = factor(rep(characteristic, each = length(sampleIDs)), levels = characteristic[rowOrder]),
                         class = rep(knownClasses, length(characteristic)),
                         Metric = unlist(metricValues))

  originalLegends <- showLegends                       
  originalmetricColours <- metricColours
  showLegends <- FALSE
  if(is.list(metricColours))
  {
    metricColours <- unlist(metricColours)
    plotData[, "Metric"] <- unlist(classedMetricValues)
  }                         
  
  classData <- data.frame(Class = knownClasses)
  classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
    ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
    ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, 0.01, 0), "npc"),
                                                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                   legend.position = ifelse(showLegends, "right", "none"),
                                                   legend.key.size = legendSize)
  
  if(!is.null(featureValues))
  {
    if(is.factor(featureValues))
    {
      featureValuesData <- data.frame(Group = featureValues)
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, ggplot2::aes(1:length(featureValues), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(name = featureName, values = groupColours) + ggplot2::geom_tile(ggplot2::aes(fill = Group)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(1, length(featureValues))) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "npc"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
    } else # Numeric data about the samples.
    {
      featureValuesData <- data.frame(Class = knownClasses, measurements = featureValues)
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, environment = environment()) +
      ggplot2::geom_point(ggplot2::aes(x = 1:length(featureValues), y = measurements, colour = Class)) +
      ggplot2::scale_colour_manual(values = classColours) +
      ggplot2::scale_x_continuous(breaks = NULL, limits = c(1, length(featureValues))) +
      ggplot2::scale_y_continuous(breaks = c(min(featureValues), max(featureValues))) +
      ggplot2::labs(x = featureName, y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 1, 0), "lines"),
                                                     axis.text = ggplot2::element_text(colour = "black"),
                                                     axis.title.x = ggplot2::element_text(size = fontSizes[4]),          
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize, panel.background = ggplot2::element_blank(),
                                                     panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
                                                     )      
    }
  }
                                                   
  metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric), show.legend = TRUE) +
    ggplot2::scale_fill_manual(values = metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                   axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                   axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                   plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                   legend.position = ifelse(showLegends, "right", "none"),
                   legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)

  classGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))
  metricGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
  if(!is.null(featureValues))
    featureValuesGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(featureValuesPlot))
    
  if(!is.null(featureValues))
    commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]], featureValuesGrob[["widths"]])
  else
    commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]])
  classGrob[["widths"]] <- commonWidth
  metricGrob[["widths"]] <- commonWidth
  if(!is.null(featureValues))
    featureValuesGrob[["widths"]] <- commonWidth
  
  if(originalLegends == TRUE)
  {
    showLegends <- TRUE
    metricColours <- originalmetricColours

    classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
      ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, 0.01, 0), "lines"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
    classGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot)) 
    
    if(!is.null(featureValues) && is.factor(featureValues))
    {
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, ggplot2::aes(1:length(featureValues), factor(1)), environment = environment()) +
        ggplot2::scale_fill_manual(name = featureName, values = groupColours) + ggplot2::geom_tile(ggplot2::aes(fill = Group)) +
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, limits = c(1, length(featureValues))) +
        ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
        ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"),
                                                       legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                       legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                       legend.position = ifelse(showLegends, "right", "none"),
                                                       legend.key.size = legendSize)
      featureValuesGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(featureValuesPlot))     
    } else {featureValuesGrobUnused <- grid::grob()}
    
    plotData[, "Metric"] <- unlist(metricValues)
    classLegend <- NULL
    if(is.list(metricColours))
      classLegend <- paste(levels(knownClasses)[1], NULL)
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric), show.legend = TRUE) +
      ggplot2::scale_fill_manual(name = paste(classLegend, metricText, sep = ''),
                                 values = if(is.list(metricColours)) metricColours[[1]] else metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 0, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
    if(!is.null(featureValues) && is.factor(featureValues))
      commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]], featureValuesGrobUnused[["widths"]])                   
    else
      commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]])                   
    metricGrobUnused[["widths"]] <- commonWidth  
    classGrobUnused[["widths"]] <- commonWidth
    if(!is.null(featureValues) && is.factor(featureValues))
      featureValuesGrobUnused[["widths"]] <- commonWidth
    classLegend <- classGrobUnused[["grobs"]][[which(sapply(classGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    if(!is.null(featureValues) && is.factor(featureValues))
      featureValuesLegend <- featureValuesGrobUnused[["grobs"]][[which(sapply(featureValuesGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    else
      featureValuesLegend <- grid::grob()
    if(showLegends == TRUE)    
      firstLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric), show.legend = TRUE) +
      ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[2], metricText),
                                 values = if(is.list(metricColours)) metricColours[[2]] else metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot)) 
    metricGrobUnused[["widths"]] <- commonWidth  
    secondLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
  }
  
  if(showLegends == TRUE)
  {
    classesHeight <- grid::unit(1 / (mapHeight + 1), "npc")
    if(!is.null(featureValues) && is.factor(featureValues))
      legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]), sum(featureValuesLegend[["widths"]]))
    else
      legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]))
    featureValuesHeight <- grid::unit(0, "cm")
    if(!is.null(featureValues))
    {
      if(is.factor(featureValues))
        featureValuesHeight <- featureValuesLegend[["heights"]][3]
      else
        featureValuesHeight <- classesHeight
    }
    
    if(is.list(metricColours))
      legendHeight <- (grid::unit(1, "npc") - featureValuesHeight) * (mapHeight / (mapHeight + 1) / 2)
    else
      legendHeight <- (grid::unit(1, "npc") - featureValuesHeight) * (mapHeight / (mapHeight + 1))
    
    widths <- grid::unit.c(unit(1, "npc") - legendWidth, legendWidth)
    heights <- grid::unit.c(unit(1 / (mapHeight + 1), "npc"), featureValuesHeight, legendHeight)
    if(is.list(metricColours))
    {
      heights <- grid::unit.c(heights, legendHeight)
    } else # Greyscale legend.
    {
      heights <- grid::unit.c(heights, unit(1, "npc") - legendHeight - featureValuesHeight - grid::unit(1 / (mapHeight + 1), "npc"))
    }
  }
  else
  {
    widths <- grid::unit(1, "npc")
    heights <- grid::unit.c(grid::unit(1 / (mapHeight + 1), "npc"), featureValuesHeight, grid::unit(1, "npc") - featureValuesHeight - grid::unit(1 / (mapHeight + 1), "npc"))
  }
  
  classLegend[["vp"]][["valid.just"]] <- c(0.7, 0.5)
  if(!is.null(featureValues) && is.factor(featureValues))
    featureValuesLegend[["vp"]][["valid.just"]] <- c(0.7, 0.33)

  grobTable <- gtable::gtable(widths, heights)
  grobTable <- gtable::gtable_add_grob(grobTable, classGrob, 1, 1)
  if(!is.null(featureValues))
    grobTable <- gtable::gtable_add_grob(grobTable, featureValuesGrob, 2, 1)
  else
    grobTable <- gtable::gtable_add_grob(grobTable, grid::grob(), 2, 1)
  if(showLegends == TRUE)
  {
    if(is.list(metricColours))
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1, 4, 1)
    else
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1, 3, 1)
  } else
  {
    grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1)
  }
  if(showLegends == TRUE)
  {  
    grobTable <- gtable::gtable_add_grob(grobTable, classLegend, 1, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, featureValuesLegend, 2, 2)
  }
  if(showLegends == TRUE && !is.list(metricColours))
  {
    firstLegend[["vp"]][["valid.just"]] <- c(0.62, 0.5)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 3, 2)
  }
  if(showLegends == TRUE && is.list(metricColours))
  {
    firstLegend[["vp"]][["valid.just"]] <- c(0.62, 0.5)
    secondLegend[["vp"]][["valid.just"]] <- c(0.62, 0.4)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 3, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, secondLegend, 4, 2)
  }
  wholePlot <- gridExtra::arrangeGrob(grobTable, top = grid::textGrob(title, vjust = 0.5, gp = grid::gpar(fontsize = fontSizes[1])))
  grid::grid.draw(wholePlot)
  
  wholePlot
})