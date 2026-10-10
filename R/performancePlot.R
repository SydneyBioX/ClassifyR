#' Plot Performance Measures for Various Classifications
#' 
#' Draws a graphical summary of a particular performance measure for a list of
#' classifications
#' 
#' If there are multiple values for a performance measure in a single result
#' object, they are plotted as a box plot or, if \code{densityStyle} is
#' \code{"violin"}, a violin plot. If a result's value of
#' \code{characteristicsList[["x"]]} is one of \code{aggregate}, its values
#' are averaged to a single number and a bar is plotted.
#' 
#' A dashed line shows the metric's expected value for predictions made by chance: 0.5 for AUC and C-index, 0 for
#' Matthews Correlation Coefficient and, with k classes, 1/k for accuracy, precision, recall and F1 metrics and 1 - 1/k for
#' error metrics. This is the value expected from predicting classes uniformly at random. No line is drawn for
#' other metrics.
#' 
#' @aliases performancePlot performancePlot,list-method
#' @param results A list of \code{\link{ClassifyResult}} objects.
#' @param aggregate A character vector of the levels of
#' \code{characteristicsList['x']} to aggregate to a single number by taking
#' the mean. This is particularly meaningful when the cross-validation is
#' leave-k-out, when k is small.
#' @param metric Default: \code{"auto"}. The name of the
#' performance measure or "auto". If the results are classification then
#' balanced accuracy will be displayed. Otherwise, the results would be survival risk
#' predictions and then C-index will be displayed. This is one of the names printed
#' in the Performance Measures field when a \code{\link{ClassifyResult}} object is
#' printed, or if none are stored, the performance metric will be calculated automatically.
#' @param characteristicsList A named list of characteristics. Each element's
#' name must be one of \code{"x"}, \code{"row"}, \code{"column"},
#' \code{"fillColour"}, or \code{"lineColour"}. The value of each element must be a
#' characteristic name, as stored in the \code{"characteristic"} column of the
#' results' characteristics table. Only \code{"x"} is mandatory. It is
#' \code{"auto"} by default, which will identify a characteristic that has a unique
#' value for each element of \code{results}. \code{"x"} represents a characteristic which will
#' form the x-axis of the plot. \code{"row"} and  \code{"column"} each specify one characteristic
#' which will form the row facet and the column facet, respectively, of a facetted plot.
#' @param coloursList A named list of plot aspects and colours for the aspects.
#' No elements are mandatory. If specified, each list element's name must be
#' either \code{"fillColours"} or \code{"lineColours"}. If a characteristic is
#' associated to fill or line by \code{characteristicsList} but this list is
#' empty, a palette of colours will be automatically chosen.
#' @param alpha Default: 1. A number between 0 and 1 specifying the transparency level of any fill.
#' @param orderingList An optional named list. Any of the variables specified
#' to \code{characteristicsList} can be the name of an element of this list and
#' the value of the element is the order in which the factors should be
#' presented in, in case alphabetical sorting is undesirable. Special values
#' \code{"performanceAscending"} and \code{"performanceDescending"} indicate that
#' the order of levels will be computed based on the median performance value of
#' the characteristic being sorted into ascending or descending order.
#' @param yLimits The minimum and maximum value of the performance metric to
#' plot. If \code{NULL}, the limits are chosen automatically.
#' @param densityStyle Default: "box". Either \code{"violin"} for violin plot or
#' \code{"box"} for box plot. If cross-validation is not repeated, then a bar chart.
#' @param fontSizes A vector of length 4. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values. The fourth number is the font size of the
#' titles of grouped plots, if any are produced. In other words, when
#' \code{rowVariable} or \code{columnVariable} are not \code{NULL}.
#' @param title An overall title for the plot.
#' @param margin The margin to have around the plot.
#' @param rotate90 Logical. IF \code{TRUE}, the plot is horizontal.
#' @param showLegend If \code{TRUE}, a legend is plotted next to the plot. If
#' FALSE, it is hidden.
#' @param ... Not used by end user.
#' @return An object of class \code{ggplot}.
#' @author Dario Strbenac
#' @examples
#' 
#'   predicted <- DataFrame(sample = sample(LETTERS[1:10], 80, replace = TRUE),
#'                          permutation = rep(1:2, each = 40),
#'                          class = factor(rep(c("Healthy", "Cancer"), 40)))
#'   actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
#'   result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Example", "t-test", "Differential Expression", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], paste("Gene", 1:100), list(paste("Gene", 1:100), paste("Gene", c(10:1, 11:100)), paste("Gene", 1:100), paste("Gene", 1:100)),
#'                             list(paste("Gene", 1:3), paste("Gene", c(2, 5, 6)), paste("Gene", 1:4), paste("Gene", 5:8)),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   result1 <- calcCVperformance(result1, "Macro F1")
#' 
#'   predicted <- DataFrame(sample = sample(LETTERS[1:10], 80, replace = TRUE),
#'                           permutation = rep(1:2, each = 40),
#'                           class = factor(rep(c("Healthy", "Cancer"), 40)))
#'                                
#'   result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Example", "Bartlett Test", "Differential Variability", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], paste("Gene", 1:100), list(paste("Gene", 1:100), paste("Gene", c(10:1, 11:100)), paste("Gene", 1:100), paste("Gene", 1:100)),
#'                             list(c(1:3), c(4:6), c(1, 6, 7, 9), c(5:8)),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   result2 <- calcCVperformance(result2, "Macro F1")
#'   
#'   performancePlot(list(result1, result2), metric = "Macro F1",
#'                   title = "Comparison")
#' 
#' @importFrom rlang sym
#' @importFrom ggupset axis_combmatrix theme_combmatrix
#' @rdname performancePlot
#' @usage NULL
#' @export
setGeneric("performancePlot", function(results, ...) standardGeneric("performancePlot"))

#' @rdname performancePlot
#' @export
setMethod("performancePlot", "ClassifyResult", function(results, ...) {
    performancePlot(list(assay = results), ...)
})

#' @rdname performancePlot
#' @export
setMethod("performancePlot", "list",
          function(results, metric = "auto",
                   characteristicsList = list(x = "auto"), aggregate = character(), coloursList = list(), alpha = 1, orderingList = list(),
                   densityStyle = c("box", "violin"), yLimits = NULL, fontSizes = c(24, 16, 12, 12), title = NULL,
                   margin = grid::unit(c(1, 1, 1, 1), "lines"), rotate90 = FALSE, showLegend = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  if(!requireNamespace("ggupset", quietly = TRUE))
    stop("The package 'ggupset' could not be found. Please install it.")              
              
  densityStyle <- match.arg(densityStyle)
  densityStyle <- ifelse(densityStyle == "box", ggplot2::geom_boxplot, ggplot2::geom_violin)
  if(characteristicsList[["x"]] == "auto")
  {
    characteristicsCounts <- table(unlist(lapply(results, function(result) result@characteristics[["characteristic"]])))
    if(max(characteristicsCounts) == length(results))
    {
      validCharacteristics <- names(characteristicsCounts)[characteristicsCounts == max(characteristicsCounts)]
      allCharacteristics <- do.call(rbind, lapply(results, function(result) result@characteristics))
      valuesPerCharacteristic <- by(allCharacteristics, allCharacteristics[, "characteristic"], function(characteristicValues) length(unique(characteristicValues[, "value"])))
      characteristicsList[["x"]] <- names(valuesPerCharacteristic)[which.max(valuesPerCharacteristic)]
    } else {
      stop("No characteristic is present for all results but must be.")
    }
  }
  if(metric == "auto")
      metric <- ifelse("risk" %in% colnames(results[[1]]@predictions), "C-index", "Balanced Accuracy")
            
  metrics <- unlist(lapply(results, function(result)
    if(!is.null(result@performance)) names(result@performance)))
  namesCounts <- table(metrics)
  commonNames <- names(namesCounts)[namesCounts == length(results)]
  if(!metric %in% commonNames)
  {
    warning(paste(metric, "not found in all elements of results. Calculating it now."))
    results <- lapply(results, function(result) calcCVperformance(result, metric))
  }
  
  baseline <- .chanceLevel(metric, results[[1]])
 
  plotData <- do.call(rbind, mapply(function(result, index)
                    {
                      if(!metric %in% names(result@performance))
                        stop(metric, " not calculated for element ", index, " of results list.")
                      row <- result@characteristics[, "characteristic"] == characteristicsList[["x"]] 
                      if(any(row) && result@characteristics[row, "value"] %in% aggregate)
                        performance <- mean(result@performance[[metric]])
                      else
                        performance <- result@performance[[metric]]
                      rows <- match(unlist(characteristicsList), result@characteristics[, "characteristic"])
                      summaryTable <- data.frame(as.list(result@characteristics[rows, "value"]), performance)
                      colnames(summaryTable) <- c(characteristicsList, metric)
                      summaryTable
                    }, results, 1:length(results), SIMPLIFY = FALSE))
  
  plotData <- plotData[, !duplicated(colnames(plotData))]
  if(length(orderingList) > 0) plotData <- .addUserLevels(plotData, orderingList, metric)

  # Fill in any missing variables needed for ggplot2 code.
  if("fillColour" %in% names(characteristicsList))
    if(!"fillColours" %in% names(coloursList)) coloursList[["fillColours"]] <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["fillColour"]]])))
  if("lineColour" %in% names(characteristicsList))
    if(!"lineColours" %in% names(coloursList)) coloursList[["lineColours"]] <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["lineColour"]]])))
  if(is.null(characteristicsList[["fillColour"]])) fillVariable <- NULL else fillVariable <- rlang::sym(characteristicsList[["fillColour"]])
  if(is.null(characteristicsList[["lineColour"]])) lineVariable <- NULL else lineVariable <- rlang::sym(characteristicsList[["lineColour"]])
  if(is.null(characteristicsList[["row"]])) rowVariable <- NULL else rowVariable <- rlang::sym(characteristicsList[["row"]])
  if(is.null(characteristicsList[["column"]])) columnVariable <- NULL else columnVariable <- rlang::sym(characteristicsList[["column"]])
  
  allCharacteristics <- unlist(characteristicsList)
  xLabel <- allCharacteristics['x']
  if(rotate90 == TRUE) plotData[, xLabel] <- factor(plotData[, xLabel], levels = rev(levels(factor(plotData[, xLabel]))))

  legendPosition <- ifelse(showLegend == TRUE, "right", "none")
  characteristicsList <- lapply(characteristicsList, rlang::sym)

  performancePlot <- ggplot2::ggplot() + ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA))
  if(!is.na(baseline)) performancePlot <- performancePlot + ggplot2::geom_hline(yintercept = baseline, linetype = 2)

  if(!is.null(yLimits) && rotate90 == FALSE) performancePlot <- performancePlot + ggplot2::coord_cartesian(ylim = yLimits)
  if("fillColour" %in% names(characteristicsList))
    performancePlot <- performancePlot + ggplot2::scale_fill_manual(values = coloursList[["fillColours"]])
  if("lineColour" %in% names(characteristicsList))
    performancePlot <- performancePlot + ggplot2::scale_colour_manual(values = coloursList[["lineColours"]])

  analysisGrouped <- split(plotData, plotData[, allCharacteristics])
  analysisGroupSizes <- sapply(analysisGrouped, nrow)
  if(any(analysisGroupSizes > 1))
  {
    multiPlotData <- do.call(rbind, analysisGrouped[analysisGroupSizes > 1])
    performancePlot <- performancePlot + densityStyle(data = multiPlotData, ggplot2::aes(x = !!characteristicsList[['x']], y = !!(rlang::sym(metric)), fill = !!fillVariable, colour = !!lineVariable), alpha = alpha)
  }
  if(any(analysisGroupSizes == 1))
  {
    singlePlotData <- do.call(rbind, analysisGrouped[analysisGroupSizes == 1])
    performancePlot <- performancePlot + ggplot2::geom_bar(data = singlePlotData, stat = "identity", ggplot2::aes(x = !!characteristicsList[['x']], y = !!(rlang::sym(metric)), fill = !!fillVariable, colour = !!lineVariable), alpha = alpha)
  }
  
  if(rotate90 == TRUE) performancePlot <- performancePlot + ggplot2::coord_flip(ylim = yLimits)
  
  performancePlot <- performancePlot + ggplot2::facet_grid(ggplot2::vars(!!rowVariable), ggplot2::vars(!!columnVariable)) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))
  performancePlot <- performancePlot + ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin)

  # Multivariate characteristic so plot upset.
  if(any(grepl(", ", plotData[, as.character(characteristicsList[['x']])])))
      performancePlot <- performancePlot + ggupset::axis_combmatrix(sep = ", ") + ggupset::theme_combmatrix(combmatrix.panel.line.size = 0, combmatrix.label.text = ggplot2::element_text(colour = "black"))
  
  performancePlot
})

# Expected value of a metric for predictions made by chance. For classes, predictions are uniformly random.
.chanceLevel <- function(metric, result)
{
  if(metric %in% c("AUC", "C-index", "Sample C-index")) return(0.5)
  if(metric == "Matthews Correlation Coefficient") return(0)
  outcome <- actualOutcome(result)
  if(!is.factor(outcome)) return(NA)
  nClasses <- length(levels(outcome))
  if(metric %in% c("Accuracy", "Balanced Accuracy", "Sample Accuracy", "Micro Precision", "Micro Recall", "Micro F1",
                   "Macro Precision", "Macro Recall", "Macro F1")) return(1 / nClasses)
  if(metric %in% c("Error", "Balanced Error", "Sample Error")) return(1 - 1 / nClasses)
  NA
}
