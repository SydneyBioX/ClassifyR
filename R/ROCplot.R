#' Plot Receiver Operating Curve Graphs for Classification Results
#' 
#' Creates one ROC plot or multiple ROC plots for a list of ClassifyResult
#' objects. One plot is created if the data set has two classes and multiple
#' plots are created if the data set has three or more classes.
#' 
#' The scores stored in the results should be higher if the sample is more
#' likely to be from the class which the score is associated with. The score
#' for each class must be in a column which has a column name equal to the
#' class name.
#' 
#' For cross-validated classification, all predictions from all iterations are
#' considered simultaneously, to calculate one curve per classification.
#' 
#' @aliases ROCplot ROCplot,list-method
#' @param results A list of \code{\link{ClassifyResult}} objects.
#' @param mode Default: \code{"merge"}. Whether to merge all predictions of all
#' iterations of cross-validation into one set or keep them separate. Keeping
#' them separate will cause separate ROC curves to be computed for each
#' iteration and confidence intervals to be drawn with the solid line being the
#' averaged ROC curve.
#' @param interval Default: 95 (percent). The percent confidence interval to
#' draw around the averaged ROC curve, if mode is \code{"average"}.
#' @param comparison Default: \code{"auto"}. The aspect of the experimental design to compare. Can be
#' any characteristic that all results share. If the data set has two classes,
#' then the slot name with factor levels to be used for colouring the lines.
#' Otherwise, it specifies the variable used for plot facetting.
#' @param lineColours Default: \code{"auto"}. A vector of colours for different levels of the
#' comparison parameter, or if there are three or more classes, the classes.
#' If \code{"auto"}, a default colour palette is automatically generated.
#' @param lineWidth A single number controlling the thickness of lines drawn.
#' @param fontSizes A vector of length 5. The first number is the size of the
#' title.  The second number is the size of the axes titles and AUC text, if it
#' is not part of the legend. The third number is the size of the axes values.
#' The fourth number is the size of the legends' titles. The fifth number is
#' the font size of the legend labels.
#' @param labelPositions Default: 0.0, 0.2, 0.4, 0.6, 0.8, 1.0. Locations where
#' to put labels on the x and y axes.
#' @param plotTitle An overall title for the plot.
#' @param legendTitle A default name is used if the value is \code{NULL}.
#' Otherwise a character name can be provided.
#' @param xLabel Label to be used for the x-axis of false positive rate.
#' @param yLabel Label to be used for the y-axis of true positive rate.
#' @param showAUC Logical. If \code{TRUE}, the AUC value of each result is
#' added to its legend text.
#' @param ... Parameters not used by the \code{ClassifyResult} method but passed to
#' the \code{list} method.
#' @return An object of class \code{ggplot} and a plot on the current graphics
#' device, if \code{plot} is \code{TRUE}.
#' @author Dario Strbenac
#' @examples
#' 
#'   predicted <- do.call(rbind, list(DataFrame(data.frame(sample = LETTERS[seq(1, 20, 2)],
#'                                Healthy = c(0.89, 0.68, 0.53, 0.76, 0.13, 0.20, 0.60, 0.25, 0.10, 0.30),
#'                                Cancer = c(0.11, 0.32, 0.47, 0.24, 0.87, 0.80, 0.40, 0.75, 0.90, 0.70),
#'                                fold = 1)),
#'                     DataFrame(sample = LETTERS[seq(2, 20, 2)],
#'                                Healthy = c(0.45, 0.56, 0.33, 0.56, 0.65, 0.33, 0.20, 0.60, 0.40, 0.80),
#'                                Cancer = c(0.55, 0.44, 0.67, 0.44, 0.35, 0.67, 0.80, 0.40, 0.60, 0.20),
#'                                fold = 2)))
#'   actual <- factor(c(rep("Healthy", 10), rep("Cancer", 10)), levels = c("Healthy", "Cancer"))
#'   result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name", "Cross-validation"),
#'                             value = c("Melanoma", "t-test", "Random Forest", "2-fold")),
#'                             LETTERS[1:20], paste("Gene", LETTERS[1:10]), list(paste("Gene", LETTERS[1:10]), paste("Gene", LETTERS[c(5:1, 6:10)])),
#'                             list(paste("Gene", LETTERS[1:3]), paste("Gene", LETTERS[1:5])),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   
#'   predicted[c(2, 6), "Healthy"] <- c(0.40, 0.60)
#'   predicted[c(2, 6), "Cancer"] <- c(0.60, 0.40)
#'   result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name", "Cross-validation"),
#'                                       value = c("Melanoma", "Bartlett Test", "Differential Variability", "2-fold")),
#'                             LETTERS[1:20], paste("Gene", LETTERS[1:10]), list(paste("Gene", LETTERS[1:10]), paste("Gene", LETTERS[c(5:1, 6:10)])),
#'                             list(paste("Gene", LETTERS[1:3]), paste("Gene", LETTERS[1:5])),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   ROCplot(list(result1, result2), plotTitle = "Cancer ROC")
#'
#' @usage NULL
#' @export
setGeneric("ROCplot", function(results, ...) standardGeneric("ROCplot"))

#' @rdname ROCplot
#' @export
setMethod("ROCplot", "ClassifyResult", function(results, ...) {
    ROCplot(list(assay = results), ...)
})

#' @rdname ROCplot
#' @export
setMethod("ROCplot", "list",
          function(results, mode = c("merge", "average"), interval = 95,
                   comparison = "auto", lineColours = "auto",
                   lineWidth = 1, fontSizes = c(24, 16, 12, 12, 12), labelPositions = seq(0.0, 1.0, 0.2),
                   plotTitle = "ROC", legendTitle = NULL, xLabel = "False Positive Rate", yLabel = "True Positive Rate", showAUC = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  mode <- match.arg(mode)
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
  resultsWithComparison <- sum(sapply(results, function(result) any(result@characteristics[, "characteristic"] == comparison)))
  if(resultsWithComparison < length(results))
    stop("Not all results have comparison characteristic ", comparison, ' but need to.')
               
  distinctClasses <- levels(actualOutcome(results[[1]]))
  numberDistinctClasses <- length(distinctClasses)
  comparisonName <- comparison
  comparisonValues <- sapply(results, function(result) result@characteristics[match(comparisonName, result@characteristics[, "characteristic"]), "value"])
  
  plotDataList <- mapply(function(result, comparisonValue)
  {
    predictions <- as.data.frame(result@predictions, optional = TRUE)
    if(mode == "average")
    {
      if("fold" %in% colnames(predictions))
      {
        if("permutation" %in% colnames(predictions))
          predictionsList <- split(predictions, paste(predictions[, "permutation"], predictions[, "fold"]))
        else # Just k folds.
          predictionsList <- split(predictions, predictions[, "fold"])
      } else if("permutation" %in% colnames(predictions))
      {
        predictionsList <- split(predictions, predictions[, "permutation"])
      }
    } else {
      predictionsList <- list(predictions)
    }

    allPRlist <- lapply(predictionsList, function(predictions)
    {
      actualClasses <- actualOutcome(result)[match(predictions[, "sample"], sampleNames(result))]
      do.call(rbind, lapply(levels(actualClasses), function(class)
      {
        rates <- .ROCrates(predictions[, class], actualClasses == class)
        rates <- data.frame(FPR = c(0, rates[["FPR"]]), TPR = c(0, rates[["TPR"]]), class = class)

        summaryTable <- data.frame(comparisonValue, rates)
        colnames(summaryTable)[1] <- comparisonName
        summaryTable
      }))
    })
  }, results, comparisonValues, SIMPLIFY = FALSE)

  if(mode == "merge") {
    plotDataList <- lapply(plotDataList, function(resultTable)
    {
      .calcArea(resultTable[[1]], distinctClasses)
    })
    plotData <- do.call(rbind, plotDataList)
  } else { # ROC curve averaging.
    # Make mean and intervals of ROC curves of each set of predictions.
    quantiles <- c((100 - interval)/100/2, 1 - (100 - interval)/100/2)
    plotData <- do.call(rbind, mapply(function(allPRtables, comparisonValue) # Process all tables for one ClassifyResult.
    {
      combinedTable <- do.call(rbind, allPRtables) # To calculate change points.
      averagedTable <- do.call(rbind, lapply(distinctClasses, function(aClass)
      {
        classTable <- subset(combinedTable, class == aClass)
        changePoints <- sort(unique(classTable[, "FPR"]))
        summaryTable <- do.call(rbind, lapply(changePoints, function(changePoint)
        {
          TPRs <- sapply(allPRtables, function(PRtable)
          {
            PRtable <- subset(PRtable, class == aClass)
            PRtable[max(which(PRtable[, "FPR"] <= changePoint)), "TPR"]
          })
          data.frame(FPR = changePoint, TPR = mean(TPRs), lower = unname(quantile(TPRs, quantiles[1])), upper = unname(quantile(TPRs, quantiles[2])), class = aClass)
        }))
        summaryTable <- rbind(data.frame(FPR = 0, TPR = 0, lower = 0, upper = 0, class = aClass), summaryTable)
        summaryTable <- data.frame(comparisonValue, summaryTable)
        colnames(summaryTable)[1] <- comparisonName
        summaryTable
      }))
      .calcArea(averagedTable, distinctClasses)
    }, plotDataList, comparisonValues, SIMPLIFY = FALSE))
    
  }
  
  if(numberDistinctClasses > 2)
    lineColour <- "class"
  else
    lineColour <- comparison
  
  if(length(lineColours) == 1 && lineColours == "auto")
      lineColours <- scales::hue_pal()(ifelse(lineColour == "class", numberDistinctClasses, max(characteristicsCounts)))
  if(is.null(legendTitle))
    legendTitle <- ifelse(lineColour == "class", "Class", comparisonName)
  
  if(showAUC == TRUE)
  {
    if(numberDistinctClasses == 2)
    {
      plotData <- subset(plotData, class == distinctClasses[2])
      plotData[, comparisonName] <- paste(plotData[, comparisonName], " (AUC ", round(plotData[, "AUC"], 2), ')', sep = '')
      plotData[, comparisonName] <- factor(plotData[, comparisonName], levels = unique(plotData[, comparisonName]))
    } else {
      plotData[, "class"] <- paste(plotData[, "class"], " (AUC ", round(plotData[, "AUC"], 2), ')', sep = '')
      plotData[, "class"] <- factor(plotData[, "class"], levels = unique(plotData[, "class"]))
    }
  }
  
  plotDataSets <- list(plotData)
  if(numberDistinctClasses > 2)
    plotDataSets <- split(plotDataSets[[1]], plotDataSets[[1]][, comparisonName])

  lineColour <- rlang::sym(lineColour)
  comparison <- rlang::sym(comparison)
  ROCplots <- lapply(plotDataSets, function(plotData)
              {
                ROCplot <- ggplot2::ggplot(plotData, ggplot2::aes(x = FPR, y = TPR, colour = !!lineColour)) +
                           ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) +
                           ggplot2::geom_line(linewidth = lineWidth) + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::labs(colour = legendTitle) + ggplot2::geom_segment(x = 0, y = 0, xend = 1, yend = 1, linewidth = lineWidth, colour = "black") + ggplot2::scale_x_continuous(breaks = labelPositions, limits = c(0, 1)) +  ggplot2::scale_y_continuous(breaks = labelPositions, limits = c(0, 1)) +
                           ggplot2::theme(axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = ggplot2::element_rect(fill = "transparent"), legend.title = ggplot2::element_text(size = fontSizes[4], hjust = 0), legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::guides(colour = ggplot2::guide_legend(title.hjust = 0.5)) + ggplot2::scale_colour_manual(values = lineColours)
                
                if(mode == "average") # Add some confidence bands.
                  ROCplot <- ROCplot + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = !!lineColour), alpha = 0.3)
                
                if(numberDistinctClasses == 2)
                  ROCplot <- ROCplot + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) + ggplot2::ggtitle(plotTitle) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5))
                
                if(length(results) == 1 && showAUC == TRUE && numberDistinctClasses == 2)
                  ROCplot <- ROCplot + ggplot2::annotate("text", x = Inf, y = 0, label = paste("AUC =", round(plotData[1, "AUC"], 2)), hjust = 1.1, size = fontSizes[2] * 5/14) + ggplot2::theme(legend.position = "none")
                
                if(length(results) > 1 && numberDistinctClasses > 2)
                  ROCplot <- ROCplot + ggplot2::facet_wrap(ggplot2::vars(!!comparison))
                
                ROCplot
              })

  if(length(ROCplots) == 1)
  {
    ROCplot <- ROCplots[[1]]
  } else {
    ROCplot <- do.call(gridExtra::arrangeGrob,
                       c(ROCplots, nrow = 1, list(top = grid::textGrob(plotTitle, gp = grid::gpar(fontsize = fontSizes[1]), vjust = 0.6),
                                                  left = grid::textGrob(yLabel, gp = grid::gpar(fontsize = fontSizes[2]), rot = 90),
                                                  bottom = grid::textGrob(xLabel, gp = grid::gpar(fontsize = fontSizes[2])))))
  }
  
  ROCplot
})
