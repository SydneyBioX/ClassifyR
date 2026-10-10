#' Plot Pair-wise Overlap of Ranked Features
#' 
#' Pair-wise overlaps can be done for two types of analyses. Firstly, each
#' cross-validation iteration can be considered within a single classification.
#' This explores the feature ranking stability. Secondly, the overlap may be
#' considered between different classification results. This approach compares
#' the feature ranking commonality between different results. Two types of
#' commonality are possible to analyse. One summary is the average pair-wise
#' overlap between all possible pairs of results. The second kind of summary is
#' the pair-wise overlap of each level of the comparison factor that is not the
#' reference level against the reference level. The overlaps are converted to
#' percentages and plotted as lineplots. The overlap of two rankings at a number of top features is the number of features
#' in both lists of top features, divided by the length of the shorter list (which is shorter than the number of top
#' features only if a ranking has fewer features than it).
#' 
#' If \code{comparison} is \code{"within"}, then the feature selection overlaps
#' are compared within a particular analysis. The result will inform how stable
#' the selections are between different iterations of cross-validation for a
#' particular analysis. Otherwise, the comparison is between different
#' cross-validation runs, and this gives an indication about how common are the
#' features being selected by different classifications.
#' 
#' Calculating all pair-wise set overlaps for a large cross-validation result
#' can be time-consuming.  This stage can be done on multiple CPUs by providing
#' the relevant options to \code{parallelParams}.
#' 
#' @aliases rankingPlot rankingPlot,list-method
#' @param results A list of \code{\link{ClassifyResult}} objects.
#' @param topRanked A sequence of thresholds of number of the best features to
#' use for overlapping.
#' @param comparison Default: \code{"within"}. The aspect of the experimental design to
#' compare. Can be any characteristic that all results share or special value
#' "within" to compared between all pairwise iterations of cross-validation.
#' @param referenceLevel The level of the comparison factor to use as the
#' reference to compare each non-reference level to. If \code{NULL}, then each
#' level has the average pairwise overlap calculated to all other levels.
#' @param characteristicsList A named list of characteristics. The name must be
#' one of \code{"lineColour"}, \code{"pointType"}, \code{"row"} or
#' \code{"column"}. The value of each element must be a characteristic name, as
#' stored in the \code{"characteristic"} column of the results' characteristics
#' table.
#' @param orderingList An optional named list. Any of the variables specified
#' to \code{characteristicsList} can be the name of an element of this list and
#' the value of the element is the order in which the factor should be
#' presented in.
#' @param sizesList Default: \code{lineWidth = 1, pointSize = 2,
#' legendLinesPointsSize = 1, fonts = c(24, 16, 12, 12, 12, 16)}. A list which
#' must contain elements named \code{lineWidth}, \code{pointSize},
#' \code{legendLinesPointsSize} and \code{fonts}. The first three specify the
#' size of lines and points in the graph, as well as in the plot legend.
#' \code{fonts} is a vector of length 6.  The first element is the size of the
#' title text. The second element is the size of the axes titles.  The third
#' element is the size of the axes values. The fourth element is the size of
#' the legends' titles.  The fifth element is the font size of the legend
#' labels. The sixth element is the font size of the titles of grouped plots,
#' if any are produced. Each list element must numeric.
#' @param lineColours A vector of colours for different levels of the line
#' colouring parameter, if one is specified by
#' \code{characteristicsList[["lineColour"]]}. If none are specified but,
#' \code{characteristicsList[["lineColour"]]} is, an automatically-generated
#' palette will be used.
#' @param xLabelPositions Locations where to put labels on the x-axis.
#' @param yMax The maximum value of the percentage to plot.
#' @param title An overall title for the plot.
#' @param yLabel Label to be used for the y-axis of overlap percentages. By default, it states that the percentage of
#' top features in common is shown.
#' @param margin The margin to have around the plot.
#' @param showLegend If \code{TRUE}, a legend is plotted next to the plot. If
#' FALSE, it is hidden.
#' @param parallelParams An object of class \code{\link{MulticoreParam}} or
#' \code{\link{SnowParam}}.
#' @param ... Not used by end user.
#' @return An object of class \code{ggplot}.
#' @author Dario Strbenac
#' @examples
#' 
#'   predicted <- DataFrame(sample = sample(10, 100, replace = TRUE),
#'                           permutation = rep(1:2, each = 50),
#'                           class = rep(c("Healthy", "Cancer"), each = 50))
#'   actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
#'   allFeatures <- sapply(1:100, function(index) paste(sample(LETTERS, 3), collapse = ''))
#'   rankList <- list(allFeatures[1:100], allFeatures[c(15:6, 1:5, 16:100)],
#'                    allFeatures[c(1:9, 11, 10, 12:100)], allFeatures[c(1:50, 61:100, 60:51)])
#'   result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name", "Cross-validation"),
#'                             value = c("Melanoma", "t-test", "Diagonal LDA", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], allFeatures, rankList,
#'                             list(rankList[[1]][1:15], rankList[[2]][1:15],
#'                                  rankList[[3]][1:10], rankList[[4]][1:10]),
#'                             list(function(oracle){}), NULL,
#'                             predicted, actual)
#'   
#'   predicted[, "class"] <- sample(predicted[, "class"])
#'   rankList <- list(allFeatures[1:100], allFeatures[c(sample(20), 21:100)],
#'   allFeatures[c(1:9, 11, 10, 12:100)], allFeatures[c(1:50, 60:51, 61:100)])
#'   result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validations"),
#'                             value = c("Melanoma", "t-test", "Random Forest", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], allFeatures, rankList,
#'                             list(rankList[[1]][1:15], rankList[[2]][1:15],
#'                                  rankList[[3]][1:10], rankList[[4]][1:10]),
#'                             list(function(oracle){}), NULL,
#'                             predicted, actual)
#'                             
#'   rankingPlot(list(result1, result2), characteristicsList = list(pointType = "Classifier Name"))
#' 
#' @usage NULL
#' @export
setGeneric("rankingPlot", function(results, ...)
standardGeneric("rankingPlot"))

#' @rdname rankingPlot
#' @export
setMethod("rankingPlot", "ClassifyResult", function(results, ...) {
    rankingPlot(list(assay = results), ...)
})

#' @rdname rankingPlot
#' @export
setMethod("rankingPlot", "list",
          function(results, topRanked = seq(10, 100, 10),
                   comparison = "within", referenceLevel = NULL,
                   characteristicsList = list(), orderingList = list(),
                   sizesList = list(lineWidth = 1, pointSize = 2, legendLinesPointsSize = 1,
                                    fonts = c(24, 16, 12, 12, 12, 16)),
                   lineColours = NULL, xLabelPositions = seq(10, 100, 10), yMax = 100,
                   title = if(comparison[1] == "within") "Feature Ranking Stability" else "Feature Ranking Commonality",
                   yLabel = if(is.null(referenceLevel)) "Average Top Features in Common (%)" else paste("Average Top Features in Common with", referenceLevel, "(%)"),
                   margin = grid::unit(c(1, 1, 1, 1), "lines"),
                   showLegend = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(comparison == "within" && !is.null(referenceLevel))
    stop("'comparison' should not be \"within\" if 'referenceLevel' is not NULL.")

  originalFeatures <- allFeatureNames(results[[1]])
  if(is.character(originalFeatures)) nFeatures <- length(originalFeatures) else nFeatures <- nrow(originalFeatures)
  error <- character()
  if(max(topRanked) > nFeatures)
    error <- paste("'topRanked' is as high as", max(topRanked))
  if(max(xLabelPositions) > nFeatures)
    error <- paste(error, if(nchar(error) > 0) "and", "'xLabelPositions' is as high as", max(xLabelPositions))
  if(length(error) > 0)
  {
    error <- paste(error, "but there are only", nFeatures, "features in the data set.")
    stop(error)
  }
  fontSizes <- sizesList[["fonts"]]

  if(!is.null(characteristicsList[["lineColour"]]))
    lineColourValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["lineColour"]], "value"])
  if(!is.null(characteristicsList[["pointType"]]))
    pointTypeValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["pointType"]], "value"])
  if(!is.null(characteristicsList[["row"]]))
    rowValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["row"]], "value"])
  if(!is.null(characteristicsList[["column"]]))
    columnValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["column"]], "value"])
  if(comparison != "within")
    referenceVar <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == comparison, "value"])
  
  allCharacteristics <- c(unlist(characteristicsList), comparison)
  allCharacteristics <- setdiff(allCharacteristics, "within")
  
  if(!is.null(referenceLevel) && comparison != "within" && !(referenceLevel %in% referenceVar))
    stop("Reference level ", referenceLevel, " is neither a level of the comparison factor (",
          paste(referenceVar, collapse = ", "), ") nor is it NULL.")
  
  allRankedList <- lapply(results, function(result)
  {
    .getFeaturesStrings(result@rankedFeatures)
  })
  
  # Fill in any missing variables needed for ggplot2 code.
  if(is.null(characteristicsList[["pointType"]])) pointTypeVariable <- NULL else pointTypeVariable <- rlang::sym(characteristicsList[["pointType"]])
  if(is.null(characteristicsList[["lineColour"]])) lineVariable <- NULL else lineVariable <- rlang::sym(characteristicsList[["lineColour"]])
  if(is.null(characteristicsList[["row"]])) rowVariable <- NULL else rowVariable <- rlang::sym(characteristicsList[["row"]])
  if(is.null(characteristicsList[["column"]])) columnVariable <- NULL else columnVariable <- rlang::sym(characteristicsList[["column"]])
  
  if(comparison == "within")
  {
    plotData <- do.call(rbind, bpmapply(function(result, rankedList)
    {
      averageOverlap <- rowMeans(do.call(cbind, mapply(function(features, index)
      {
        otherFeatures <- rankedList[(index + 1):length(rankedList)]
        sapply(otherFeatures, function(other)
        {
          sapply(topRanked, function(top)
          {
            .topOverlap(features, other, top)
          })
        })
      }, rankedList[1:(length(rankedList) - 1)], 1:(length(rankedList) - 1), SIMPLIFY = FALSE)))

      if(length(allCharacteristics) > 0)
      {
        characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
        characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
        summaryTable <- data.frame(characteristicsList, top = topRanked, overlap = averageOverlap, check.names = FALSE)
        colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
      } else {
        summaryTable <- data.frame(top = topRanked, overlap = averageOverlap)
      }
      summaryTable
    }, results, allRankedList, BPPARAM = parallelParams, SIMPLIFY = FALSE))
  } else { # Commonality analysis.
    groupingFactor <- paste(if("lineColour" %in% names(characteristicsList) && characteristicsList[["lineColour"]] != comparison) lineColourValues,
                            if("pointType" %in% names(characteristicsList) && characteristicsList[["pointType"]] != comparison) pointTypeValues,
                            if("row" %in% names(characteristicsList) && characteristicsList[["row"]] != comparison) rowValues,
                            if("column" %in% names(characteristicsList) && characteristicsList[["column"]] != comparison) columnValues, sep = " ")
    if(length(groupingFactor) == 0) groupingFactor <- rep("None", length(results))
    compareIndices <- split(1:length(results), groupingFactor)
    plotData <- do.call(rbind, unname(unlist(bplapply(compareIndices, function(indicesSet)
    {
      if(is.null(referenceLevel))
      {
        indiciesCombinations <- lapply(1:length(indicesSet),
                                       function(index) c(indicesSet[index], indicesSet[-index]))
      } else { # Compare each factor level other than the reference level to the reference level.
        indiciesCombinations <- list(c(indicesSet[match(referenceLevel, referenceVar[indicesSet])],
                                       indicesSet[setdiff(1:length(indicesSet), match(referenceLevel, referenceVar[indicesSet]))]))
      }

      unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        rankedList <- allRankedList[[indiciesCombination[1]]]
        otherDatasetIndices <- indiciesCombination[-1]
        otherDatasets <- results[otherDatasetIndices]
        
        overlapToOther <- do.call(cbind, unname(lapply(otherDatasetIndices, function(otherIndex) # Other data sets to compare to.
        { 
          rowMeans(do.call(cbind, lapply(rankedList, function(rankings) # List of ranked features of a data set.
          {
            otherRankedList <- allRankedList[[otherIndex]]
            sapply(otherRankedList, function(otherRanked) # List of ranked features of another data set.
            {
              sapply(topRanked, function(top)
              {
                .topOverlap(rankings, otherRanked, top)
              })          
            })
          })))
        })))

        if(is.null(referenceLevel))
        {
          overlapToOther <- rowMeans(overlapToOther)
          
          characteristicsOrder <- match(allCharacteristics, aDataset@characteristics[["characteristic"]])
          characteristicsList <- as.list(aDataset@characteristics[["value"]])[characteristicsOrder]
          summaryTable <- data.frame(characteristicsList, top = topRanked, overlap = overlapToOther, check.names = FALSE)
          colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
          summaryTable
        } else { # Each other level has been compared to the reference level of the factor.
          overlapToOther <- as.numeric(overlapToOther) # Convert matrix of overlaps to a vector.
          topRankedAll <- rep(topRanked, length.out = length(overlapToOther))
          
          if(length(allCharacteristics) > 0)
            {
            summaryTable <- do.call(rbind, lapply(otherDatasets, function(otherDataset)
            {
              characteristicsOrder <- match(allCharacteristics, otherDataset@characteristics[["characteristic"]])
              characteristicsList <- as.list(otherDataset@characteristics[["value"]])[characteristicsOrder]
            }))
            colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
            summaryTable <- summaryTable[rep(1:nrow(summaryTable), each = length(topRanked)), , drop = FALSE]
            data.frame(summaryTable, top = topRankedAll, overlap = overlapToOther,
                       check.names = FALSE)
          } else {
            data.frame(top = topRankedAll, overlap = overlapToOther)
          }
        }
      }))
    }, BPPARAM = parallelParams), recursive = FALSE)))
  }
  rownames(plotData) <- NULL # Easier for viewing during maintenance.
  if(is.null(lineColours) && "lineColour" %in% names(characteristicsList))
    lineColours <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["lineColour"]]])))
  legendPosition <- ifelse(showLegend == TRUE, "right", "none")
  
  if(length(orderingList) > 0) plotData <- .addUserLevels(plotData, orderingList, "overlap")
  if(length(characteristicsList) > 0) characteristicsList <- lapply(characteristicsList, rlang::sym)
  
  overlapPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = top, y = overlap, colour = !!lineVariable, shape = !!pointTypeVariable)) +
                          ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) +
                          ggplot2::geom_line(linewidth = sizesList[["lineWidth"]]) + ggplot2::geom_point(size = sizesList[["pointSize"]]) + ggplot2::scale_x_continuous(breaks = xLabelPositions, limits = range(xLabelPositions)) + ggplot2::coord_cartesian(ylim = c(0, yMax)) +
                          ggplot2::xlab("Top Features") + ggplot2::ylab(yLabel) + ggplot2::ggtitle(title) + ggplot2::scale_colour_manual(values = lineColours) +
                          ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.position = legendPosition, legend.title = ggplot2::element_text(size = fontSizes[4]), legend.text = ggplot2::element_text(size = fontSizes[5]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin) +
                          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = sizesList[["legendLinesPointsSize"]])),
                                          shape = ggplot2::guide_legend(override.aes = list(size = sizesList[["legendLinesPointsSize"]])))
  
  overlapPlot <- overlapPlot + ggplot2::facet_grid(ggplot2::vars(!!rowVariable), ggplot2::vars(!!columnVariable)) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[6]))
  
  overlapPlot
})

# Percentage of the top features of two rankings which are in common. The denominator is the length of the shorter
# list of top features, so that a ranking with fewer features than top is not penalised.
.topOverlap <- function(ranking, otherRanking, top)
{
  ranking <- head(ranking, top)
  otherRanking <- head(otherRanking, top)
  length(intersect(ranking, otherRanking)) / min(length(ranking), length(otherRanking)) * 100
}
