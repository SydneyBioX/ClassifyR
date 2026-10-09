#' Plot Pair-wise Overlap, Variable Importance or Selection Size Distribution of Selected Features
#' 
#' Pair-wise overlaps can be done for two types of analyses. Firstly, each
#' cross-validation iteration can be considered within a single classification.
#' This explores the feature selection stability. Secondly, the overlap may be
#' considered between different classification results. This approach compares
#' the feature selection commonality between different selection methods. Two
#' types of commonality are possible to analyse. One summary is the average
#' pair-wise overlap between all levels of the comparison factor and the other
#' summary is the pair-wise overlap of each level of the comparison factor that
#' is not the reference level against the reference level. The overlaps are
#' converted to percentages and plotted as lineplots.
#' 
#' Additionally, a heatmap of selection size frequencies can be made by
#' specifying size as the comparison to make.
#' 
#' Lastly, a plot showing the distribution of performance metric changes when
#' features are excluded from training can be made if variable importance
#' calculation was turned on during cross-validation.
#' 
#' If \code{comparison} is \code{"within"}, then the feature selection overlaps
#' are compared within a particular analysis. The result will inform how stable
#' the selections are between different iterations of cross-validation for a
#' particular analysis. Otherwise, the comparison is between different
#' cross-validation runs, and this gives an indication about how common are the
#' features being selected by different classifications.
#' 
#' Calculating all pair-wise set overlaps can be time-consuming. This stage can
#' be done on multiple CPUs by providing the relevant options to
#' \code{parallelParams}. The percentage is calculated as the intersection of
#' two sets of features divided by the union of the sets, multiplied by 100.
#' 
#' For the feature selection size mode, \code{binsList} is used to create bins
#' which include the lowest value for the first bin, and the highest value for
#' the last bin using \code{\link{cut}}.
#' 
#' @aliases selectionPlot selectionPlot,list-method
#' @param results A list of \code{\link{ClassifyResult}} objects.
#' @param comparison Default: \code{"within"}. The aspect of the experimental design to
#' compare. Can be any characteristic that all results share or either one of
#' the special values \code{"within"} to compare between all pairwise
#' iterations of cross-validation. or \code{"size"}, to draw a bar chart of the
#' frequency of selected set sizes, or \code{"importance"} to plot the variable importance
#' scores of selected variables. \code{"importance"} only usable if \code{doImportance} was
#' \code{TRUE} during cross-validation.
#' @param referenceLevel The level of the comparison factor to use as the
#' reference to compare each non-reference level to. If \code{NULL}, then each
#' level has the average pairwise overlap calculated to all other levels.
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
#' @param binsList Used only if \code{comparison} is \code{"size"}. A list with
#' elements named \code{"setSizes"} and \code{"frequencies"} Both elements are
#' mandatory. \code{"setSizes"} specifies the bin boundaries for bins of
#' interest of feature selection sizes (e.g. 0, 10, 20, 30).
#' \code{"frequencies"} specifies the bin boundaries for the relative frequency
#' percentages to plot (e.g. 0, 20, 40, 60, 80, 100).
#' @param yMax Used only if \code{comparison} is not \code{"size"}. The maximum
#' value of the percentage overlap to plot.
#' @param densityStyle Default: "box". Either \code{"violin"} for violin plot or
#' \code{"box"} for box plot. If cross-validation is not repeated, then a bar chart.
#' @param fontSizes A vector of length 4. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values. The fourth number is the font size of the
#' titles of grouped plots, if any are produced. In other words, when
#' \code{rowVariable} or \code{columnVariable} are not \code{NULL}.
#' @param title An overall title for the plot. By default, specifies whether
#' stability or commonality is shown.
#' @param yLabel Label to be used for the y-axis of overlap percentages. By
#' default, specifies whether stability or commonality is shown.
#' @param margin The margin to have around the plot.
#' @param rotate90 Logical. If \code{TRUE}, the boxplot is horizontal.
#' @param showLegend If \code{TRUE}, a legend is plotted next to the plot. If
#' FALSE, it is hidden.
#' @param parallelParams An object of class \code{\link{MulticoreParam}} or
#' \code{\link{SnowParam}}.
#' @param ... Not used by end user.
#' @return An object of class \code{ggplot}. If \code{comparison} is \code{"importance"} and
#' \code{characteristicsList} has \code{"row"} or \code{"column"}, a grob of one plot for each
#' group, which can be drawn with \code{grid::grid.draw}.
#' @author Dario Strbenac
#' @examples
#' 
#'   predicted <- DataFrame(sample = sample(10, 100, replace = TRUE),
#'                           class = rep(c("Healthy", "Cancer"), each = 50))
#'   actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
#'   allFeatures <- sapply(1:100, function(index) paste(sample(LETTERS, 3), collapse = ''))
#'   rankList <- list(allFeatures[1:100], allFeatures[c(5:1, 6:100)],
#'                    allFeatures[c(1:9, 11, 10, 12:100)], allFeatures[c(1:50, 60:51, 61:100)])
#'   result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validations"),
#'                             value = c("Melanoma", "t-test", "Random Forest", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], allFeatures, rankList,
#'                             list(rankList[[1]][1:15], rankList[[2]][1:15],
#'                                  rankList[[3]][1:10], rankList[[4]][1:10]),
#'                             list(function(oracle){}), NULL,
#'                             predicted, actual)
#'   
#'   predicted[, "class"] <- sample(predicted[, "class"])
#'   rankList <- list(allFeatures[1:100], allFeatures[c(sample(20), 21:100)],
#'                    allFeatures[c(1:9, 11, 10, 12:100)], allFeatures[c(1:50, 60:51, 61:100)])
#'   result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Melanoma", "t-test", "Diagonal LDA", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], allFeatures, rankList,
#'                             list(rankList[[1]][1:15], rankList[[2]][1:25],
#'                                  rankList[[3]][1:10], rankList[[4]][1:10]),
#'                             list(function(oracle){}), NULL,
#'                             predicted, actual)
#'   cList <- list(x = "Classifier Name", fillColour = "Classifier Name")
#'   selectionPlot(list(result1, result2), characteristicsList = cList)
#'   
#'   cList <- list(x = "Classifier Name", fillColour = "size")
#'   selectionPlot(list(result1, result2), comparison = "size",
#'                 characteristicsList = cList,
#'                 binsList = list(frequencies = seq(0, 100, 10), setSizes = seq(0, 25, 5))
#'                 )
#' @import grid               
#' @usage NULL
#' @export
setGeneric("selectionPlot", function(results, ...)
standardGeneric("selectionPlot"))

#' @rdname selectionPlot
#' @export
setMethod("selectionPlot", "ClassifyResult", function(results, ...) {
    selectionPlot(list(assay = results), ...)
})

#' @rdname selectionPlot
#' @export
setMethod("selectionPlot", "list", 
          function(results,
                   comparison = "within", referenceLevel = NULL,
                   characteristicsList = list(x = "auto"), coloursList = list(), alpha = 1, orderingList = list(), binsList = list(),
                   yMax = 100, densityStyle = c("box", "violin"), fontSizes = c(24, 16, 12, 16), title = if(comparison == "within") "Feature Selection Stability" else if(comparison == "size") "Feature Selection Size" else if(comparison == "importance") "Variable Importance" else "Feature Selection Commonality",
                   yLabel = if(is.null(referenceLevel) && !comparison %in% c("size", "importance")) "Common Features (%)" else if(comparison == "size") "Set Size" else if(comparison == "importance") tail(names(results[[1]]@importance), 1) else paste("Common Features with", referenceLevel, "(%)"),
                   margin = grid::unit(c(1, 1, 1, 1), "lines"), rotate90 = FALSE, showLegend = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  if(!requireNamespace("ggupset", quietly = TRUE))
    stop("The package 'ggupset' could not be found. Please install it.")              
  if(comparison == "within" && !is.null(referenceLevel))
    stop("'comparison' should not be \"within\" if 'referenceLevel' is not NULL.")
              
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
  
  allFeaturesList <- lapply(results, function(result)
  {
    .getFeaturesStrings(chosenFeatureNames(result))
  })

  if(!is.null(characteristicsList[['x']]))
    xValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[['x']], "value"])
  if(!is.null(characteristicsList[["row"]]))
    rowValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["row"]], "value"])
  if(!is.null(characteristicsList[["column"]]))
    columnValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["column"]], "value"])
  if(!comparison %in% c("within", "importance"))
    referenceVar <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == comparison, "value"])
  
  allCharacteristics <- c(unlist(characteristicsList), comparison)
  allCharacteristics <- setdiff(allCharacteristics, c("within", "size", "importance"))

  if(!is.null(referenceLevel) && !(referenceLevel %in% referenceVar))
    stop("Reference level is neither a level of the comparison factor nor is it NULL.")
  
  # Fill in any missing variables needed for ggplot2 code.
  if(is.null(characteristicsList[["fillColour"]])) fillVariable <- NULL else fillVariable <- rlang::sym(characteristicsList[["fillColour"]])
  if(is.null(characteristicsList[["lineColour"]])) lineVariable <- NULL else lineVariable <- rlang::sym(characteristicsList[["lineColour"]])
  if(is.null(characteristicsList[["row"]])) rowVariable <- NULL else rowVariable <- rlang::sym(characteristicsList[["row"]])
  if(is.null(characteristicsList[["column"]])) columnVariable <- NULL else columnVariable <- rlang::sym(characteristicsList[["column"]])
  
  if(comparison == "within")
  {
    plotData <- do.call(rbind, bpmapply(function(result, featuresList)
    {
      percentOverlaps <- unlist(mapply(function(features, index)
      {
        otherFeatures <- featuresList[(index + 1):length(featuresList)]
        sapply(otherFeatures, function(other)
        {
          length(intersect(features, other)) / length(union(features, other)) * 100
        })
      }, featuresList[1:(length(featuresList) - 1)], 1:(length(featuresList) - 1), SIMPLIFY = FALSE))      

      characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
      characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
      summaryTable <- data.frame(characteristicsList, overlap = percentOverlaps, check.names = FALSE)
      colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
      summaryTable
    }, results, allFeaturesList, BPPARAM = parallelParams, SIMPLIFY = FALSE))
  } else if(comparison == "size") {
    plotData <- do.call(rbind, mapply(function(result, featuresList)
                {
                  setSizes <- sapply(featuresList, length)
                  characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
                  characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
                  summaryTable <- data.frame(characteristicsList, size = setSizes, check.names = FALSE)
                  colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
                  summaryTable
                }, results, allFeaturesList, SIMPLIFY = FALSE))

    plotData[, "size"] <- cut(plotData[, "size"], breaks = binsList[["setSizes"]])
    selectionSizes <- as.data.frame(table(plotData)) # Spaces will be destroyed.
    colnames(selectionSizes)[1:length(characteristicsList)] <- unlist(characteristicsList)
    groupingIDs <- apply(selectionSizes[, -match(c("size", "Freq"), colnames(selectionSizes)), drop = FALSE], 1, paste, collapse = '.')
    plotData <- do.call(rbind, by(selectionSizes, groupingIDs, function(dataSubset) {
                                  dataSubset[, "Freq"] <- dataSubset[, "Freq"] / sum(dataSubset[, "Freq"]) * 100
                                  dataSubset                                                   
                                  })
                        )
    
    plotData[, "Freq"] <- cut(plotData[, "Freq"], breaks = binsList[["frequencies"]])
    if(min(binsList[["frequencies"]]) == 0)
    {
      levelsRestore <- levels(plotData[, "Freq"])
      plotData[, "Freq"] <- as.character(plotData[, "Freq"])
      plotData[is.na(plotData[, "Freq"]), "Freq"] <- 0
      plotData[, "Freq"] <- factor(plotData[, "Freq"], levels = c(0, levelsRestore))
    }
  } else if(comparison == "importance") {
    plotData <- do.call(rbind, lapply(results, function(result)
                {
                  if(length(allCharacteristics) > 0)
                  {
                  characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
                  characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
                  summaryTable <- data.frame(characteristicsList, result@importance, check.names = FALSE)
                  colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
                  
                  } else {
                      summaryTable <- result@importance
                  }
                 if("assay" %in% colnames(summaryTable))
                 {
                    summaryTable[, "feature"] <- paste(summaryTable[, "assay"], summaryTable[, "feature"])
                    summaryTable <- summaryTable[, -match("assay", colnames(summaryTable))]
                  }
                  summaryTable
                }))
    plotData <- as.data.frame(plotData, optional = TRUE)
  } else { # Commonality analysis.
    groupingFactor <- paste(if('x' %in% names(characteristicsList) && characteristicsList[['x']] != comparison) xValues,
                            if("row" %in% names(characteristicsList) && characteristicsList[["row"]] != comparison) rowValues,
                            if("column" %in% names(characteristicsList) && characteristicsList[["column"]] != comparison) columnValues, sep = " ")
    if(length(groupingFactor) == 0) groupingFactor <- rep("None", length(results))
    compareIndices <- split(1:length(results), groupingFactor)

    plotData <- do.call(rbind, bplapply(compareIndices, function(indicesSet)
    {
      if(is.null(referenceLevel))
      {
        indiciesCombinations <- lapply(1:length(indicesSet),
                                       function(index) c(indicesSet[index], indicesSet[-index]))
      } else { # Compare each factor level other than the reference level to the reference level.
        indiciesCombinations <- list(c(indicesSet[match(referenceLevel, referenceVar[indicesSet])],
                                       indicesSet[setdiff(1:length(indicesSet), match(referenceLevel, referenceVar[indicesSet]))]))
      }

      do.call(rbind, unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        featuresList <- allFeaturesList[[indiciesCombination[1]]]
        otherResultIndices <- indiciesCombination[-1]
        otherResults <- results[otherResultIndices]
        
        overlapToOther <- lapply(otherResultIndices, function(otherIndex) # Other data sets to compare to.
        {
          unlist(lapply(featuresList, function(features) # List of features of a data set.
          {
            otherFeaturesList <- allFeaturesList[[otherIndex]]
            sapply(otherFeaturesList, function(otherFeatures) # List of features of another data set.
            {
              length(intersect(features, otherFeatures)) / length(union(features, otherFeatures)) * 100
            })
          }))
        })
        
        overlapToOther <- unlist(overlapToOther) # Convert all overlaps to a vector.
        if(is.null(referenceLevel))
        {
          characteristicsOrder <- match(allCharacteristics, aDataset@characteristics[["characteristic"]])
          characteristicsList <- as.list(aDataset@characteristics[["value"]])[characteristicsOrder]
          summaryTable <- data.frame(characteristicsList, check.names = FALSE)
          colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
        } else { # Each other level has been compared to the reference level of the factor.
          otherSelections <- sapply(otherResults, length)
          
          summaryTable <- do.call(rbind, mapply(function(otherResult, otherIndex)
          {
            selectTimes <- length(featuresList) * length(allFeaturesList[[otherIndex]]) # Number of overlaps calculated.
            do.call(cbind, lapply(allCharacteristics, function(characteristic)
              rep(otherResult@characteristics[otherResult@characteristics[, "characteristic"] == characteristic, "value"], selectTimes)))
          }, otherResults, otherResultIndices, SIMPLIFY = FALSE))
          colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
        }      
        
        data.frame(summaryTable, overlap = overlapToOther, check.names = FALSE)
      })))
    }, BPPARAM = parallelParams))
  }
  rownames(plotData) <- NULL # Easier for viewing during maintenance.
  if(!"fillColours" %in% names(coloursList) && "fillColour" %in% names(characteristicsList))
  {
    if(characteristicsList[["fillColour"]] == "size")
    { # Automatically grey for zero.
      colours <- character()
      if(any(plotData[, "Freq"] == '0'))
        colours <- c("grey", scales::hue_pal()(length(unique(plotData[, "Freq"])) - 1))
      else
        colours <- scales::hue_pal()(length(unique(plotData[, "Freq"])))
      coloursList[["fillColours"]] <- colours
    } else {
      coloursList[["fillColours"]] <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["fillColour"]]])))
    }
  }
  if(!"lineColours" %in% names(coloursList) && "lineColour" %in% names(characteristicsList))
    coloursList[["lineColours"]] <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["lineColour"]]])))
  
  if(comparison == "importance")
  {
      xLabel <- "Feature"
      xColumn <- "feature"
  } else {
      xLabel <- characteristicsList[['x']]
      xColumn <- xLabel
  }
  xData <- plotData[, xColumn]
  
  if(rotate90 == TRUE) plotData[, xColumn] <- factor(plotData[, xColumn], levels = rev(levels(factor(plotData[, xColumn]))))
  if(length(orderingList) > 0) plotData <- .addUserLevels(plotData, orderingList, "overlap")
  
  if(!comparison %in% c("size", "importance"))
  {
    characteristicsList <- lapply(characteristicsList, rlang::sym)
    legendPosition <- ifelse(showLegend == TRUE, "right", "none")
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = !!characteristicsList[['x']], y = overlap, fill = !!fillVariable, colour = !!lineVariable)) +
                            ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) +
                            ggplot2::coord_cartesian(ylim = c(0, yMax)) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
                            ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin)
    if(max(table(xData)) == 1) selectionPlot <- selectionPlot + ggplot2::geom_bar(stat = "identity", alpha = alpha) else selectionPlot <- selectionPlot + densityStyle(alpha = alpha)
  } else if(comparison == "importance") {
    changeName <- tail(colnames(plotData), 1)
    performanceName <- gsub("Change in ", '', changeName)
    better <- .ClassifyRenvir[["performanceInfoTable"]][.ClassifyRenvir[["performanceInfoTable"]][, "type"] == performanceName, "better"]

    if(any(c("row", "column") %in% names(characteristicsList)))
    {
      grouping <- unique(plotData[, c(characteristicsList[["row"]], characteristicsList[["column"]]), drop = FALSE])
      plotDataList <- list()
      for(groupIndex in 1:nrow(grouping))
      {
        inGroup <- rep(TRUE, nrow(plotData)) # Rows with every characteristic of the group.
        for(characteristic in colnames(grouping))
          inGroup <- inGroup & plotData[, characteristic] == grouping[groupIndex, characteristic]
        plotDataGroup <- plotData[inGroup, ]
        featureCounts <- table(plotDataGroup[, "feature"])
        keepFeatures <- names(featureCounts)[featureCounts >= 3]
        plotDataGroup <- plotDataGroup[plotDataGroup[, "feature"] %in% keepFeatures, ]
        featuresRanked <- sort(by(plotDataGroup[, ncol(plotDataGroup)], plotDataGroup[, "feature"], median), decreasing = ifelse(better == "lower", TRUE, FALSE))
        featuresTop <- names(featuresRanked)[1:min(length(featuresRanked), 10)]
        plotDataGroup <- plotDataGroup[plotDataGroup[, "feature"] %in% featuresTop, ]
        plotDataGroup[, "feature"] <- factor(plotDataGroup[, "feature"], levels = featuresTop)
        plotDataList <- append(plotDataList, list(plotDataGroup))
      }
      
      characteristicsListSym <- lapply(characteristicsList, rlang::sym)
      plotList <- lapply(plotDataList, function(plotDataGroup)
      {
        aPlot <- ggplot2::ggplot(plotDataGroup, ggplot2::aes(x = feature, y = !!rlang::sym(changeName), fill = !!fillVariable, colour = !!lineVariable)) + ggplot2::labs(x = NULL, y = NULL) +
                 ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) +
                 ggplot2::theme(axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        aPlot <- aPlot + ggplot2::geom_hline(yintercept = 0, colour = "red", linetype = "dashed") + densityStyle(alpha = alpha)
        if("row" %in% names(characteristicsList))
          aPlot <- aPlot + ggplot2::facet_grid(rows = ggplot2::vars(!!rowVariable))
        if("column" %in% names(characteristicsList))
          aPlot <- aPlot + ggplot2::facet_grid(cols = ggplot2::vars(!!columnVariable))
        aPlot
      })
      nRows <- ifelse("row" %in% names(characteristicsList), length(unique(plotData[, characteristicsList[["row"]]])), 1)
      nColumns <- ifelse("column" %in% names(characteristicsList), length(unique(plotData[, characteristicsList[["column"]]])), 1)
      g <- gridExtra::arrangeGrob(grobs = plotList, layout_matrix = matrix(seq_along(plotList), nrow = nRows, ncol = nColumns, byrow = TRUE),
                                  top = textGrob("Variable Importance", gp = grid::gpar(vjust = 0.6)),
                                  left = textGrob(yLabel, rot = 90),
                                  bottom = textGrob("Feature"))
      
      return(g)
    } else {
    legendPosition <- ifelse(showLegend == TRUE, "right", "none")
    featureCounts <- table(plotData[, "feature"])
    keepFeatures <- names(featureCounts)[featureCounts >= 2]
    plotData <- plotData[plotData[, "feature"] %in% keepFeatures, ]
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = feature, y = !!rlang::sym(changeName), fill = !!fillVariable, colour = !!lineVariable)) +
      ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) +
      ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
      ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin) + densityStyle(alpha = alpha)
    }
  } else {
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = !!rlang::sym(characteristicsList[['x']]), y = size)) +
                     ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) +
                     ggplot2::geom_tile(ggplot2::aes(fill = Freq), alpha = alpha) + ggplot2::ggtitle(title) + ggplot2::labs(x = xLabel, y = yLabel) + ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5)) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Frequency (%)"))
  }
  
  if(length(coloursList[["fillColours"]]) > 0)
    selectionPlot <- selectionPlot + ggplot2::scale_fill_manual(values = coloursList[["fillColours"]])
  if(length(coloursList[["lineColours"]]) > 0)
    selectionPlot <- selectionPlot + ggplot2::scale_colour_manual(values = coloursList[["lineColours"]])
  
  if(rotate90 == TRUE)
    selectionPlot <- selectionPlot + ggplot2::coord_flip(ylim = c(0, yMax))
  
  if(any(c("row", "column") %in% names(characteristicsList)))
    selectionPlot <- selectionPlot + ggplot2::facet_grid(ggplot2::vars(!!characteristicsList[["row"]]), ggplot2::vars(!!characteristicsList[["column"]])) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))

  # Multivariate characteristic so plot upset.
  if(any(grepl(", ", plotData[, as.character(characteristicsList[['x']])])))
      selectionPlot <- selectionPlot + ggupset::axis_combmatrix(sep = ", ") + ggupset::theme_combmatrix(combmatrix.panel.line.size = 0, combmatrix.label.text = ggplot2::element_text(colour = "black"))
  
  selectionPlot
})
