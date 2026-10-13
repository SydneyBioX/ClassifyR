#' Plot Density, Scatterplot, Parallel Plot or Bar Chart for Features By Class
#' 
#' Allows the visualisation of measurements in the data set. If \code{useFeatures}
#' is of type \code{\link{Pairs}}, then a parallel plot is automatically drawn.
#' If it's a single categorical variable, then a bar chart is automatically
#' drawn.
#' 
#' 
#' @aliases plotFeatureClasses plotFeatureClasses,matrix-method
#' plotFeatureClasses,DataFrame-method
#' plotFeatureClasses,MultiAssayExperiment-method
#' @param measurements A \code{\link{matrix}}, \code{\link{DataFrame}} or a
#' \code{\link{MultiAssayExperiment}} object containing the data.  For a
#' matrix, the rows are for features and the columns are for samples.  A column
#' with name \code{"class"} must be present in the \code{DataFrame} stored in
#' the \code{colData} slot.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' or if the measurements are of class \code{DataFrame} a character vector of
#' length 1 containing the column name in \code{measurement} is also permitted.
#' Not used if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param useFeatures If \code{measurements} is a \code{matrix} or
#' \code{DataFrame}, then a vector of numeric or character indices or the
#' feature identifiers corresponding to the feature(s) to be plotted. If
#' \code{measurements} is a \code{MultiAssayExperiment}, then a
#' \code{DataFrame} of 2 columns must be specified. The first column contains
#' the names of the assays and the second contains the names of the variables,
#' thus each row unambiguously specifies a variable to be plotted.
#' @param classesColumn If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the class column in the table extracted by \code{colData(multiAssayExperiment)}
#' that contains each sample's outcome to use for prediction.
#' @param groupBy If \code{measurements} is a \code{DataFrame}, then a
#' character vector of length 1, which contains the name of a categorical
#' feature, may be specified.  If \code{measurements} is a
#' \code{MultiAssayExperiment}, then a character vector of length 2, which
#' contains the name of a data table as the first element and the name of a
#' categorical feature as the second element, may be specified.  Additionally,
#' the value \code{"clinical"} may be used to refer to the column annotation
#' stored in the \code{colData} slot of the of the \code{MultiAssayExperiment}
#' object. A density plot will have additional lines of different line types
#' for each category. A strip chart plot will have a separate strip chart
#' created for each category and the charts will be drawn in a single column on
#' the graphics device. A parallel plot and bar chart plot will similarly be
#' laid out.
#' @param groupingName A label for the grouping variable to be used in plots.
#' @param ... Unused variables by the three top-level methods passed to the
#' internal method which generates the plot(s).
#' @param whichNumericFeaturePlots If the feature is a single feature and has
#' numeric measurements, this option specifies which types of plot(s) to draw.
#' The default value is \code{"both"}, which draws a density plot and also a
#' stip chart below the density plot. Other options are \code{"density"} for
#' drawing only a density plot and \code{"stripchart"} for drawing only a strip
#' chart.
#' @param measurementLimits The minimum and maximum expression values to plot.
#' Default: \code{NULL}.  By default, the limits are automatically computed
#' from the data values.
#' @param lineWidth Numeric value that alters the line thickness for density
#' plots. Default: 1.
#' @param dotBinWidth Numeric value that alters the diameter of dots in the
#' strip chart. Default: 1.
#' @param xAxisLabel The axis label for the plot's horizontal axis. Default:
#' \code{NULL}.
#' @param yAxisLabels A character vector of length 1 or 2. If the feature's
#' measurements are numeric an \code{whichNumericFeaturePlots} has the value
#' \code{"both"}, the first value is the y-axis label for the density plot and
#' the second value is the y-axis label for the strip chart. Otherwise, if the
#' feature's measurements are numeric and only one plot is drawn, then a
#' character vector of length 1 specifies the y-axis label for that particular
#' plot. Ignored if the feature's measurements are categorical.
#' @param showXtickLabels Logical. Default: \code{TRUE}. If set to
#' \code{FALSE}, the x-axis labels are hidden.
#' @param showYtickLabels Logical. Default: \code{TRUE}. If set to
#' \code{FALSE}, the y-axis labels are hidden.
#' @param xLabelPositions Either \code{"auto"} or a vector of values. The positions of
#' labels on the x-axis.  If \code{"auto"}, the placement of labels is automatically
#' calculated.
#' @param yLabelPositions Either \code{"auto"} or a vector of values. The positions of
#' labels on the y-axis.  If \code{"auto"}, the placement of labels is automatically
#' calculated.
#' @param fontSizes A vector of length 5. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values. The fourth number is the size of the
#' legends' titles. The fifth number is the font size of the legend labels.
#' @param colours The colours to plot data of each class in. The length of this
#' vector must be as long as the distinct number of classes in the data set.
#' @param showAssayName Logical. Default: \code{TRUE}. If \code{TRUE} and the
#' data is in a \code{MultiAssayExperiment} object, the the name of the table
#' in which the feature is stored in is added to the plot title.
#' @return Plots are created on the current graphics device and a list of plot
#' objects is invisibly returned. The classes of the plot object are determined
#' based on the type of data plotted and the number of plots per feature
#' generated. If the plotted variable is discrete or if the variable is numeric
#' and one plot type was specified, the list element is an object of class
#' \code{ggplot}. Otherwise, if the variable is numeric and both the density
#' and stripchart plot types were made, the list element is an object of class
#' \code{TableGrob}.
#' 
#' Settling \code{lineWidth} and \code{dotBinWidth} to the same value doesn't
#' result in the density plot and the strip chart having elements of the same
#' size. Some manual experimentation is required to get similarly sized plot
#' elements.
#' @author Dario Strbenac
#' @examples
#' 
#'   # First 25 samples and first 5 genes are mixtures of two normals. Last 25 samples are
#'   # one normal.
#'   genesMatrix <- sapply(1:15, function(geneColumn) c(rnorm(5, 5, 1)))
#'   genesMatrix <- cbind(genesMatrix, sapply(1:10, function(geneColumn) c(rnorm(5, 15, 1))))
#'   genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) c(rnorm(5, 9, 2))))
#'   genesMatrix <- rbind(genesMatrix, sapply(1:50, function(geneColumn) rnorm(95, 9, 3)))
#'   genesMatrix <- t(genesMatrix)
#'   rownames(genesMatrix) <- paste("Sample", 1:50)
#'   colnames(genesMatrix) <- paste("Gene", 1:100)
#'   classes <- factor(rep(c("Poor", "Good"), each = 25), levels = c("Good", "Poor"))
#'   plotFeatureClasses(genesMatrix, classes, useFeatures = "Gene 4",
#'                      xAxisLabel = bquote(log[2]*'(expression)'), dotBinWidth = 0.5)
#'                      
#'                      
#'   
#'   infectionResults <- c(rep(c("No", "Yes"), c(20, 5)), rep(c("No", "Yes"), c(5, 20)))
#'   genders <- factor(rep(c("Male", "Female"), each = 10, length.out = 50))
#'   clinicalData <- DataFrame(Gender = genders, Sugar = runif(50, 4, 10),
#'                               Infection = factor(infectionResults, levels = c("No", "Yes")),
#'                             row.names = rownames(genesMatrix))
#'   plotFeatureClasses(clinicalData, classes, useFeatures = "Infection")
#'   plotFeatureClasses(clinicalData, classes, useFeatures = "Infection", groupBy = "Gender")
#'   
#'   genesMatrix <- t(genesMatrix) # MultiAssayExperiment needs features in rows.
#'   dataContainer <- MultiAssayExperiment(list(RNA = genesMatrix),
#'                                         colData = cbind(clinicalData, class = classes))
#'   targetFeatures <- DataFrame(assay = "RNA", feature = "Gene 50")                                     
#'   plotFeatureClasses(dataContainer, useFeatures = targetFeatures, classesColumn = "class",
#'                      groupBy = c("clinical", "Gender"), # Table name, feature name.
#'                      xAxisLabel = bquote(log[2]*'(expression)'), dotBinWidth = 0.5)
#' 
#' @importFrom dplyr mutate n
#' @importFrom tidyr gather
#' @usage NULL
#' @export
setGeneric("plotFeatureClasses", function(measurements, ...)
  standardGeneric("plotFeatureClasses"))

#' @rdname plotFeatureClasses
#' @export
setMethod("plotFeatureClasses", "matrix", function(measurements, ...)
{
  plotFeatureClasses(S4Vectors::DataFrame(measurements, check.names = FALSE), ...)
})

#' @rdname plotFeatureClasses
#' @export
setMethod("plotFeatureClasses", "DataFrame", function(measurements, classes, useFeatures, groupBy = NULL,
                                                      groupingName = NULL, whichNumericFeaturePlots = c("both", "density", "stripchart"),
                                                      measurementLimits = NULL, lineWidth = 1, dotBinWidth = 1,
                                                      xAxisLabel = NULL, yAxisLabels = c("Density", "Classes"),
                                                      showXtickLabels = TRUE, showYtickLabels = TRUE,
                                                      xLabelPositions = "auto", yLabelPositions = "auto",
                                                      fontSizes = c(24, 16, 12, 12, 12),
                                                      colours = c("#3F48CC", "#880015"),
                                                      showAssayName = TRUE)
{
  if(missing(useFeatures))
    stop("'useFeatures' must be specified.")

  if(is.character(groupBy))
  {
    groupingName <- groupBy
    groupBy <- measurements[, groupBy]
    levelsOrder <- levels(groupBy)
    groupBy <- list(legends = factor(groupBy, levels = levelsOrder),
                    facets = factor(paste(groupingName, "is", groupBy), levels = paste(groupingName, "is", levelsOrder)))
  }
  if(is.character(classes) || is.integer(classes))
    classes <- measurements[, classes] # Otherwise an independent factor.
  if(is.character(classes)) classes <- factor(classes)

  if(length(colours) != length(levels(classes)))
    colours <- scales::hue_pal()(length(levels(classes)))
  
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(!requireNamespace("cowplot", quietly = TRUE))
    stop("The package 'cowplot' could not be found. Please install it.")  
  if(!requireNamespace("grid", quietly = TRUE))
    stop("The package 'grid' could not be found. Please install it.")
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")
  
  plotTheme <- ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA), panel.grid.major = ggplot2::element_line(colour = "grey", linetype = "dashed"))
  whichNumericFeaturePlots <- match.arg(whichNumericFeaturePlots)
  if(xLabelPositions[1] == "auto")
    xLabelPositions <- ggplot2::waiver()
  if(yLabelPositions[1] == "auto")
    yLabelPositions <- ggplot2::waiver()
  
  # Subsetting of measurements to the features of interest.
  if(!is(useFeatures, "data.frame") && !is(useFeatures, "DataFrame"))
  {
    if(!is(useFeatures, "Pairs")) # A simple vector.
      measurements <- tryCatch(measurements[useFeatures], error = function(error) stop("Parameter 'useFeatures' not in measurements, subscript contains out-of-bounds indices.", call. = FALSE))
    else # Pairs object.
      measurements <- tryCatch(measurements[union(S4Vectors::first(useFeatures), S4Vectors::second(useFeatures))], error = function(error) stop("Parameter 'useFeatures' not in measurements, subscript contains out-of-bounds indices.", call. = FALSE))
  }
  
  if(!is(useFeatures, "Pairs"))
  {
    invisible(lapply(seq_along(measurements), function(columnIndex) # Plot a single feature at a time.
    {
      plotData <- data.frame(measurement = measurements[, columnIndex], class = classes)
      if(!is.null(groupBy))
      {
        plotData[, "legends grouping"] <- groupBy[["legends"]]
        plotData[, "facets grouping"] <- groupBy[["facets"]]
      }
      
      if(is.null(S4Vectors::mcols(measurements)))
      {
        featureText <- colnames(measurements)[columnIndex]
      } else {
        featureText <- S4Vectors::mcols(measurements)[columnIndex, "feature"]
        if(showAssayName == TRUE && !all(S4Vectors::mcols(measurements)[, "assay"] == "assay"))
        {
          featureText <- paste(featureText, paste('(', S4Vectors::mcols(measurements)[columnIndex, "assay"], ')', sep = ''))
        }
      }
      
      if(is.numeric(plotData[, "measurement"]))
      {
        if(whichNumericFeaturePlots %in% c("both", "density"))
        {
          densPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = measurement, colour = class)) + plotTheme +
            ggplot2::stat_density(ggplot2::aes(y = ggplot2::after_stat(density)), geom = "path", position = "identity", linewidth = lineWidth) +
            ggplot2::scale_colour_manual("Class", values = colours) + ggplot2::coord_cartesian(xlim = measurementLimits) +
            ggplot2::scale_x_continuous(breaks = xLabelPositions) + ggplot2::scale_y_continuous(breaks = yLabelPositions)
          
          if(!is.null(groupBy))
          {
            densPlot <- densPlot + ggplot2::aes(linetype = `legends grouping`) + ggplot2::scale_linetype_discrete(name = groupingName)
          }
        }
        
        if(whichNumericFeaturePlots %in% c("both", "stripchart"))
        {
          yLabel <- ifelse(whichNumericFeaturePlots == "both", yAxisLabels[2], yAxisLabels[1])
          stripPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = class, y = measurement)) + plotTheme +
            ggplot2::geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", ggplot2::aes(colour = class), binwidth = dotBinWidth) +
            ggplot2::scale_colour_manual("Class", values = colours) + ggplot2::xlab(yLabel) + ggplot2::ylab(xAxisLabel) + ggplot2::scale_y_continuous(limits = measurementLimits) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                           axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                           axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                           legend.text = ggplot2::element_text(size = fontSizes[5]))
          stripPlot <- stripPlot + ggplot2::coord_flip()
          
          if(!is.null(groupBy))
          {
            stripPlot <- stripPlot + ggplot2::facet_wrap(~ `facets grouping`, ncol = 1, strip.position = "left")
          }
        }
        
        if(whichNumericFeaturePlots == "both")
        {
          densPlot <- densPlot + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                                                axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                                axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                                axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                                axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                                legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                legend.text = ggplot2::element_text(size = fontSizes[5])) +
            ggplot2::labs(x = NULL, y = yAxisLabels[1])
          
          alignedPlots <- cowplot::align_plots(densPlot, stripPlot, align = 'v', axis = "lr")
          bothGraphics <- gridExtra::arrangeGrob(alignedPlots[[1]], alignedPlots[[2]], nrow = 2,
                                                 top = grid::textGrob(featureText, gp = grid::gpar(fontsize = fontSizes[1]), vjust = 1))
          bothGraphics
        } else if(whichNumericFeaturePlots == "density")
        {
          densPlot <- densPlot + ggplot2::xlab(xAxisLabel) +  ggplot2::ylab(yAxisLabels[1]) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                           axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                           legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::ggtitle(featureText)
          
          densPlot
        } else {
          stripPlot <- stripPlot + ggplot2::ggtitle(featureText)
          stripPlot
        }
      } else { # Plotting variable is a factor.
        barPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = measurement, fill = class)) + plotTheme +
          ggplot2::geom_bar() + ggplot2::xlab(xAxisLabel) + ggplot2::ylab("Number of Samples") +
          ggplot2::scale_fill_manual("Class", values = colours) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                         axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                         axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                         axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                         axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                         axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                         legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::ggtitle(featureText)
        if(!is.null(groupBy))
        {
          barPlot <- barPlot + ggplot2::facet_wrap(~ `facets grouping`, ncol = 1, strip.position = "left")
        }
        barPlot
      }
    }))
  } else { # Plot parallel plots.
    invisible(lapply(useFeatures, function(featurePair) # Plot a single feature at a time.
    {
      plotData <- data.frame(first = measurements[, S4Vectors::first(featurePair)],
                             second = measurements[, S4Vectors::second(featurePair)], class = classes)
      if(!is.null(groupBy))
        plotData[, "facets grouping"] <- groupBy[["facets"]]
      
      groupedPlotData <- tidyr::gather(dplyr::mutate(plotData, ID = 1:dplyr::n()), key, value, 1:2)
      
      pairsPlot <- ggplot2::ggplot(groupedPlotData, ggplot2::aes(key, value, group = ID, colour = class)) + plotTheme +
        ggplot2::geom_line(linewidth = lineWidth) +
        ggplot2::scale_colour_manual("Class", values = colours) + ggplot2::coord_cartesian(ylim = measurementLimits) +
        ggplot2::scale_x_discrete(expand = c(0.05, 0.05), labels = c(S4Vectors::first(featurePair), S4Vectors::second(featurePair))) +
        ggplot2::scale_y_continuous(breaks = yLabelPositions) + ggplot2::ggtitle("Pairs Plot") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                       axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                       axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                       axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                       axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                       legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::labs(x = xAxisLabel, y = yAxisLabels[1])
      
      if(!is.null(groupBy))
      {
        pairsPlot <- pairsPlot + ggplot2::facet_wrap(~ `facets grouping`, ncol = 1, strip.position = "left")
      }
      pairsPlot
    }))
  }
})

#' @rdname plotFeatureClasses
#' @export
setMethod("plotFeatureClasses", "MultiAssayExperiment",
          function(measurements, useFeatures, classesColumn, groupBy = NULL, groupingName = NULL, showAssayName = TRUE, ...)
          {
            if(missing(useFeatures))
              stop("'useFeatures' must be specified by the user.")
            if(!all(useFeatures[, 1] %in% c(names(measurements), "clinical")))
              stop("Some table names in 'useFeatures' are not assay names in 'measurements' or \"clinical\".")  
            
            assaysUseFeatures <- useFeatures[useFeatures[, 1] != "clinical", ]
            clinicalUseFeatures <- useFeatures[useFeatures[, 1] == "clinical", ]
            measurements <- measurements[assaysUseFeatures[, 2], , assaysUseFeatures[, 1]]
            classes <- MultiAssayExperiment::colData(measurements)[, classesColumn]
            
            if(!is.null(groupBy))
            {
              if(is.null(groupingName))
                groupingName <- groupBy[2]
              groupingTable <- groupBy[1]
              if(groupingTable == "clinical")
              {
                groupBy <- MultiAssayExperiment::colData(measurements)[, groupBy[2]]
              } else { # One of the omics tables.
                groupBy <- measurements[groupBy[2], , groupingTable]
                if(showAssayName == TRUE)
                  groupingName <- paste(groupingName, groupingTable)
              }
              levelsOrder <- levels(groupBy)
              groupBy <- list(legends = factor(groupBy, levels = levelsOrder),
                              facets = {groupText <- paste(groupingName, "is", groupBy)
                              factor(groupText, levels = paste(groupingName, "is", levelsOrder))}
              )
            }
            
            MultiAssayExperiment::colData(measurements) <- MultiAssayExperiment::colData(measurements)[colnames(MultiAssayExperiment::colData(measurements)) %in% clinicalUseFeatures[, 2]]
            measurements <- MultiAssayExperiment::wideFormat(measurements, colDataCols = seq_along(MultiAssayExperiment::colData(measurements)), check.names = FALSE, collapse = ':')
            measurements <- measurements[, -1, drop = FALSE] # Remove sample IDs.
            S4Vectors::mcols(measurements)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(measurements)[, "sourceName"])
            colnames(S4Vectors::mcols(measurements))[1] <- "assay"
            S4Vectors::mcols(measurements)[, "feature"] <- S4Vectors::mcols(measurements)[, "rowname"]
            missingIndices <- is.na(S4Vectors::mcols(measurements)[, "feature"])
            S4Vectors::mcols(measurements)[missingIndices, "feature"] <- colnames(measurements)[missingIndices]
            S4Vectors::mcols(measurements) <- S4Vectors::mcols(measurements)[, c("assay", "feature")]
            
            plotFeatureClasses(measurements, classes, S4Vectors::mcols(measurements), groupBy, groupingName, showAssayName = showAssayName, ...)
          })