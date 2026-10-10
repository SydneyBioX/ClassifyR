#' Convert Individual Features into Differences Between Binary Interactors
#' Based on Known Sub-networks
#' 
#' This conversion is useful for creating a meta-feature table for classifier
#' training and prediction based on sub-networks that were selected based on
#' their differential correlation between classes.
#' 
#' The pairs of features known to interact with each other are specified by
#' \code{networkSets}.
#' 
#' @aliases interactorDifferences interactorDifferences,matrix-method
#' interactorDifferences,DataFrame-method
#' interactorDifferences,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are samples, and the columns are features.
#' @param featurePairs A object of type \code{\link{Pairs}}.
#' @param absolute If TRUE, then the absolute values of the differences are
#' returned.
#' @param useFeatures If \code{measurements} is a \code{MultiAssayExperiment},
#' \code{"all"} or a two-column table of features to use. If a table, the first column must have
#' assay names and the second column must have feature names found for that assay.
#' \code{"clinical"} is also a valid assay name and refers to the clinical data table.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return An object of class \code{\link{DataFrame}} with one column for each
#' interactor pair difference and one row for each sample. Additionally,
#' \code{mcols(resultTable)} prodvides a \code{\link{DataFrame}} with a column
#' named "original" containing the name of the sub-network each meta-feature
#' belongs to.
#' @author Dario Strbenac
#' @references Dynamic modularity in protein interaction networks predicts
#' breast cancer outcome, Ian W Taylor, Rune Linding, David Warde-Farley,
#' Yongmei Liu, Catia Pesquita, Daniel Faria, Shelley Bull, Tony Pawson, Quaid
#' Morris and Jeffrey L Wrana, 2009, \emph{Nature Biotechnology}, Volume 27
#' Issue 2, \url{https://www.nature.com/articles/nbt.1522}.
#' @examples
#' 
#'   pairs <- S4Vectors::Pairs(rep(c('A', 'G'), each = 3), c('B', 'C', 'D', 'H', 'I', 'J'))
#'                            
#'   # Consistent differences for interactors of A.                                           
#'   measurements <- matrix(c(5.7, 10.1, 6.9, 7.7, 8.8, 9.1, 11.2, 6.4, 7.0, 5.5,
#'                            3.6, 7.6, 4.0, 4.4, 5.8, 6.2, 8.1, 3.7, 4.4, 2.1,
#'                            8.5, 13.0, 9.9, 10.0, 10.3, 11.9, 13.8, 9.9, 10.7, 8.5,
#'                            8.1, 10.6, 7.4, 10.7, 10.8, 11.1, 13.3, 9.7, 11.0, 9.1,
#'                            round(rnorm(60, 8, 0.3), 1)), nrow = 10)
#'                          
#'   rownames(measurements) <- paste("Patient", 1:10)
#'   colnames(measurements) <- LETTERS[1:10]
#'   
#'   interactorDifferences(measurements, pairs)
#'
#' @export
#' @usage NULL
setGeneric("interactorDifferences", function(measurements, ...)
           standardGeneric("interactorDifferences"))

#' @rdname interactorDifferences
#' @export
setMethod("interactorDifferences", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  interactorDifferences(S4Vectors::DataFrame(measurements, check.names = FALSE), ...)
})

#' @rdname interactorDifferences
#' @export
setMethod("interactorDifferences", "DataFrame", # Possibly mixed data types.
          function(measurements, featurePairs = NULL, absolute = FALSE, verbose = 3)
{
  if(is.null(featurePairs))
    stop("'featurePairs' is NULL but must be provided.")

  if(verbose == 3)
    message("Calculating differences between the specified interactors.")
            
  keep <- S4Vectors::first(featurePairs) %in% colnames(measurements) & S4Vectors::second(featurePairs) %in% colnames(measurements)
  featurePairs <- featurePairs[keep]
  interactorTable <- as(measurements[, S4Vectors::first(featurePairs)], "matrix") # Coerce to basic matrix for calculation speed.
  otherInteractorTable <- as(measurements[, S4Vectors::second(featurePairs)], "matrix")
  differences <- otherInteractorTable - interactorTable
  if(absolute == TRUE)
    differences <- abs(differences)
  differences <- S4Vectors::DataFrame(differences)
  colnames(differences) <- paste(S4Vectors::second(featurePairs), '-', S4Vectors::first(featurePairs))
  differences
})

#' @rdname interactorDifferences
#' @export
setMethod("interactorDifferences", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, useFeatures = "all", ...)
{
  measurementsDF <- prepareData(measurements, useFeatures = useFeatures)[["measurements"]]
  interactorDifferences(measurementsDF, ...)
})