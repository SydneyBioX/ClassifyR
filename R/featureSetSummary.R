#' Transform a Table of Feature Abundances into a Table of Feature Set
#' Abundances.
#' 
#' Represents a feature set by the mean or median feature measurement of a
#' feature set for all features belonging to a feature set.
#' 
#' This feature transformation method is unusual because the mean or median
#' feature of a feature set for one sample may be different to another sample,
#' whereas most other feature transformation methods do not result in different
#' features being compared between samples during classification.
#' 
#' @aliases featureSetSummary featureSetSummary,matrix-method
#' featureSetSummary,DataFrame-method
#' featureSetSummary,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param target If the input is a \code{\link{MultiAssayExperiment}}, this
#' specifies which data set will be transformed. Can either be an integer index or a
#' character string specifying the name of the table. Must have length 1.
#' @param location Default: The median. The type of location to summarise a set
#' of features belonging to a feature set by.
#' @param featureSets An object of type \code{\link{FeatureSetCollection}}
#' which defines the feature sets.
#' @param minimumOverlapPercent The minimum percentage of overlapping features
#' between the data set and a feature set defined in \code{featureSets} for
#' that feature set to not be discarded from the anaylsis.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return The same class of variable as the input variable \code{measurements}
#' is, with the individual features summarised to feature sets. The number of
#' samples remains unchanged, so only one dimension of \code{measurements} is
#' altered.
#' @author Dario Strbenac
#' @references Network-based biomarkers enhance classical approaches to
#' prognostic gene expression signatures, Rebecca L Barter, Sarah-Jane Schramm,
#' Graham J Mann and Yee Hwa Yang, 2014, \emph{BMC Systems Biology}, Volume 8
#' Supplement 4 Article S5,
#' \url{https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-8-S4-S5}.
#' @examples
#' 
#'   sets <- list(Adhesion = c("Gene 1", "Gene 2", "Gene 3"),
#'                `Cell Cycle` = c("Gene 8", "Gene 9", "Gene 10"))
#'   featureSets <- FeatureSetCollection(sets)
#'   
#'   # Adhesion genes have a median gene difference between classes.
#'   genesMatrix <- matrix(c(rnorm(5, 9, 0.3), rnorm(5, 7, 0.3), rnorm(5, 8, 0.3),
#'                         rnorm(5, 6, 0.3), rnorm(10, 7, 0.3), rnorm(70, 5, 0.1)),
#'                         nrow = 10)
#'   rownames(genesMatrix) <- paste("Patient", 1:10)
#'   colnames(genesMatrix) <- paste("Gene", 1:10)
#'   classes <- factor(rep(c("Poor", "Good"), each = 5)) # But not used for transformation.
#'   
#'   featureSetSummary(genesMatrix, featureSets = featureSets)
#'
#' @export
#' @usage NULL
setGeneric("featureSetSummary", function(measurements, ...)
           standardGeneric("featureSetSummary"))

#' @rdname featureSetSummary
#' @export
setMethod("featureSetSummary", "matrix", # Matrix of numeric measurements.
          function(measurements, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
{
  if(class(featureSets) != "FeatureSetCollection")
    stop("'featureSets' is not of type FeatureSetCollection but must be.")

  assayedFeatures <- colnames(measurements)
  featureSets <- featureSets@sets
  keepSets <- sapply(featureSets, function(featureSet)
    length(intersect(featureSet, assayedFeatures)) / length(featureSet) * 100 > minimumOverlapPercent)
  if(all(keepSets == FALSE))
    stop("No feature sets had an overlap of at least ", minimumOverlapPercent,
         "% with the data set's feature identifiers.")
 
  if(any(keepSets == FALSE)) # Filter out those sets without adequate identifier overlap.
  {
    if(verbose == 3)
      message("Based on ", paste(minimumOverlapPercent, "% overlap rule, reducing ", sep = ''), length(featureSets), " feature sets to ", sum(keepSets), " feature sets.")
    featureSets <- featureSets[keepSets]
  }
  
  # Reduce set representations to only those features which were assayed.
  featureSets <- lapply(featureSets, function(featureSet) intersect(featureSet, assayedFeatures))
  
  location <- match.arg(location)
  if(location == "mean")
    locationFunction <- mean
  else
    locationFunction <- median
  
    if(verbose == 3)
      message("Summarising features to feature sets.")

  # Transform measurements into one feature per set.
  t(apply(measurements, 1, function(sampleMeasurements)
  {
    sapply(featureSets, function(featureSet) locationFunction(sampleMeasurements[featureSet]))
  })) # Columns are sets, rows are samples.
})

#' @rdname featureSetSummary
#' @export
setMethod("featureSetSummary", "DataFrame", # Possibly mixed data types.
          function(measurements, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
{
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  if(class(featureSets) != "FeatureSetCollection")
    stop("'featureSets' is not of type FeatureSetCollection but must be.")

  assayedFeatures <- colnames(measurements)
  featureSets <- featureSets@sets
  keepSets <- sapply(featureSets, function(featureSet)
    length(intersect(featureSet, assayedFeatures)) / length(featureSet) * 100 > minimumOverlapPercent)
  if(all(keepSets == FALSE))
    stop("No feature sets had an overlap of at least ", minimumOverlapPercent,
         "% with the data set's feature identifiers.")
  
  if(any(keepSets == FALSE)) # Filter out those sets without adequate identifier overlap.
  {
    if(verbose == 3)
      message("Based on ", paste(minimumOverlapPercent, "% overlap rule, reducing ", sep = ''), length(featureSets), " feature sets to ", sum(keepSets), " feature sets.")
    featureSets <- featureSets[keepSets]
  }
  
  # Reduce set representations to only those features which were assayed.
  featureSets <- lapply(featureSets, function(featureSet) intersect(featureSet, assayedFeatures))
  
  location <- match.arg(location)
  if(location == "mean")
    locationFunction <- mean
  else
    locationFunction <- median
  
    if(verbose == 3)
      message("Summarising features to feature sets.")

  # Transform measurements into one feature per set.
  measurements <- as.matrix(measurements)
  measurementsCollapsed <- t(apply(measurements, 1, function(sampleMeasurements)
  {
    sapply(featureSets, function(featureSet) locationFunction(sampleMeasurements[featureSet]))
  }))

  S4Vectors::DataFrame(measurementsCollapsed, check.names = FALSE)
})

#' @rdname featureSetSummary
#' @export
setMethod("featureSetSummary", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
{
  if(is.null(target))
    stop("'target' is NULL but must specify one of the data sets in 'measurements'.")
  if(class(featureSets) != "FeatureSetCollection")
    stop("'featureSets' is not of type FeatureSetCollection but must be.")

  assayUsed <- measurements[[target]]                                   
  assayedFeatures <- rownames(assayUsed)
  featureSets <- featureSets@sets
  keepSets <- sapply(featureSets, function(featureSet)
    length(intersect(featureSet, assayedFeatures)) / length(featureSet) * 100 > minimumOverlapPercent)
  if(all(keepSets == FALSE))
    stop("No feature sets had an overlap of at least ", minimumOverlapPercent,
         "% with the data set's feature identifiers.")
  
  if(any(keepSets == FALSE)) # Filter out those sets without adequate identifier overlap.
  {
    if(verbose == 3)
      message("Based on ", paste(minimumOverlapPercent, "% overlap rule, reducing ", sep = ''), length(featureSets), " feature sets to ", sum(keepSets), " feature sets.")
    featureSets <- featureSets[keepSets]
  }
  
  # Reduce set representations to only those features which were assayed.
  featureSets <- lapply(featureSets, function(featureSet) intersect(featureSet, assayedFeatures))
  
  location <- match.arg(location)
  if(location == "mean")
    locationFunction <- mean
  else
    locationFunction <- median
  
    if(verbose == 3)
      message("Summarising features to feature sets.")

  # Transform measurements into one feature per set.
  transformed <- apply(assayUsed, 2, function(sampleMeasurements)
  {
    sapply(featureSets, function(featureSet) locationFunction(sampleMeasurements[featureSet]))
  })
  
  measurements[[target]] <- transformed
  measurements
})