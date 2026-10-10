# Automated Selection of Previously Selected Features
previousSelection <- function(measurementsTrain, classesTrain, classifyResult, minimumOverlapPercent = 80,
                              .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Choosing previous features.")

  previousIDs <- chosenFeatureNames(classifyResult)[[.iteration]]
  # Match on the original feature names, which prepareData keeps when it makes the column names syntactic.
  featuresInfo <- S4Vectors::mcols(measurementsTrain)
  if(is.character(previousIDs))
  {
    wantedIDs <- unique(previousIDs)
    if(!is.null(featuresInfo) && "feature" %in% colnames(featuresInfo)) featuresIDs <- featuresInfo[, "feature"] else featuresIDs <- colnames(measurementsTrain)
  } else { # A data frame describing the assay and variable name of the chosen feature.
    wantedIDs <- unique(paste(previousIDs[, "assay"], previousIDs[, "feature"], sep = '\r'))
    featuresIDs <- paste(featuresInfo[, "assay"], featuresInfo[, "feature"], sep = '\r')
  }

  # Percentage of the previously selected features which are in the current data set.
  commonFeatures <- intersect(wantedIDs, featuresIDs)
  overlapPercent <- length(commonFeatures) / length(wantedIDs) * 100
  if(overlapPercent < minimumOverlapPercent)
    warning("Only ", round(overlapPercent), "% of the previously selected features are in the current data set, ",
            "fewer than the minimum of ", minimumOverlapPercent, "%.")

  match(commonFeatures, featuresIDs) # Return indices, not identifiers.
}
attr(previousSelection, "name") <- "previousSelection"