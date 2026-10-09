# Automated Selection of Previously Selected Features
previousSelection <- function(measurementsTrain, classesTrain, classifyResult, minimumOverlapPercent = 80,
                              .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Choosing previous features.")

  previousIDs <- chosenFeatureNames(classifyResult)[[.iteration]]
  featuresIDs <- colnames(measurementsTrain)
  if(is.character(previousIDs))
  {
    safeIDs <- unique(make.names(previousIDs))
  } else { # A data frame describing the assay and variable name of the chosen feature.
    oldSafeIDs <- rownames(previousIDs)
    safeIDs <- unique(gsub("clinical_", '', oldSafeIDs)) # wideFormat doesn't prefix the clinical data, unlike all assays.
  }

  # Percentage of the previously selected features which are in the current data set.
  commonFeatures <- intersect(safeIDs, featuresIDs)
  overlapPercent <- length(commonFeatures) / length(safeIDs) * 100
  if(overlapPercent < minimumOverlapPercent)
    warning("Only ", round(overlapPercent), "% of the previously selected features are in the current data set, ",
            "fewer than the minimum of ", minimumOverlapPercent, "%.")

  match(commonFeatures, featuresIDs) # Return indices, not identifiers.
}
attr(previousSelection, "name") <- "previousSelection"