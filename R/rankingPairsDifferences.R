# Ranking of Pairs of Features that are Different Between Classes
# featurePairs must be of type Pairs from S4Vectors package.
pairsDifferencesRanking <- function(measurementsTrain, classesTrain, featurePairs = NULL, verbose = 3)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")
  
  suppliedPairs <- length(featurePairs)
  keepPairs <- S4Vectors::first(featurePairs) %in% colnames(measurementsTrain) & S4Vectors::second(featurePairs) %in% colnames(measurementsTrain)
  featurePairs <- featurePairs[keepPairs]
  
  if(verbose == 3)
    message(Sys.time(), ": ", suppliedPairs, " pairs input and ", length(featurePairs), " pairs remain after filtering based on data set row names.")
  
  if(verbose == 3)
    message(Sys.time(), ": Selecting pairs of features with consistent differences.")

  oneClassTraining <- which(classesTrain == levels(classesTrain)[1])
  otherClassTraining <- which(classesTrain == levels(classesTrain)[2])
  oneClassMeasurements <- measurementsTrain[oneClassTraining, ]
  otherClassMeasurements <- measurementsTrain[otherClassTraining, ]

  numerator <- as.matrix(oneClassMeasurements[, S4Vectors::first(featurePairs)])
  denominator <- as.matrix(oneClassMeasurements[, S4Vectors::second(featurePairs)])
  oneClassDifferences <- colMeans(numerator - denominator)
  
  numerator <- as.matrix(otherClassMeasurements[, S4Vectors::first(featurePairs)])
  denominator <- as.matrix(otherClassMeasurements[, S4Vectors::second(featurePairs)])
  otherClassDifferences <- colMeans(numerator - denominator)
  
  pairsClassDifferences <- otherClassDifferences - oneClassDifferences
  
  order(abs(pairsClassDifferences), decreasing = TRUE)
}
attr(pairsDifferencesRanking, "name") <- "pairsDifferencesRanking"