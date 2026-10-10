# Ranking of Differential Distributions with Likelihood Ratio Statistic (normal distribution)
likelihoodRatioRanking <- function(measurementsTrain, classesTrain, alternative = c(location = "different", scale = "different"),
                   ..., verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Ranking features by likelihood ratio test statistic.")

  measurementsTrain <- as.matrix(measurementsTrain)
  allDistribution <- getLocationsAndScales(measurementsTrain, ...)
  # Log-likelihood of each feature (column) under a normal distribution, in closed form for all features at once.
  normalLogLikelihoods <- function(measurements, locations, scales)
    -nrow(measurements) * (log(2 * pi) / 2 + log(scales)) - colSums(sweep(measurements, 2, locations)^2) / (2 * scales^2)

  classesLogLikelihoods <- lapply(levels(classesTrain), function(class)
  {
    classMeasurements <- measurementsTrain[which(classesTrain == class), , drop = FALSE]
    if(nrow(classMeasurements) == 0) return(0)
    classDistribution <- getLocationsAndScales(classMeasurements, ...)
    normalLogLikelihoods(classMeasurements,
                         switch(alternative[["location"]], same = allDistribution[[1]], different = classDistribution[[1]]),
                         switch(alternative[["scale"]], same = allDistribution[[2]], different = classDistribution[[2]]))
  })
  logLikelihoodRatios <- normalLogLikelihoods(measurementsTrain, allDistribution[[1]], allDistribution[[2]]) -
                         Reduce(`+`, classesLogLikelihoods)
  
  order(logLikelihoodRatios) # From smallest to largest.
}
attr(likelihoodRatioRanking, "name") <- "likelihoodRatioRanking"
