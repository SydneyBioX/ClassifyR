# Ranking of Differential Distributions with Likelihood Ratio Statistic (normal distribution)
likelihoodRatioRanking <- function(measurementsTrain, classesTrain, alternative = c(location = "different", scale = "different"),
                   ..., verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Ranking features by likelihood ratio test statistic.")

  allDistribution <- getLocationsAndScales(measurementsTrain, ...)
  logLikelihoodRatios <- unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  measurementsTrain, allDistribution[[1]], allDistribution[[2]])) -
  rowSums(sapply(levels(classesTrain), function(class)
  {
    classMeasurements <- measurementsTrain[which(classesTrain == class), ]
    classDistribution <- getLocationsAndScales(classMeasurements, ...)
    
    unlist(mapply(function(featureMeasurements, scale, location)
    sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
    classMeasurements,
    switch(alternative[["location"]], same = allDistribution[[1]], different = classDistribution[[1]]),
    switch(alternative[["scale"]], same = allDistribution[[2]], different = classDistribution[[2]])))    
  }))
  
  order(logLikelihoodRatios) # From smallest to largest.
}
attr(likelihoodRatioRanking, "name") <- "likelihoodRatioRanking"
