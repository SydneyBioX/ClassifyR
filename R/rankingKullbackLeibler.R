# Ranking of Differential Distributions with Kullback-Leibler Distance
KullbackLeiblerRanking <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Selecting features by Kullback-Leibler divergence.")

  oneClassMeasurements <- measurementsTrain[classesTrain == levels(classesTrain)[1], ]
  otherClassMeasurements <- measurementsTrain[classesTrain == levels(classesTrain)[2], ]
  oneClassDistribution <- getLocationsAndScales(oneClassMeasurements, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassMeasurements, ...)
  locationDifference <- oneClassDistribution[[1]] - otherClassDistribution[[1]]
  divergence <- 1/2 * (locationDifference^2 / ((oneClassDistribution[[2]])^2) +
                         locationDifference^2 / ((otherClassDistribution[[2]])^2) +
                         ((oneClassDistribution[[2]])^2) / ((otherClassDistribution[[2]])^2) +
                         ((otherClassDistribution[[2]])^2) / ((oneClassDistribution[[2]])^2))

  order(divergence, decreasing = TRUE)
}
attr(KullbackLeiblerRanking, "name") <- "KullbackLeiblerRanking"