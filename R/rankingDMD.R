# Differences in Mean or Median and/or Three Kinds of Deviations.
DMDranking <- function(measurementsTrain, classesTrain, differences = c("both", "location", "scale"),
                   ..., verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Selecting features by DMD.")
  differences <- match.arg(differences)
  
  allClassesLocationsScales <- lapply(levels(classesTrain), function(class)
  {
    aClassMeasurements <- measurementsTrain[which(classesTrain == class), ]
    getLocationsAndScales(aClassMeasurements, ...)
  })
  allClassesLocations <- sapply(allClassesLocationsScales, "[[", 1)
  allClassesScales <- sapply(allClassesLocationsScales, "[[", 2)
  locationsDifferences <- apply(allClassesLocations, 1, function(locations) sum(abs(c(dist(locations)))))
  scalesDifferences <- apply(allClassesScales, 1, function(scales) sum(abs(c(dist(scales)))))
  
  divergence <- 0
  if(differences %in% c("both", "location"))
    divergence <- divergence + locationsDifferences
  if(differences %in% c("both", "scale"))
    divergence <- divergence + scalesDifferences

  order(divergence, decreasing = TRUE)
}
attr(DMDranking, "name") <- "DMDranking"