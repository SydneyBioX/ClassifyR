# Selection of Differential Variability with Levene Statistic
leveneRanking <- function(measurementsTrain, classesTrain, verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Calculating Levene statistic.")

  # Levene's test centred on the median, as car::leveneTest does: an F-test of the absolute deviations
  # of each value from its class median. Calculated for all features at once.
  classesTrain <- droplevels(classesTrain)
  measurementsMatrix <- as.matrix(measurementsTrain)
  classesMedians <- do.call(cbind, lapply(levels(classesTrain), function(class)
                      apply(measurementsMatrix[classesTrain == class, , drop = FALSE], 2, median)))
  deviations <- abs(measurementsMatrix - t(classesMedians)[as.integer(classesTrain), , drop = FALSE])
  pValues <- genefilter::rowFtests(t(deviations), classesTrain)[, "p.value"]
  
  order(pValues) # From smallest to largest.
}
attr(leveneRanking, "name") <- "leveneRanking"
