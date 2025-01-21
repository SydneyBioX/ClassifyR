# Ranking of Differentially Abundant Features
differentMeansRanking <- function(measurementsTrain, classesTrain, verbose = 3)
{
  if(!requireNamespace("genefilter", quietly = TRUE))
    stop("The package 'genefilter' could not be found. Please install it.")

  pValues <- NULL
  categorical <- sapply(measurementsTrain, class) %in% c("character", "factor")
  if(any(categorical)) {
     pValues[categorical] <- sapply(which(categorical), function(featureIndex){
                             pval <- 1
                             if(length(unique(measurementsTrain[, featureIndex])) > 1)
                               pval <- chisq.test(measurementsTrain[, featureIndex], classesTrain)$p.value
                             pval
    })
  }

  measurementsMatrix <- t(as.matrix(measurementsTrain[, !categorical, drop = FALSE]))
  if(any(!categorical))
  {
    if(length(levels(classesTrain)) == 2)
    {
      if(verbose == 3)
        message(Sys.time(), ": Ranking features based on t-statistic.")
      pValues[!categorical] <- genefilter::rowttests(measurementsMatrix, classesTrain)[, "p.value"]
    } else {
      if(verbose == 3)
        message(Sys.time(), ": Ranking features based on F-statistic.")
      pValues[!categorical]  <- genefilter::rowFtests(measurementsMatrix, classesTrain)[, "p.value"]
    }
  }
  
  order(pValues) # From smallest to largest.
}
attr(differentMeansRanking, "name") <- "differentMeansRanking"
