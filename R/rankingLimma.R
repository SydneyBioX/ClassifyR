# Ranking of Differentially Abundant Features (moderated F-test)
limmaRanking <- function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(!requireNamespace("limma", quietly = TRUE))
    stop("The package 'limma' could not be found. Please install it.")

  fitParams <- list(t(as.matrix(measurementsTrain)), model.matrix(~ classesTrain))
  if(!missing(...))
    fitParams <- c(fitParams, list(...))
  linearModel <- do.call(limma::lmFit, fitParams)
  linearModel <- limma::eBayes(linearModel)
  linearModel <- linearModel[, -1] # Get rid of intercept.
  
  order(linearModel[["F.p.value"]]) # From smallest to largest.
}
attr(limmaRanking, "name") <- "limmaRanking"