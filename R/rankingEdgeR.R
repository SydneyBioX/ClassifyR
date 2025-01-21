# Feature Ranking Based on Likelihood Ratio Test for Count Data
edgeRranking <- function(countsTrain, classesTrain, normFactorsOptions = NULL, dispOptions = NULL, fitOptions = NULL, verbose = 3)
{
  if(verbose == 3)
    message(Sys.time(), ": Doing edgeR LRT feature ranking")
  if(!requireNamespace("edgeR", quietly = TRUE))
    stop("The package 'edgeR' could not be found. Please install it.")
  
  # DGEList stores features as rows and samples as columns.          
  countsList <- edgeR::DGEList(t(as.matrix(countsTrain)), group = classesTrain)
  paramList <- list(countsList)
  if(!is.null(normFactorsOptions))
    paramList <- append(paramList, normFactorsOptions)
  if(verbose == 3)
    message(Sys.time(), ": Calculating scaling factors.")
  countsList <- do.call(edgeR::calcNormFactors, paramList)
  paramList <- list(countsList, model.matrix(~ classesTrain))
  if(!is.null(dispOptions))
    paramList <- append(paramList, dispOptions)
  if(verbose == 3)
    message(Sys.time(), ": Estimating dispersion.")
  countsList <- do.call(edgeR::estimateDisp, paramList)
  paramList <- list(countsList, model.matrix(~ classesTrain))
  if(!is.null(fitOptions))
    paramList <- append(paramList, fitOptions)
  if(verbose == 3)
    message(Sys.time(), ": Fitting linear model.")
  fit <- do.call(edgeR::glmFit, paramList)
  test <- edgeR::glmLRT(fit, coef = 2:length(levels(classesTrain)))[["table"]]
  
  order(test[, "PValue"]) # From smallest to largest.
}
attr(edgeRranking, "name") <- "edgeRranking"