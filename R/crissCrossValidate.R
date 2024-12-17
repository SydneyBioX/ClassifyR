#' A function to perform pairwise cross validation
#'
#' This function has been designed to perform cross-validation and model prediction on datasets in a pairwise manner.
#'
#' @param measurements A \code{list} of either \code{\link{DataFrame}}, \code{\link{data.frame}} or \code{\link{matrix}} class measurements.
#' @param outcomes A \code{list} of vectors that respectively correspond to outcomes of the samples in \code{measurements} list.
#' @param nFeatures The number of features to be used for modelling.
#' @param selectionMethod Default: \code{"auto"}. A character keyword of the feature algorithm to be used. If \code{"auto"}, t-test (two categories) /
#' F-test (three or more categories) ranking and top \code{nFeatures} optimisation is done. Otherwise, the ranking method is per-feature Cox proportional
#' hazards p-value.
#' @param selectionOptimisation A character of "Resubstitution", "Nested CV" or "none" specifying the approach used to optimise nFeatures.
#' @param trainType Default: \code{"modelTrain"}. A keyword specifying whether a fully trained model is used to make predictions on the test
#' set or if only the feature identifiers are chosen using the training data set and a number of training-predictions are made by cross-validation
#' in the test set.
#' @param classifier Default: \code{"auto"}. A character keyword of the modelling algorithm to be used. If \code{"auto"}, then a random forest is used
#' for a classification task or Cox proportional hazards model for a survival task.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation. 
#' @param performanceType Default: \code{"auto"}. If \code{"auto"}, then balanced accuracy for classification or C-index for survival. Otherwise, any one of the
#' options described in \code{\link{calcPerformance}} may otherwise be specified.
#' @param doRandomFeatures Default: \code{FALSE}. Whether to perform random feature selection to establish a baseline performance. Either \code{FALSE} or \code{TRUE}
#' are permitted values.
#' @param verbose Default: 0. A number between 0 and 3 for the amount of progress messages to give.  A higher number will produce more messages as
#' more lower-level functions print messages.
#' @return A list with elements \code{"real"} for the matrix of pairwise performance metrics using real
#' feature selection, \code{"random"} if \code{doRandomFeatures} is \code{TRUE} for metrics of random selection and
#' \code{"params"} for a list of parameters used during the execution of this function.
#' @author Harry Robertson
#'
#' @export

crissCrossValidate <- function(measurements, outcomes, 
                               nFeatures = 20, selectionMethod = "auto",
                               selectionOptimisation = "Resubstitution",
                               trainType = c("modelTrain", "modelTest"),
                               performanceType = "auto",
                               doRandomFeatures = FALSE,
                               classifier = "auto",
                               nFolds = 5, nRepeats = 20, nCores = 1, verbose = 0)
{
  trainType <- match.arg(trainType)
  if(!is.list(measurements)) stop("'measurements' is not of type list but is of type", class(measurements))
  if(is.null(names(measurements))) stop("Each element of 'measurements' must be named by the name of the data set.")
  if(!is.list(outcomes)) stop("'outcomes' is not of type list but is of type", class(outcomes))
  isCategorical <- is.character(outcomes[[1]]) && (length(outcomes[[1]]) == 1 || length(outcomes[[1]]) == nrow(measurements[[1]])) || is.factor(outcomes[[1]])
  if(performanceType == "auto")
    if(isCategorical) performanceType <- "Balanced Accuracy" else performanceType <- "C-index"
  if(length(selectionMethod) == 1 && selectionMethod == "auto")
    if(isCategorical) selectionMethod <- "t-test" else selectionMethod <- "CoxPH"
  if(length(classifier) == 1 && classifier == "auto")
    if(isCategorical) classifier <- "randomForest" else classifier <- "CoxPH"
  
  dataCleaned <- mapply(function(measurementsOne, outcomesOne)
  {
    prepareData(measurementsOne, outcomesOne)
  }, measurements, outcomes, SIMPLIFY = FALSE)
  measurements <- lapply(dataCleaned, "[[", 1)
  outcomes <- lapply(dataCleaned, "[[", 2)
  
  # If trainType is modelTrain, then build a model on a data set and test it on every data set.
  if(trainType == "modelTrain")
  {
    if(verbose > 0) message("Using built training models on all test data sets.")
    # Build a model for each dataset.
    trainedModels <- mapply(function(measurementsOne, outcomesOne)
    {
      train(measurementsOne, outcomesOne,
            nFeatures = nFeatures,
            selectionMethod = selectionMethod, selectionOptimisation = selectionOptimisation,
            classifier = classifier, multiViewMethod = "none", verbose = verbose)
     }, measurements, outcomes, SIMPLIFY = FALSE)

    # Perform pair-wise model assessment.
    performanceAllPairs <- lapply(trainedModels, function(trainedModel)
    {
      mapply(function(testData, testOutcomes)
      {
        predictions <- predict(trainedModel, testData, verbose = verbose)
        if(performanceType != "AUC") # Then expect a pair of parallel factor vectors.
            predictions <- predictions[, na.omit(match(c("class", "risk"), colnames(predictions)))]
        # Else, keep the table of class probabilities as it is.
        calcExternalPerformance(testOutcomes, predictions, performanceType)
      }, measurements, outcomes)
    })

    realPerformance <- matrix(unlist(performanceAllPairs), ncol = length(measurements), byrow = TRUE,
                              dimnames = list(paste("Select and Train", names(measurements)), paste("Predict", names(measurements))))
    realPerformance <- round(realPerformance, 2)
  } else { # trainType is "modelTest".
    if(verbose > 0) message("Using features chonsen in training to do cross-validation in the test data sets.")
    trainedModels <- mapply(function(measurementsOne, outcomesOne)
    {
      crossValidate(measurementsOne, outcomesOne,
                    nFeatures = nFeatures,
                    selectionMethod = selectionMethod,
                    selectionOptimisation = selectionOptimisation,
                    classifier = classifier,
                    multiViewMethod = "none",
                    nFolds = nFolds,
                    nCores = nCores,
                    nRepeats = nRepeats, verbose = verbose)
     }, measurements, outcomes, SIMPLIFY = FALSE)
    
    # Make it for runTests, which allows existing results to be passed into selection process.
    crossValParams <- generateCrossValParams(nRepeats, nFolds, nCores, selectionOptimisation)

    performanceAllPairs <- lapply(trainedModels, function(trainedModel)
    {
      mapply(function(measurementsOne, outcomesOne)
      {
        classifierParams <- .classifierKeywordToParams(classifier, NULL)
        modellingParams <- ModellingParams(selectParams = SelectParams("previousSelection", intermediate = ".iteration", classifyResult = trainedModel),
                                           trainParams = classifierParams$trainParams,
                                           predictParams = classifierParams$predictParams)
        
        result <- runTests(measurementsOne, outcomesOne, crossValParams, modellingParams)
        mean(performance(calcCVperformance(result, performanceType))[[performanceType]])
      }, measurements, outcomes, SIMPLIFY = FALSE)
     })
    
    realPerformance <- matrix(unlist(performanceAllPairs), ncol = length(measurements), byrow = TRUE,
                              dimnames = list(paste("Select", names(measurements)), paste("Cross-validate", names(measurements))))
    realPerformance <- round(realPerformance, 2)
  }

  # Return matrix of pair-wise model accuracy. 
  # I've made this a list so that I can add things to it later on. 
  result <- list(real = realPerformance)
  # We want to include a set of nFeatures to compare between our feature selection method.
  if(doRandomFeatures == TRUE){
    message("Starting random feature selection procedure.")
    # Sample nFeatures randomly from each dataset.
    randomFeatures <- lapply(measurements, function(dataset) sample(colnames(dataset), nFeatures))
    performanceAllPairs <- lapply(randomFeatures, function(randomFeaturesSet)
    {
      mapply(function(testData, testOutcomes)
      {
        result <- crossValidate(testData[, randomFeaturesSet], testOutcomes,
                    nFeatures = nFeatures,
                    selectionMethod = "none",
                    classifier = classifier,
                    multiViewMethod = "none",
                    nFolds = nFolds,
                    nCores = nCores,
                    nRepeats = nRepeats)
        mean(performance(calcCVperformance(result, performanceType))[[performanceType]])
      }, measurements, outcomes)
    })

    randomPerformance <- matrix(unlist(performanceAllPairs), ncol = length(measurements), byrow = TRUE,
                              dimnames = list(paste("Random Select", names(measurements)), paste("Cross-validate", names(measurements))))
    randomPerformance <- round(randomPerformance, 2)
    
    result$random <- randomPerformance
  }
  
# Add information about the params to the output.
  result$params <- list(nFeatures = nFeatures, selectionMethod = selectionMethod,
                     selectionOptimisation = selectionOptimisation,
                     classifier = classifier, nFolds = nFolds, nRepeats = nRepeats, nCores = nCores,
                     trainType = trainType, performanceType = performanceType,
                     doRandomFeatures = doRandomFeatures)
  
  result
}

#' A function to plot the output of the crissCrossValidate function.
#'
#' This function has been designed to give a heatmap output of the crissCrossValidate function.
#'
#' @param crissCrossResult The output of the crissCrossValidate function.
#' @param includeValues If TRUE, then the values of the matrix will be included in the plot.
#' @author Harry Robertson
#' 
#' @import ggplot2
#' @import reshape2
#' @import ggpubr
#' 
#' @export

crissCrossPlot <- function(crissCrossResult, includeValues = FALSE){

  attach(crissCrossResult)
  scalebar_title <- params$performanceType

  # If the user does not want to compare features.
  if(params$trainType == "modelTrain"){
    melted_cormat <- reshape2::melt(real, na.rm = TRUE)
  
    ggheatmap <- ggplot(melted_cormat, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(high = "#e25563ff", mid = "white", low = "#094bacff", 
                          midpoint = 0.5, limit = c(0,1), space = "Lab", 
                          name=as.character(scalebar_title)) +
      theme_bw() + xlab("Training Dataset") + ylab("Testing Dataset") +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) + 
      theme(axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)) +
      coord_fixed() 
    
    if(includeValues == TRUE) ggheatmap <- ggheatmap + geom_text(aes(label = value), color = "black", size = 3)
  }
  
  else if(params$trainType == "modelTest"){
    melted_cormat_1 <- melt(real, na.rm = TRUE)
    ggheatmap_1 <- ggplot(melted_cormat_1, aes(Var1, Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(high = "#e25563ff", mid = "white", low = "#094bacff", 
                          midpoint = 0.5, limit = c(0,1), space = "Lab", 
                          name=as.character(scalebar_title)) +
      theme_bw() + xlab("Features Extracted") + ylab("Dataset Tested") +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) + 
      theme(axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)) +
      coord_fixed()
    if(includeValues == TRUE) ggheatmap_1 <- ggheatmap_1 + geom_text(aes(label = value), color = "black", size = 3)
    
    if(params$doRandomFeatures == TRUE){
      melted_cormat_2 <- melt(random, na.rm = TRUE)
      ggheatmap_2 <- ggplot(melted_cormat_2, aes(Var1, Var2, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(high = "#e25563ff", mid = "white", low = "#094bacff", 
                            midpoint = 0.5, limit = c(0,1), space = "Lab", 
                            name=as.character(scalebar_title)) +
        theme_bw() + xlab("Features Extracted") + ylab("Dataset Tested") +
        theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) + 
        theme(axis.text.y = element_text(vjust = 1, size = 8, hjust = 1)) +
        coord_fixed()
      if(includeValues == TRUE) ggheatmap_2 <- ggheatmap_2 + geom_text(aes(label = value), color = "black", size = 3)

        ggheatmap <- ggarrange(ggheatmap_1, ggheatmap_2, labels = c("A - Feature Selection", "B - Random Features"), 
                           ncol = 2, common.legend = TRUE, legend = "right")
    } else {
      ggheatmap <- ggheatmap_1
    }
  }
  print(ggheatmap)
}