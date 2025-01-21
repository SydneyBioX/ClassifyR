# An Interface for ranger Package's randomForest Function
randomForestTrainInterface <- function(measurementsTrain, outcomeTrain, mTryProportion = 0.5, ..., verbose = 3)
{
  if(!requireNamespace("ranger", quietly = TRUE))
    stop("The package 'ranger' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting random forest classifier to training data.")
  mtry <- round(mTryProportion * ncol(measurementsTrain)) # Number of features to try.
  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  fittedModel <- ranger::ranger(x = as(measurementsTrain, "data.frame"), y = outcomeTrain, mtry = mtry, ...)
  forImportance <- ranger::ranger(x = as(measurementsTrain, "data.frame"), y = outcomeTrain, mtry = mtry, importance = "impurity_corrected", ...)
  attr(fittedModel, "forImportance") <- forImportance
  fittedModel
}
attr(randomForestTrainInterface, "name") <- "randomForestTrainInterface"
    
# forest is of class ranger
randomForestPredictInterface <- function(forest, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  classes <- forest$forest$levels
  if(verbose == 3)
    message("Predicting using random forest.")  
  measurementsTest <- as.data.frame(measurementsTest)
  
  predictions <- predict(forest, measurementsTest)
  if(predictions$treetype == "Classification")
  {
    classPredictions <- predictions$predictions
    classScores <- predict(forest, measurementsTest, predict.all = TRUE)[[1]]
    classScores <- t(apply(classScores, 1, function(sampleRow) table(factor(classes[sampleRow], levels = classes)) / forest$forest$num.trees))
    rownames(classScores) <- names(classPredictions) <- rownames(measurementsTest)
    switch(returnType, class = classPredictions,
           score = classScores,
           both = data.frame(class = classPredictions, classScores, check.names = FALSE))
  } else { # It is "Survival".
      -1 * rowSums(predictions$survival) # Make it a risk score.
  }
}

################################################################################
#
# Get selected features
#
################################################################################

forestFeatures <- function(forest)
                  {
                    forImportance <- attr(forest, "forImportance")
                    rankedFeaturesIndices <- order(ranger::importance(forImportance), decreasing = TRUE)
                    selectedFeaturesIndices <- which(ranger::importance(forImportance) > 0)
                    list(rankedFeaturesIndices, selectedFeaturesIndices)
                  }