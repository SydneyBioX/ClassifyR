# An Interface for ranger Package's randomForest Function
randomForestTrainInterface <- function(measurementsTrain, outcomeTrain, mTryProportion = NULL, ..., verbose = 3)
{
  if(!requireNamespace("ranger", quietly = TRUE))
    stop("The package 'ranger' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting random forest classifier to training data.")
  # Number of features to try at each split. By default, ranger's own: the square root of the number of features.
  mtry <- if(!is.null(mTryProportion)) round(mTryProportion * ncol(measurementsTrain))
  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  measurementsTrain <- .asDataFrame(measurementsTrain) # ranger needs a data.frame.
  # Features are ranked by the impurity importance of the forest itself, unless another importance mode is given.
  if(is.null(list(...)[["importance"]]))
    ranger::ranger(x = measurementsTrain, y = outcomeTrain, mtry = mtry, importance = "impurity", ...)
  else
    ranger::ranger(x = measurementsTrain, y = outcomeTrain, mtry = mtry, ...)
}
attr(randomForestTrainInterface, "name") <- "randomForestTrainInterface"
    
# forest is of class ranger
randomForestPredictInterface <- function(forest, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("ranger", quietly = TRUE))
    stop("The package 'ranger' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  classes <- forest$forest$levels
  if(verbose == 3)
    message("Predicting using random forest.")  
  measurementsTest <- .asDataFrame(measurementsTest)
  
  predictions <- predict(forest, measurementsTest)
  if(predictions$treetype == "Classification")
  {
    classPredictions <- predictions$predictions
    classScores <- predict(forest, measurementsTest, predict.all = TRUE)[[1]]
    # Share of trees voting for each class.
    classScores <- vapply(seq_along(classes), function(classIndex) rowSums(classScores == classIndex), numeric(nrow(classScores))) / forest$forest$num.trees
    if(!is.matrix(classScores)) classScores <- matrix(classScores, nrow = 1)
    colnames(classScores) <- classes
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
                    # Models made by earlier versions keep a second forest, grown only to rank features.
                    forImportance <- attr(forest, "forImportance")
                    if(is.null(forImportance)) forImportance <- forest
                    rankedFeaturesIndices <- order(ranger::importance(forImportance), decreasing = TRUE)
                    selectedFeaturesIndices <- which(ranger::importance(forImportance) > 0)
                    list(rankedFeaturesIndices, selectedFeaturesIndices)
                  }