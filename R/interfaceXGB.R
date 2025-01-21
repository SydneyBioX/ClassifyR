# An Interface for xgboost Package's xgboost Function
extremeGradientBoostingTrainInterface <- function(measurementsTrain, outcomeTrain, mTryProportion = 0.5, nrounds = 10, ..., verbose = 3)
{
  if(!requireNamespace("xgboost", quietly = TRUE))
    stop("The package 'xgboost' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting extreme gradient boosting classifier to training data and making predictions on test
            data.")
  measurementsTrain <- as(measurementsTrain, "data.frame")
  # Convert to one-hot encoding as xgboost doesn't understand factors. Need to get rid of intercept afterwards.
  measurementsTrain <- MatrixModels::model.Matrix(~ 0 + ., data = measurementsTrain, sparse = TRUE)
  
  isClassification <- FALSE
  numClasses <- NULL
  if(is(outcomeTrain, "Surv")) # xgboost only knows about numeric vectors.
  {
    time <- outcomeTrain[, "time"]
    event <- as.numeric(outcomeTrain[, "status"])
    if(max(event) == 2) event <- event - 1
    outcomeTrain <- time * ifelse(event == 1, 1, -1) # Negative for censoring.
    objective <- "survival:cox"
    trained <- xgboost::xgboost(measurementsTrain, outcomeTrain, objective = objective, nrounds = nrounds,
                                colsample_bynode = mTryProportion, verbose = 0, nthread = 1, ...)
  } else { # Classification task.
    isClassification <- TRUE
    classes <- levels(outcomeTrain)
    numClasses <- length(classes)
    objective <- "multi:softprob"
    outcomeTrain <- as.numeric(outcomeTrain) - 1 # Classes are represented as 0, 1, 2, ...
    trained <- xgboost::xgboost(measurementsTrain, outcomeTrain, objective = objective, nrounds = nrounds,
                              num_class = numClasses, colsample_bynode = mTryProportion, nthread = 1, verbose = 0, ...)
  }

  if(isClassification)
    attr(trained, "classes") <- classes # Useful for factor predictions in predict method.
  attr(trained, "featureNames") <- colnames(measurementsTrain)
  attr(trained, "featureGroups") <- measurementsTrain@assign

  trained
}
attr(extremeGradientBoostingTrainInterface, "name") <- "extremeGradientBoostingTrainInterface"
    
# booster is of class xgb.Booster
extremeGradientBoostingPredictInterface <- function(booster, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using boosted random forest.")  
  measurementsTest <- as(measurementsTest, "data.frame")
  # Convert to one-hot encoding as xgboost doesn't understand factors. Need to get rid of intercept afterwards.
  measurementsTest <- MatrixModels::model.Matrix(~ 0 + ., data = measurementsTest, sparse = TRUE)
  measurementsTest <- measurementsTest[, attr(booster, "featureNames")]
  
  scores <- predict(booster, measurementsTest, reshape = TRUE)
  colnames(scores) <- attr(booster, "classes")
  if(!is.null(attr(booster, "classes"))) # It is a classification task.
  {
    classPredictions <- attr(booster, "classes")[apply(scores, 1, function(sampleRow) which.max(sampleRow)[1])]
    classPredictions <- factor(classPredictions, levels = attr(booster, "classes"))
    rownames(scores) <- names(classPredictions) <- rownames(measurementsTest)
    result <- switch(returnType, class = classPredictions,
                     score = scores,
                     both = data.frame(class = classPredictions, scores, check.names = FALSE))
  } else { # A survival task.
     result <- scores
  }
  result
}

################################################################################
#
# Get selected features
#
################################################################################

XGBfeatures <- function(booster)
                  {
                    importanceGains <- xgboost::xgb.importance(model = booster)[["Gain"]]
                    gains <- rep(0, length(unique(attr(booster, "featureGroups"))))
                    featureGroups <- attr(booster, "featureGroups")[match(xgboost::xgb.importance(model = booster)[["Feature"]], attr(booster, "featureNames"))]
                    maxGains <- by(importanceGains, featureGroups, max)
                    indicesUsed <- as.numeric(names(maxGains))
                    gains[indicesUsed]  <- maxGains # Put into particular indexes.
                    rankedFeaturesIndices <- order(gains, decreasing = TRUE)
                    selectedFeaturesIndices <- indicesUsed
                    list(rankedFeaturesIndices, selectedFeaturesIndices)
                  }
