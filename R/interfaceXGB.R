# An Interface for xgboost Package's xgb.train Function
extremeGradientBoostingTrainInterface <- function(measurementsTrain, outcomeTrain, mTryProportion = 0.5, nrounds = 100, ..., verbose = 3)
{
  if(!requireNamespace("xgboost", quietly = TRUE))
    stop("The package 'xgboost' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting extreme gradient boosting classifier to training data and making predictions on test
            data.")
  # Convert to one-hot encoding as xgboost doesn't understand factors. A sparse matrix, so that zeros are
  # treated as missing values by xgboost.
  measurementsTrain <- .encodeTrain(measurementsTrain)
  trainMatrix <- as(measurementsTrain, "CsparseMatrix")
  
  # Other arguments are booster parameters, unless they are arguments of xgb.train or xgb.DMatrix.
  extras <- list(...)
  trainArguments <- intersect(names(extras), setdiff(names(formals(xgboost::xgb.train)), c("params", "data", "nrounds", "verbose", "...")))
  dataArguments <- intersect(names(extras), c("missing", "weight"))
  params <- list(colsample_bynode = mTryProportion, nthread = 1)
  isParameter <- setdiff(names(extras), c(trainArguments, dataArguments))
  params[isParameter] <- extras[isParameter] # User-specified parameters, such as nthread, take precedence.
  
  isClassification <- FALSE
  if(is(outcomeTrain, "Surv")) # xgboost only knows about numeric vectors.
  {
    time <- outcomeTrain[, "time"]
    event <- as.numeric(outcomeTrain[, "status"])
    if(max(event) == 2) event <- event - 1
    outcomeTrain <- time * ifelse(event == 1, 1, -1) # Negative for censoring.
    params[["objective"]] <- "survival:cox"
  } else { # Classification task.
    isClassification <- TRUE
    classes <- levels(outcomeTrain)
    params[["objective"]] <- "multi:softprob"
    params[["num_class"]] <- length(classes)
    outcomeTrain <- as.numeric(outcomeTrain) - 1 # Classes are represented as 0, 1, 2, ...
  }
  trainData <- do.call(xgboost::xgb.DMatrix, c(list(trainMatrix, label = outcomeTrain, nthread = params[["nthread"]]),
                                               extras[dataArguments]))
  trained <- do.call(xgboost::xgb.train, c(list(params = params, data = trainData, nrounds = nrounds, verbose = 0),
                                           extras[trainArguments]))

  if(isClassification)
    attr(trained, "classes") <- classes # Useful for factor predictions in predict method.
  attr(trained, "featureNames") <- colnames(measurementsTrain)
  attr(trained, "featureGroups") <- attr(measurementsTrain, "assign")
  attr(trained, "encoding") <- attr(measurementsTrain, "encoding")

  trained
}
attr(extremeGradientBoostingTrainInterface, "name") <- "extremeGradientBoostingTrainInterface"
    
# booster is of class xgb.Booster
extremeGradientBoostingPredictInterface <- function(booster, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using boosted random forest.")  
  # Same one-hot encoding, columns and column order as the training data.
  measurementsTest <- .encodeTest(measurementsTest, booster)
  testMatrix <- as(measurementsTest, "CsparseMatrix")
  
  scores <- predict(booster, testMatrix)
  if(!is.null(attr(booster, "classes"))) # It is a classification task.
  {
    scores <- matrix(scores, nrow = nrow(measurementsTest), dimnames = list(rownames(measurementsTest), attr(booster, "classes")))
    classPredictions <- attr(booster, "classes")[apply(scores, 1, function(sampleRow) which.max(sampleRow)[1])]
    classPredictions <- factor(classPredictions, levels = attr(booster, "classes"))
    names(classPredictions) <- rownames(measurementsTest)
    result <- switch(returnType, class = classPredictions,
                     score = scores,
                     both = data.frame(class = classPredictions, scores, check.names = FALSE))
  } else { # A survival task.
     result <- setNames(as.numeric(scores), rownames(measurementsTest))
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
