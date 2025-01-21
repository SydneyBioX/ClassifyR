# An Interface for randomForestSRC Package's rfsrc random forest survival Function
rfsrcTrainInterface <- function(measurementsTrain, survivalTrain, mTryProportion = 0.5, ..., verbose = 3)
{
  if(!requireNamespace("randomForestSRC", quietly = TRUE))
    stop("The package 'randomForestSRC' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting rfsrc classifier to training data and making predictions on test data.")

  # Surv objects store survival information as a two-column table, time and event, in that order.    
  bindedMeasurements <- cbind(measurementsTrain, time = survivalTrain[, 1], event = survivalTrain[, 2])
  mtry <- round(mTryProportion * ncol(measurementsTrain)) # Number of features to try.
  randomForestSRC::rfsrc(Surv(time, event) ~ ., data = as.data.frame(bindedMeasurements), mtry = mtry,
                          var.used = "all.trees", importance = TRUE, ...)
}
attr(rfsrcTrainInterface, "name") <- "rfsrcTrainInterface"

# model is of class rfsrc
rfsrcPredictInterface <- function(model, measurementsTest, ..., verbose = 3)
{
  predictedOutcome = predict(model, as.data.frame(measurementsTest), ...)$predicted
  names(predictedOutcome) = rownames(measurementsTest)
  predictedOutcome
}

rfsrcFeatures <- function(forest)
                  {
                    rankedFeaturesIndices <- order(forest[["importance"]], decreasing = TRUE)
                    selectedFeaturesIndices <- which(forest[["var.used"]] > 0)
                    list(rankedFeaturesIndices, selectedFeaturesIndices)
                  }