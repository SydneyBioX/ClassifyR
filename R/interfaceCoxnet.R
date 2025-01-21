# An Interface for glmnet Package's coxnet Function. Survival modelling with sparsity.

coxnetTrainInterface <- function(measurementsTrain, survivalTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message(Sys.time(), ": Fitting coxnet model to data.")
    
  measurementsTrain <- data.frame(measurementsTrain, check.names = FALSE)
  measurementsMatrix <- MatrixModels::model.Matrix(~ 0 + ., data = measurementsTrain)
  
  # The response variable is a Surv class of object.
  fit <- glmnet::cv.glmnet(measurementsMatrix, survivalTrain, family = "cox", type = "C", ...)
  fitted <- fit$glmnet.fit
  
  offset <- -mean(predict(fitted, measurementsMatrix, s = fit$lambda.min, type = "link"))
  attr(fitted, "tune") <- list(lambda = fit$lambda.min, offset = offset)
  attr(fitted, "featureNames") <- colnames(measurementsMatrix)
  attr(fitted, "featureGroups") <- measurementsMatrix@assign
  
  class(fitted) <- class(fitted)[1] # Get rid of glmnet which messes with dispatch. 
  fitted
}
attr(coxnetTrainInterface, "name") <- "coxnetTrainInterface"

# model is of class coxnet.
coxnetPredictInterface <- function(model, measurementsTest, survivalTest = NULL, lambda, ..., verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using cox model.")
  
  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.
  
  testMatrix <- glmnet::makeX(as(measurementsTest, "data.frame"))
  testMatrix <- testMatrix[, rownames(model[["beta"]])]
  
  offset <- attr(model, "tune")[["offset"]]
  model$offset <- TRUE
  
  survScores <- predict(model, testMatrix, s = lambda, type = "response", newoffset = offset)
  rownames(survScores) <- rownames(measurementsTest)
  survScores[, 1]
}